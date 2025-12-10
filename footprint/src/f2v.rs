// keeps 500 bp from center  of motif region (sums up every postiion for each cell)
// do the same but for bias 
// then we can do the scatter plot in R or in python to see the distribution across cells 
// next we add in the weight, so we can see the distribution of peak_scores and understand why it is not changing 
use std::{
    io,
    fs,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
    io::Write,
};

use ndarray::{Array1};
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use flate2::Compression;
use log::error;
use log::info;
use log::warn;
use rustc_hash::FxHashMap;
use gzp::{
    deflate::Gzip,
    ZWriter,
    par::compress::{ParCompress, ParCompressBuilder},
};

use crate::bias::{read_bias_factors, cal_sum_bias, cal_mean_bias, cal_max_bias, calculate_mean, calculate_std}; 

pub fn f2v(matches: &clap::ArgMatches) -> Result<(), Box<dyn Error>> {

    let frag_file = Path::new(matches.get_one::<String>("fragments").unwrap())
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    let bed_file = Path::new(matches.get_one::<String>("bed").unwrap())
        .canonicalize()
        .expect("Can't find path to input BED file");
    info!("Received BED file: {:?}", bed_file);

    let cell_file = Path::new(matches.get_one::<String>("cells").unwrap())
        .canonicalize()
        .expect("Can't find path to input cell file");
    info!("Received cell file: {:?}", cell_file);

    let output_directory = matches.get_one::<String>("outdir").unwrap();
    info!("Received output directory: {:?}", output_directory);

    let group = matches.get_one::<usize>("group").copied();
    info!("Grouping peaks: {:?}", group);
    
    let output_path = Path::new(output_directory);

    let num_threads = *matches.get_one::<usize>("threads").unwrap();

    // Create the directory if it does not exist
    if !output_path.exists() {
        if let Err(e) = fs::create_dir_all(output_path) {
            eprintln!("Failed to create output directory: {}", e);
            std::process::exit(1);
        }
    }
    
    // make sure output is a directory
    match fs::metadata(output_path) {
        Ok(metadata) => {
            if metadata.is_dir() {
                info!("{:?} is a directory.", output_path);
            } else {
                eprintln!("Provided output is not a directory: {}", output_path.display());
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Failed to get metadata for {:?}: {}", output_path, e);
            std::process::exit(1);
        }
    }

    fcount(&frag_file, &bed_file, &cell_file, output_path, group, num_threads)?;
    
    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    output: &Path,
    group: Option<usize>,
    num_threads: usize,

) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );
  
    // let h5_path = "/home/users/astar/gis/limtw/scratch/proj_TF/motif/Tn5Bias/hg38Tn5Bias.h5";
    let h5_path = "/home/users/astar/gis/limtw/scratch/proj_tfa/hg38Tn5Bias.h5";
    let bias_factors = read_bias_factors(h5_path);

    let (tf_vector, mut lapper_map) = match feature_intervals(bed_file, group, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };
    
    println!("number of tf:{}", tf_vector.len());
    println!("lapper_map: {}", lapper_map.len());
    // create hashmap for cell barcodes {cell barcode: cell index}
    let cellreader = File::open(cell_file)
        .map(BufReader::new)?;
    // cell barcode: cell index 
    let mut cells: FxHashMap<String, u32> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        let index_u32 = index as u32;
        cells.insert(line, index_u32);
    }

    // create 2 files: 1 for actual Tn5 insertions, one for bias insertion 
    // column 1: cell barcode, column 2: TF feature, column 3-503: positional insertions 
    let mut aggregated_counts: FxHashMap<(String, String), Vec<usize>> = FxHashMap::default();
    let mut aggregated_bias: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    // column 1: cell barcode, column 2: TF feature,  column 3: total number of insertions, column 3-500+: culmulative positional bias -> to be divided by total number of insertions?? 
    // 1D array for fragment counts per barcode 
    let mut frag_per_cell: Array1<f32> = Array1::zeros(cells.keys().len());

    // Additional 3 files: 1. Upstream & Downstream Bias Distribution (Mean & Std) 
                        // 2. Bias Corrected Count Vector 
    let mut bias_fraction_upstream: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    let mut bias_fraction_downstream: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    // 1 * (1 / positional_bias) 
    let mut corrected_counts_one: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    // 1 * (avg_bias / positional_bias) 
    let mut corrected_counts_two: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    
    // frag file reading
    let frag_file = File::open(frag_file)?;
    let mut reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(frag_file));

    let mut line_count: u64 = 0;
    let update_interval = 1_000_000;
    let mut line_str = String::new();
    let mut startpos: usize;
    let mut endpos: usize;

    let mut current_chrom = String::new();
    let mut current_lapper: Option<&mut Lapper<usize, usize>> = None;
    let mut cursor: usize = 0;
    let mut check_end: bool;


    loop { // fragment file 
        match reader.read_line(&mut line_str) {
            Ok(0) => break,
            Ok(_) => {},
            Err(e) => {
                error!("Error reading fragment file: {}", e);
                return Err(e);
            }
        }
        let line = &line_str[..line_str.len() - 1];

        // # Skip header lines that start with #
        if line.starts_with('#') {
            line_str.clear();
            continue;
        }

        line_count += 1;
        if line_count % update_interval == 0 {
            print!("\rProcessed {} M fragments", line_count / 1_000_000);
            std::io::stdout().flush().expect("Can't flush output");
        }

        // Parse BED entry
        let fields: Vec<&str> = line.split('\t').collect();

        // Check if cell is to be included
        let cell_barcode: &str = fields[3];
        
        if let Some(&cell_index) = cells.get(cell_barcode) {
            // println!("{}", cell_barcode);
            check_end = true;
            
            frag_per_cell[cell_index as usize] += 1.0;
            
            let seqname: &str = fields[0];

            if seqname != current_chrom {
                current_chrom = seqname.to_string();
                current_lapper = lapper_map.get_mut(&current_chrom);
                cursor = 0;
            }

            // try to parse the coordinates, skip the line if parsing fails
            startpos = match fields[1].trim().parse() {
                Ok(num) => num,
                Err(e) => {
                    warn!("Failed to parse start position: {:?}. Error: {}", line_count, e);
                    line_str.clear();
                    continue;
                }
            };
            
            endpos = match fields[2].trim().parse() {
                Ok(num) => num,
                Err(e) => {
                    warn!("Failed to parse end position: {:?}. Error: {}", line_count, e);
                    line_str.clear();
                    continue;
                }
            };
            
        

            if let Some(lapper) = &mut current_lapper {
                // seems to be a problem with seek if lapper has one element
                // set cursor to 0
                if lapper.intervals.len() == 1 {
                    cursor = 0;
                }
            
                for interval in lapper.seek(startpos as usize, startpos as usize + 1, &mut cursor) {
                    
                    let interval_start = interval.start;
                    let interval_end = interval.stop;
                    // 
                    let (tf_id, upstream_frac, downstream_frac, strand) = &tf_vector[interval.val];

                    
                    let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.clone()))
                        .or_insert_with(|| vec![0;2000]);
                    let insertion_corrected_one = corrected_counts_one.entry((cell_barcode.to_string(), tf_id.clone()))
                        .or_insert_with(|| vec![0.0;2000]);
                    let insertion_corrected_two = corrected_counts_two.entry((cell_barcode.to_string(), tf_id.clone()))
                        .or_insert_with(|| vec![0.0;2000]);
                    
                    let relative_position: usize  = startpos - interval_start; 

                    let mut positional_bias: f32 = 0.167; 

                    if let Some(bias) = &bias_factors.get(&current_chrom) {
                        if let Some(&value) = bias.get(startpos) {
                            positional_bias = value;
                        }
                    }
                    
                    let avg_bias: f32 = cal_mean_bias(&current_chrom, startpos-50, startpos + 50, &bias_factors);

                    
                    
                    if strand == "+" {
                        insertion_counts[relative_position] += 1; 
                        insertion_corrected_one[relative_position] += (1.0/ positional_bias) ; 
                        insertion_corrected_two[relative_position] += (avg_bias / positional_bias) ;
                    } else if strand == "-" {
                        insertion_counts[1999 - relative_position] += 1; 
                        insertion_corrected_one[1999 - relative_position] += (1.0/ positional_bias) ; 
                        insertion_corrected_two[1999 - relative_position] += (avg_bias / positional_bias) ;
                    }
                     
                    
                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    // distal 0-50 , flanking 150 - 250, flanking 250 - 350, 
                    // only insert into bias vector if relative position is between 150 - 350 
                    let bias_landscape = aggregated_bias.entry((cell_barcode.to_string(), tf_id.clone()))
                            .or_insert_with(|| vec![0.0;1001]);
                    let bias_frac_up = bias_fraction_upstream.entry((cell_barcode.to_string(), tf_id.clone())).or_insert_with(Vec::new);
                    let bias_frac_dw = bias_fraction_downstream.entry((cell_barcode.to_string(), tf_id.clone())).or_insert_with(Vec::new);

                    // also only add to the motif bias fraction 
                    // use if let which is similar to "match 
                    if let 900..=1100 = relative_position {
                       
                        // Update Bias Vector 
                        if let Some(bias) = bias_factors.get(&current_chrom) {
                            bias_landscape[0] += 1.0;
                            if strand == "+" {
                                bias_landscape[1..=1000].iter_mut()
                                    .zip(bias[interval_start+500..interval_start+1500].iter())
                                    .for_each(|(a, b)| *a += b);
                            } else if strand == "-" {
                                bias_landscape[1..=1000].iter_mut()
                                    .zip(bias[interval_start+500..interval_start+1500].iter().rev()) // reverse iterator 
                                    .for_each(|(a, b)| *a += b);  
                            }
                        } else {
                            println!("No bias data for {}", current_chrom);
                        }
                        // increment ocunt of insertions at bias_landscape[0]
                        // efficient modification of bias_landscape in-place with no extra memory allocation
                        // Push Bias_fraction_up and down
                        if relative_position < 1000{
                            if strand == "+" {
                                bias_frac_up.push(*upstream_frac);
                            } else {
                                bias_frac_dw.push(*downstream_frac);
                            }
                        } else {
                            if strand == "-" {
                                bias_frac_up.push(*upstream_frac);
                            } else {
                                bias_frac_dw.push(*downstream_frac);
                            }
                        }
                    }
                    // Check if fragment end also within interval 
                    if endpos < interval_end {
                        check_end = false;

                        let relative_position: usize  = endpos - interval_start; 
                        //////////////////////////////////////////////////////////////////
                        let positional_bias:f32 = cal_max_bias(&current_chrom, endpos-1, endpos+2, &bias_factors);
                        let avg_bias: f32 = cal_mean_bias(&current_chrom, endpos-50, endpos + 50, &bias_factors);
                    
                        if strand == "+" {
                            insertion_counts[relative_position] += 1; 
                            insertion_corrected_one[relative_position] += (1.0/ positional_bias) ; 
                            insertion_corrected_two[relative_position] += (avg_bias / positional_bias) ;
                        } else if strand == "-" {
                            insertion_counts[1999 - relative_position] += 1; 
                            insertion_corrected_one[1999 - relative_position] += (1.0/ positional_bias) ; 
                            insertion_corrected_two[1999 - relative_position] += (avg_bias / positional_bias) ;
                        }
                        if let 400..=600 = relative_position {
                        // Update Bias Vector 
                            if let Some(bias) = bias_factors.get(&current_chrom) {
                                bias_landscape[0] += 1.0;
                                if strand == "+" {
                                    bias_landscape[1..=1000].iter_mut()
                                        .zip(bias[interval_start+500..interval_start+1500].iter())
                                        .for_each(|(a, b)| *a += b);
                                }else if strand == "-" {
                                    bias_landscape[1..=1000].iter_mut()
                                        .zip(bias[interval_start+500..interval_start+1500].iter().rev()) // reverse iterator 
                                        .for_each(|(a, b)| *a += b);  
                                }
                            } else {
                                println!("No bias data for {}", current_chrom);
                            }

	                    // Push Bias_fraction_up and down
                            if relative_position < 1000{
                            	if strand == "+" {
                                   bias_frac_up.push(*upstream_frac);
                           	} else {
                                   bias_frac_dw.push(*downstream_frac);
                            	}
                            } else {
                            	if strand == "-" {
                                   bias_frac_up.push(*upstream_frac);
                            	} else {
                               	   bias_frac_dw.push(*downstream_frac);
                            	}
                            }

                        }
                    }
                }
                    
                if check_end {
                    for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                        let interval_start = interval.start;
                        let interval_end = interval.stop;
                        let (tf_id, upstream_frac, downstream_frac, strand) = &tf_vector[interval.val];

                        let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.clone()))
                            .or_insert_with(|| vec![0;2000]);
                        let insertion_corrected_one = corrected_counts_one.entry((cell_barcode.to_string(), tf_id.clone()))
                            .or_insert_with(|| vec![0.0;2000]);
                        let insertion_corrected_two = corrected_counts_two.entry((cell_barcode.to_string(), tf_id.clone()))
                            .or_insert_with(|| vec![0.0;2000]);

                        
                        let relative_position: usize = endpos - interval_start; 
                        //let positional_bias:f32 = cal_max_bias(&current_chrom, endpos-1, endpos+2, &bias_factors);
                        let mut positional_bias: f32 = 0.167; 

                        if let Some(bias) = &bias_factors.get(&current_chrom) {
                            if let Some(&value) = bias.get(endpos) {
                                positional_bias = value;
                            }
                        }
                        let avg_bias: f32 = cal_mean_bias(&current_chrom, endpos-50, endpos + 50, &bias_factors);

                        if strand == "+" {
                            insertion_counts[relative_position] += 1; 
                            insertion_corrected_one[relative_position] += (1.0/ positional_bias) ; 
                            insertion_corrected_two[relative_position] += (avg_bias / positional_bias) ;
                        } else if strand == "-" {
                            insertion_counts[1999 - relative_position] += 1; 
                            insertion_corrected_one[1999 - relative_position] += (1.0/ positional_bias) ; 
                            insertion_corrected_two[1999 - relative_position] += (avg_bias / positional_bias) ;
                        }
                        
                        let bias_landscape = aggregated_bias.entry((cell_barcode.to_string(), tf_id.clone()))
                            .or_insert_with(|| vec![0.0;1001]);
                        let bias_frac_up = bias_fraction_upstream.entry((cell_barcode.to_string(), tf_id.clone())).or_insert_with(Vec::new);
                        let bias_frac_dw = bias_fraction_downstream.entry((cell_barcode.to_string(), tf_id.clone())).or_insert_with(Vec::new);
                        
                        if let 900..=1100 = relative_position {
                            // Update Bias Vector 
                            if let Some(bias) = bias_factors.get(&current_chrom) {
                                bias_landscape[0] += 1.0;
                                if strand == "+" {
                                    bias_landscape[1..=1000].iter_mut()
                                        .zip(bias[interval_start+500..interval_start+1500].iter())
                                        .for_each(|(a, b)| *a += b);
                                }else if strand == "-" {
                                    bias_landscape[1..=1000].iter_mut()
                                        .zip(bias[interval_start+500..interval_start+1500].iter().rev()) // reverse iterator 
                                        .for_each(|(a, b)| *a += b);  
                                }
                            } else {
                                println!("No bias data for {}", current_chrom);
                            }
                            if relative_position < 1000{
                            	if strand == "+" {
                                   bias_frac_up.push(*upstream_frac);
                           	} else {
                                   bias_frac_dw.push(*downstream_frac);
                            	}
                            } else {
                            	if strand == "-" {
                                   bias_frac_up.push(*upstream_frac);
                            	} else {
                               	   bias_frac_dw.push(*downstream_frac);
                            	}
                            }
                        }
                    }
                }       
            }
        }
        line_str.clear();
    }
    eprintln!();
    //println!("Aggregated counts length: {}", aggregated_counts.len());
    let output_path = output.join("aggregated_counts.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in aggregated_counts {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}", cell, tf, counts_str).expect("Failed to write to file");
    }

    //println!("Aggregated bias length: {}", aggregated_bias.len());
    let output_path = output.join("aggregated_bias.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in aggregated_bias {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}\n", cell, tf, counts_str).expect("Failed to write to file");
    }

    //  let mut bias_fraction_upstream: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    let output_path = output.join("bias_fraction.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), frac_up_vec) in bias_fraction_upstream {
        let mean_up = calculate_mean(&frac_up_vec);
        let std_up = calculate_std(&frac_up_vec);
        
        let default_vec = vec![0.0;2]; // Vec::new();
        let frac_down_vec = bias_fraction_downstream.get(&(cell.clone(), tf.clone())).unwrap_or(&default_vec);
        let mean_down = calculate_mean(frac_down_vec);
        let std_down = calculate_std(frac_down_vec);
        
        
    
        writeln!(output_file, "{}\t{}\t{}\t{}\t{}\t{}\n", cell, tf, mean_up,std_up, mean_down,std_down).expect("Failed to write to file");
    }
    
    // 1 * (1 / positional_bias) 
    let output_path = output.join("corrected_counts_1.tsv");
    let mut output_file = File::create(output_path).expect("Failed  to create output file");
    for ((cell, tf), counts) in corrected_counts_one {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}", cell, tf, counts_str).expect("Failed to write to file");
    }
    
     // 1 * (avg_bias / positional_bias) 
    let output_path = output.join("corrected_counts_2.tsv");
    let mut output_file = File::create(output_path).expect("Failed  to create output file");
    for ((cell, tf), counts) in corrected_counts_two {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}", cell, tf, counts_str).expect("Failed to write to file");
    }

    let mut writer = File::create(output.join("fragment_counts.tsv"))?;
    let mut foutput = String::new();
    for (cell_barcode, &cell_index) in &cells {
        foutput.push_str(&format!("{}\t{}\n", cell_barcode, frag_per_cell[cell_index as usize]));
    }
    writer.write_all(foutput.as_bytes())?;

    Ok(())
}


// Reads BED file and organizes intervals :-> 500bp bin by chromosome 
fn feature_intervals(
    motif_match_file: &Path, 
    group: Option<usize>,
    num_threads: usize,
) -> io::Result<(Vec<(String, f32, f32, String)>, FxHashMap<String, Lapper<usize, usize>>)> {
    // Bias Factor 
    let h5_path = "/home/users/astar/gis/limtw/scratch/proj_tfa/hg38Tn5Bias.h5";
    let bias_factors = read_bias_factors(h5_path);

    // bed file reader
    let file = File::open(motif_match_file)?;
    let reader = BufReader::new(file);
    // hashmap of features for each chromosome  chr: [Interval(start, end , tf_id}] 
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<usize, usize>>> = FxHashMap::default();

    let mut tf_num: usize = 0; 

    // define the tf_vector[tf_index]: (motif_id, upstream bias fraction, downstream bias fraction, strand) 
    
    let mut tf_vector: Vec<(String, f32, f32, String)> = Vec::new();
    let mut vec_index: usize = 0; 
    
    for (index, line) in reader.lines().enumerate() {
        
        let line = line.expect("Error reading line");
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue; 
        }
        let chr = fields[0].to_string();
        let start: usize = fields[1].parse().expect(&format!("Invalid start position {}", index + 1));
        let end: usize = fields[2].parse().expect(&format!("Invalid end position {}", index + 1)); // non inclusive 
        
        let tf_id = fields[group.unwrap()].to_string();
    
        let strand = fields[5].to_string();
        
       // MOTIF BIAS 
        let motif_bias = cal_sum_bias(&chr, start, end, &bias_factors);
        // UPSTREAM BIAS 
        let upstream_bias = cal_sum_bias(&chr, start.saturating_sub(50), start, &bias_factors);
        let upstream_frac = motif_bias / (motif_bias + upstream_bias);
        // DOWNSTREAM BIAS 
        let downstream_bias =cal_sum_bias(&chr, end, end+50, &bias_factors);
        let downstream_frac = motif_bias /(motif_bias + downstream_bias); 
        
        let mid: usize = start + ((end-start)/2); // floor()

        let intervals = chromosome_trees.entry(chr.clone()).or_insert_with(Vec::new);

        intervals.push(Interval {start: mid - 1000 , stop: mid + 1000 , val: vec_index});
        
        if strand == "+" {
            tf_vector.push((tf_id, upstream_frac, downstream_frac, strand));
        } else if strand == "-" {
            tf_vector.push((tf_id, downstream_frac, upstream_frac, strand));
        }
        
        vec_index += 1;
        
    }

    let lapper_map = chromosome_trees.into_iter()
        .map(|(chr, intervals)| {
            println!("Chromosome: {}, Number of intervals: {}", chr, intervals.len());
            (chr, Lapper::new(intervals))
        })
        .collect();
    
    Ok((tf_vector lapper_map))
}
