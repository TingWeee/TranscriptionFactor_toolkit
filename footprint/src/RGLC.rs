// Rapid Genome landscape Characterisation at Cell Level 
// bed file format: chr, start, end, strand 

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


pub fn rglc(matches: &clap::ArgMatches) -> Result<(), Box<dyn Error>> {

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

    fcount(&frag_file, &bed_file, &cell_file, output_path, num_threads)?;
    
    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    output: &Path,
    num_threads: usize,

) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );
    
    // bed intervals 
    let mut lapper_map = match feature_intervals(bed_file, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };
    
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
    // result with strand correction 
    let mut aggregated_counts: FxHashMap<(String, String), Vec<usize>> = FxHashMap::default();
    // result without strand correction 
    let mut counts_nostrand: FxHashMap<(String, String), Vec<usize>> = FxHashMap::default();
    // 1D array for fragment counts per barcode 
    let mut frag_per_cell: Array1<f32> = Array1::zeros(cells.keys().len());
    // frag file reading
    let frag_file = File::open(frag_file)?;
    let mut reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(frag_file));

    let mut line_count: u64 = 0;
    let update_interval = 1_000_000;
    let mut line_str = String::new();
    let mut startpos: usize;
    let mut endpos: usize;

    let mut current_chrom = String::new();
    let mut current_lapper: Option<&mut Lapper<usize, String>> = None;
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
                    let strand = &interval.val; 
                    let tf_id = "t"; // lazy change data structure 
                    
                    let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![0;2000]);
                    let ns_counts = counts_nostrand.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![0;2000]);
                    let relative_position: usize  = startpos - interval_start; 

                    if strand == "+" {
                        insertion_counts[relative_position] += 1; 
                    } else if strand == "-" {
                        insertion_counts[2000 - 1 - relative_position] += 1; 
                    }
                     
                    ns_counts[relative_position] += 1;
                    
                    // Check if fragment end also within interval 
                    if endpos < interval_end {
                        check_end = false;
                        let relative_position: usize  = endpos - interval_start; 
                        //////////////////////////////////////////////////////////////////
    
                        if strand == "+" {
                            insertion_counts[relative_position] += 1; 
                        } else if strand == "-" {
                            insertion_counts[2000 - 1 - relative_position] += 1; 
                        }
                        ns_counts[relative_position] += 1;
                    }
                }
                    
                if check_end {
                    
                    for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                        let interval_start = interval.start;
                        let interval_end = interval.stop;
                        let strand = &interval.val; 
                        let tf_id = "t";
                        let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.to_string()))
                            .or_insert_with(|| vec![0;2000]);
                        let ns_counts = counts_nostrand.entry((cell_barcode.to_string(), tf_id.to_string()))
                            .or_insert_with(|| vec![0;2000]);
                        
                        let relative_position: usize = endpos - interval_start; 
                       
                        if strand == "+" {
                            insertion_counts[relative_position] += 1; 
                        } else if strand == "-" {
                            insertion_counts[2000 - 1 - relative_position] += 1; 
                        }
                        ns_counts[relative_position] += 1;
                    }
                }       
            }
        }
        line_str.clear();
    }
    eprintln!();
    //println!("Aggregated counts length: {}", aggregated_counts.len());
    let output_path = output.join("aggregated_counts_strand.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in aggregated_counts {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}", cell, tf, counts_str).expect("Failed to write to file");
    }

    //println!("Aggregated bias length: {}", aggregated_bias.len());
    let output_path = output.join("aggregated_counts_nostrand.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in counts_nostrand {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}\n", cell, tf, counts_str).expect("Failed to write to file");
    }

    let mut writer = File::create(output.join("fragment_counts.tsv"))?;
    let mut foutput = String::new();
    for (cell_barcode, &cell_index) in &cells {
        foutput.push_str(&format!("{}\t{}\n", cell_barcode, frag_per_cell[cell_index as usize]));
    }
    writer.write_all(foutput.as_bytes())?;

    Ok(())
}


// Reads BED file and organizes intervals 
fn feature_intervals(
    motif_match_file: &Path, 
    num_threads: usize,
) -> io::Result<FxHashMap<String, Lapper<usize, String >>> {

    // bed file reader
    let file = File::open(motif_match_file)?;
    let reader = BufReader::new(file);
    // hashmap of features for each chromosome  chr: [Interval(start, end, strand}]
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<usize, String>>> = FxHashMap::default();
    
    for (index, line) in reader.lines().enumerate() {
        
        let line = line.expect("Error reading line");
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue; 
        }
        let chr = fields[0].to_string();
        let start: usize = fields[1].parse().expect(&format!("Invalid start position {}", index + 1));
        let end: usize = fields[2].parse().expect(&format!("Invalid end position {}", index + 1)); // non inclusive 
        let strand = fields[3].to_string();
        
        let mid: usize = start + ((end-start)/2); // floor()
        let intervals = chromosome_trees.entry(chr.clone()).or_insert_with(Vec::new);
        intervals.push(Interval {start: mid - 1000 , stop: mid + 1000 , val: strand}); // interval is 2001 in length 
    }

    let lapper_map = chromosome_trees.into_iter()
        .map(|(chr, intervals)| {
            println!("Chromosome: {}, Number of intervals: {}", chr, intervals.len());
            (chr, Lapper::new(intervals))
        })
        .collect();
    
    Ok(lapper_map)
}
