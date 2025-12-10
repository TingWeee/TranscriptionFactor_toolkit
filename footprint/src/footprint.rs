// not sure why there is f2v.rs and RGLC.rs 
// collapse into 1 file and upload to github 
use std::{
    error::Error,fs, fs::File, io, io::BufRead, io::BufReader, io::Write, path::Path, sync::mpsc,
    thread, collections::HashSet,
};

use bio::io::fasta::IndexedReader;
use serde_json::json;
use ndarray::{Array1};
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use flate2::Compression;
use log::error;
use log::info;
use log::warn;
use rustc_hash::{FxHashMap,FxHashSet};
use gzp::{
    deflate::Gzip,
    ZWriter,
    par::compress::{ParCompress, ParCompressBuilder},
};
use smallvec::SmallVec;
use crate::bias::{read_bias_factors}; 

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

    let feature_list_file = Path::new(matches.get_one::<String>("features").unwrap())
        .canonicalize()
        .expect("Can't find path to input TF_id list file");
    info!("Received features (TF) file: {:?}", feature_list_file);

    let output_directory = matches.get_one::<String>("outdir").unwrap();
    info!("Received output directory: {:?}", output_directory);
    let output_path = Path::new(output_directory);

    let num_threads = *matches.get_one::<usize>("threads").unwrap();

    let track_bias = matches.get_flag("bias");

    let track_seq_consensus = matches.get_flag("seq");

    let track_control = matches.get_flag("control");


    // create directory if it does not exist
    if !output_path.exists() {
        if let Err(e) = fs::create_dir_all(output_path) {
            eprintln!("Failed to create output directory: {}", e);
            std::process::exit(1);
        }
    }

    // make sure output is a directory 
    match fs::metadata(output_path){
        Ok(metadata) => {
            if metadata.is_dir() {
                info!("{:?} is a directory.", output_path);
            } else {
                eprintln!("Provided output is not a directory: {}", output_path.display());
                std::process::exit(1);
            }
        } 
        Err(e) => {
            eprintln!("Failed to get metadata from {:?}: {}", output_path, e);
            std::process::exit(1);
        }
    }

    fcount(&frag_file, &bed_file, &cell_file, &feature_list_file, output_path, num_threads, 
           track_bias, track_seq_consensus, track_control)?;
    Ok(())
}

fn fcount(
    frag_file: &Path, 
    bed_file: &Path, 
    cell_file: &Path, 
    feature_list_file: &Path, 
    output: &Path, 
    num_threads: usize, 
    track_bias: bool, 
    track_seq_consensus: bool,
    track_control: bool,
    
) -> io::Result<()> {
    info!("Processing fragment file: {:?}, BED file: {:?}, Cell file:{:?}", frag_file, bed_file, cell_file);

    // get list of sites in bedfile as a 
    let mut lapper_map = match parse_file_to_lapper(bed_file, feature_list_file, track_control) {
        Ok(trees) => trees, 
        Err(e) => {
            error!("Failed to read Index or feature list file: {}", e);
            return Err(e);
        }
    };

    // create hashmap for cell barcodes {cell barcodes: cell index}
    let cellreader = File::open(cell_file).map(BufReader::new)?;
    let mut cells: FxHashMap<String, u32> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate(){
        let line = line?;
        let index_u32 = index as u32;
        cells.insert(line, index_u32);
    }

    // intialize the empty files 
    let mut aggregated_counts: FxHashMap<(String, String), Vec<usize>> = FxHashMap::default();
    let mut aggregated_bias: FxHashMap<(String, String), Vec<f32>> = FxHashMap::default();
    let mut consensus_seq: FxHashMap<(String, String), Vec<[u32;4]>> = FxHashMap::default(); // [A_count, C_count, G_count, T_count] 

    // read in the sequence fasta file and the bias vector file
    let h5_path = "/home/users/astar/gis/limtw/scratch/proj_TF/motif/Tn5Bias/hg38Tn5Bias.h5";
    let bias_factors = read_bias_factors(h5_path);

    let ref_path = "/home/users/astar/gis/limtw/scratch/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa";
    let mut faidx = load_fasta(ref_path).expect("Could not load ref fasta");
    let mut seq = Vec::new();
    // usage : 
    // let chr_trimmed = &chrom[3..];
    // faidx.fetch(&chr_trimmed, motif_start as u64, motif_end as u64).expect("Couldn't fetch interval 


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

    loop {
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
            //check_end = true;
            
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
            
                for interval in lapper.find(startpos as usize, startpos as usize + 1) {
                    
                    let interval_start = interval.start;
                    let interval_end = interval.stop;
                    let tf_id = interval.val.clone(); 

                    
                    let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![0;500]);
                    let relative_position: usize  = startpos - interval_start; 
                    insertion_counts[relative_position] += 1;
                    
                    let bias_landscape = aggregated_bias.entry((cell_barcode.to_string(), tf_id.to_string()))
                            .or_insert_with(|| vec![0.0;501]);
                    if let Some(bias) = bias_factors.get(&current_chrom) {
                        bias_landscape[0] += 1.0; // first position in bias vector is number of insertions
                        bias_landscape[1..=500].iter_mut()
                            .zip(bias[interval_start..interval_start+500].iter())
                            .for_each(|(a, b)| *a += b);
                    } else {
                        println!("No bias data for {}", current_chrom);
                    }

                    // fetch the sequence 
                    seq.clear();
                    let chr_trimmed = &current_chrom[3..];
                    faidx.fetch(&chr_trimmed, interval_start as u64, interval_end as u64).expect("Couldn't fetch sequence of interval");
                    faidx.read(&mut seq).expect("Couldn't read the sequence interval");
                    let pwm = consensus_seq.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![[0u32;4]; 500]); // interval_start to end should be 500
                    for (i, base) in seq.iter().enumerate() {
                        match base { 
                            b'A' | b'a' => pwm[i][0] += 1,
                            b'C' | b'c' => pwm[i][1] += 1,
                            b'G' | b'g' => pwm[i][2] += 1,
                            b'T' | b't' => pwm[i][3] += 1,
                            _ => {} // don't count Ns
                            }
                    }
                }
                for interval in lapper.find(endpos, endpos +1){
                    let interval_start = interval.start;
                    let interval_end = interval.stop;
                    let tf_id = interval.val.clone(); 

                    let insertion_counts = aggregated_counts.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![0;500]);
                    let relative_position: usize  = endpos - interval_start; 
                    insertion_counts[relative_position] += 1;
                    
                    let bias_landscape = aggregated_bias.entry((cell_barcode.to_string(), tf_id.to_string()))
                            .or_insert_with(|| vec![0.0;501]);
                    if let Some(bias) = bias_factors.get(&current_chrom) {
                        bias_landscape[0] += 1.0; // first position in bias vector is number of insertions
                        bias_landscape[1..=500].iter_mut()
                            .zip(bias[interval_start..interval_start+500].iter())
                            .for_each(|(a, b)| *a += b);
                    } else {
                        println!("No bias data for {}", current_chrom);
                    }

                    // fetch the sequence 
                    seq.clear();
                    let chr_trimmed = &current_chrom[3..];
                    faidx.fetch(&chr_trimmed, interval_start as u64, interval_end as u64).expect("Couldn't fetch sequence of interval");
                    faidx.read(&mut seq).expect("Couldn't read the sequence interval");
                    let pwm = consensus_seq.entry((cell_barcode.to_string(), tf_id.to_string()))
                        .or_insert_with(|| vec![[0u32;4]; 500]); // interval_start to end should be 500
                    for (i, base) in seq.iter().enumerate() {
                        match base { 
                            b'A' | b'a' => pwm[i][0] += 1,
                            b'C' | b'c' => pwm[i][1] += 1,
                            b'G' | b'g' => pwm[i][2] += 1,
                            b'T' | b't' => pwm[i][3] += 1,
                            _ => {} // don't count Ns
                        }
                    }
                }
                    
            }
        
        }
        line_str.clear();
    }
    eprintln!();

    // Write file 
    let output_path = output.join("aggregated_counts.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in aggregated_counts {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}", cell, tf, counts_str).expect("Failed to write to file");
    }
    let output_path = output.join("aggregated_bias.tsv");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), counts) in aggregated_bias {
        let counts_str = counts.iter().map(|v| v.to_string()).collect::<Vec<String>>().join("\t");
        writeln!(output_file, "{}\t{}\t{}\n", cell, tf, counts_str).expect("Failed to write to file");
    }
    // save consensus sequence as json file
    let output_path = output.join("consensus_pwm.jsonl");
    let mut output_file = File::create(output_path).expect("Failed to create output file");
    for ((cell, tf), pwm) in consensus_seq {
        // serialize to a JSON object
        let json_obj = json!({
            "cell":cell, 
            "tf":tf,
            "pwm": pwm,
        });
        writeln!(output_file, "{}", json_obj.to_string()).expect("Failed to write PWM JSON");
    }
    let mut writer = File::create(output.join("fragment_counts.tsv"))?;
    let mut foutput = String::new();
    for (cell_barcode, &cell_index) in &cells {
        foutput.push_str(&format!("{}\t{}\n", cell_barcode, frag_per_cell[cell_index as usize]));
    }
    writer.write_all(foutput.as_bytes())?;

    Ok(())      
}

// Reads in the special index file and organizes intervals
fn parse_file_to_lapper(
    index_file: &Path, 
    feature_list_file: &Path, 
    track_control: bool, 
) -> io::Result<FxHashMap<String, Lapper<usize, String >>>{

    // Read in list of tfs to track
    let reader_feature = File::open(feature_list_file)
        .map(BufReader::new)?;
    
    //let file_feature = File::open(feature_list_file);
    //let reader_feature = BufReader::new(file_feature);
    // hashset O(1) average-time lookup while Vec<32> gives O(n) lookup
    let mut keep_tf_names: FxHashSet<String> = FxHashSet::default();
    for line in reader_feature.lines() {
        if let Ok(l) = line {
            let t = l.trim();
            if !t.is_empty() {
                keep_tf_names.insert(t.to_string());
            }
        }
    }
    let mut keep_tf_indices: FxHashSet<usize> = FxHashSet::default();
    
    // Channel for decompressed lines 
    let (tx, rx) = mpsc::sync_channel::<Vec<String>>(50);
    let index_file = index_file.to_path_buf();

    // Spawn decompression thread 
    let reader_handle = thread::spawn(move || {
        let file = match File::open(&index_file) {
            Ok(f) => f, 
            Err(e) => {
                error!("Failed to open index file: {}", e);
                return;
            }
        };
        let reader = BufReader::with_capacity(1024*1024, MultiGzDecoder::new(file));
        const CHUNK_SIZE: usize = 10_000;
        let mut lines = Vec::with_capacity(CHUNK_SIZE);

        for line_result in reader.lines() {
            match line_result {
                Ok(line) => {
                    lines.push(line);
                    if lines.len() >= CHUNK_SIZE {
                        let chunk = std::mem::replace(&mut lines, Vec::with_capacity(CHUNK_SIZE));
                        if tx.send(chunk).is_err(){
                            break;
                        }
                    }
                }
                Err(e) => {
                    error!("Error reading index file: {}", e);
                    break;
                }
            }
        }

        // Send remaining lines
        if !lines.is_empty() {
            let _ = tx.send(lines);
        }
        info!("Finished reading index file");
    });

    // Initialize data structures
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<usize, String>>> = FxHashMap::default();

    let mut in_sites: bool = false;
    let mut in_controls: bool = false;
    let mut in_site_positions: bool = false;
    let mut in_control_positions: bool = false;

    let mut tf_index_to_name: FxHashMap<usize, String> = FxHashMap::default();

    // Process chunks from reader thread
    for chunk in rx {
        for line in chunk {
            let trimmed = line.trim();
            if trimmed.is_empty() { 
                continue;
            }

            // Handle section markers
            if trimmed.starts_with("##[") {
                info!("Section: {}", trimmed);
                if track_control {
                    in_sites = trimmed == "##[controls]";
                    in_control_positions = trimmed == "##[control positions]";
                } else {
                    in_sites = trimmed == "##[sites]";
                    in_site_positions = trimmed == "##[site positions]";
                }
                // info!("Section: {}", trimmed);
                // in_sites = trimmed == "##[sites]";
                // in_controls = trimmed == "##[controls]";
                // in_site_positions = trimmed == "##[site positions]";
                // in_control_positions = trimmed == "##[control positions]";
                continue;
            }
// Skip comments
            if trimmed.starts_with('#') {
                continue;
            }
            
            let fields: SmallVec<[&str; 8]> = trimmed.split('\t').collect();
            
            // parse [sites]
            if in_sites {
                if fields.len() >= 3 {
                    let tf_name = fields[0].to_string();
                    let tf_index: usize = fields[1].parse().expect("Invalid TF index");
                    if !keep_tf_names.contains(&tf_name){
                        continue;
                    }
                    keep_tf_indices.insert(tf_index);
                    tf_index_to_name.insert(tf_index, tf_name);
                }
            }
            
            // parse [controls]
            if in_controls {
               continue; 
            }
            
            // parse [site positions] - each interval has exactly 1 site index
            if in_site_positions {
                if fields.len() < 4 {
                    continue;
                }
                
                let chr = fields[0].to_string();
                let motif_start: usize = fields[1]
                    .parse()
                    .expect("Invalid start position");
                let motif_end: usize = fields[2]
                    .parse()
                    .expect("Invalid end position");
                let tf_index: usize = fields[3]
                    .parse()
                    .expect("Invalid TF index");

                if !keep_tf_indices.contains(&tf_index) {
                    continue;
                }
                let mid: usize = motif_start + (motif_end - motif_start)/2;
                
                let tf_name: String = tf_index_to_name.get(&tf_index).expect("tf_index not found").clone();
                
                let intervals = chromosome_trees.entry(chr).or_insert_with(Vec::new);
                intervals.push(Interval {
                    start: mid - 250,
                    stop: mid + 250,
                    val: tf_name,
                });
            }
            
            if in_control_positions {
                if fields.len() < 4 {
                    continue;
                }
                let chr = fields[0].to_string();
                let motif_start: usize = fields[1]
                    .parse()
                    .expect("Invalid start position");
                let motif_end: usize = fields[2]
                    .parse()
                    .expect("Invalid end position");

                let mid: usize = motif_start + (motif_end - motif_start)/2;

                let tf_indices_str = fields[3];
                let tf_indices: Vec<usize> = tf_indices_str
                    .split(',')
                    .map(|x| x.trim().parse::<usize>().expect("Invalid index"))
                    .collect();
                
                for tf_index in tf_indices {
                    if !keep_tf_indices.contains(&tf_index) {
                        continue;
                    }
                    let tf_name: String = tf_index_to_name.get(&tf_index).expect("tf_index not found").clone();
                    let intervals = chromosome_trees.entry(chr.clone()).or_insert_with(Vec::new);
                    intervals.push(Interval {
                        start: mid - 250,
                        stop: mid + 250,
                        val: tf_name,
                    });
                }
            }
        }
    }
    
    // Wait for reader thread to finish
    reader_handle.join().expect("Reader thread panicked");
    
    // Build lapper structure
    let lapper_map = chromosome_trees.into_iter()
        .map(|(chr, intervals)| {
            println!("Chromosome: {}, Number of intervals: {}", chr, intervals.len());
            (chr, Lapper::new(intervals))
        })
        .collect();
    
    Ok(lapper_map)
}
            
pub fn load_fasta(fasta_path: &str) -> Result<IndexedReader<File>, std::io::Error>{
    let fasta_path = Path::new(fasta_path);
    let fai_path = fasta_path.with_file_name(
    format!("{}.fai", fasta_path.file_name().unwrap().to_string_lossy())
);
    println!("fai_path:{}", fai_path.display());

    let fasta_file = File::open(fasta_path)?;
    println!("opened fasta file");
    let fai_file = File::open(fai_path)?;
    println!("opened fai file");
    let faidx = IndexedReader::new(fasta_file, fai_file)?;
    // if dont want result<>, just replace "?" with .unwrap()
    
    Ok(faidx)
    //https://docs.rs/bio/latest/bio/io/fasta/
} 
    
    
    
    
