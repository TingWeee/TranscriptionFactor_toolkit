#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: std::alloc::System = std::alloc::System;

use clap::{Command, Arg};
use std::error::Error;
use clap::ArgAction;
mod footprint;
mod bias; 

fn main() -> Result<(), Box<dyn Error>> {
    
    let matches = Command::new("footprint")
        .about("Genomic Landscape Characterisation")
        .arg_required_else_help(true)
        .subcommand(
            Command::new("vector")
                .about("Create a position vector x cell matrix from a fragment file")
                .arg(
                    Arg::new("fragments")
                        .short('f')
                        .long("fragments")
                        .help("Path to the fragment file")
                        .required(true),
                )
                .arg(
                    Arg::new("bed")
                        .short('b')
                        .long("bed")
                        .help("BED file containing non-overlapping genomic regions to quantify")
                        .required(true),
                )
                .arg(
                    Arg::new("cells")
                        .short('c')
                        .long("cells")
                        .help("File containing cell barcodes to include")
                        .required(true),
                )
                .arg(
                    Arg::new("features")
                        .long("features")
                        .help("Txt file containing TF ids to include")
                        .required(true),
                )
                .arg(
                    Arg::new("outdir")
                        .short('o')
                        .long("outdir")
                        .help("Output directory name")
                        .long_help("Output directory name. Directory will be created if it does not exist. \
                               The output directory will contain matrix.mtx.gz, features.tsv, barcodes.tsv")
                        .required(true),
                )
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .help("Number of compression threads to use")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("4")
                        .required(false),
                )
                .arg(
                    Arg::new("bias")
                        .long("bias")
                        .help("Whether to keep track of the bias")
                        .action(ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("seq")
                        .long("seq")
                        .help("Whether to keep track of sequence consensus")
                        .action(ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("control")
                    .long("control")
                    .help("Default is to get motif sites, set control to aggregate control regions")
                    .action(ArgAction::SetTrue)
                )
                //.arg(
                 //   Arg::new("length")
                  //  .long("length")
                   // .help("length of region to keep track")
                   // .default_value("500")
                   // .required(false)
               // )
        )   
        .get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        Some(("vector", sub_matches)) => footprint::f2v(sub_matches)?,
        _ => {

        }
    }

    Ok(())
}