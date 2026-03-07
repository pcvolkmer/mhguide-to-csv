use crate::export_record::Record;
use crate::files::read_file;
use clap::Parser;
use rayon::prelude::*;

mod cli;
mod export_record;
pub mod files;
mod hgnc;
mod mhguide;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();
    let mhguide = read_file(&cli.input_file)?;

    let records = if cli.all_variants {
        mhguide.all_variants()
    } else if cli.oncogenic {
        mhguide.oncogenic_variants()
    } else {
        mhguide.relevant_variants(cli.no_artifacts)
    }
    .par_iter()
    .filter(|variant| variant.gene_symbol.is_some())
    .map(|variant| {
        Record::from_variant(
            &mhguide.general.patient_identifier.h_number,
            &mhguide.general.ref_genome_version,
            variant,
        )
    })
    .collect::<Vec<_>>();

    if cli.xlsx {
        return files::write_xlsx_file(&cli.input_file, &records);
    }

    files::write_csv_file(&cli.input_file, &records)
}
