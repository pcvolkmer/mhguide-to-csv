use crate::export_record::Record;
use clap::Parser;
use rayon::prelude::*;
use std::fs;

mod cli;
mod export_record;
mod hgnc;
mod mhguide;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();
    let json = fs::read_to_string(cli.input_file.clone())?;
    let mhguide = serde_json::from_str::<mhguide::MhGuide>(&json)?;

    let records = if cli.all_variants {
        mhguide.all_variants()
    } else if cli.oncogenic {
        mhguide.oncogenic_variants()
    } else {
        mhguide.relevant_variants()
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

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .escape(b'"')
        .delimiter(b';')
        .from_writer(vec![]);

    for record in records {
        let _ = writer.serialize(record);
    }

    let mut output_file = cli.input_file.clone();
    output_file.set_extension("csv");

    fs::write(output_file, writer.into_inner()?)?;

    Ok(())
}
