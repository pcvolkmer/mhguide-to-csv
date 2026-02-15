use crate::mhguide::{three_letter_protein_modification, RefGenomeVersion};
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use std::fs;

mod cli;
mod mhguide;

#[derive(Debug, Serialize)]
struct Csv {
    #[serde(rename = "H-Nummer")]
    h_nummer: String,
    #[serde(rename = "Referenz-Genom")]
    ref_genome: String,
    #[serde(rename = "Ergebnis")]
    variantenart: String,
    #[serde(rename = "Gen")]
    gene: String,
    #[serde(rename = "Genomposition (g.)")]
    genomic_position: String,
    #[serde(rename = "cDNA Nomenklatur (c.)")]
    cdna: String,
    #[serde(rename = "Proteinebene Nomenklatur (p.)")]
    protein: String,
    #[serde(rename = "Chromosom")]
    chromosome: String,
    #[serde(rename = "Start")]
    start: String,
    #[serde(rename = "Ende")]
    end: String,
    #[serde(rename = "Alternative Nucleotide")]
    alt_allele: String,
    #[serde(rename = "Reference Nucleotide")]
    ref_allele: String,
    #[serde(rename = "Allelfrequenz (%)")]
    allelic_frequency: String,
    #[serde(rename = "dbSNP ID")]
    dbsnp: String,
    #[serde(rename = "Type")]
    cnv_type: String,
    #[serde(rename = "Total CN")]
    total_copy_number: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();
    let json = std::fs::read_to_string(cli.input_file.clone())?;
    let mhguide = serde_json::from_str::<mhguide::MhGuide>(&json)?;

    let records = mhguide
        .variants
        .par_iter()
        .filter(|variant| variant.gene_symbol.is_some())
        .map(|variant| Csv {
            h_nummer: mhguide.general.patient_identifier.h_number.clone(),
            ref_genome: match &mhguide.general.ref_genome_version {
                RefGenomeVersion::Hg19 => "HG19",
                RefGenomeVersion::Hg38 => "HG38",
            }
            .to_string(),
            variantenart: if variant
                .protein_modification
                .clone()
                .unwrap_or_default()
                .to_ascii_lowercase()
                .starts_with("copy number")
            {
                "Copy Number Variation"
            } else {
                "Einfache Variante (?)"
            }
            .to_string(),
            gene: variant.gene_symbol.clone().unwrap_or_default(),
            chromosome: variant.chromosome.clone().unwrap_or_default(),
            cdna: if variant
                .transcript_hgvs_modified_object
                .clone()
                .unwrap_or_default()
                .starts_with("c.")
            {
                variant
                    .transcript_hgvs_modified_object
                    .clone()
                    .unwrap_or_default()
            } else {
                String::new()
            },
            protein: if variant
                .protein_modification
                .clone()
                .unwrap_or_default()
                .starts_with("p.")
            {
                three_letter_protein_modification(&variant.protein_modification.clone().unwrap_or_default())
            } else {
                String::new()
            },
            genomic_position: match variant
                .chromosome_modification
                .clone()
                .unwrap_or_default()
                .as_str()
            {
                "TMB" | "MSS" | "MSI" | "HRD-positive" => "",
                value => {
                    if value.starts_with("Chr") {
                        ""
                    } else {
                        value
                    }
                }
            }
            .to_string(),
            start: variant.start(),
            end: variant.end(),
            ref_allele: variant.ref_allele(),
            alt_allele: variant.alt_allele(),
            allelic_frequency: match variant.variant_allele_frequency_in_tumor {
                Some(value) => format!("{value:.2}"),
                None => String::new(),
            },
            dbsnp: variant.db_snp.clone().unwrap_or_default(),
            total_copy_number: match variant.copy_number {
                Some(value) => format!("{value:.2}"),
                None => String::new(),
            },
            cnv_type: match variant
                .protein_modification
                .clone()
                .unwrap_or_default()
                .to_ascii_lowercase()
                .as_str()
            {
                "copy number loss" => "loss",
                "copy number gain" => match variant.copy_number {
                    Some(value) => {
                        if value >= 3.0 {
                            "high level gain"
                        } else {
                            "low level gain"
                        }
                    }
                    None => "",
                },
                _ => "",
            }
            .to_string(),
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
