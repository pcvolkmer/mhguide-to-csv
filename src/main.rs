use crate::export_record::{BiomarkerRecord, CopyNumberRecord, FusionRecord, SimpleVariantRecord};
use crate::files::read_file;
use crate::mhguide::ResultType;
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

    let variants = if cli.all_variants {
        mhguide.all_variants()
    } else if cli.oncogenic {
        mhguide.oncogenic_variants()
    } else {
        mhguide.relevant_variants(cli.no_artifacts)
    };

    let simple_variant_records = variants
        .par_iter()
        .filter(|variant| variant.gene_symbol.is_some())
        .filter(|variant| match &variant.display_variant_type {
            Some(ResultType::SimpleVariant(_)) => true,
            Some(_) => false,
            None => matches!(
                &variant.protein_variant_type,
                Some(ResultType::SimpleVariant(_))
            ),
        })
        .map(|variant| {
            SimpleVariantRecord::from_variant(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                variant,
            )
        })
        .collect::<Vec<_>>();

    let copy_number_records = variants
        .par_iter()
        .filter(|variant| variant.gene_symbol.is_some())
        .filter(|variant| match &variant.display_variant_type {
            Some(ResultType::CopyNumberVariant) => true,
            Some(_) => false,
            None => matches!(
                &variant.protein_variant_type,
                Some(ResultType::CopyNumberVariant)
            ),
        })
        .map(|variant| {
            CopyNumberRecord::from_variant(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                variant,
            )
        })
        .collect::<Vec<_>>();

    let fusion_records = mhguide
        .fusions()
        .par_iter()
        .map(|fusion| {
            FusionRecord::from_fusion(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                fusion,
            )
        })
        .collect::<Vec<_>>();

    let mut biomarker_records = vec![];
    if let Some(value) = mhguide.hrd_score() {
        biomarker_records.push(BiomarkerRecord::from_hrd(
            &mhguide.general.patient_identifier.h_number,
            &mhguide.general.ref_genome_version,
            value,
        ));
    }
    if let Some(value) = mhguide.msi_score() {
        biomarker_records.push(BiomarkerRecord::from_msi(
            &mhguide.general.patient_identifier.h_number,
            &mhguide.general.ref_genome_version,
            value,
        ));
    }
    if let Some(value) = mhguide.tmb_value() {
        biomarker_records.push(BiomarkerRecord::from_tmb(
            &mhguide.general.patient_identifier.h_number,
            &mhguide.general.ref_genome_version,
            value,
        ));
    }

    if cli.xlsx {
        return files::write_xlsx_file(
            &cli.input_file,
            &simple_variant_records,
            &copy_number_records,
            &fusion_records,
            &biomarker_records,
        );
    }

    files::write_csv_file(
        &cli.input_file,
        &simple_variant_records,
        &copy_number_records,
        &fusion_records,
        &biomarker_records,
    )
}
