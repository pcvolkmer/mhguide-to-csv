use crate::cli::{Cli, SubCommand};
use crate::export::json::OsMolekulargenetik;
use crate::files::read_file;
use clap::Parser;
use dialoguer::{Input, Password};
use export::table::{BiomarkerRecord, CopyNumberRecord, FusionRecord, SimpleVariantRecord};
use mhguide_umr::{MhGuide, ResultType, Variant};
use rayon::prelude::*;
use schemars::schema_for;

mod cli;
mod export;
mod files;
mod hgnc;

fn simple_variant_records(mhguide: &MhGuide, cli: &Cli) -> Vec<SimpleVariantRecord> {
    let variants = selected_variants(mhguide, cli);
    variants
        .par_iter()
        .filter(|&&variant| variant.gene_symbol.is_some())
        .filter(|&&variant| match &variant.display_variant_type {
            Some(ResultType::SimpleVariant(_)) => true,
            Some(_) => false,
            None => matches!(
                &variant.protein_variant_type,
                Some(ResultType::SimpleVariant(_))
            ),
        })
        .map(|&variant| {
            SimpleVariantRecord::from_variant(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                variant,
            )
        })
        .collect::<Vec<_>>()
}

fn copy_number_records(mhguide: &MhGuide, cli: &Cli) -> Vec<CopyNumberRecord> {
    let variants = selected_variants(mhguide, cli);
    variants
        .par_iter()
        .filter(|&&variant| variant.gene_symbol.is_some())
        .filter(|&&variant| match &variant.display_variant_type {
            Some(ResultType::CopyNumberVariant) => true,
            Some(_) => false,
            None => matches!(
                &variant.protein_variant_type,
                Some(ResultType::CopyNumberVariant)
            ),
        })
        .map(|&variant| {
            CopyNumberRecord::from_variant(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                variant,
            )
        })
        .collect::<Vec<_>>()
}

fn fusion_records(mhguide: &MhGuide) -> Vec<FusionRecord> {
    mhguide
        .fusions()
        .par_iter()
        .map(|fusion| {
            FusionRecord::from_fusion(
                &mhguide.general.patient_identifier.h_number,
                &mhguide.general.ref_genome_version,
                fusion,
            )
        })
        .collect::<Vec<_>>()
}

fn biomarker_records(mhguide: &MhGuide) -> Vec<BiomarkerRecord> {
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

    biomarker_records
}

fn selected_variants<'a>(mhguide: &'a MhGuide, cli: &Cli) -> Vec<&'a Variant> {
    match &cli.command {
        SubCommand::Convert { convert_args, .. } | SubCommand::Push { convert_args, .. } => {
            if convert_args.all_variants {
                mhguide.all_variants()
            } else if convert_args.oncogenic {
                mhguide.oncogenic_variants()
            } else {
                mhguide.relevant_variants(convert_args.no_artifacts)
            }
        }
        _ => vec![],
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        SubCommand::Convert { convert_args, .. } => {
            let mhguide = read_file(&convert_args.input_file)?;

            let simple_variant_records = simple_variant_records(&mhguide, &cli);
            let copy_number_records = copy_number_records(&mhguide, &cli);
            let fusion_records = fusion_records(&mhguide);
            let biomarker_records = biomarker_records(&mhguide);

            if convert_args.xlsx {
                return files::write_xlsx_file(
                    &convert_args.input_file,
                    &simple_variant_records,
                    &copy_number_records,
                    &fusion_records,
                    &biomarker_records,
                );
            }

            if convert_args.json {
                return files::write_json_file(
                    &convert_args.input_file,
                    &mhguide,
                    &simple_variant_records,
                    &copy_number_records,
                    &fusion_records,
                    &biomarker_records,
                );
            }

            files::write_csv_file(
                &convert_args.input_file,
                &simple_variant_records,
                &copy_number_records,
                &fusion_records,
                &biomarker_records,
            )
        }
        SubCommand::Push { convert_args, url } => {
            #[derive(serde::Deserialize)]
            struct Response {
                success: bool,
                #[serde(rename = "addedVariants")]
                added_variants: usize,
                #[serde(rename = "removedVariants")]
                removed_variants: usize,
                #[serde(rename = "addedBiomarkers")]
                added_biomarkers: usize,
                #[serde(rename = "removedBiomarkers")]
                removed_biomarkers: usize,
            }

            let mhguide = read_file(&convert_args.input_file)?;

            let simple_variant_records = simple_variant_records(&mhguide, &cli);
            let copy_number_records = copy_number_records(&mhguide, &cli);
            let fusion_records = fusion_records(&mhguide);
            let biomarker_records = biomarker_records(&mhguide);

            let json = files::get_json_data(
                &convert_args.input_file,
                &mhguide,
                &simple_variant_records,
                &copy_number_records,
                &fusion_records,
                &biomarker_records,
            )?;

            let username: String = Input::new().with_prompt("Benutzername").interact_text()?;
            let password: String = Password::new().with_prompt("Passwort").interact()?;

            let client = reqwest::blocking::Client::builder()
                .user_agent("mhguide-to-csv")
                .build()?;

            let response = client
                .put(
                    format!(
                        "{}/x-api/patient/{}/molgen/{}",
                        url.strip_suffix("/").unwrap_or(url),
                        mhguide.general.patient_identifier.pid,
                        mhguide.general.patient_identifier.h_number
                    )
                    .as_str(),
                )
                .basic_auth(username, Some(password))
                .header("Content-Type", "application/json")
                .body(json)
                .send();

            if let Ok(response) = response {
                if response.status().is_success() {
                    println!("Daten erfolgreich in Onkostar veröffentlicht");
                    if let Ok(response) = serde_json::from_str::<Response>(&response.text()?) {
                        println!(
                            "Varianten - hinzugefügt: {}, gelöscht: {}",
                            response.added_variants, response.removed_variants
                        );
                        println!(
                            "Biomarker - hinzugefügt: {}, gelöscht: {}",
                            response.added_biomarkers, response.removed_biomarkers
                        );
                    }
                } else {
                    println!(
                        "Fehler beim Veröffentlichen der Daten in Onkostar: {}",
                        response.status()
                    );
                }
            } else {
                println!("Fehler beim Senden der Daten an Onkostar");
            }
            Ok(())
        }
        SubCommand::JsonSchema => {
            let schema = schema_for!(OsMolekulargenetik);
            println!("{}", serde_json::to_string_pretty(&schema)?);
            Ok(())
        }
    }
}
