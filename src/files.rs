use crate::export_record::{BiomarkerRecord, CopyNumberRecord, FusionRecord, SimpleVariantRecord};
use crate::mhguide::MhGuide;
use itertools::Itertools;
use mv64e_mtb_dto::{
    Chromosome, Cnv, CnvCoding, CnvCodingCode, Coding, NgsReportResults, Position, Reference, Snv,
    TranscriptId, TranscriptIdSystem,
};
use rust_xlsxwriter::{Format, Workbook};
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::str::FromStr;

fn map_chromosome(s: &str) -> Result<Chromosome, ()> {
    match s {
        "chr1" => Ok(Chromosome::Chr1),
        "chr2" => Ok(Chromosome::Chr2),
        "chr3" => Ok(Chromosome::Chr3),
        "chr4" => Ok(Chromosome::Chr4),
        "chr5" => Ok(Chromosome::Chr5),
        "chr6" => Ok(Chromosome::Chr6),
        "chr7" => Ok(Chromosome::Chr7),
        "chr8" => Ok(Chromosome::Chr8),
        "chr9" => Ok(Chromosome::Chr9),
        "chr10" => Ok(Chromosome::Chr10),
        "chr11" => Ok(Chromosome::Chr11),
        "chr12" => Ok(Chromosome::Chr12),
        "chr13" => Ok(Chromosome::Chr13),
        "chr14" => Ok(Chromosome::Chr14),
        "chr15" => Ok(Chromosome::Chr15),
        "chr16" => Ok(Chromosome::Chr16),
        "chr17" => Ok(Chromosome::Chr17),
        "chr18" => Ok(Chromosome::Chr18),
        "chr19" => Ok(Chromosome::Chr19),
        "chr20" => Ok(Chromosome::Chr20),
        "chr21" => Ok(Chromosome::Chr21),
        "chr22" => Ok(Chromosome::Chr22),
        "chrX" => Ok(Chromosome::ChrX),
        "chrY" => Ok(Chromosome::ChrY),
        _ => Err(()),
    }
}

fn read_json_content(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    match path.extension() {
        Some(ext) if ext == "json" => Ok(fs::read_to_string(path)?),
        Some(ext) if ext == "zip" => {
            let mut archive = zip::ZipArchive::new(fs::File::open(path)?)?;
            if archive.len() != 1
                || !Path::new(archive.by_index(0)?.name())
                    .extension()
                    .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
            {
                return Err(
                    "ZIP archive does not contain a single JSON file. Only JSON files and ZIP compressed JSON files are supported."
                        .into(),
                );
            }
            let mut file = archive.by_index(0)?;
            let mut result = String::new();
            file.read_to_string(&mut result)?;
            Ok(result)
        }
        _ => Err(
            "Unsupported file format. Only JSON files and ZIP compressed JSON files are supported."
                .into(),
        ),
    }
}

pub(crate) fn read_file(path: &Path) -> Result<MhGuide, Box<dyn std::error::Error>> {
    let json = read_json_content(path)?;
    Ok(serde_json::from_str::<MhGuide>(&json)?)
}

pub(crate) fn write_csv_file(
    path: &Path,
    simple_variant_records: &[SimpleVariantRecord],
    copy_number_records: &[CopyNumberRecord],
    fusion_records: &[FusionRecord],
    biomarker_records: &[BiomarkerRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .flexible(true)
        .escape(b'"')
        .delimiter(b';')
        .from_writer(vec![]);

    if !simple_variant_records.is_empty() {
        let _ = writer.serialize(SimpleVariantRecord::csv_headlines());
        for record in simple_variant_records {
            let _ = writer.serialize(record);
        }
        let _ = writer.serialize(vec![""]);
    }

    if !copy_number_records.is_empty() {
        let _ = writer.serialize(CopyNumberRecord::csv_headlines());
        for record in copy_number_records {
            let _ = writer.serialize(record);
        }
        let _ = writer.serialize(vec![""]);
    }

    if !fusion_records.is_empty() {
        let _ = writer.serialize(FusionRecord::csv_headlines());
        for record in copy_number_records {
            let _ = writer.serialize(record);
        }
        let _ = writer.serialize(vec![""]);
    }

    if !biomarker_records.is_empty() {
        let _ = writer.serialize(BiomarkerRecord::csv_headlines());
        for record in biomarker_records {
            let _ = writer.serialize(record);
        }
        let _ = writer.serialize(vec![""]);
    }

    let mut output_file = path.to_path_buf();
    output_file.set_extension("csv");

    fs::write(output_file, writer.into_inner()?).map_err(Into::into)
}

pub(crate) fn write_xlsx_file(
    path: &Path,
    simple_variant_records: &[SimpleVariantRecord],
    copy_number_records: &[CopyNumberRecord],
    fusion_records: &[FusionRecord],
    biomarker_records: &[BiomarkerRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    fn write_worksheet<T>(
        workbook: &mut Workbook,
        name: &str,
        records: &[T],
    ) -> Result<(), Box<dyn std::error::Error>>
    where
        T: Serialize + for<'de> Deserialize<'de>,
    {
        let worksheet = workbook.add_worksheet();
        worksheet.set_name(name)?;

        worksheet.deserialize_headers_with_format::<T>(0, 0, &Format::new().set_bold())?;
        worksheet.serialize(&records)?;

        worksheet.autofit();

        Ok(())
    }

    let mut workbook = Workbook::new();

    if !simple_variant_records.is_empty() {
        write_worksheet(&mut workbook, "Einfache Varianten", simple_variant_records)?;
    }

    if !copy_number_records.is_empty() {
        write_worksheet(&mut workbook, "Copy Number Varianten", copy_number_records)?;
    }

    if !fusion_records.is_empty() {
        write_worksheet(&mut workbook, "Fusionen", fusion_records)?;
    }

    if !biomarker_records.is_empty() {
        write_worksheet(&mut workbook, "Biomarker", biomarker_records)?;
    }

    let mut output_file = path.to_path_buf();
    output_file.set_extension("xlsx");
    workbook.save(output_file).map_err(Into::into)
}

pub(crate) fn write_json_file(
    path: &Path,
    simple_variant_records: &[SimpleVariantRecord],
    copy_number_records: &[CopyNumberRecord],
    fusion_records: &[FusionRecord],
    biomarker_records: &[BiomarkerRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let simple_variants = simple_variant_records
        .iter()
        .map(|record| Snv {
            allelic_frequency: f64::from_str(&record.allelic_frequency.replace(',', "."))
                .unwrap_or(0.0),
            alt_allele: if record.alt_allele.is_empty() {
                "-".to_string()
            } else {
                record.alt_allele.clone()
            },
            // Use Chromosome::ChrMt as placeholder for non present value
            chromosome: map_chromosome(&record.chromosome).unwrap_or(Chromosome::ChrMt),
            dna_change: record.cdna.clone(),
            exon_id: None,
            external_ids: None,
            gene: Coding {
                code: record.gene.clone(),
                display: Some(record.hgnc_name.clone()),
                system: None,
                version: None,
            },
            id: String::new(),
            interpretation: None,
            localization: None,
            patient: Reference {
                display: None,
                id: String::new(),
                reference_type: None,
                system: None,
            },
            position: Position {
                start: record.start.parse().unwrap_or(0.0),
                end: record.end.parse().ok(),
            },
            protein_change: if record.protein.clone().is_empty() {
                None
            } else {
                Some(record.protein.clone())
            },
            read_depth: record.read_depth.parse().unwrap_or(0), // To be interpreted as "not present"
            ref_allele: if record.ref_allele.is_empty() {
                "-".to_string()
            } else {
                record.ref_allele.clone()
            },
            transcript_id: TranscriptId {
                system: TranscriptIdSystem::EnsemblOrg,
                value: record.ensembl_id.clone(),
            },
        })
        .collect_vec();

    let copy_number_variants = copy_number_records
        .iter()
        .map(|record| Cnv {
            // Use Chromosome::ChrMt as placeholder for non present value
            chromosome: map_chromosome(&record.chromosome).unwrap_or(Chromosome::ChrMt),
            cn_a: None,
            cn_b: None,
            cnv_type: CnvCoding {
                code: if record.cnv_type.contains("loss") {
                    CnvCodingCode::Loss
                } else if record
                    .total_copy_number
                    .replace(',', ".")
                    .parse()
                    .unwrap_or(0.0)
                    < 3.0
                {
                    CnvCodingCode::LowLevelGain
                } else {
                    CnvCodingCode::HighLevelGain
                },
                display: Some(record.cnv_type.clone()),
                system: None,
                version: None,
            },
            copy_number_neutral_lo_h: None,
            end_range: None,
            external_ids: None,
            id: String::new(),
            localization: None,
            patient: Reference {
                display: None,
                id: String::new(),
                reference_type: None,
                system: None,
            },
            relative_copy_number: None,
            reported_affected_genes: None,
            reported_focality: None,
            start_range: None,
            total_copy_number: record
                .total_copy_number
                .replace(',', ".")
                .parse::<i64>()
                .ok(),
        })
        .collect_vec();

    let ngs_report_results = NgsReportResults {
        brcaness: None,
        copy_number_variants: if copy_number_variants.is_empty() {
            None
        } else {
            Some(copy_number_variants)
        },
        dna_fusions: None,
        hrd_score: None,
        rna_fusions: None,
        rna_seqs: None,
        simple_variants: if simple_variants.is_empty() {
            None
        } else {
            Some(simple_variants)
        },
        tmb: None,
        tumor_cell_content: None,
    };
    let json_content = serde_json::to_string_pretty(&ngs_report_results)?;
    let mut output_file = path.to_path_buf();
    output_file.set_extension("dnpm.json");
    fs::write(output_file, json_content).map_err(Into::into)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::files::read_json_content;
    use std::path::PathBuf;
    use std::str::FromStr;

    const TEST_CONTENT: &str = include_str!("../testfiles/sv-mhguide.json");

    #[test]
    fn test_should_read_json_content() {
        let actual =
            read_json_content(&PathBuf::from_str("./testfiles/sv-mhguide.json").unwrap()).unwrap();
        assert_eq!(actual, TEST_CONTENT);
    }

    #[test]
    fn test_should_read_zip_content() {
        let actual =
            read_json_content(&PathBuf::from_str("./testfiles/sv-mhguide.json.zip").unwrap())
                .unwrap();
        assert_eq!(actual, TEST_CONTENT);
    }
}
