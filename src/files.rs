use crate::export::json;
use crate::export::json::{OsMolekulargenUntersuchung, OsMolekulargenetik};
use crate::export::table::{BiomarkerRecord, CopyNumberRecord, FusionRecord, SimpleVariantRecord};
use itertools::Itertools;
use mhguide_umr::{MhGuide, PathogenicClassification};
use rust_xlsxwriter::{Format, Workbook};
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::str::FromStr;

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

#[allow(unused)]
pub(crate) fn write_json_file(
    path: &Path,
    mh_guide: &MhGuide,
    simple_variant_records: &[SimpleVariantRecord],
    copy_number_records: &[CopyNumberRecord],
    fusion_records: &[FusionRecord],
    biomarker_records: &[BiomarkerRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let json_content = get_json_data(
        path,
        mh_guide,
        simple_variant_records,
        copy_number_records,
        fusion_records,
        biomarker_records,
    )?;
    let mut output_file = path.to_path_buf();
    output_file.set_extension("os_molgen.json");
    fs::write(output_file, json_content).map_err(Into::into)
}

#[allow(unused)]
pub(crate) fn get_json_data(
    path: &Path,
    mh_guide: &MhGuide,
    simple_variant_records: &[SimpleVariantRecord],
    copy_number_records: &[CopyNumberRecord],
    fusion_records: &[FusionRecord],
    biomarker_records: &[BiomarkerRecord],
) -> Result<String, Box<dyn std::error::Error>> {
    fn map_pathogenic_classification(s: &str) -> Option<String> {
        let result = match PathogenicClassification::from_str(s) {
            Ok(pathogenic_classification) => match pathogenic_classification {
                PathogenicClassification::Benign => "1",
                PathogenicClassification::LikelyBenign => "2",
                PathogenicClassification::Vus => "3",
                PathogenicClassification::LikelyPathogenic => "4",
                PathogenicClassification::Pathogenic => "5",
            },
            Err(()) => return None,
        };

        Some(result.to_string())
    }

    fn empty_option(s: &str) -> Option<String> {
        if !s.is_empty() {
            return Some(s.to_string());
        }

        None
    }

    let mut variants = simple_variant_records
        .iter()
        .map(|record| OsMolekulargenUntersuchung::SimpleVariant {
            untersucht: record.gene.clone(),
            genomposition: record.genomic_position.clone(),
            cdnanomenklatur: record.cdna.clone(),
            proteinebenenomenklatur: empty_option(&record.protein),
            evchromosom: record.chromosome.clone(),
            evensemblid: record.ensembl_id.clone(),
            evhgncid: record.hgnc_id.clone(),
            evhgncsymbol: empty_option(&record.gene),
            evhgncname: empty_option(&record.hgnc_name),
            evstart: u64::from_str(&record.start).ok(),
            evende: u64::from_str(&record.end).ok(),
            evaltnucleotide: if record.alt_allele.is_empty() {
                "*".to_string()
            } else {
                record.alt_allele.clone()
            },
            evrefnucleotide: if record.ref_allele.is_empty() {
                "*".to_string()
            } else {
                record.ref_allele.clone()
            },
            evreaddepth: u64::from_str(&record.read_depth).ok(),
            allelfrequenz: f64::from_str(&record.allelic_frequency.replace(',', ".")).ok(),
            evdbsnpid: empty_option(&record.dbsnp),
            pathogenitaetsklasse: map_pathogenic_classification(&record.classification),
        })
        .collect_vec();

    let mut copy_number_variants = copy_number_records
        .iter()
        .map(|record| OsMolekulargenUntersuchung::CopyNumberVariant {
            untersucht: record.gene.clone(),
            copynumbervariation: if record.cnv_type.to_ascii_lowercase().contains("loss") {
                "L"
            } else if record.cnv_type.to_ascii_lowercase().contains("gain") {
                "GAIN"
            } else {
                ""
            }
            .to_string(),
            cnvchromosom: record.chromosome.clone(),
            cnvensemblid: empty_option(&record.ensembl_id),
            cnvhgncid: empty_option(&record.hgnc_id),
            cnvhgncsymbol: empty_option(&record.gene),
            cnvhgncname: empty_option(&record.hgnc_name),
            cnvtotalcndouble: f64::from_str(&record.total_copy_number.replace(',', ".")).ok(),
            pathogenitaetsklasse: map_pathogenic_classification(&record.classification),
        })
        .collect_vec();

    variants.append(&mut copy_number_variants);

    let mut fusion_variants = fusion_records
        .iter()
        .map(|record| OsMolekulargenUntersuchung::RnaFusion {
            untersucht: record.gene.clone(),
            fusioniertesgen: record.fusion_gene.clone(),
            fusionrna5ensemblid: empty_option(&record.ensembl_id_5),
            fusionrna5hgncid: record.hgnc_id_5.clone(),
            fusionrna5hgncsymbol: empty_option(&record.hgnc_name_5),
            fusionrna5hgncname: empty_option(&record.hgnc_name_5),
            fusionrna5transcriptid: record.transcript_id_5.clone(),
            fusionrna5exonid: record.exon_id_5.clone(),
            fusionrna5transposition: u64::from_str(&record.transcript_position_5.clone()).ok(),
            fusionrna5strand: record.strand_5.clone(),
            fusionrna3ensemblid: empty_option(&record.ensembl_id_3),
            fusionrna3hgncid: record.hgnc_id_3.clone(),
            fusionrna3hgncsymbol: empty_option(&record.hgnc_name_3),
            fusionrna3hgncname: empty_option(&record.hgnc_name_3),
            fusionrna3transcriptid: record.transcript_id_3.clone(),
            fusionrna3exonid: record.exon_id_3.clone(),
            fusionrna3transposition: u64::from_str(&record.transcript_position_3.clone()).ok(),
            fusionrna3strand: record.strand_3.clone(),
            fusionrnareportednumread: u64::from_str(&record.number_reported_reads.clone()).ok(),
            pathogenitaetsklasse: map_pathogenic_classification(&record.classification),
        })
        .collect_vec();

    variants.append(&mut fusion_variants);

    let biomarker = biomarker_records
        .iter()
        .filter_map(|record| {
            let hrd = f64::from_str(&record.hrd.replace(',', ".")).ok();
            if hrd.is_some() {
                return Some(json::Biomarker::Hrd {
                    score: hrd,
                    lst: None,
                    loh: None,
                    tai: None,
                    bewertung: None,
                });
            }

            let msi = f64::from_str(&record.msi.replace(',', ".")).ok();
            if msi.is_some() {
                return Some(json::Biomarker::Msi {
                    seqprozentwert: msi,
                    pcrergebnis: None,
                    immunergebnismsi: None,
                });
            }

            let tmb = f64::from_str(&record.tmb.replace(',', ".")).ok();
            if tmb.is_some() {
                return Some(json::Biomarker::Tmb {
                    tumormutationalburden: tmb,
                    bewertung: None,
                });
            }

            None
        })
        .collect_vec();

    let result = OsMolekulargenetik {
        patient_id: mh_guide.general.patient_identifier.pid.clone(),
        datum: chrono::NaiveDate::parse_from_str(&mh_guide.report_signed_formatted, "%d %b %Y")
            .unwrap_or_default(),
        einsendenummer: mh_guide.general.patient_identifier.h_number.clone(),
        referenzgenom: mh_guide.general.ref_genome_version.to_string(),
        molekulargenuntersuchung: if variants.is_empty() {
            None
        } else {
            Some(variants)
        },
        biomarker: if biomarker.is_empty() {
            None
        } else {
            Some(biomarker)
        },
    };

    serde_json::to_string_pretty(&result).map_err(Into::into)
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
