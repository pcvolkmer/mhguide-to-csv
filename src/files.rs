use crate::export_record::{BiomarkerRecord, CopyNumberRecord, FusionRecord, SimpleVariantRecord};
use crate::mhguide::MhGuide;
use rust_xlsxwriter::{Format, Workbook};
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Read;
use std::path::Path;

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
    simple_variant_records: &Vec<SimpleVariantRecord>,
    copy_number_records: &Vec<CopyNumberRecord>,
    fusion_records: &Vec<FusionRecord>,
    biomarker_records: &Vec<BiomarkerRecord>,
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
    simple_variant_records: &Vec<SimpleVariantRecord>,
    copy_number_records: &Vec<CopyNumberRecord>,
    fusion_records: &Vec<FusionRecord>,
    biomarker_records: &Vec<BiomarkerRecord>,
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
