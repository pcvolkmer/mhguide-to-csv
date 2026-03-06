use crate::export_record::Record;
use crate::mhguide::MhGuide;
use rust_xlsxwriter::{Format, Workbook};
use std::fs;
use std::io::Read;
use std::path::PathBuf;

fn read_json_content(path: &PathBuf) -> Result<String, Box<dyn std::error::Error>> {
    match path.extension() {
        Some(ext) if ext == "json" => Ok(fs::read_to_string(path)?),
        Some(ext) if ext == "zip" => {
            let mut archive = zip::ZipArchive::new(fs::File::open(path)?)?;
            if archive.len() != 1
                || !std::path::Path::new(archive.by_index(0)?.name())
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

pub(crate) fn read_file(path: &PathBuf) -> Result<MhGuide, Box<dyn std::error::Error>> {
    let json = read_json_content(path)?;
    Ok(serde_json::from_str::<MhGuide>(&json)?)
}

pub(crate) fn write_csv_file(
    path: &PathBuf,
    records: &Vec<Record>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .escape(b'"')
        .delimiter(b';')
        .from_writer(vec![]);

    for record in records {
        let _ = writer.serialize(record);
    }

    let mut output_file = path.clone();
    output_file.set_extension("csv");

    fs::write(output_file, writer.into_inner()?).map_err(Into::into)
}

pub(crate) fn write_xlsx_file(
    path: &PathBuf,
    records: &Vec<Record>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();
    worksheet.set_name("Variants")?;

    worksheet.deserialize_headers_with_format::<Record>(0, 0, &Format::new().set_bold())?;
    worksheet.serialize(&records)?;

    worksheet.autofit();

    let mut output_file = path.clone();
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
