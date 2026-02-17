use regex::Regex;
use serde::{Deserialize, Deserializer};
use std::str::FromStr;

#[derive(Debug, Deserialize, PartialEq)]
pub(crate) struct MhGuide {
    #[serde(rename = "GENERAL")]
    pub(crate) general: General,
    #[serde(rename = "VARIANT_LONG_LIST")]
    pub(crate) variants: Vec<Variant>,
}

#[derive(Debug, PartialEq)]
pub(crate) enum RefGenomeVersion {
    Hg19,
    Hg38,
}

impl<'de> Deserialize<'de> for RefGenomeVersion {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value: u8 = Deserialize::deserialize(deserializer)?;
        match value {
            37 => Ok(RefGenomeVersion::Hg19),
            38 => Ok(RefGenomeVersion::Hg38),
            _ => Err(serde::de::Error::custom(format!(
                "Invalid RefGenomeVersion: {value}"
            ))),
        }
    }
}

#[derive(Debug, PartialEq)]
pub(crate) struct PatientIdentifier {
    pub(crate) h_number: String,
    pub(crate) pid: String,
}

impl<'de> Deserialize<'de> for PatientIdentifier {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value: String = Deserialize::deserialize(deserializer)?;
        let parts: Vec<&str> = value.split('_').collect();
        if parts.len() != 2 {
            return Err(serde::de::Error::custom("Invalid PatientIdentifier format"));
        }
        let h_number = parts[0].to_string();
        let pid = parts[1].to_string();
        Ok(PatientIdentifier { h_number, pid })
    }
}

#[derive(Debug, Deserialize, PartialEq)]
pub(crate) struct General {
    #[serde(rename = "ORDER_DATE")]
    pub(crate) order_date: String,
    #[serde(rename = "REF_GENOME_VERSION")]
    pub(crate) ref_genome_version: RefGenomeVersion,
    #[serde(rename = "PATIENT_IDENTIFIER")]
    pub(crate) patient_identifier: PatientIdentifier,
}

#[derive(Debug, Deserialize, PartialEq)]
pub(crate) struct Variant {
    #[serde(rename = "GENE_SYMBOL")]
    pub(crate) gene_symbol: Option<String>,
    #[serde(rename = "PROTEIN_MODIFICATION")]
    pub(crate) protein_modification: Option<String>,
    #[serde(rename = "PROTEIN_VARIANT_TYPE")]
    pub(crate) protein_variant_type: Option<String>,
    #[serde(rename = "CHROMOSOMAL_MODIFIED_OBJECT")]
    pub(crate) chromosome: Option<String>,
    #[serde(rename = "CHROMOSOMAL_MODIFICATION")]
    pub(crate) chromosome_modification: Option<String>,
    #[serde(rename = "TRANSCRIPT_HGVS_MODIFIED_OBJECT")]
    pub(crate) transcript_hgvs_modified_object: Option<String>,
    #[serde(rename = "VARIANT_ALLELE_FREQUENCY_IN_TUMOR")]
    pub(crate) variant_allele_frequency_in_tumor: Option<f32>,
    #[serde(rename = "DBSNP")]
    pub(crate) db_snp: Option<String>,
    #[serde(rename = "COPY_NUMBER")]
    pub(crate) copy_number: Option<f32>,
}

impl Variant {
    pub(crate) fn dna_change(&self) -> DnaChange {
        DnaChange::from_str(
            self.chromosome_modification
                .as_ref()
                .unwrap_or(&String::new()),
        )
        .unwrap_or_default()
    }
}

#[allow(clippy::expect_used)]
pub(crate) fn three_letter_protein_modification(short: &str) -> String {
    fn map_value(value: &str) -> String {
        match value {
            "*" => "*",
            "=" => "=",
            "fs" => "fs",
            "F" => "Phe",
            "L" => "Leu",
            "S" => "Ser",
            "Y" => "Tyr",
            "C" => "Cys",
            "W" => "Trp",
            "P" => "Pro",
            "H" => "His",
            "Q" => "Gln",
            "R" => "Arg",
            "I" => "Ile",
            "M" => "Met",
            "T" => "Thr",
            "N" => "Asn",
            "K" => "Lys",
            "V" => "Val",
            "A" => "Ala",
            "D" => "Asp",
            "E" => "Glu",
            "G" => "Gly",
            _ => value,
        }
        .to_string()
    }

    let regex = Regex::new(r"^p\.(?<refA>[*FLSYCWPHQRIMTNKVADEG])?(?<posA>\d+)?(?<sep>_)?(?<refB>[*FLSYCWPHQRIMTNKVADEG])(?<posB>\d+)(?<type>del|ins|delins|dup)?(?<alt>[*=FLSYCWPHQRIMTNKVADEG]+|fs)?$")
        .expect("Invalid regex");

    if let Some(captures) = regex.captures(short) {
        let refa_capture = match captures.name("refA") {
            Some(m) => m.as_str(),
            None => "",
        };
        let posa_capture = match captures.name("posA") {
            Some(m) => m.as_str(),
            None => "",
        };
        let sep_capture = match captures.name("sep") {
            Some(m) => m.as_str(),
            None => "",
        };
        let refb_capture = match captures.name("refB") {
            Some(m) => m.as_str(),
            None => "",
        };
        let posb_capture = match captures.name("posB") {
            Some(m) => m.as_str(),
            None => "",
        };
        let type_capture = match captures.name("type") {
            Some(m) => m.as_str(),
            None => "",
        };
        let alt_capture = match captures.name("alt") {
            Some(m) => m
                .as_str()
                .chars()
                .collect::<Vec<_>>()
                .iter()
                .map(|c| map_value(&c.to_string()))
                .collect::<String>(),
            None => String::new(),
        };
        return format!(
            "p.{}{}{}{}{}{}{}",
            map_value(refa_capture),
            posa_capture,
            map_value(sep_capture),
            map_value(refb_capture),
            posb_capture,
            type_capture,
            alt_capture
        );
    }

    short.to_string()
}

#[derive(Debug, Default, PartialEq)]
pub(crate) struct DnaChange {
    pub(crate) start: String,
    pub(crate) end: String,

    pub(crate) ref_allele: String,
    pub(crate) alt_allele: String,
}

impl FromStr for DnaChange {
    type Err = String;

    #[allow(clippy::expect_used)]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let regexes = [
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?del$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?dup$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)_-?(?P<end>-?\d+)ins(?P<alt>[ACGT]+)$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)_-?(?P<end>-?\d+)delins(?P<alt>[ACGT]+)$")
                .expect("Invalid regex"),
        ];

        for regex in &regexes {
            if let Some(captures) = regex.captures(s) {
                let start = captures["start"].parse::<i128>().unwrap_or_default();
                let end = captures
                    .name("end")
                    .map_or(0, |m| m.as_str().parse::<i128>().unwrap_or_default());
                let ref_allele = captures
                    .name("ref")
                    .map(|m| m.as_str().into())
                    .unwrap_or_default();
                let alt_allele = captures
                    .name("alt")
                    .map(|m| m.as_str().into())
                    .unwrap_or_default();

                let start = if start == 0 {
                    String::new()
                } else {
                    start.to_string()
                };
                let end = if end == 0 {
                    String::new()
                } else {
                    end.to_string()
                };
                return Ok(DnaChange {
                    start,
                    end,
                    ref_allele,
                    alt_allele,
                });
            }
        }
        Err("Invalid DNA change format".to_string())
    }
}

#[cfg(test)]
mod tests {
    use crate::mhguide::*;
    use rstest::rstest;

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_sv_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: Some("p.A123V".to_string()),
                    protein_variant_type: Some("SNV".to_string()),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345678G>A".to_string()),
                    transcript_hgvs_modified_object: Some("c.123C>T".to_string()),
                    variant_allele_frequency_in_tumor: Some(42.42),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None
                }]
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_sv_del_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv_del-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: None,
                    protein_variant_type: None,
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345670_12345678del".to_string()),
                    transcript_hgvs_modified_object: Some("c.120-1_128_1del".to_string()),
                    variant_allele_frequency_in_tumor: Some(42.42),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None
                }]
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_cnv_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/cnv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: None,
                    protein_variant_type: None,
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34)
                }]
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_ref_allele() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(mhguide.variants.len(), 1);
        assert_eq!(
            mhguide.variants.first().unwrap().dna_change().ref_allele,
            "G"
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_alt_allele() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(mhguide.variants.len(), 1);
        assert_eq!(
            mhguide.variants.first().unwrap().dna_change().alt_allele,
            "A"
        );
    }

    #[rstest]
    #[case("c.123C>T",
        DnaChange{ start: "123".to_string(), end: String::new(), ref_allele: "C".to_string(), alt_allele: "T".to_string() }
    )]
    #[case("c.-123C>T",
        DnaChange{ start: "-123".to_string(), end: String::new(), ref_allele: "C".to_string(), alt_allele: "T".to_string() }
    )]
    #[case("c.123_124insA",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: "A".to_string() }
    )]
    #[case("c.123_124del",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    #[case("c.-123_123del",
        DnaChange{ start: "-123".to_string(), end: "123".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    #[case("c.123_124delinsCTGA",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: "CTGA".to_string() }
    )]
    #[case("g.41149933A>G",
        DnaChange{ start: "41149933".to_string(), end: String::new(), ref_allele: "A".to_string(), alt_allele: "G".to_string() }
    )]
    #[case("g.41149933_41150000dup",
        DnaChange{ start: "41149933".to_string(), end: "41150000".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    fn test_dna_change_parsing(#[case] case: &str, #[case] expected: DnaChange) {
        let actual = DnaChange::from_str(case);
        assert_eq!(actual, Ok(expected));
    }

    #[rstest]
    #[case("p.F123G", "p.Phe123Gly")]
    #[case("p.L123F", "p.Leu123Phe")]
    #[case("p.S123L", "p.Ser123Leu")]
    #[case("p.Y123S", "p.Tyr123Ser")]
    #[case("p.C123Y", "p.Cys123Tyr")]
    #[case("p.W123C", "p.Trp123Cys")]
    #[case("p.P123W", "p.Pro123Trp")]
    #[case("p.H123P", "p.His123Pro")]
    #[case("p.Q123H", "p.Gln123His")]
    #[case("p.R123Q", "p.Arg123Gln")]
    #[case("p.I123R", "p.Ile123Arg")]
    #[case("p.M123I", "p.Met123Ile")]
    #[case("p.T123M", "p.Thr123Met")]
    #[case("p.N123T", "p.Asn123Thr")]
    #[case("p.K123N", "p.Lys123Asn")]
    #[case("p.V123K", "p.Val123Lys")]
    #[case("p.A123V", "p.Ala123Val")]
    #[case("p.D123A", "p.Asp123Ala")]
    #[case("p.E123D", "p.Glu123Asp")]
    #[case("p.G123E", "p.Gly123Glu")]
    #[case("p.Y123=", "p.Tyr123=")]
    #[case("p.Y123fs", "p.Tyr123fs")]
    #[case("p.S123_I125delinsF", "p.Ser123_Ile125delinsPhe")]
    #[case("p.S123_I125delinsFE", "p.Ser123_Ile125delinsPheGlu")]
    #[case("p.S123_I125del", "p.Ser123_Ile125del")]
    #[case("p.Y123dup", "p.Tyr123dup")]
    // Examples from Onkostar Notices
    #[case("p.L858R", "p.Leu858Arg")]
    #[case("p.*del*", "p.*del*")]
    #[case("p.V600*", "p.Val600*")]
    // Not mappable - keep as is
    #[case("p.X123X", "p.X123X")]
    #[case("c.123A>C", "c.123A>C")]
    fn test_three_letter_protein_modification(#[case] short: &str, #[case] long: &str) {
        assert_eq!(three_letter_protein_modification(short), long);
    }
}
