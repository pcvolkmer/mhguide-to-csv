use crate::hgnc::Genes;
use crate::mhguide;
use crate::mhguide::{RefGenomeVersion, VariantType, three_letter_protein_modification};
use serde::{Deserialize, Serialize};
use std::sync::LazyLock;

static GENES: LazyLock<Genes> = LazyLock::new(Genes::new);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct SimpleVariantRecord {
    #[serde(rename = "H-Nummer")]
    h_nummer: String,
    #[serde(rename = "Referenz-Genom")]
    ref_genome: String,
    #[serde(rename = "Ergebnis")]
    ergebnis: String,
    #[serde(rename = "Gen")]
    gene: String,
    #[serde(rename = "Genomposition (g.)")]
    genomic_position: String,
    #[serde(rename = "cDNA Nomenklatur (c.)")]
    cdna: String,
    #[serde(rename = "Proteinebene (original)")]
    protein_orig: String,
    #[serde(rename = "Proteinebene Nomenklatur (p.)")]
    protein: String,
    #[serde(rename = "Chromosom")]
    chromosome: String,
    #[serde(rename = "EnsemblID")]
    ensembl_id: String,
    #[serde(rename = "HGNC ID")]
    hgnc_id: String,
    #[serde(rename = "HGNC Name")]
    hgnc_name: String,
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
    #[serde(rename = "Pathogenitätsklasse")]
    classification: String,
}

impl SimpleVariantRecord {
    /// Constructs a `SimpleVariantRecord` from the input variant details.
    ///
    /// This function generates a `SimpleVariantRecord` by extracting and transforming information
    /// from the provided variant, genomic reference version, and H-number (identifier).
    ///
    /// # Arguments
    ///
    /// * `h_number` - A reference to the H-number string, which serves as an identifier.
    /// * `ref_genome_version` - The reference genome version.
    /// * `variant` - A reference to a `mhguide::Variant` object which provides variant information.
    ///
    /// # Returns
    ///
    /// * `SimpleVariantRecord` - An instance of the `SimpleVariantRecord` struct populated with various variant details,
    /// including gene information, genomic changes, protein modifications, and various metadata.
    ///
    /// # Behavior
    ///
    /// * Identifies the corresponding gene based on the variant's gene symbol, falling back
    ///   to previous symbols if necessary.
    /// * Converts raw DNA and protein modifications into standardized representations.
    /// * Extracts genomic position, allelic information, and variant-specific frequencies.
    /// * Retrieves additional metadata such as gene identifiers and database references
    ///   from `hgnc::Genes`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let variant = mhguide::Variant::new();
    /// let record = SimpleVariantRecord::from_variant("H12345", &RefGenomeVersion::HG38, &variant);
    /// println!("{:?}", record);
    /// ```
    pub(crate) fn from_variant(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        variant: &mhguide::Variant,
    ) -> SimpleVariantRecord {
        let gene = match GENES.find_by_symbol(&variant.gene_symbol.clone().unwrap_or_default()) {
            Some(gene) => gene,
            None => GENES
                .find_by_previous_symbol(&variant.gene_symbol.clone().unwrap_or_default())
                .unwrap_or_default(),
        };

        let dna_change = variant.dna_change();

        SimpleVariantRecord {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
            ergebnis: match &variant.display_variant_type {
                Some(variant_type) => variant_type.to_string(),
                None => match &variant.protein_variant_type {
                    Some(variant_type) => variant_type.to_string(),
                    None => VariantType::default().to_string(),
                },
            },
            gene: variant.gene_symbol.clone().unwrap_or_default(),
            chromosome: variant.chromosome.clone().unwrap_or_default(),
            ensembl_id: gene.ensembl_id.unwrap_or_default(),
            hgnc_id: gene.hgnc_id,
            hgnc_name: gene.name,
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
            protein_orig: if variant
                .protein_modification
                .clone()
                .unwrap_or_default()
                .starts_with("p.")
            {
                variant.protein_modification.clone().unwrap_or_default()
            } else {
                String::new()
            },
            protein: if variant
                .protein_modification
                .clone()
                .unwrap_or_default()
                .starts_with("p.")
            {
                three_letter_protein_modification(
                    &variant.protein_modification.clone().unwrap_or_default(),
                )
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
            start: dna_change.start,
            end: dna_change.end,
            ref_allele: dna_change.ref_allele,
            alt_allele: dna_change.alt_allele,
            allelic_frequency: match variant.variant_allele_frequency_in_tumor {
                Some(value) => format!("{value:.2}").replace('.', ","),
                None => String::new(),
            },
            dbsnp: variant.db_snp.clone().unwrap_or_default(),
            classification: variant.classification_name.clone().unwrap_or_default(),
        }
    }

    pub(crate) fn csv_headlines() -> Vec<String> {
        vec![
            "H-Nummer".to_string(),
            "Referenz-Genom".to_string(),
            "Ergebnis".to_string(),
            "Gen".to_string(),
            "Genomposition (g.)".to_string(),
            "cDNA Nomenklatur (c.)".to_string(),
            "Proteinebene (original)".to_string(),
            "Proteinebene Nomenklatur (p.)".to_string(),
            "Chromosom".to_string(),
            "EnsemblID".to_string(),
            "HGNC ID".to_string(),
            "HGNC Name".to_string(),
            "Start".to_string(),
            "Ende".to_string(),
            "Alternative Nucleotide".to_string(),
            "Reference Nucleotide".to_string(),
            "Allelfrequenz (%)".to_string(),
            "dbSNP ID".to_string(),
            "Pathogenitätsklasse".to_string(),
        ]
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct CopyNumberRecord {
    #[serde(rename = "H-Nummer")]
    h_nummer: String,
    #[serde(rename = "Referenz-Genom")]
    ref_genome: String,
    #[serde(rename = "Ergebnis")]
    ergebnis: String,
    #[serde(rename = "Type")]
    cnv_type: String,
    #[serde(rename = "Gen")]
    gene: String,
    #[serde(rename = "Chromosom")]
    chromosome: String,
    #[serde(rename = "EnsemblID")]
    ensembl_id: String,
    #[serde(rename = "HGNC ID")]
    hgnc_id: String,
    #[serde(rename = "HGNC Name")]
    hgnc_name: String,
    #[serde(rename = "Total CN")]
    total_copy_number: String,
    #[serde(rename = "Pathogenitätsklasse")]
    classification: String,
}

impl CopyNumberRecord {
    /// Constructs a `CopyNumberRecord` from the input variant details.
    ///
    /// This function generates a `CopyNumberRecord` by extracting and transforming information
    /// from the provided variant, genomic reference version, and H-number (identifier).
    ///
    /// # Arguments
    ///
    /// * `h_number` - A reference to the H-number string, which serves as an identifier.
    /// * `ref_genome_version` - The reference genome version.
    /// * `variant` - A reference to a `mhguide::Variant` object which provides variant information.
    ///
    /// # Returns
    ///
    /// * `CopyNumberRecord` - An instance of the `CopyNumberRecord` struct populated with related details.
    ///
    /// # Behavior
    ///
    /// * Identifies the corresponding gene based on the variant's gene symbol, falling back
    ///   to previous symbols if necessary.
    /// * Classifies copy number variations (CNVs) as 'loss', 'low level gain', or 'high level gain',
    ///   depending on variant data.
    /// * Retrieves additional metadata such as gene identifiers and database references
    ///   from `hgnc::Genes`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let variant = mhguide::Variant::new();
    /// let record = CopyNumberRecord::from_variant("H12345", &RefGenomeVersion::HG38, &variant);
    /// println!("{:?}", record);
    /// ```
    pub(crate) fn from_variant(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        variant: &mhguide::Variant,
    ) -> CopyNumberRecord {
        let gene = match GENES.find_by_symbol(&variant.gene_symbol.clone().unwrap_or_default()) {
            Some(gene) => gene,
            None => GENES
                .find_by_previous_symbol(&variant.gene_symbol.clone().unwrap_or_default())
                .unwrap_or_default(),
        };

        CopyNumberRecord {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
            ergebnis: match &variant.display_variant_type {
                Some(variant_type) => variant_type.to_string(),
                None => match &variant.protein_variant_type {
                    Some(variant_type) => variant_type.to_string(),
                    None => VariantType::default().to_string(),
                },
            },
            gene: variant.gene_symbol.clone().unwrap_or_default(),
            chromosome: variant.chromosome.clone().unwrap_or_default(),
            ensembl_id: gene.ensembl_id.unwrap_or_default(),
            hgnc_id: gene.hgnc_id,
            hgnc_name: gene.name,
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
            classification: variant.classification_name.clone().unwrap_or_default(),
        }
    }

    pub(crate) fn csv_headlines() -> Vec<String> {
        vec![
            "H-Nummer".to_string(),
            "Referenz-Genom".to_string(),
            "Ergebnis".to_string(),
            "Type".to_string(),
            "Gen".to_string(),
            "Chromosom".to_string(),
            "EnsemblID".to_string(),
            "HGNC ID".to_string(),
            "HGNC Name".to_string(),
            "Total CN".to_string(),
            "Pathogenitätsklasse".to_string(),
        ]
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct BiomarkerRecord {
    #[serde(rename = "H-Nummer")]
    h_nummer: String,
    #[serde(rename = "Referenz-Genom")]
    ref_genome: String,
    #[serde(rename = "Ergebnis")]
    ergebnis: String,
    #[serde(rename = "HRD - Score/Ergebnis")]
    hrd: String,
    #[serde(rename = "MSI - Prozentwert")]
    msi: String,
    #[serde(rename = "TMB - Tumor Mutational Burden")]
    tmb: String,
}

impl BiomarkerRecord {
    pub(crate) fn from_hrd(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        value: f32,
    ) -> BiomarkerRecord {
        BiomarkerRecord {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
            ergebnis: VariantType::HRD.to_string(),
            hrd: format!("{value:.2}").replace('.', ","),
            msi: String::new(),
            tmb: String::new(),
        }
    }

    pub(crate) fn from_msi(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        value: f32,
    ) -> BiomarkerRecord {
        BiomarkerRecord {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
            ergebnis: VariantType::MSI.to_string(),
            hrd: String::new(),
            msi: format!("{value:.2}").replace('.', ","),
            tmb: String::new(),
        }
    }

    pub(crate) fn from_tmb(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        value: f32,
    ) -> BiomarkerRecord {
        BiomarkerRecord {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
            ergebnis: VariantType::TMB.to_string(),
            hrd: String::new(),
            msi: String::new(),
            tmb: format!("{value:.2}").replace('.', ","),
        }
    }

    pub(crate) fn csv_headlines() -> Vec<String> {
        vec![
            "H-Nummer".to_string(),
            "Referenz-Genom".to_string(),
            "Ergebnis".to_string(),
            "HRD - Score/Ergebnis".to_string(),
            "MSI - Prozentwert".to_string(),
            "TMB - Tumor Mutational Burden".to_string(),
        ]
    }
}
