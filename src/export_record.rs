use crate::hgnc::Genes;
use crate::mhguide;
use crate::mhguide::{RefGenomeVersion, three_letter_protein_modification};
use serde::{Deserialize, Serialize};
use std::sync::LazyLock;

static GENES: LazyLock<Genes> = LazyLock::new(Genes::new);

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct Record {
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
    #[serde(rename = "Type")]
    cnv_type: String,
    #[serde(rename = "Total CN")]
    total_copy_number: String,
    #[serde(rename = "Pathogenitätsklasse")]
    classification: String,
}

impl Record {
    pub(crate) fn from_variant(
        h_number: &str,
        ref_genome_version: &RefGenomeVersion,
        variant: &mhguide::Variant,
    ) -> Record {
        let gene = match GENES.find_by_symbol(&variant.gene_symbol.clone().unwrap_or_default()) {
            Some(gene) => gene,
            None => GENES
                .find_by_previous_symbol(&variant.gene_symbol.clone().unwrap_or_default())
                .unwrap_or_default(),
        };

        let dna_change = variant.dna_change();

        Record {
            h_nummer: h_number.to_string(),
            ref_genome: ref_genome_version.to_string(),
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
            classification: variant.classification_name.clone().unwrap_or_default(),
        }
    }
}
