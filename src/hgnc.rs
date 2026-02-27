use csv::ReaderBuilder;
use serde::Deserialize;

#[derive(Clone, Debug, Default, Deserialize)]
pub(crate) struct Gene {
    #[serde(rename = "HGNC ID")]
    pub(crate) hgnc_id: String,
    #[serde(rename = "Approved symbol")]
    pub(crate) symbol: String,
    #[serde(rename = "Previous symbols")]
    pub(crate) previous_symbols: Vec<String>,
    #[serde(rename = "Approved name")]
    pub(crate) name: String,
    #[serde(rename = "Ensembl ID(supplied by Ensembl)")]
    pub(crate) ensembl_id: Option<String>,
}

pub(crate) struct Genes {
    genes: Vec<Gene>,
}

impl Genes {
    pub(crate) fn new() -> Genes {
        static GENE_LIST: &str = include_str!("../resources/genes.csv");
        let mut reader = ReaderBuilder::default()
            .delimiter(b'\t')
            .from_reader(GENE_LIST.as_bytes());

        let genes = reader
            .records()
            .flatten()
            .map(|record| Gene {
                hgnc_id: record.get(0).unwrap_or_default().to_string(),
                symbol: record
                    .get(1)
                    .unwrap_or_default()
                    .to_string()
                    .replace(' ', ""),
                previous_symbols: record
                    .get(2)
                    .unwrap_or_default()
                    .split(',')
                    .map(ToString::to_string)
                    .map(|s| s.trim().to_string())
                    .collect(),
                name: record.get(3).unwrap_or_default().to_string(),
                ensembl_id: record.get(5).map(ToString::to_string),
            })
            .collect::<Vec<_>>();

        Genes { genes }
    }

    pub(crate) fn find_by_symbol(&self, symbol: &str) -> Option<Gene> {
        self.genes
            .iter()
            .find(|gene| gene.symbol == symbol)
            .cloned()
    }

    pub(crate) fn find_by_previous_symbol(&self, symbol: &str) -> Option<Gene> {
        self.genes
            .iter()
            .find(|&gene| gene.previous_symbols.contains(&symbol.to_string()))
            .cloned()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::hgnc::Genes;

    #[test]
    fn test_should_find_gene_by_symbol() {
        let genes = Genes::new();
        let gene = genes.find_by_symbol("BRCA1");
        assert!(gene.is_some());
        let gene = gene.unwrap();
        assert_eq!(gene.symbol, "BRCA1");
        assert_eq!(gene.hgnc_id, "HGNC:1100");
        assert_eq!(gene.name, "BRCA1 DNA repair associated");
        assert_eq!(gene.ensembl_id, Some("ENSG00000012048".to_string()));
    }

    #[test]
    fn test_should_find_gene_by_previous_symbol() {
        let genes = Genes::new();
        let gene = genes.find_by_previous_symbol("FAM83H");
        assert!(gene.is_some());
        let gene = gene.unwrap();
        assert_eq!(gene.symbol, "SACK1H");
        assert_eq!(gene.hgnc_id, "HGNC:24797");
        assert_eq!(gene.name, "scaffolding CK1 anchoring protein H");
        assert_eq!(gene.ensembl_id, Some("ENSG00000180921".to_string()));
    }
}
