use csv::ReaderBuilder;
use serde::Deserialize;

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct Gene {
    #[serde(rename = "HGNC ID")]
    pub(crate) hgnc_id: String,
    #[serde(rename = "Approved symbol")]
    pub(crate) symbol: String,
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
                name: record.get(2).unwrap_or_default().to_string(),
                ensembl_id: record.get(4).map(ToString::to_string),
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
}
