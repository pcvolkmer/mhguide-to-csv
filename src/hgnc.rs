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
    /// Creates a new instance of the `Genes` struct by loading and parsing gene data from a CSV file.
    ///
    /// This function reads a tab-delimited CSV file located in the `../resources/genes.csv` file path
    /// at compile time (via the `include_str!` macro). The file is expected to have the following format:
    /// - Column 0: HGNC ID
    /// - Column 1: Gene symbol
    /// - Column 2: Previous symbols (comma-separated list)
    /// - Column 3: Gene name
    /// - Column 5: Ensembl ID
    ///
    /// Each record from the CSV is mapped to a `Gene` struct, which includes:
    /// - `hgnc_id` (String): The HGNC ID of the gene.
    /// - `symbol` (String): The gene's apprived symbol, with spaces removed.
    /// - `previous_symbols` (Vec<String>): A list of previous symbols for the gene.
    /// - `name` (String): The name of the gene.
    /// - `ensembl_id` (Option<String>): The Ensembl ID, if provided.
    ///
    /// # Returns
    /// A `Genes` instance containing a vector of all parsed `Gene` records.
    ///
    /// # Panics
    /// This function does not explicitly handle panics, but internal usage of `unwrap_or_default` ensures
    /// that malformed or missing CSV data is replaced with defaults when possible. The function panics
    /// if there are issues with the compile-time inclusion of the resource file.
    ///
    /// # Example
    /// ```rust
    /// let genes = Genes::new();
    /// println!("{:?}", genes);
    /// ```
    ///
    /// # Notes
    /// - The delimiter for the CSV file is a tab (`\t`).
    /// - Fields that are missing or empty in the CSV file are replaced with default values.
    /// - File `../resources/genes.csv` will be compiled into the binary.
    ///
    /// # Dependencies
    /// - The `csv` crate is used to parse the CSV file.
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

    /// Finds a `Gene` in the collection by its symbol.
    ///
    /// This method searches the internal `genes` collection for a `Gene`
    /// that has a matching `symbol`. If a gene with the provided symbol
    /// is found, it returns a clone of the `Gene` wrapped in an `Option`.
    /// If no matching gene is found, it returns `None`.
    ///
    /// # Arguments
    ///
    /// * `symbol` - A string slice representing the symbol of the gene to find.
    ///
    /// # Returns
    ///
    /// * `Option<Gene>` - Returns `Some(Gene)` if a gene with the given symbol
    ///   exists in the collection, or `None` if no such gene is found.
    ///
    /// # Example
    ///
    /// ```rust
    /// let genes = Genes::new();
    ///
    /// if let Some(gene) = genes.find_by_symbol("BRCA1") {
    ///     println!("Found gene: {}", gene.symbol);
    /// } else {
    ///     println!("Gene not found");
    /// }
    /// ```
    pub(crate) fn find_by_symbol(&self, symbol: &str) -> Option<Gene> {
        self.genes
            .iter()
            .find(|gene| gene.symbol == symbol)
            .cloned()
    }

    /// Searches for a `Gene` in the collection that has the specified `symbol`
    /// in its list of previous symbols.
    ///
    /// # Arguments
    ///
    /// * `symbol` - A string slice representing the symbol to search for
    /// within the previous symbols of a gene.
    ///
    /// # Returns
    ///
    /// * `Option<Gene>` - Returns `Some(Gene)` if a matching gene is found
    /// whose `previous_symbols` contains the specified `symbol`.
    /// Returns `None` if no matching gene is found.
    ///
    /// # Example
    ///
    /// ```rust
    /// let genes = Genes::new();
    /// if let Some(gene) = genes.find_by_previous_symbol("OLD1") {
    ///     println!("Found gene: {:?}", gene);
    /// } else {
    ///     println!("No gene found with the given previous symbol.");
    /// }
    /// ```
    ///
    /// # Notes
    ///
    /// * This function performs a linear search through the collection of genes.
    /// * The search is case-sensitive.
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
