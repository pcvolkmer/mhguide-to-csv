use schemars::JsonSchema;
use serde::Serialize;

#[derive(Debug, Clone, Serialize, JsonSchema)]
#[schemars(
    deny_unknown_fields,
    title = "mhguide-to-csv and xapi JSON schema for OS.Molekulargenetik"
)]
pub(crate) struct OsMolekulargenetik {
    pub(crate) patient_id: String,
    pub(crate) datum: chrono::NaiveDate,
    pub(crate) einsendenummer: String,
    pub(crate) referenzgenom: String,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) molekulargenuntersuchung: Option<Vec<OsMolekulargenUntersuchung>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) biomarker: Option<Vec<Biomarker>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
#[serde(tag = "ergebnis")]
pub(crate) enum OsMolekulargenUntersuchung {
    #[serde(rename = "P")]
    SimpleVariant {
        untersucht: String,
        genomposition: String,
        cdnanomenklatur: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        proteinebenenomenklatur: Option<String>,
        evchromosom: String,
        evensemblid: String,
        evhgncid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        evhgncsymbol: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        evhgncname: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        evstart: Option<u64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        evende: Option<u64>,
        evaltnucleotide: String,
        evrefnucleotide: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        evreaddepth: Option<u64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        allelfrequenz: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        evdbsnpid: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pathogenitaetsklasse: Option<String>,
    },
    #[serde(rename = "CNV")]
    CopyNumberVariant {
        untersucht: String,
        copynumbervariation: String,
        cnvchromosom: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvensemblid: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvhgncid: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvhgncsymbol: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvhgncname: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvtotalcndouble: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pathogenitaetsklasse: Option<String>,
    },
    #[serde(rename = "F")]
    RnaFusion {
        untersucht: String,
        fusioniertesgen: String,

        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna5ensemblid: Option<String>,
        fusionrna5hgncid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna5hgncsymbol: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna5hgncname: Option<String>,
        fusionrna5transcriptid: String,
        fusionrna5exonid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna5transposition: Option<u64>,
        fusionrna5strand: String,

        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna3ensemblid: Option<String>,
        fusionrna3hgncid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna3hgncsymbol: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna3hgncname: Option<String>,
        fusionrna3transcriptid: String,
        fusionrna3exonid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna3transposition: Option<u64>,
        fusionrna3strand: String,

        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrnareportednumread: Option<u64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pathogenitaetsklasse: Option<String>,
    },
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
#[serde(tag = "komplexerbiomarker")]
pub(crate) enum Biomarker {
    #[serde(rename = "HRD")]
    Hrd {
        #[serde(skip_serializing_if = "Option::is_none")]
        score: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        lst: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        loh: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        tai: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        bewertung: Option<MolMarkerBewertung>,
    },
    #[serde(rename = "MSI")]
    Msi {
        #[serde(skip_serializing_if = "Option::is_none")]
        seqprozentwert: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pcrergebnis: Option<MsiPcrErgebnis>,
        #[serde(skip_serializing_if = "Option::is_none")]
        immunergebnismsi: Option<MsiImmunergebnis>,
    },
    #[serde(rename = "TMB")]
    Tmb {
        #[serde(skip_serializing_if = "Option::is_none")]
        tumormutationalburden: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        bewertung: Option<MolMarkerBewertung>,
    },
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub(crate) enum MsiPcrErgebnis {
    #[serde(rename = "MSS")]
    Stable,
    #[serde(rename = "H")]
    MsiHigh,
    #[serde(rename = "L")]
    MsiLow,
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub(crate) enum MsiImmunergebnis {
    #[serde(rename = "MMP")]
    MmrProficient,
    #[serde(rename = "MMD")]
    MmrDeficient,
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub(crate) enum MolMarkerBewertung {
    #[serde(rename = "H")]
    High,
    #[serde(rename = "M")]
    Intermediate,
    #[serde(rename = "L")]
    Low,
    #[serde(rename = "X")]
    Unknown,
}
