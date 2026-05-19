use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub(crate) struct OsMolekulargenetik {
    pub(crate) patient_id: String,
    pub(crate) datum: String,
    pub(crate) einsendenummer: String,
    pub(crate) referenzgenom: String,

    pub(crate) molekulargenuntersuchung: Vec<OsMolekulargenUntersuchung>,
    pub(crate) biomarker: Vec<Biomarker>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "ergebnis")]
pub(crate) enum OsMolekulargenUntersuchung {
    #[serde(rename = "P")]
    SimpleVariant {
        untersucht: String,
        genomposition: String,
        cdnanomenklatur: String,
        proteinebenenomenklatur: String,
        evchromosom: String,
        evensemblid: String,
        evhgncid: String,
        evhgncsymbol: String,
        evhgncname: String,
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
        evdbsnpid: String,
        pathogenitaetsklasse: String,
    },
    #[serde(rename = "CNV")]
    CopyNumberVariant {
        untersucht: String,
        copynumbervariation: String,
        cnvchromosom: String,
        cnvensemblid: String,
        cnvhgncid: String,
        cnvhgncsymbol: String,
        cnvhgncname: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvtotalcndouble: Option<f64>,
        pathogenitaetsklasse: String,
    },
    #[serde(rename = "F")]
    RnaFusion {
        untersucht: String,
        fusioniertesgen: String,
        fusionrna5ensemblid: String,
        fusionrna5hgncid: String,
        fusionrna5hgncsymbol: String,
        fusionrna5hgncname: String,
        fusionrna5transcriptid: String,
        fusionrna5exonid: String,
        fusionrna5transposition: Option<u64>,
        fusionrna5strand: String,
        fusionrna3ensemblid: String,
        fusionrna3hgncid: String,
        fusionrna3hgncsymbol: String,
        fusionrna3hgncname: String,
        fusionrna3transcriptid: String,
        fusionrna3exonid: String,
        fusionrna3transposition: Option<u64>,
        fusionrna3strand: String,

        fusionrnareportednumread: Option<u64>,
        pathogenitaetsklasse: String,
    },
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "komplexerbiomarker")]
pub(crate) enum Biomarker {
    #[serde(rename = "HRD")]
    Hrd {
        score: Option<f64>,
        lst: Option<f64>,
        loh: Option<f64>,
        tai: Option<f64>,
        bewertung: Option<MolMarkerBewertung>,
    },
    #[serde(rename = "MSI")]
    Msi {
        seqprozentwert: Option<f64>,
        pcrergebnis: Option<MsiPcrErgebnis>,
        immunergebnismsi: Option<MsiImmunergebnis>,
    },
    #[serde(rename = "TMB")]
    Tmb {
        tumormutationalburden: Option<f64>,
        bewertung: Option<MolMarkerBewertung>,
    },
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize)]
pub(crate) enum MsiPcrErgebnis {
    #[serde(rename = "MSS")]
    Stable,
    #[serde(rename = "H")]
    MsiHigh,
    #[serde(rename = "L")]
    MsiLow,
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize)]
pub(crate) enum MsiImmunergebnis {
    #[serde(rename = "MMP")]
    MmrProficient,
    #[serde(rename = "MMD")]
    MmrDeficient,
}

#[allow(unused)]
#[derive(Debug, Clone, Serialize)]
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
