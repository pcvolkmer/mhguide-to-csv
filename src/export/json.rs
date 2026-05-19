use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub(crate) struct OsMolekulargenetik {
    #[serde(skip_serializing_if = "String::is_empty")]
    pub(crate) patient_id: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub(crate) datum: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub(crate) einsendenummer: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub(crate) referenzgenom: String,

    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) molekulargenuntersuchung: Vec<OsMolekulargenUntersuchung>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) biomarker: Vec<Biomarker>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "ergebnis")]
pub(crate) enum OsMolekulargenUntersuchung {
    #[serde(rename = "P")]
    SimpleVariant {
        #[serde(skip_serializing_if = "String::is_empty")]
        untersucht: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        genomposition: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cdnanomenklatur: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        proteinebenenomenklatur: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evchromosom: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evensemblid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evhgncid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evhgncsymbol: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evhgncname: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        evstart: Option<u64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        evende: Option<u64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        evaltnucleotide: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        evrefnucleotide: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        evreaddepth: Option<u64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        allelfrequenz: Option<f64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        evdbsnpid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        pathogenitaetsklasse: String,
    },
    #[serde(rename = "CNV")]
    CopyNumberVariant {
        #[serde(skip_serializing_if = "String::is_empty")]
        untersucht: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        copynumbervariation: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cnvchromosom: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cnvensemblid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cnvhgncid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cnvhgncsymbol: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        cnvhgncname: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        cnvtotalcndouble: Option<f64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        pathogenitaetsklasse: String,
    },
    #[serde(rename = "F")]
    RnaFusion {
        #[serde(skip_serializing_if = "String::is_empty")]
        untersucht: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusioniertesgen: String,

        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5ensemblid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5hgncid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5hgncsymbol: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5hgncname: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5transcriptid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5exonid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna5transposition: Option<u64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna5strand: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3ensemblid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3hgncid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3hgncsymbol: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3hgncname: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3transcriptid: String,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3exonid: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrna3transposition: Option<u64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        fusionrna3strand: String,

        #[serde(skip_serializing_if = "Option::is_none")]
        fusionrnareportednumread: Option<u64>,
        #[serde(skip_serializing_if = "String::is_empty")]
        pathogenitaetsklasse: String,
    },
}

#[derive(Debug, Clone, Serialize)]
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
