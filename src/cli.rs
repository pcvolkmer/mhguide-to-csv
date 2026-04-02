use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about)]
#[command(arg_required_else_help(true))]
pub(crate) struct Cli {
    #[arg(help = "Zu lesende JSON-Datei")]
    pub(crate) input_file: PathBuf,

    #[arg(
        long,
        conflicts_with = "oncogenic",
        help = "Alle Varianten verwenden, nicht nur '(Likely) oncogenic' oder aus 'REPORT_NARRATIVE'"
    )]
    pub(crate) all_variants: bool,

    #[arg(
        long,
        conflicts_with = "all_variants",
        help = "Nur Varianten mit '(Likely) oncogenic' verwenden, keine aus 'REPORT_NARRATIVE'"
    )]
    pub(crate) oncogenic: bool,

    #[arg(
        long,
        conflicts_with = "all_variants",
        conflicts_with = "oncogenic",
        help = "Entferne Artefakte aus 'REPORT_NARRATIVE'"
    )]
    pub(crate) no_artifacts: bool,

    #[arg(
        long,
        conflicts_with = "json",
        help = "Exportiere im XLSX-Format (Excel 2007-365)"
    )]
    pub(crate) xlsx: bool,

    #[arg(
        long,
        conflicts_with = "xlsx",
        help = "Exportiere JSON gemäß DNPM-Datenmodell 2.1"
    )]
    pub(crate) json: bool,
}
