use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about)]
#[command(arg_required_else_help(true))]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) command: SubCommand,
}

#[derive(Subcommand)]
pub(crate) enum SubCommand {
    #[command(name = "convert", about = "Konvertiere Daten in angegebenes Format")]
    Convert {
        #[command(flatten)]
        convert_args: ConvertArgs,
    },
    #[command(
        name = "push",
        about = "Veröffentliche und aktualisiere Daten in Onkostar"
    )]
    Push {
        #[command(flatten)]
        convert_args: ConvertArgs,
        #[arg(
            long,
            help = "Onkostar-URL",
            default_value = "http://localhost:8080/onkostar"
        )]
        url: String,
    },
    #[command(
        name = "json-schema",
        about = "Schreibe JSON-Schema für JSON-Export in Standardausgabe"
    )]
    JsonSchema,
    #[command(
        name = "check-report-narrative",
        about = "Prüfe auf Probleme in 'REPORT_NARRATIVE'"
    )]
    CheckReportNarrative {
        #[command(flatten)]
        convert_args: ConvertArgs,
    },
}

#[derive(Args)]
pub(crate) struct ConvertArgs {
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
        help = "Exportiere in 'OS.Molekulargenetik' ähnlichem JSON-Format"
    )]
    pub(crate) json: bool,
}
