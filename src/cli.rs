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
        help = "Alle Varianten verwenden, nicht nur '(Likely) oncogenic' oder '(Likely) benign'"
    )]
    pub(crate) all_variants: bool,

    #[arg(
        long,
        conflicts_with = "all_variants",
        help = "Nur Varianten mit '(Likely) oncogenic' verwenden, nicht '(Likely) benign'"
    )]
    pub(crate) oncogenic: bool,
}
