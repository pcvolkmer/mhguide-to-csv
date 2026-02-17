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
        help = "Alle Varianten verwenden, nicht nur '(Likely) oncogenic' oder '(Likely) benign'"
    )]
    pub(crate) all_variants: bool,
}
