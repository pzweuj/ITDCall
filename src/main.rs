use anyhow::Result;
use clap::Parser;

use itdcall::cli::Cli;
use itdcall::pipeline::run_cli;

fn main() -> Result<()> {
    let cli = Cli::parse();
    run_cli(&cli)
}
