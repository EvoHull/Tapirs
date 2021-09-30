# Reproducible Analyses

Reproducibility is important professionally, but the main person it is helping is future-you.

## Exit Strategy

When an experiment is finished you should have an 'exit strategy' checklist to make sure your work is as reproducible as possible. We hope that we have made this achievable in Tapirs

- data provenance
  - data archive is possible
- list of all software, sources, and versions
  - conda export of environment, all software conda installable
- workflow
- human readable reports
- easy to archive

## Reports

You can generate an overall snakemake report on what was run and the provinance of the data for each results with the command:
`snakemake --report reports/snakemake_report.html`

Other reports are written in subdirectories in `reports/` by the analysis programs

## Software list and versions

The full list of software, their dependencies and version numbers is written to `envs/archived_envs/environment.yaml` at the end of the run. This file can be used to reproduce the experimental software conditions.

## Archiving

Snakemake can be asked to [make an archive](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#sustainable-and-reproducible-archiving) of code, config, and input files: `snakemake --archive my-workflow.tar.gz`
