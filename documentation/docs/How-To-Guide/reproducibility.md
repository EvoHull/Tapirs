# Reproducible Analyses

Reproducibility is important. The main person it is helping is you.

# Exit Strategy
When an experiment is finished you should have an 'exit strategy' checklist to make sure your work is as reproducible as possible. We hope that we have made this achievable in Tapirs
- data provinance
  + data archive is possible
- list of all software, sources, and versions
  + conda export of environment, all software conda installable
- workflow
- human readable reports
- easy to archive

# Software list and versions

The full list of software, their dependencies and version numbers called `environment.yaml` is written to envs/ directory at the end of the run. This file can be used to reproduce the experimental software conditions.


Snakemake can be asked to archive it all.

Transparent workflow
