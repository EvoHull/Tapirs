# Reproducible Analyses

Reproducibility is important professionally, but the main person it is helping is future-you.

# Exit Strategy
When an experiment is finished you should have an 'exit strategy' checklist to make sure your work is as reproducible as possible. We hope that we have made this achievable in Tapirs

- data provenance
    - data archive is possible
- list of all software, sources, and versions
    - conda export of environment, all software conda installable
- workflow
- human readable reports
- easy to archive

# Software list and versions

The full list of software, their dependencies and version numbers called `environment.yaml` is written to envs/archived_envs directory at the end of the run. This file can be used to reproduce the experimental software conditions.

# Archiving
Snakemake can be asked to [make an archive](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#sustainable-and-reproducible-archiving) of code, config, and input files:
```
snakemake --archive my-workflow.tar.gz
```

# Transparent workflow