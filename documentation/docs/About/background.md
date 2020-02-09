# Background

We developed Tapirs primarily for our own work in the EvoHull group at the University of Hull, UK. It has been common for us to carry out short (Illumina) metabarcoding on environmental or eDNA community samples.

# Motivation
We wanted a transparent workflow that was not dependent on the person who developed the software, could be understood easily by new people in the lab, and enabled us to take reproducibility seriously. We wanted to be able to vary conditions and reanalyse data experimentally, and that meant that changing parameters and re-running should be easy.

We had previously used our own metabarcoding software [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT), written as a large python script, and while very productive we found that supporting it was increasingly difficult, and adding new analyses to this pipeline is not straightforward. We learned a great deal from designing and using metaBEAT, and it greatly increased our interest in workflow management software.

# Other software
In addition to metaBEAT there are many other wonderful software packages for sophisticated analysis of metabarcoding data. QIIME2 is particularly powerful, as is Mothur, and too many others to mention. We felt that other metabarcoding software did not allow us to easily change the analysis software, but rather presented a fixed 'solution'.

We saw a lot value for our research in having a standard workflow manager, with design simplicity, rather than the bigger more powerful but more opaque and fixed packages designed primarily for 16S bacterial sequencing.

# Using the Snakemake workflow management system
We decided that it would be best practice to use a well-designed workflow manager rather than link together our software components in a more ad hoc manner with our own scripts. We chose the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/index.html). This is well-designed software, heavily used and tested across a large community of bioinformatics researchers. We explicitly decided to use a specialist workflow manager, rather than designing a pipeline system ourselves, as we valued the ability to modify and experiment with our workflow, and this is by definition what workflow managers do.

# Reproducibility
We take reproducible science very seriously, and so should you. Not only is it the highest quality science, but it is also the easiest science, saving you time and effort. We hope that the way we have designed Tapirs will give others a very good chance to exactly repeat what you did without significant suffering, and also to use their own data with your exact approach to build upon your work. Please see the [Reproducibility page](../How-To-Guide/reproducibility.md) for details of maximising the reproducibility of your experiment. We would greatly appreciate your thoughts on better reproducibility in metabarcoding.

# Contributions
Very many people have contributed to Tapirs. EvoHull: Many scientists have worked in EvoHull, discussed metabarcoding software and approaches, made suggestions about good practice, explored data analysis, critiqued this work, made suggestions, and taught us how to proceed. 

The code here was physically written by Dave Lunt, Mike Winter, Graham Sellers, and Marco Benucci.
