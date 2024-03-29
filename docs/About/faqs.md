![tapirs_logo](../images/faq.png)

# Frequently Asked Questions

## Why is it called Tapirs? What does it stand for?

Tapirs is not an acronym, we just called it Tapirs because an easy name helps everyone, and tapirs are cool.

## Why did you make this and not use X

Nothing met our needs, see the [background](background.md) page for discussion. We wanted an extensible workflow, built by a specialist workflow manager, that would allow real reproducibility.

## Can you add in my favourite software for me?

No, sorry. We really hope however that choosing a workflow manager (Snakemake) means that this is much easier than for any other metabarcoding package. Look at the [extending Tapirs How-To guide](../Running-Tapirs/extending.md).

## I've added some tools, do you want my changes?

Yes please, thank you. Do this via a pull request on github.

## There is an error/hole in the documentation

You should just be able to edit and fix it on Github. Open the page and click the pencil icon, top right, to edit. Thanks for your help. If you think that parts of the documentation are lacking, and you would like to significantly extend them, we also welcome that. If you prefer to discuss it before you start then please do get in touch. We build the documentation with [mkdocs](https://www.mkdocs.org/), and you might find it useful to install this program.

## It works on one machine but not on another

Run-anywhere portability is a real challenge in bioinformatics, not just for Tapirs. Tapirs has been tested successfully on OSX, and Linux Ubuntu, but problems may still arise. Like most bioinformatics software it is harder but still possible to run on MS Windows machines. We think that a Linux virtual machine or the new Windows Subsystem for Linux are most promising.

Conda environments have proved really useful for us, but they aren't perfect. As a start just check that both machines are running the most recent version of Conda and you have the most recent version of Tapirs. Make sure Conda is activated to the correct environment - `tapirs` not `base`. You could also take the very detailed exported environment (reports/archived_envs) from the successful machine and use it to create a new environment on the not-running machine (`conda env create --file environment.yaml`). Activate this environment. Both machines should then have exactly the same software installed.

## It is incredibly slow, it's been 10 mins and hasn't started yet

This is normal the first time that you run the software as it will need to find, download, and install all the software and their dependencies. This will only happen one time however. In order to avoid this frustration you could build a conda environment of all the software before first run as described in the installation instructions and problem-solving section.

## Why can't I point it at data somewhere else than the data dir?

You can do this with minor modifications to the config.yaml file. We don't generally recommend this approach however as it makes it difficult to keep together the resources for the entire experiment, and to make a suitable archive at the end.

## Is vsearch creating cluster OTUs?

The default configuration of Tapirs (at time of writing) uses **exact sequence variants** (ESVs) not cluster OTUs, though this cluster based approach is also available.

## Who do I contact about the software?

The best way is to flag this on GitHub. If you need to get in touch personally then contact Dave Lunt on Gmail.

## Is it licensed?

Yes and no, all parts that we have written are free under CC0 public domain and you may do as you wish. You may consider this a licence if you wish. Yes, you may include whatever you wish into your own software, we encourage this. All software components that we have included (rather than written, see environment.yaml) have been chosen because they have permissive open source licences, but you may wish to check the details for yourself.

## Is there a citation?

We greatly appreciate you citing us wherever you are able. A suitable citation might be:
```
Title: Tapirs, an extensible workflow for reproducible metabarcoding
Authors:
URL: https://github.com/EvoHull/Tapirs
```

Please also cite the analysis software that you have used. An appropriate way to do this would be to include references to the original software: "A reproducible metabarcoding workflow was implemented in Tapirs [1] using blast [2] and Kraken2 [3]."
