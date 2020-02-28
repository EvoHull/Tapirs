![tapirs_logo](../images/faq.png)

Frequently Asked Questions

# Why is it called Tapirs? What does it stand for?
Tapirs is not an acronym, we just called it Tapirs because an easy name helps everyone, and tapirs are cool.

# Why did you make this and not use X
Nothing met our needs, see the [background](background.md) page for discussion. We wanted an extensible workflow, built by a specialist workflow manager, that would allow real reproducibility.

# Can you add in my favourite software for me?
No, sorry. We really hope however that choosing a workflow manager (Snakemake) means that this is much easier than for any other metabarcoding package. Look at the [extending Tapirs How-To guide](../How-To-Guide/extending.md).

# I've added some tools, do you want my changes?
Yes please. Do this via a pull request on github.

# There is an error/hole in the documentation
You should just be able to edit and fix it on Github. Open the page and click the pencil icon, top right, to edit. Thanks for your help. If you think that parts of the documentation are lacking, and you would like to significantly extend them, we also welcome that. If you prefer to discuss it before you start then please do get in touch.

# It works on one machine but not on another
Run-anywhere portability is a real challenge in bioinformatics, not just for Tapirs. Conda environments have proved really useful for us, but they aren't perfect. As a start just check that both machines are running the tapirs environment, and have not returned to the base environment. You could also take the very detailed exported environment (reports/archived_envs) from the successful machine and use it to create a new environment on the not-running machine. They should then have exactly the same software installed. Also it is worth double-checking both have the exact same version of Tapirs.

# Who do I contact about the software?
The best way is to flag this on GitHub. If you need to get in touch personally then contact Dave Lunt on Gmail.

# Is it licensed?
No, all parts that we have written are free under CC0 public domain and you may do as you wish. Yes, you may include whatever you wish into your own software, we encourage this. All software components that we have included (rather than written, see environment.yaml) have been chosen because they have permissive open source licences, but you may wish to check the details for yourself.

# Is there a citation
We greatly appreciate you citing us wherever you are able. A suitable citation might be:
```
Title: Tapirs, an extensible workflow for reproducible metabarcoding
Authors:
URL: https://github.com/davelunt/Tapirs
```

An appropriate way to do this would be to include references to the original analysis software: "A reproducible metabarcoding workflow was implemented in Tapirs [1] using blast [2], Kraken2 [3], SINTAX [4] and Krona [5]."
