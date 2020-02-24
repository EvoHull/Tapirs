Advice on extending the workflow by adding tools, and other approaches.

## Rationale
What is the best approach for your data analysis?  Whatever you think is the answer you will need to gather evidence. Workflow systems like snakemake allow relatively easy addition and extension of the analyses. This is very useful in comparative analyses. If a new paper claims taxonomic assignment method Y is better than method X then you may wish to carry out both and compare the results, while standardising the other parts of you workflow.

## Adding tools
Although familiarity with Snakemake is helpful, in many situations you may be able to add a tool by modifying an existing rule.

First make sure that your new method produces the results you expect when run at the command line. If it doesn't work at the command line it won't work in Snakemake.

Snakemake rules have 3 integral parts; input, output, and the command to turn the first into the second.

## Common problems
Teaching and trouble-shooting Snakemake and bioinformatics are beyond the scope of this document. A couple of pointers however can save a lot of time.

Most problems are because you have a typo

If you have multiple lines of input or output each line except the last must finish in a comma.

```
rule test:
  input:
    file1: "firstfile.fasta",
    file2: "secondfile.fasta"
  output:
    "allseqs.fasta"
  shell:
    "cat {input.file1} {input.file2} >> {output}"
```

## Contributing your improvements to Tapirs
We would love to hear from you about the improvements you've made. A pull-request for your git branch would probably be best.

Why not help improve the documentation? A new and slightly different tutorial is always welcome.
