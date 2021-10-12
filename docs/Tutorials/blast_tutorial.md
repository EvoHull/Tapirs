# Basic blast tutorial

## Reference data

Download the sequences you wish to include in your reference database.

## Build a blast database

There are many online guides to building a blast database, including detailed ones at the NCBI.

## blast parameters

**Minimum percentage identity:** `BLAST_min_perc_ident` Default = 75. Blast hits with an identity of 75% are kept. This value should be equal to or lower than `MLCA_identity` (see below).

**e-value:** `BLAST_min_evalue` Default = 1e-20. This may work for you, again you can explore how stringent to be with this parameter. Using the NCBI web blast page for a few sequences can sometimes give you an idea of the sort of e-values you would expect for your data.

**number of hits:** `BLAST_max_target_seqs:` Default = 50. Returning a maximum of 50 hits per query is a good default, though you may also change if you wish.

## mlca

Investigating the correct parameters for LCA analysis is an overlooked part of many experimental designs. LCA however is responsible for turning your 50 blast hits into an assigned taxonomy, so it can be exceptionally important. MLCA is short for Majority Lowest Common Ancestor, and our script allows you to modify several parameters of the LCA estimation, including whether to accept the majority taxonomy returned, or require ALL hits to have the same taxonomy.

**bitscore threshold:** `MLCA_bitscore` This will restrict the hits to a certain percentage of the top bitscore. Default = 2 (within 2% of top hit) which may be a little too strict if you do not have close blast hits. Experiment with this value.

**identity:** `MLCA_identity` Minimum percentage identity value accepted for the hits. Default = 98. "100" restricts to perfect matches only. This value should be equal to or higher than `BLAST_min_perc_ident` (see above).

**coverage:** `MLCA_identity` Minimum percentage coverage value accepted for the hits. Default = 90.

**majority taxonomic concensus:** `MLCA_majority` Percentage match requirement across hit taxonomies. Default = 80. "80" indicates that 2 out of 10 hits may disagree at a taxonomic rank without the assigned taxonomy moving back to a higher shared taxonomic rank.

**minimum number of hits:** `MLCA_hits` Minimum hits required to perform LCA. Default = 1. "1" is not a true LCA but rather takes the single top hit and assigns it to the query - useful with small reference databases. Increasing this vlue will ensure all assignemnts are via true LCA.
