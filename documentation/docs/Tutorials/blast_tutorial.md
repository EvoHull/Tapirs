# Basic blast tutorial

## Reference data
Download the sequences you wish to include in your reference database.

## Build a blast database
There are many online guides to building a blast database, including detailed ones at the NCBI.

## blast parameters
**Minimum percentage identity** `min_perc_ident = "100"` identifies perfect matches but may be a little restrictive unless you know that you have all the potential hits in you database. Lower the value on a few test cases until you feel you are getting too much noise and not enough signal.

**e-value** `min_evalue = "1e-20"` may work for you, again you can explore. Using the NCBI web blast page for a few sequences can sometimes give you an idea of the sort of e-values you would expect for your data.

**number of hits** `descriptions = "50"` ie return maximum of 50 hits is a good default, though you may also change if you wish.

## mlca
Investigating the correct parameters for LCA analysis is an overlooked part of many experimental designs. LCA however is responsible for turning your 50 blast hits into an assigned taxonomy, so it can be exceptionally important.

**bitscore threshold** This will restrict the hits to a certain percentage of the to bitscore. Default=10 (within 10%) which may be a little too strict if you do not have close blast hits.

**identity** Minimum percentage identity value accepted for the hits. "100" restricts to perfect matches only.

**coverage** Percentage coverage of query on top hit. Default = "60"

**majority taxonomic concensus** Percentage match requirement across hit taxonomies. "90" indicates that 1 out of 10 may disagree in its taxonomy without the assigned taxonomy moving back to a deeper common ancestor.

**minimum number of hits** Useful with small reference databases. Default = 2, ie requires at least 2 hits in the blast database. hits = "1" is not a true LCA but rather takes the top hit.

## outputs
