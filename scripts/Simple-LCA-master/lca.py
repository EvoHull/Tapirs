#!/usr/bin/python
"""
Simple-LCA   V1.0    martenhoogeveen@naturalis.nl
This script finds the lowest common ancestor of multiple BLAST hits is a simple way.
For every query only a top percentage is used for finding the lca. These top hits are determined by the top hits
bitscore treshold. The blast hits need to have a bitscore higher then the treshold percentage of the hit with the
highest bitscore. So in other words, if the highest bitscore is 200 and your treshold is 8% only hits with a
bitscore higher then 184 are used. The other tresholds are identity and coverage. This script finds the lca of the
the BLAST hits that passes all tree tresholds. There is no specific order of checking the tresholds, they are all
checked at the same time. There is also an option to not find the lca but output a species level identification if
the best hit passes the tophit tresholds. So if the best blast hit is above the top hit tresholds it will be determined
with a species name in the output. If it is lower then the top hit tresholds the highest possible taxon in genus.
The hits can also be filtered by taxonomy, this option in introduced because of wrong determinations. Sequences can be
from the species 'environmental sample' or 'uncultered sample'. The filter hits (-fh) parameter is a comma
separated parameter. If the string in the parameter occurs in the taxonomy column the hit will be removed and not used
in later steps. The filter lca hits (-flh) parameter is a comma separated parameter if the parameter occurs in a taxon
rank it will be ignored while finding the lca.
"""
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input_file', metavar='galaxy blast output', dest='input', type=str,
			help='input data in galaxy blast format', default='', required=True)
parser.add_argument('-o', '--output_file', metavar='output file', dest='output', type=str,
			help='results file in tabular', required=True)
parser.add_argument('-b', '--bitscore', metavar='bitscore top percentage treshold', dest='top', type=str,
			help='top hits to find the lowest common ancestor', required=True)
parser.add_argument('-id', metavar='identity', dest='id', type=str,
			help='identity treshold', required=True)
parser.add_argument('-cov', metavar='coverage', dest='cov', type=str,
			help='coverage treshold', required=True)
parser.add_argument('-t','--tophit', metavar='tophit', dest='tophit', type=str,
			help='Check de best hit first, if it is above the gives treshold the tophit will become the output', required=False, choices=['no', 'yes'], nargs='?', default='no')
parser.add_argument('-tid', metavar='top_hit_identity', dest='topid', type=str,
			help='identity treshold for the tophit', required=False, default='100')
parser.add_argument('-tcov', metavar='top_hit_coverage', dest='topcoverage', type=str,
			help='query coverage treshold for the tophit', required=False,  default='100')
parser.add_argument('-fh', metavar='filter hits', dest='filterHitsParam', type=str,
			help='filter out hits that contain unwanted taxonomy', required=False, default="",nargs='?')
parser.add_argument('-flh', metavar='filter lca hits', dest='filterLcaHits', type=str,
			help='ignore this in de lca determination', required=False, default="",nargs='?')
args = parser.parse_args()

def filter_check(filterParam, line):
    filterHitsParam = filterParam.strip()
    if filterParam[-1] == ",":
        filterHitsParam = filterParam[:-1]
    filter = filterHitsParam.split(",")
    for x in filter:
        if x.lower() in line.lower():
            a = False
            break
        else:
            a = True
    return a

def remove_hits(otu):
        filteredOtu = []
        for line in otu:
            if filter_check(args.filterHitsParam, line[-1]):
                filteredOtu.append(line)
        return filteredOtu

def remove_taxon(zippedTaxonomy):
    filteredZipper = []
    for level in zippedTaxonomy:
        filteredRank = []
        for item in level:
            if filter_check(args.filterLcaHits, item):
                filteredRank.append(item)
        if not filteredRank:
            filteredRank.append("no identification")
        filteredZipper.append(filteredRank)
    return filteredZipper

def check_best_hit(otu):
    if float(otu[0][4]) >= float(args.topid) and float(otu[0][5]) >= float(args.topcoverage):
        taxonomy = map(str.strip, otu[0][-1].split(" / "))
        return otu[0][0] + "\tspecies\t" + taxonomy[-1] + "\t" + "\t".join(taxonomy) + "\tbest hit\n"
    else:
        return False

def calculate_bitscore_treshold(otu):
    highestScore = 0
    for x in otu:
        bitscore = float(x[7])
        highestScore = bitscore if bitscore > highestScore else highestScore
    topTreshold = float(highestScore) * (1 - (float(args.top)/100))
    return topTreshold

def zip_taxonomy_column(otu, topTreshold):
    taxons = []
    for tax in otu:
        if float(tax[7]) >= topTreshold and float(tax[4]) >= float(args.id) and float(tax[5]) >= float(args.cov):
            taxons.append(map(str.strip, tax[-1].split(" / ")))
    #use zip function for all* taxon lists
    zippedTaxonomy = zip(*taxons)
    return zippedTaxonomy

def find_lca(zippedTaxonomy):
    count = 0
    taxonomy = []
    if zippedTaxonomy:
        zippedTaxonomy[-1] = "no identification"
    for y in zippedTaxonomy:
        if len(set(y)) > 1:
            break
        count += 1
        taxonomy.append(list(set(y))[0])
    t = True
    newTaxonomy = []
    for rev in reversed(taxonomy):
        if rev == "no identification" and t:
            pass
        else:
            t = False
            newTaxonomy.append(rev)
    return list(reversed(newTaxonomy))

def generate_output_line(taxonomy, otu):
    taxonLevels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxonLevel = len(taxonomy) - 1
    while len(taxonomy) < 7:
        taxonomy.append("no identification")

    if all(rank == "no identification" for rank in taxonomy):
        return otu[0][0] + "\tno identification\tno identification\t" + "\t".join(taxonomy).strip() + "\tno identification\n"
    else:
        if taxonLevel == 6:
            taxonomy[-1] = ""
            return otu[0][0] + "\t" + taxonLevels[taxonLevel-1] + "\t" + taxonomy[taxonLevel-1] + "\t" + "\t".join(taxonomy).strip()+"\tno identification\tlca\n"
        else:
            return otu[0][0] + "\t" + taxonLevels[taxonLevel] + "\t" + taxonomy[taxonLevel] + "\t" + "\t".join(taxonomy) + "\tlca\n"

def get_lca(otu):
    #find highest bitscore and calculate lowest bitscore treshold
    topTreshold = calculate_bitscore_treshold(otu)

    #place the taxon column in lists for the zip function
    zippedTaxonomy = zip_taxonomy_column(otu, topTreshold)

    #filter the taxonomy levels, taxons that match 'filter lca hits' parameters are removed
    if args.filterLcaHits:
        zippedTaxonomy = remove_taxon(zippedTaxonomy)
    #find the lca and make the line for the output file
    outputLine = generate_output_line(find_lca(zippedTaxonomy), otu)

    return outputLine

def determine_taxonomy(otu):
    if args.filterHitsParam.strip():
        otu_filtered = remove_hits(otu)
    with open(args.output, "a") as output:
        if otu_filtered:
            if args.tophit == "yes":
                bestHit = check_best_hit(otu_filtered)
                if bestHit:
                    pass
                    output.write(bestHit)
                else:
                    resultingTaxonomy = get_lca(otu_filtered)
                    output.write(resultingTaxonomy)
            else:
                resultingTaxonomy = get_lca(otu_filtered)
                output.write(resultingTaxonomy)
        else:
            taxonomy = ["no identification", "no identification", "no identification", "no identification",
                        "no identification", "no identification", "no identification"]
            output.write(otu[0][0] + "\tno identification\tno identification\t" + "\t".join(taxonomy).strip() + "\tfiltered out\n")

def linecount():
    i = 0
    with open(args.input, "r") as f:
        for i, l in enumerate(f):
            pass
    return i

def write_header():
    """
    Write a header line to the output file
    """
    with open(args.output, "a") as output:
        output.write("#Query\t#lca rank\t#lca taxon\t#kingdom\t#phylum\t#class\t#order\t#family\t#genus\t#species\t#method\n")

def lca():
    """
    This method loops trough the BLAST output and all the hits per otu will be the input for the get_lca method.
    Every hit is stored in the otuLines variable.
    """
    write_header()
    lastLineCount = linecount()
    with open(args.input, "r") as input:
        otuList = []
        otuLines = []
        for num, line in enumerate(input):
            if line.split("\t")[0] not in otuList and num != 0 or num == lastLineCount:
                if num == lastLineCount:
                    otuList.append(line.split("\t")[0])
                    otuLines.append(line.split("\t"))
                # find the lca for the query
                determine_taxonomy(otuLines)
                otuList = []
                otuLines = []
                otuList.append(line.split("\t")[0])
                otuLines.append(line.split("\t"))
            else:
                otuList.append(line.split("\t")[0])
                otuLines.append(line.split("\t"))

def main():
    lca()

if __name__ == "__main__":
    main()
