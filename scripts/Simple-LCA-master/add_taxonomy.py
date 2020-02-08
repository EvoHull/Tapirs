#!/usr/bin/python
"""
Simple-LCA   V1.0    martenhoogeveen@naturalis.nl
This script adds the taxonomy to the BLAST output.
"""
import json, sys, argparse, os
# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description='Add taxonomy to BLAST output')
parser.add_argument('-i', '--blast_input', metavar='BLAST custom outfmt 6 output', dest='blastinput', type=str,help='blast inputfile', required=True)
parser.add_argument('-t', '--taxonomy_reference', metavar='taxonomy reference', dest='rankedlineage', type=str, help='reference json taxonomy file', required=False, nargs='?', default="taxonomy_reference.json")
parser.add_argument('-m', '--merged', metavar='merged taxonids', dest='merged', type=str, help='merged taxon id json', required=False, nargs='?', default="merged_taxonomy.json")
parser.add_argument('-o', '--output', metavar='output', dest='output', type=str, help='output file, BLAST hits with taxonomy', required=False, nargs='?', default="")
args = parser.parse_args()

def merged_taxonomy():
    mergedDict = {}
    with open(args.merged) as merged:
        for taxid in merged:
            a = map(str.strip, taxid.split("|"))
            mergedDict[a[0]]=a[1]
    return mergedDict

def reference_taxonomy():
    taxonomyDict = {}
    with open(args.rankedlineage) as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split("|")
            taxonid = tax[0]
            species = tax[1].strip() if tax[1].strip() else "unknown species"
            genus = tax[3].strip() if tax[3].strip() else "unknown genus"
            family = tax[4].strip() if tax[4].strip() else "unknown family"
            order = tax[5].strip() if tax[5].strip() else "unknown order"
            classe = tax[6].strip() if tax[6].strip() else "unknown class"
            phylum = tax[7].strip() if tax[7].strip() else "unknown phylum"
            kingdom = tax[8].strip() if tax[8].strip() else "unknown kingdom"
            superkingdom = tax[9].strip() if tax[9].strip() else "unknown superkingdom"
            taxonomyDict[str(tax[0].strip())] = {"species":species, "genus":genus, "family":family, "order":order, "class":classe, "phylum":phylum, "kingdom":kingdom,"superkingdom":superkingdom}
    return taxonomyDict

def check_merged_taxonomy(taxid, mergedTaxonDict):
    try:
        a = mergedTaxonDict[taxid]
        return a
    except:
        return taxid

def add_taxonomy(taxonomyDict, mergedTaxonDict):
    outputFile = args.output if args.output else "added_taxonomy_"+str(os.path.basename(args.blastinput))
    with open(args.blastinput) as blasthits, open(outputFile, "a") as output:
        for hit in blasthits:
            taxid = hit.split("\t")[3]
            if taxid == "N/A":
                output.write(hit.strip()+"\t"+"unknown kingdom / unknown phylum / unknown class / unknown order / family / genus / species\n")
            else:
                taxonomydb = taxonomyDict
                try:
                    kingdom = taxonomydb[taxid]["kingdom"]
                    superkingdom = taxonomydb[taxid]["superkingdom"]
                except KeyError:
                    taxid = check_merged_taxonomy(taxid, mergedTaxonDict)
                    kingdom = taxonomydb[taxid]["kingdom"]
                    superkingdom = taxonomydb[taxid]["superkingdom"]

                if kingdom and kingdom != "unknown kingdom":
                    output.write(hit.strip()+"\t"+taxonomydb[taxid]["kingdom"]+" / "+taxonomydb[taxid]["phylum"]+ " / " +taxonomydb[taxid]["class"]+" / "+taxonomydb[taxid]["order"]+" / "+taxonomydb[taxid]["family"]+" / "+taxonomydb[taxid]["genus"]+" / "+taxonomydb[taxid]["species"]+"\n")
                elif superkingdom and superkingdom != "unknown superkingdom":
                    output.write(hit.strip()+"\t"+taxonomydb[taxid]["superkingdom"]+" / "+taxonomydb[taxid]["phylum"]+ " / " +taxonomydb[taxid]["class"]+" / "+taxonomydb[taxid]["order"]+" / "+taxonomydb[taxid]["family"]+" / "+taxonomydb[taxid]["genus"]+" / "+taxonomydb[taxid]["species"]+"\n")
                else:
                    output.write(hit.strip()+"\t"+taxonomydb[taxid]["kingdom"]+" / "+taxonomydb[taxid]["phylum"]+ " / " +taxonomydb[taxid]["class"]+" / "+taxonomydb[taxid]["order"]+" / "+taxonomydb[taxid]["family"]+" / "+taxonomydb[taxid]["genus"]+" / "+taxonomydb[taxid]["species"]+"\n")

def main():
    taxonomyDict = reference_taxonomy()
    mergedTaxonDict = merged_taxonomy()
    add_taxonomy(taxonomyDict, mergedTaxonDict)


if __name__ == "__main__":
    main()