"""
File: step_0.py
Author: Maya Varma
Usage: python step_0.py <input folder>
Input: Input folder contains a group of text files with variants and corresponding coefficient scores for each cross fold. 
Description: Step 0 of the PLINK-LD filtering procedure. Outputs a file called variant_ranges.tsv with all variants from the cross fold analysis.
Output: variant_ranges.tsv, vcf_input.txt
variant_ranges format: CHR START_BP_POS END_BP_POS SOMETHING1 SOMETHING2
"""

import sys
import numpy as np

NUM_FOLDS = 5

def readFolder(folder):
    """
    Reads all cross fold text files in folder and returns a union of variants. 
    Input: folder name
    Output: list of all variants in text files
    """
    allVariants = set()
    for i in range(1, NUM_FOLDS+1):
        with open(folder + str(i) + '.txt', 'r') as f:
            for line in f.readlines():
                allVariants.add(line.split()[0])
    return sorted(list(allVariants))

def writeOutput(folder, allVariants):
    """
    Writes all variants to .tsv file in the format that plink requires.
    """
    result_plink = ""
    result_vcf = ""
    for variant in allVariants:
        chrom = variant.split('-')[0]
        pos = variant.split('-')[1]
        result_plink += chrom + "\t" + pos + "\t" + pos + "\t" + "filler" + "\n";
        result_vcf += chrom + "\t" + pos + "\t" + pos + "\n";
    with open(folder + 'variant_ranges.tsv', 'w') as f:
        f.write(result_plink)
    with open(folder + 'vcf_input.txt', 'w') as f:
        f.write(result_vcf)

def main():
    #check that the file is being properly used
    if (len(sys.argv) != 2):
        print("Please specify an input folder.")
        return
    folder = sys.argv[1]
    allVariants = readFolder(folder)
    writeOutput(folder, allVariants)
    return

if __name__=="__main__":
    main()
