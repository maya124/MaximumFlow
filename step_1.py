"""
File: step_1.py
Author: Maya Varma
Usage: python3 step_1.py <prefix>
Input: experiment name is the name of the feature selection procedure
Description: Step 1 of the PLINK-LD filtering procedure. Read all bim files and assign each variant a unique identifier
"""
import sys
import numpy as np

def main():
    #check that the file is being properly used
    if (len(sys.argv) != 2):
        print("Please specify an input folder.")
        return
    prefix = sys.argv[1]
    
    chromNums = [str(x) for x in range(1, 23)]
    chromNums.append('X')
    for i in chromNums: 
        result = ""
        counter = 0
        with open(prefix+'.'+i+'.bim', 'r') as f:
            for line in f.readlines():
                counter+=1
                elems = line.split()
                elems[1] = 'var%d' % counter
                result+='\t'.join(elems)+"\n"
        #rewrite bim file with variant names
        with open(prefix+'.'+i+'.bim', 'w') as f:
            f.write(result)
        #create count file with number of variants
        with open(prefix+'.'+i+'.count.txt', 'w') as f:
            f.write(str(counter))

    return

if __name__=="__main__":
    main()
