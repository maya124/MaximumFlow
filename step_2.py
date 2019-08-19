"""
File: step_2.py
Author: Maya Varma
Usage: python step_2.py <input folder> <prefix> <feature set 1> <feature set 2>
Input: experiment name is the name of the feature selection procedure, feature set numbers are the ids of the feature set (1-5)
Description: Calculates pairwise LD values for variants in pairs of feature sets
Output: Python map objects mapping pairs of variants to R^2 and D' values
"""
import itertools
import pickle
import sys
#Load in all of the variants in the cross-fold feature set and map each variant position to its assigned variant ID bim file
def loadAndMapFeatures(set_num, experiment_name, prefix):
    chromToPos = {}
    with open(str(experiment_name)+"/"+str(set_num)+".txt", "r") as f:
        for line in f.readlines():
            variantName = line.split()[0]
            chrom = variantName.split('-')[0]
            pos = variantName.split('-')[1]
            if(chrom not in chromToPos): chromToPos[chrom] = []
            chromToPos[chrom].append(pos)

    #Map each of the variant positions to their assigned variant IDs in the bim file
    chromToId = {}
    for key in chromToPos:
        if(key not in chromToId): chromToId[key] = []
        with open(prefix+"."+str(key)+".bim", 'r') as f:
            for line in f.readlines():
                elems = line.split()
                if(elems[3] in chromToPos[key]): chromToId[key].append(elems[1])
    
    return chromToId

def computeLDScores(experiment_name, prefix, chromToId1, chromToId2, set_num1, set_num2):
    chromToLD_D={}
    pairToLD_D={}
    chromToLD_R2={}
    pairToLD_R2={}
    for key in chromToId1:
        print("\tProcessing Chromosome %s" % key)
        if(key not in chromToLD_D): chromToLD_D[key]=[]
        if(key not in chromToLD_R2): chromToLD_R2[key]=[]
        for pair in list(itertools.product(chromToId1[key], chromToId2[key])):
            filename = prefix+"."+key
            variant0, variant1 = pair
            outfilename=prefix+"."+str(set_num1)+"."+str(set_num2)
            output = get_ipython().getoutput(u'! plink2 --bfile $filename --ld $variant0 $variant1 --out $outfilename')
            if('Error' in output[-1]): continue
            pairToLD_D[(key, variant0, variant1)]=[]
            pairToLD_R2[(key, variant0, variant1)]=[]
            if(len(output)>24 and output[-24]=="Solution #1:"): #two solutions
                chromToLD_D[key].append(float(output[-23].split()[-1]))
                chromToLD_D[key].append(float(output[-11].split()[-1]))
                pairToLD_D[(key, variant0, variant1)].append(float(output[-23].split()[-1]))
                pairToLD_D[(key, variant0, variant1)].append(float(output[-11].split()[-1]))
                chromToLD_R2[key].append(float(output[-23].split()[2]))
                chromToLD_R2[key].append(float(output[-11].split()[2]))
                pairToLD_R2[(key, variant0, variant1)].append(float(output[-23].split()[2]))
                pairToLD_R2[(key, variant0, variant1)].append(float(output[-11].split()[2]))
            else: #single solution
                chromToLD_D[key].append(float(output[-11].split()[-1]))
                pairToLD_D[(key, variant0, variant1)].append(float(output[-11].split()[-1]))
                chromToLD_R2[key].append(float(output[-11].split()[2]))
                pairToLD_R2[(key, variant0, variant1)].append(float(output[-11].split()[2]))
            
    return(chromToLD_D, pairToLD_D, chromToLD_R2, pairToLD_R2)

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def main():
    input_folder= sys.argv[1]
    prefix = sys.argv[2]
    set_num1 = sys.argv[3]
    set_num2 = sys.argv[4]
    print('Currently Running Feature Set Pair %s %s' % (set_num1, set_num2))
    chromToId1 = loadAndMapFeatures(set_num1, input_folder, prefix)
    chromToId2 = loadAndMapFeatures(set_num2, input_folder, prefix)
    chromToLD_D, pairToLD_D, chromToLD_R2, pairToLD_R2 = computeLDScores(input_folder, prefix, chromToId1, chromToId2, set_num1, set_num2)
    save_object(chromToLD_D, "chromToLD_D_"+str(set_num1)+"_"+str(set_num2))
    save_object(pairToLD_D, "pairToLD_D_"+str(set_num1)+"_"+str(set_num2))
    save_object(chromToLD_R2, "chromToLD_R2_"+str(set_num1)+"_"+str(set_num2))
    save_object(pairToLD_R2, "pairToLD_R2_"+str(set_num1)+"_"+str(set_num2))

if __name__=="__main__":
    main()
