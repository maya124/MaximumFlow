import numpy as np
import os
from sklearn import linear_model
from sklearn.model_selection import train_test_split, cross_val_score, KFold, GroupKFold, GroupShuffleSplit, ShuffleSplit
from sklearn.metrics import f1_score, make_scorer, roc_curve
from sklearn.metrics import roc_auc_score
import sklearn.metrics as skm
import random
import csv
import scipy
import scipy.sparse as sp
from scipy.stats import chi2_contingency

class featureMatrix(object):
    def __init__(self, mat, samples, features):
        self.mat = mat
        self.samples = np.array(samples)
        self.features = features
        
    def getMat(self):
        return(self.mat)
    
    def getSamples(self):
        return(self.samples)
    
    def getFeatures(self):
        return(self.features)
    
    def filterSamplesByName(self, samplesToInclude):
        include = [(x in samplesToInclude) for x in self.samples]
        self.mat = self.mat[include, :]
        self.samples = self.samples[include]
        
    def nSamples(self):
        return(len(self.samples))
    
    def getSampleSet(self, indices):
        return(featureMatrix(self.mat[indices, :], self.samples[indices,], self.features))

    def getFeatureSet(self, feats):
        indices = []
        matfeats = self.getFeatures()
        for i in range(len(feats)):
            indices.append(matfeats.index(feats[i]))
        return(featureMatrix(self.mat[:, indices], self.samples, np.take(self.features,indices)))

    
#Returns: list of binarized phenotypes for each sample (1=autism, 0=control)
def getPhenotypeData(featMat):
    samples = list(featMat.getSamples())
    y = np.zeros(len(samples))    
    with open('phenotype.csv', 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            sample = row[0].split(',')
            if(sample[1] in samples):
                if(sample[5]=='2'): diag = 1
                else: diag = 0
                y[samples.index(sample[1])] = diag
    print ("Phenotype Data:", y)
    return(y)

#Input: featType1, featType2 (type of noncoding features (PSP and Autism respectively), matches folder names in data/)
#       allASD (If TRUE, the final matrix contains ~4000 ASD patients, with control marked 0 and case marked 1. 
#               If FALSE (default), the final matrix contains ~2000 ASD case patients (1) and ~300 PSP control patients (0)
#Returns: binary feature matrix, with samples as rows and variants as columns
def getFeatureIntersectionDataSparse_PSP(featType1, featType2, allASD=False):
        total = 0
        features_psp = [] #stores all variants present in PSP patients
        features_autism = [] #stores all variants presnet in Autism patients
        samples = []
        
        #Iterate through all saved .npz files and find intersection of variant features (in order to initialize matrix)
        for file in os.listdir('../'+featType2+'_vcf/feature_matrix/'):
            if file.endswith(".npz"):
                file_path = os.path.join('../'+featType2+'_vcf/feature_matrix/', file)
                chrom = np.load(file_path)
                features_autism.extend([elem for elem in chrom['col_labels']])
                if(total==0):
                    samples = [x.decode("ascii") for x in chrom['row_labels']] 
                    total = 1
        divIndex = len(samples)-1
        print ('Matrix Division Between PSP and iHart occurs at index: ', divIndex)
        
        total = 0
        for file in os.listdir('../'+featType1+'_vcf/feature_matrix/'):
            if file.endswith(".npz"):
                file_path = os.path.join('../'+featType1+'_vcf/feature_matrix/', file)
                chrom = np.load(file_path)
                features_psp.extend([elem for elem in chrom['col_labels']])
                if(total==0):
                    samples.extend([x.decode("ascii") for x in chrom['row_labels']])
                    total = 1

        print ('Number of iHart Variants (including Autism and Control): ', len(features_autism))
        print ('Number of PSP Variants: ', len(features_psp))
        features = list(set(features_autism).intersection(set(features_psp)))

        index = 0
        featuresIndexMap = {}
        allFeatures = []
        for elem in set(features):
            featuresIndexMap[elem] = index
            allFeatures.append(elem)
            index += 1
        features = allFeatures

        print ('Number of Variants in the Intersection: ', len(features))
        matrix = np.zeros((len(samples),len(features)), dtype=np.int8)

        #Now that we have the correct size for the binary matrix, we will load in all of the data
        for file in os.listdir('../'+featType2+'_vcf/feature_matrix/'):
            if file.endswith(".npz"):
                file_path = os.path.join('../'+featType2+'_vcf/feature_matrix/', file)
                print ('Reading file:', file_path)
                chrom = np.load(file_path)
                binMatrix = sp.csr_matrix((chrom['data'], chrom['indices'], chrom['indptr']), shape=chrom['shape'])
                binVariants = [elem for elem in chrom['col_labels']]
                matrixIndices = []
                binVariantIndices = []
                for i in range(len(binVariants)):
                    if(binVariants[i] in featuresIndexMap):
                        matrixIndices.append(featuresIndexMap[binVariants[i]])
                        binVariantIndices.append(i)
                matrix[0:divIndex+1,matrixIndices] = binMatrix.todense()[:,binVariantIndices]
        
        for file in os.listdir('../'+featType1+'_vcf/feature_matrix/'):
            if file.endswith(".npz"):
                file_path = os.path.join('../'+featType1+'_vcf/feature_matrix/', file)
                print ('Reading file:', file_path)
                chrom = np.load(file_path)
                binMatrix = sp.csr_matrix((chrom['data'], chrom['indices'], chrom['indptr']), shape=chrom['shape'])
                binVariants = [elem for elem in chrom['col_labels']]
                matrixIndices = []
                binVariantIndices = []
                for i in range(len(binVariants)):
                    if(binVariants[i] in featuresIndexMap):
                        matrixIndices.append(featuresIndexMap[binVariants[i]])
                        binVariantIndices.append(i)
                matrix[divIndex+1:,matrixIndices] = binMatrix.todense()[:,binVariantIndices]
                
        featMat = matrix
        X = featureMatrix(featMat, samples, features)
        if(allASD):
            #Remove elements that are PSP (we only want iHART autism + iHART control)
            y = getPhenotypeData(X)
            delete = [index for index in range(y.shape[0]) if y[index]==0 and index>4607]
            y = np.delete(y, delete)
            samples = X.getSamples()
            features = X.getFeatures()
            X_mat = X.getMat()
            
            index = indices(matrix[0:divIndex+1, :], matrix[divIndex+1:,:])
            print (index)
            X_mat = X_mat[:, [j for j in range(X_mat.shape[1]) if j not in index]]
            features = np.delete(features, index)
            samples = np.delete(samples, delete)
        else:
            #Remove elements that are iHART control (we only want iHART autism + PSP control)
            y = getPhenotypeData(X)
            delete = [index for index in range(y.shape[0]) if y[index]==0 and index<=4607]
            y = np.delete(y, delete)
            samples = X.getSamples()
            features = X.getFeatures()
            samples = np.delete(samples, delete)
            X_mat = X.getMat()
            
            index = indices(matrix[0:divIndex+1, :], matrix[divIndex+1:,:])
            print (index)
            X_mat = X_mat[:, [j for j in range(X_mat.shape[1]) if j not in index]]
            features = np.delete(features, index)
        
        print ("Feature Vector Shape:", features.shape)
        print ("Sample Vector Shape:", samples.shape)
        
        X_mat = np.delete(X_mat, delete, axis=0)

        X = featureMatrix(X_mat, samples, list(features))
        print ('Array Dimensions: ', X_mat.shape)
        print ('Size of Feature Matrix: ', np.prod(X_mat.shape))
        print ('Number of Zeroes: ', (np.prod(X_mat.shape)-np.count_nonzero(X_mat)))
        print ('Number of Ones: ', np.count_nonzero(X_mat))
        percentage = ((float)(np.count_nonzero(X_mat))/ (np.prod(X_mat.shape)-np.count_nonzero(X_mat)))
        print ('Percentage of Ones: ', percentage)
        
        return(X, index)
    
# Load in .npz matrices and run the batch effect filter. Returns list of indices of variants that should be filtered out
def indices(matrix1, matrix2):
    ihart = matrix1
    psp = matrix2

    m, n = ihart.shape

    pvalues = np.ones((n,))
    missing_pvalues = np.ones((n,))
    gen_pvalues = np.ones((n,))
    for i in range(n):
        if(i%1000==0): print (i)
        dc = np.asarray([[np.sum(ihart[:, i]<0), np.sum(psp[:, i]<0)],
                         [np.sum(ihart[:, i]==0), np.sum(psp[:, i]==0)],
                         [np.sum(ihart[:, i]==1), np.sum(psp[:, i]==1)],
                         [np.sum(ihart[:, i]==2), np.sum(psp[:, i]==2)],
                        ])
        # get rid of 0 rows/columns
        dc = dc[np.any(dc>0, axis=1), :]

        # calculate pvalue
        if dc.shape[0] > 1 and dc.shape[1] > 1:
            pvalues[i] = chi2_contingency(dc, correction=True)[1]

        dc = np.asarray([[np.sum(ihart[:, i]<0), np.sum(psp[:, i]<0)],
                         [np.sum(ihart[:, i]>=0), np.sum(psp[:, i]>=0)]
                        ])
        # get rid of 0 rows/columns
        dc = dc[np.any(dc>0, axis=1), :]

        # calculate pvalue
        if dc.shape[0] > 1 and dc.shape[1] > 1:
            missing_pvalues[i] = chi2_contingency(dc, correction=True)[1]

        dc = np.asarray([[np.sum(ihart[:, i]==0), np.sum(psp[:, i]==0)],
                         [np.sum(ihart[:, i]==1), np.sum(psp[:, i]==1)],
                         [np.sum(ihart[:, i]==2), np.sum(psp[:, i]==2)],
                        ])
        # get rid of 0 rows/columns
        dc = dc[np.any(dc>0, axis=1), :]

        # calculate pvalue
        if dc.shape[0] > 1 and dc.shape[1] > 1:
            gen_pvalues[i] = chi2_contingency(dc, correction=True)[1]
            
    print(np.sum((-np.log10(pvalues) < -np.log10(0.05/n)) & (-np.log10(missing_pvalues) < -np.log10(0.05/n)) & (-np.log10(gen_pvalues) < -np.log10(0.05/n))))
    boolMat = (((-np.log10(pvalues) < -np.log10(0.05/n)) & (-np.log10(missing_pvalues) < -np.log10(0.05/n)) & (-np.log10(gen_pvalues) < -np.log10(0.05/n))))
    indices = np.where(boolMat == 0)[0]
    return indices

