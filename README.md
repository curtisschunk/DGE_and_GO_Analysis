# Differential Gene Expression Analysis and GO Analysis

#### by Curtis Schunk

## Step 1 - Data Preprocessing

For my Data Preprocessing, I ran the gene expression data through CCST's data preprocessing code, data_generation_merfish.py (https://github.com/xiaoyeye/CCST/blob/main/data_generation_merfish.py)

There is info how to run the preprocessing on the CCST Github. In short, you need a file of gene expression data (columns are the cells, rows are the genes) and a file of the spatial location of the cells (column 1 is the x-coords, column 2 is the y-coords, rows are the cells)

**Steps of CCST Data Preprocessing**
1. Remove Genes with too low average values
    - Default threshold = 1
2. Normalization
    - Divide each gene value by the average of all values for that gene
3. Remove Low Variance Genes
    - Default threshold = 0.4

## Step 2 - Data Structure

1. Take the preprocessed gene expression data and break it into two seperate datasets: Expression data of the cluster to be analysed AND Expression data of all other clusters
    - *example*: if you have 3 clusters and want to find the differential gene expression of cluster 1, you would make a matrix of Gene Expression Data containing only cluster 1 and a matrix of Gene Expression Data containing clusters 2 and 3

2. Make a list of the genes, in the same order that they appear in the datasets

3. After this, you should have 2 datasets of gene expression data (one dataset that you want to test, another to use as a control) and an array of the names of genes

**Gene Expresssion Data Format**
- in a numpy array, with shape = (#OfGenes, #OfCells)
- *example*: an array of 1889 genes and 301 cells 
    <img width="1075" alt="Screen Shot 2022-07-20 at 1 43 43 PM" src="https://user-images.githubusercontent.com/108155131/180058768-ef9922ea-f3ee-4b02-8728-7139b80a53b9.png">

**List of Genes Data Format**
- in a numpy array
- *example*: an array of 1889 genes
    <img width="1069" alt="Screen Shot 2022-07-20 at 1 42 46 PM" src="https://user-images.githubusercontent.com/108155131/180058604-a0b3a427-2f6d-47c3-9325-1b42333a59a9.png">

## Step 3 - Differential Gene Expression Analysis

#### 1. Packages

```python
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy import stats
```

#### 2. Functions

- Mann-Whitney U Test (*test*)
  - Takes in two data arrays, conducts a Mann-Whitney U Test, and outputs the p-value between the two arrays

```python
def test(data1, data2):
    # mannwhitneyu
    stat, p = stats.mannwhitneyu(data1, data2)
    return p
```

- Differential Gene Analysis for Individual Genes (*p_val_finder*)
  - Takes in the gene expression array for the test group, the gene expression array for the control group, and the index number of the gene to test. This function uses the previous function to find the p-value of the selected gene between the two groups. This function also finds whether the gene is overexpressed (cur=1) or underexpressed (cur=0) in the control group in comparison to the test group.

```python
#for specified gene, finds p-val and whether it is over or under expressed

def p_val_finder(test_array,ctrl_array,gene_test_num):
    
    test_vals=test_array[gene_test_num,:]
    ctrl_vals=ctrl_array[gene_test_num,:]

    p_val=test(test_vals,ctrl_vals)
    log_p_val=-np.log10(p_val)

    if np.mean(test_vals)>np.mean(ctrl_vals): #overexpressed gene in test group
        cur=True
    else: #underexpressed gene in test group
        cur=False

    return p_val,log_p_val,cur
    
# cur -> if = 1, expression is higher in the group being compared to the other groups
#        if = 0, expression is lower in the group being comared to the other groups
```

- Differential Gene Analysis for All Genes (*dif_gene_analysis*)
  - Takes in the gene expression array for the test group, the gene expression array for the control group, and the array of gene names. This function uses the previous function to find the p-value and expression level (over/under) of each gene in the test group in comparison to the control group. The function then combines the results into a DataFrame, and sorts the DataFrame based off of p-values and wheter the gene is overexpressed (only overexpressed genes are kept for use in GO Analysis). The function outputs a list of the top 200 Differentially Expressed Genes that are overexpressed in the test group.

```python
def dif_gene_analysis(test_array,ctrl_array,gene_list): 
    gene_len=len(gene_list)
    dif_gene_array=np.zeros([gene_len,3])
    dif_gene_list=list()
    for i in range(gene_len):

        dif_gene_list.append(gene_list[i])
        dif_gene_array[i,0],dif_gene_array[i,1],dif_gene_array[i,2]=p_val_finder(test_array,ctrl_array,i)
        
    tot_df=pd.DataFrame([dif_gene_list,dif_gene_array[:,0],dif_gene_array[:,1],dif_gene_array[:,2]]).transpose()
    tot_df=tot_df.set_axis(['Gene','p-val','-log10(p-val)','cur'], axis=1, inplace=False)
    
    sorted_df=tot_df.sort_values(by=['cur', '-log10(p-val)'],ascending=False)
    
    top200_genes=sorted_df['Gene'][0:200]
    
    #return sorted_df
    return top200_genes
```

#### 3. Running

```python
dif_genes=dif_gene_analysis(test_array,ctrl_array,gene_list)
```

#### 4. Export

```python
dif_genes.to_csv('200_Dif_Genes.csv')
```

## Step 4 - GO Analysis
1) Copy the 200 Differentially Expressed Genes
2) Go to one of the GO Analysis Tools
    - http://geneontology.org/
    - https://biit.cs.ut.ee/gprofiler/gost
    - https://go.princeton.edu/cgi-bin/GOTermFinder
3) In the GO Tool make sure to select the correct Gene Annotation (Homo sapiens) and the correct Gene Ontology Terms (for identifying cell cycles, I used the 'GO Biological Process' terms)
4) Paste Differentially Expressed Genes into GO Tool and run (should not take long)
5) Download the results and sort GO Terms by highest significance (lowest p-values)
