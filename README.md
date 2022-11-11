#### by Curtis Schunk

# Differential Gene Expression and Gene Ontology (GO) Analysis (Python)

## Step 1 - Data Preprocessing

For my Data Preprocessing, I ran the gene expression data through CCST's data preprocessing code, data_generation_merfish.py (https://github.com/xiaoyeye/CCST/blob/main/data_generation_merfish.py)

This data preprocessing code acts to normalize data and correct batch effects in a dataset. The batch effects are corrected using the function *scanorama.correct* from the python tool Scanorama (https://github.com/brianhie/scanorama). This tool was published in the 2019 *Nature Biotechnology* article titled *Efficient integration of heterogeneous single-cell transcriptomes using Scanorama* (https://www.nature.com/articles/s41587-019-0113-3).

This is info on how to run the preprocessing code from the CCST Github. In short, you need a file of gene expression data (columns are the cells, rows are the genes) and a file of the spatial location of the cells (column 1 is the x-coords, column 2 is the y-coords, rows are the cells)

If using the MERFISH dataset from the CCST GitHub, the preprocessed dataset is data which has gone through this step.

**Steps of CCST Data Preprocessing**
1. Batch Correction (if needed)
2. Remove Genes with too low average values
    - Default threshold = 1
3. Normalization
    - Divide each gene value by the average of all values for that gene
4. Remove Low Variance Genes
    - Default threshold = 0.4

***How Scanorama batch correction is run in CCST Preprocessing:***

```python
 # Batch correction
    import scanorama
    
    corrected, _ = scanorama.correct(datasets, genes_lists)

    features_corrected = []
    for i, corrected_each_batch in enumerate(corrected):
        if i == 0:
            features_corrected = corrected_each_batch.A
        else:
            features_corrected = np.vstack((features_corrected, corrected_each_batch.A))
    features_corrected = np.array(features_corrected)
    np.save( generated_data_path + 'features.npy', features_corrected)
    print('corrected size: ', features_corrected.shape)
    return features_corrected
```

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

## Step 3 - Differential Gene Expression Analysis (DGE.py in DGE_Python_Code directory)

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
    
    dif_genes=np.array(top200_genes)
    fixed_label_dif_genes=list()
    for i in range(len(dif_genes)):
        current_gene_working_with=str(dif_genes[i])
        current_gene_working_with=current_gene_working_with.replace("['","")
        current_gene_working_with=current_gene_working_with.replace("']","")
        fixed_label_dif_genes.append(current_gene_working_with)
    
    return  fixed_label_dif_genes
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
4) (Optional) Enter in the Background Genes (more information below)
5) Paste Differentially Expressed Genes into GO Tool and run (should not take long)
6) Download the results and sort GO Terms by highest significance (lowest p-values)

### Background Genes
The background genes is a list of all of the genes that were used when performing differential gene expression anaysis.
The default for GO analysis is to use all known genes in your selected genome as background.
However, if you differential gene expression analysis was performed using less genes, it will increase the accuracy of GO anaysis to enter your list of background genes.

#### Steps
1) Get your list of background genes (all genes used in differential gene expression analysis)
2) Enter the list of background genes into the correct section of the tool. This is going to be in a different spot depending on what tool you are using.

##### Tool: http://geneontology.org/
a) Enter your differentially expressed genes and run the GO tool <br>
b) When you get the the *Panther Classification System* site, click the 'change' box by 'Reference List' <br>
c) Upload a file containing your background genes

##### Tool: https://biit.cs.ut.ee/gprofiler/gost
a) Prior to running the GO tool, open the 'Advanced Options' section <br>
b) Select 'Custom' for 'Statistical domain scope' and then enter your list of background genes <br>

##### Tool: https://go.princeton.edu/cgi-bin/GOTermFinder
a) Prior to running the GO tool, go to the 'Optional Advanced Input Options' section <br>
b) Click the 'Choose File' button next to the text that says 'provide a list of genes for the background population' <br>
c) Upload a file containg your background genes

# Dataset used for Differential Gene Expression Anaysis

To run Differential Gene Expression Analysis, I used the MERFISH dataset published on the CCST GitHub. Paper citation and link:

<br>

**CCST Paper:** <br>

**Citation:** Li, J., Chen, S., Pan, X., Yuan, Y., & Shen, H. B. (2022). Cell clustering for spatial transcriptomics data with graph neural networks. Nature Computational Science, 2(6), 399-408. <br>
**Link to Paper:** https://www.nature.com/articles/s43588-022-00266-5 <br>
**Link to GitHub:** https://github.com/xiaoyeye/CCST 

**Paper that MERFISH dataset is from:** <br>

**Citation:** Xia C., Fan J., Emanuel G., Hao J., & Zhuang X. (2019). Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. Proceedings of the National Academy of Sciences 116(39):19490-19499. <br>
**Link to Paper:** https://www.pnas.org/doi/10.1073/pnas.1912459116

<br>

### Downloading Gene Expression and Spatial Location Data

***note:*** the MERFISH dataset published on CCST's GitHub contains data that was measured in three seperate batches. The **raw** gene expression data does not provide batch correction for this data. The different batches are denoted in the name of each cell. 'B1' refers to cells in batch 1, 'B2' to cells in batch 2, 'B3' to cells in batch 3. The **preprocessed** dataset does provide some batch correction in the form of a python tool, Scanorama (https://github.com/brianhie/scanorama). This tool was published in the 2019 *Nature Biotechnology* article titled, "*Efficient integration of heterogeneous single-cell transcriptomes using Scanorama*" (https://www.nature.com/articles/s41587-019-0113-3).

1) Go to https://github.com/xiaoyeye/CCST 

##### Downloading the Spatial Location Data

2) Go to the ‘dataset’ folder
3) download the file titled ‘merfish.zip’
4) Unzip the file
5) The spatial location data is contained in the file titled **‘pnas.1912459116.sd15.xlsx’**

##### Downloading the Gene Expression Data

6) Decide if you would like to download the **Raw Dataset** or the **Dataset Preprocessed using CCST**
    1) **Raw Dataset:**
        1) Go to the ‘dataset’ folder
        2) download the file titled ‘merfish.zip’
        3) Unzip the file
        4) Gene Expression data is contained in the file titled ‘pnas.1912459116.sd12.csv’
    2) **Preprocessed Dataset:**
        1) Go to the ‘generated_data’ folder
        2) Go to the ‘MERFISH’ folder
        3) Download the files titled ‘features.npy’ and ‘gene_names.txt’
            - **‘features.npy’** - the gene expression data
            - **‘gene_names.txt’** - the gene names for the data in ‘features.npy’
