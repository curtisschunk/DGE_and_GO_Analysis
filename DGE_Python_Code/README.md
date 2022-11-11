# Running .py code

This python code runs differential gene expression analysis following the workflow outlined in the *README.md* file in the *Differential_Gene_Expression_Analysis* directory.

## Inputs when Running Code

**Required** <br>
-c : File for Control Gene Expression <br>
-t : File for Test Gene Expression <br>
-g : File for List of Genes <br>

**Not Required** <br>
-n : Number of Differential Genes to Export, default=200 <br>
-s : Save Path and Save File Name, default='Differential_Genes' <br> <br>
*the default will export a .csv file titled 'Differential_Genes' containing the top 200 differential genes to the current directory*

**Example Running Code** <br>
```bash
python DGE.py -c EXAMPLE_control_gene_expression.csv -t EXAMPLE_test_gene_expression.csv -g EXAMPLE_gene_list.csv
```

## Data Structure

1. Take preprocessed gene expression data and break it into two seperate datasets: Expression data of the cluster to be analysed (test gene expression data) AND Expression data of all other clusters (control gene expression data)
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
