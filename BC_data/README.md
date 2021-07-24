### Input data information

#### File: 77cancerproteomesCPTACitraq.csv
* RefSeq protein ID (each protein has a unique ID in a RefSeq database)
* gene_symbol: a symbol unique to each gene (every protein is encoded by some gene)
* gene_name: a full name of that gene, remaining columns: log2 iTRAQ ratios for each sample (protein expression data, most important), three last columns are from healthy individuals

#### File: clinicaldatabreast_cancer.csv
First column "Complete TCGA ID" is used to match the sample IDs in the main cancer proteomes file. All other columns have self-explanatory names, contain data about the cancer classification of a given sample using different methods. 'PAM50 mRNA' classification is being used in the example script.

#### File: PAM50_proteins.csv
Contains the list of genes and proteins used by the PAM50 classification system. The column RefSeqProteinID contains the protein IDs that can be matched with the IDs in the main protein expression data set.

Source: https://www.kaggle.com/piotrgrabo/breastcancerproteomes
