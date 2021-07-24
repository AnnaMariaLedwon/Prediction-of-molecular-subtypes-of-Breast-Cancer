# Prediction of molecular subtypes of Breast Cancer

### Data information
This data set contains published iTRAQ proteome profiling of 77 breast cancer samples generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). It contains expression values for ~12.000 proteins for each sample, with missing values present when a given protein could not be quantified in a given sample.

![image](https://user-images.githubusercontent.com/45319617/126872780-8c5bfe4b-8ea2-47c2-b372-18af1a3b1ed7.png)

The data set consists of:
* **77cancerproteomesCPTACitraq.csv** : The decription of quantitative mass-spectrometry-based proteomic and phosphoproteomic analyses of 105 genomically annotated breast cancers, of which 77 provided high-quality data. More information here: https://www.nature.com/articles/nature18003
* **clinicaldatabreast_cancer.csv** :  The data about the cancer classification of a given sample using different methods. For me it was important "mRNA.PAM50" method, which I used as a training and testing data set in the Machine Learning models.
* **PAM50_proteins.csv** : The list of genes and proteins used by the PAM50 classification system.

Source of dataset: https://www.kaggle.com/piotrgrabo/breastcancerproteomes 

### About project ...
Using lasso regression, 50 proteins were selected to train and test machine learning models to classify and predict the molecular subtype of breast cancer. Three models were built, namely the support vector machine, random forest and artificial neural networks. In addition, the prediction using PAM50 proteins was checked on the same classification models. It turned out that the models based on protein generated by lasso regressions performed better during the evaluation.

All *Prediction of molecular subtypes of Breast Cancer* code you can find in **BC_analysis.R** file!

#### 1.   Data modification:
* *BC_missingData.pdf* : plot with information about missing data;
* *BC_subtypes_barplot.pdf* : barplot showing the amount of concrete cancer subtype in our set;

#### 2.   PCA and T-SNE Analysis:
* *BC_subtypes_PCA.pdf* : PCA results;
* *BC_subtypes_TSNE.pdf* : T-SNE results;

#### 3.   Feature selection:
* *BC_LassoProteins_lassoRegression.pdf* : first 50 proteins was obtaind by lasso regression - showing times for each protein;
* *Lasso_proteins_boxplot.pdf* : visualisation of the lasso proteins using boxplot;
* *Lasso_proteins_volinplot.pdf* : visualisation of the lasso proteins using violin plot;
* *Lasso_proteins.csv* : table of protein ID, gene ID, description and times from lasso regression were created using biomaRt package;

#### 4. Models buildings for Lasso proteins:
* *Lasso_train&testData.pdf* : information about the number of patients of breast cancer subtype in the training and testing set
* *Lasso_SVM_results.txt* : confusion matrix for Support Vector Machine Model
* *Lasso_RandomForest_ErrorsTrees.pdf* plot for Random Forest Model
* *Lasso_RandomForest_results.txt* : confusion matrix for Random Forest Model
* *Lasso_NeuralNetworkClassifier.pdf* : visualization of Neural Network
* *Lasso_NeuralNetworkClassifier_results.txt* :confusion matrix for Neural Network Classifier

#### 5. Models buildings for PAM50 proteins:
* *PAM50_train&testData.pdf* : information about the number of patients of breast cancer subtype in the training and testing set
* *PAM50_proteins_boxplot.pdf* : visualisation of the PAM50 proteins using boxplot;
* *PAM50_proteins_volinplot.pdf* : visualisation of the PAM50 proteins using violin plot;
* *PAM50_SVM_results.txt* : confusion matrix for Support Vector Machine Model
* *PAM50_RandomForest_ErrorsTrees.pdf* plot for Random Forest Model
* *PAM50_RandomForest_results.txt* : confusion matrix for Random Forest Model
* *PAM50_NeuralNetworkClassifier.pdf* : visualization of Neural Network
* *PAM50_NeuralNetworkClassifier_results.txt* :confusion matrix for Neural Network Classifier

### The interesting scientific articles about Breast Cancer...
* *Molecular Subtypes of Breast Cancer:* https://www.breastcancer.org/symptoms/types/molecular-subtypes
* *Proteomic maps of breast cancer subtypes:* https://www.nature.com/articles/ncomms10259
* *Proteogenomics connects somatic mutations to signalling in breast cancer:* https://www.nature.com/articles/nature18003


