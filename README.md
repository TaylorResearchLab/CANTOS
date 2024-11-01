# CT-Embedding-Paper
The results of this study is used to standardize the tumor names in CT database, so they can be integrated with other biomedical databases for further downstream analysis and understanding the therapeutic agents and drug-target landscape for a given tumor.   </br>


## Embeddings Data Download 
Following are the steps for running the pipeline: </br>
1. Clone this Github repository to your local machine </br> 
2. Navigate to the following DOI:10.17605/OSF.IO/DBGWN and download the zip file labeled Embeddings.zip </br>
3. Unzip the file and store the OPEN AI embeddings files in the data directory of the cloned GitHub repository: </br>

| File Name             | Directory    | 
| :---------------------|:-------------| 
| CT_Embeddings_ADA2.csv| `CANTOS/data`  | 
| CT_Embeddings_V3.csv	| `CANTOS/data`  | 
| NCIT_Embeddings_V3.csv| `CANTOS/data`  | 
| WHO_Aggregate_ADA2.csv| `CANTOS/data`  | 
| WHO_Terms_All_V3.csv	| `CANTOS/data`  | 

</br> 
Please note that the ADA002 embeddings file for NCIT is contained in the following directory: </br>

`CANTOS/data/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv`</br>

## Run Instructions for CANTOS
We ran CANTOS on RStudio Version 2023.09.1+494 (2023.09.1+494) using R version 4.4.0 (2024-04-24). Users can also run CANTOS from the command line from the following directory
`CANTOS/analysis` using the following command: <br/>

`bash CANTOS.sh`

Please ensure your machine has R installed on it. It can be downloaded from the following website: https://www.r-project.org/ <br/>

Before running CANTOS, please ensure the libraries listed in the Library section below. 


## Description
This repository contains the code, tables, and plots associated with the CT Embedding paper. The pipeline built in this repository does the following task: <br/>

1. Extract tumor names from the CT database if they are associated with an NCT ID and have an associated drug belonging to the categories  of Drug, Biological,Combination Product,Genetic. A total 50,410 condition names are extracted.<br/>

2. These 50410 condition names are flagged as tumors and non tumors by the pipeline, which are then further manually annotated pediatric and adult tumors. A total of 13,230 tumors are identified from the 50,410 conditions and out of the 13,230 tumors,  6,324 were classified as pediatric tumors. <br/> 

3. We compute the distance of each  13,230 clinical trials tumors, 4720 WHO tumors, and 1395 NCIT tumors.  Distance metrics used are Levenshtein, Cosine , and Jarro-Winkler. </br>

4. We find the closest matching WHO term for each tumor for each distance metric and then also group the top 0.05% closest matching group of tumors. Each of the closest match terms are standardized to their closest matching WHO Term. </br>

5. We also use the distance matrices computed to perform 3 levels of nested affinity clustering and group the tumors. After grouping the tumors they are standardized to their closest matching WHO Term. </br>

6. We generate embeddings for each tumor terms (CT, WHO, NCIT) using Open AI's ADA2.0 and V-3 Large text-embedding models. We then compute the closest matching (Euclidean Distance) WHO terms for each tumor.  </br>

7. We also perform PCA on each of the embedding types and then run K-means and Affinity Clustering to group the tumors together. We refine the clusters by filtering outliers using isolation forest and local outlier factor.  </br>

8. After cluster refinement, each cluster is standardized to the WHO term that matches a majority of the members of that cluster.  
   




## Scripts

**00-generate-ct-disease-file.R**:  </br> 
This script loads data from clinical trials and select only the diseases with NCT ID , and associated with Intervention types of Drug, Biological,Combination Product,Genetic. Totally 50410 diseases are extracted. </br> 


**01-generate-disease-annotation-for-manual-review.R** </br>
This script annotates the 50K diseases automatically as cancer or not. </br>

**02-calculate-edit-distance.R** </br>
This script loads the manually annotated disease file with pediatric and adult cancer annotation and computes the edit distance matrices. WHO database all editions was used in this script <br/>
**02-calculate-edit-distance-5thed.R**</br>
This script loads the manually annotated disease file with pediatric and adult cancer annotation and computes the edit distance matrices. WHO database 5th editions was used in this script. <br/>

**03-edit-distance-clustering.R** </br>
This Script performs affinity propagation clustering using edit distances. WHO database all editions was used in this script.  </br>
**03-edit-distance-clustering-5thed.R** </br>
This Script performs affinity propagation clustering using edit distances. WHO database 5th editions was used in this script. </br>

**04A-preprocess-embedding-pca.R** </br>
These script loads ADA002 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database all editions was used in this script. </br>
**04A-preprocess-embedding-pca-ADA2-5thed.R**</br>
These script loads ADA002 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database 5th editions was used in this script. </br>


**04B-preprocess-embedding-pca-v3.R** </br>
These script loads V3 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database all editions was used in this script. </br>
**04B-preprocess-embedding-pca-v3-5thed.R** </br>
These script loads V3 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database 5th editions was used in this script. </br>

**05A-cluster-on-ADA2-embedding-Kmeans.R**</br>
This script computes Kmeans cluster using ADA002 embeddings  and also computes silhouette index.WHO database all editions was used in this script. </br>
**05A-ADA2-embedding-Kmeans-5thed.R** </br>
This script computes Kmeans cluster using ADA002 embeddings  and also computes silhouette index.WHO database 5th editions was used in this script. </br>

**05B-v3-embedding-Kmeans.R** </br>
This script computes Kmeans cluster using V3 embeddings  and also computes silhouette index.WHO database all editions was used in this script. </br>
**05B-v3-embedding-Kmeans-5thed.R** </br>
This script computes Kmeans cluster using V3 embeddings  and also computes silhouette index.WHO database 5th editions was used in this script. </br>

**06A-cluster-on-ADA-embedding-affinity.R** </br>
This script computes affinity propagation clustering using ADA002 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database all editions was used in this script. </br>
**06A-cluster-on-ADA-embedding-affinity-5thed.R** </br>
This script computes affinity propagation clustering using ADA002 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database 5th editions was used in this script.</br>

**06B-cluster-on-V3-embedding-affinity.R** </br>
This script computes affinity propagation cluster using V3 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database all editions was used in this script.</br>
**06B-cluster-on-V3-embedding-affinity-5thed.R** </br>
This script computes affinity propagation cluster using V3 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database 5th editions was used in this script. </br>


**07A-annotate-cluster-result-NCIT-WHO.R**</br>
This script annotates Affinity propagation cluster results of ADA002 embeddings. WHO database all editions was used in this script. </br>
**07A-annotate-cluster-result-NCIT-WHO-5thed.R**</br>
This script annotates Affinity propagation cluster results of ADA002 embeddings. WHO database 5th editions was used in this script. </br>

**07B-annotate-cluster-result-V3-NCIT-WHO.R** </br>
This script annotates Affinity propagation cluster results of V3 embeddings. WHO database all editions was used in this script. </br>
**07B-annotate-cluster-result-V3-NCIT-WHO-5thed.R** </br>
This script annotates Affinity propagation cluster results of V3 embeddings. WHO database 5th editions was used in this script. </br>
**08-outlier-detection-embeddings.R** </br>
This script is used to detect if Affinity propagation cluster members are outliers using LOF and Isolation Forest.We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database all editions was used in this script </br>
**08-outlier-detection-embeddings-5thed.R** </br>
This script is used to detect if Affinity propagation cluster members are outliers using LOF and Isolation Forest.We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database 5th editions was used in this script </br>

**09-cluster-reassignment-outlier.R** </br>
This script performs reannotates Affinity cluster  after outlier detection. We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database all editions was used in this script. </br>
**09-cluster-reassignment-outlier-5thed.R** </br>
This script performs reannotates Affinity cluster  after outlier detection. We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database 5th editions was used in this script. </br>

**10-assign-who-ncit-outlier-kmeans-editdistance-clustering.R **</br>
This script to detect outliers for embedding-based-Kmeans and editdistance based standardization. WHO database all editions was used in this script. </br>
**10-assign-who-ncit-outlier-kmeans-editdistance-clustering-5thed.R** </br>
This script to detect outliers for embedding-based-Kmeans and editdistance based standardization. WHO database 5th editions was used in this script. </br>

**11-generate-records-annotation.R**</br>
This script is used to annotate the types of ground truth found for each of the 1600 tumors sampled. </br>

**12-summarize-results.R**</br>
Prints the accuracy of each standardization method. </br>

**13-plot-script.R** </br>
Script is used for generating silhouette plots. </br>

## Libraries 
1. apcluster </br>
2. biomaRt </br>
3. cluster </br>
4. data.table </br>
5. dbscan </br>
6. DescTools </br>
7. doParallel </br>
8. dplyr </br>
9. factoextra</br>
10. foreach</br>
11. ggplot2 </br>
12. ggpubr </br>
13. ghql </br>
14. httr </br>
15. isotree </br>
16. jsonlite </br>
17. magrittr </br>
18. qdapRegex </br>
19. readxl </br>
20. stringdist </br>
21. stringi </br>
22. stringr </br>
23. tidyverse </br>
