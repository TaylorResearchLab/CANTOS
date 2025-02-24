# CT-Embedding-Paper
The results of this study is used to standardize the tumor names in CT database, so they can be integrated with other biomedical databases for further downstream analysis and understanding the therapeutic agents and drug-target landscape for a given tumor.   </br>


## Embeddings Data Download 
Following are the steps for running the pipeline: </br>
1. Clone this Github repository to your local machine </br> 
2. Navigate to the following Open Science Foundation website to find the embeddings data files needed to run CANTOS: https://doi.org/10.17605/OSF.IO/DBGWN </br>
3. Under the files section, locate the file titled "Embeddings.zip" and click it. </br>
4. On the next page click the menu with three vertical dots, and select the download button, a zip file will be downloaded. </br>
5. Unzip the file and store the embeddings files in the data directory of the cloned GitHub repository: </br>

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
Before running CANTOS, please ensure <br/>
1. Please ensure your machine has R installed on it. It can be downloaded from the following website: https://www.r-project.org/ <br/>
2. The libraries listed in the Library section below are installed. <br/>
3. In the following scripts make sure to edit the makeCluster argument with the number of cores available on your machine. The number of available cores cab found using the R command `detectCores()` in the `parallel` library <br/>
  **02-calculate-edit-distance-5thed.R** <br/>
  **02-calculate-edit-distance.R** <br/>
  **07A-annotate-cluster-result-NCIT-WHO-5thed.R**<br/>
  **07A-annotate-cluster-result-NCIT-WHO.R**<br/>
  **07B-annotate-cluster-result-V3-NCIT-WHO-5thed.R** <br/>
  **07B-annotate-cluster-result-V3-NCIT-WHO.R**<br/>
  
We ran CANTOS on RStudio Version 2023.09.1+494 (2023.09.1+494) using R version 4.4.0 (2024-04-24). Users can also run CANTOS from the command line from the following directory
`CANTOS/analysis` using the following command: <br/>

`bash CANTOS.sh`


## Description
This repository contains the code, tables, and plots associated with the CT Embedding paper. The pipeline built in this repository does the following task: <br/>

1. Extract tumor names from the CT database if they are associated with an NCT ID and have an associated drug belonging to the categories  of Drug, Biological,Combination Product,Genetic. A total 50,410 condition names are extracted.<br/>

2. These 50410 condition names are flagged as tumors and non tumors by the pipeline, which are then further manually annotated pediatric and adult tumors. A total of 13,230 tumors are identified from the 50,410 conditions and out of the 13,230 tumors,  6,324 were classified as pediatric tumors. <br/> 

3. Compute the distance of each  13,230 clinical trials tumors, 4720 WHO tumors, and 1395 NCIT tumors.  Distance metrics used are Levenshtein, Cosine , and Jarro-Winkler. </br>

4. Find the closest matching WHO term for each tumor for each distance metric and then also group the top 0.05% closest matching group of tumors. Each of the closest match terms are standardized to their closest matching WHO Term. </br>

5. Use the distance matrices computed to perform 3 levels of nested affinity clustering and group the tumors. After grouping the tumors they are standardized to their closest matching WHO Term. </br>

6. Generate embeddings for each tumor terms (CT, WHO, NCIT) using Open AI's models  text-embedding-3-large (LTE-3) 
 and text-embedding-ada-002 (ADA-002). We then identify the closest matching (Euclidean Distance) WHO terms for each tumor name from the CTR.  </br>

7. Perform PCA on each of the embedding types and then run K-means and Affinity Clustering to group the tumors together. We refine the clusters by filtering outliers using isolation forest and local outlier factor.  </br>

8. After cluster refinement, each cluster is standardized to the WHO term that matches a majority of the members of that cluster.  

9. We generate embeddings for each tumor terms (CT, WHO, NCIT) using embeddings obtained from LLama 3.3, Llama 3.2,LLama 3.0  
,MedLLama2,MedLLama13B,BioBERT,PubMedBERT-Abstract (PubMedBERT),ModernBERT-Large (ModernBERT) ,Phi-4,,e5-large,e5-large-v2, all-mpnet-base-v2,gtr-t5-large,all-MiniLM-L12-v2,all-MiniLM-L6-v2,all-roberta-large-v1 ,SapBERT ,ClinicalBERT,LaBSE,BioGPT,
DeepSeek_8B,SciBERT,nomic-embed-text, and Cohere: embed-english-v2.0.

10. These embeddings are provided as inputs to CANTOS which then identifies the closest matching (Euclidean Distance) WHO terms for each tumor name from the CTR using these embeddings .  </br>
   
11. We randomly sampled 1600 tumor names from the CTR and manually annotated their ground truths obtained from the WHO System. We observed the methods LTE-3+Euclidean Dist and all-MiniLM-L12-v2+Euclidean Dist had the highest standardization accuracy against WHO all editions and WHO 5th edition respectively. We filtered any CTR tumor names (from the 1600 sampled) that did not have a ground truth to evaluate the standardization accuracy.  We then plotted the distributions of the Euclidean distances of the remaining CTR terms to their respective WHO terms as identified by these two methods and segregated the distribution based on correct and incorrect standardization.</br>



## Scripts

**00-generate-ct-disease-file.R**:  </br> 
This script loads data from clinical trials and select only the diseases with NCT ID , and associated with Intervention types of Drug, Biological,Combination Product,Genetic. Totally 50410 diseases are extracted. </br> 


**01-generate-disease-annotation-for-manual-review.R** </br>
This script annotates the 50K diseases automatically as cancer or not. </br>

**02-calculate-edit-distance.R** </br>
This script loads the manually annotated disease file with pediatric and adult cancer annotation and computes the edit distance matrices. WHO  all editions was used in this script <br/>
**02-calculate-edit-distance-5thed.R**</br>
This script loads the manually annotated disease file with pediatric and adult cancer annotation and computes the edit distance matrices. WHO  5th editions was used in this script. <br/>

**03-edit-distance-clustering.R** </br>
This Script performs affinity propagation clustering using edit distances. WHO all editions was used in this script.  </br>
**03-edit-distance-clustering-5thed.R** </br>
This Script performs affinity propagation clustering using edit distances. WHO 5th editions was used in this script. </br>

**04A-preprocess-embedding-pca.R** </br>
These script loads AD-A002 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO all editions was used in this script. </br>
**04A-preprocess-embedding-pca-ADA2-5thed.R**</br>
These script loads ADA-002 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO 5th editions was used in this script. </br>


**04B-preprocess-embedding-pca-v3.R** </br>
These script loads LTE-3 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database all editions was used in this script. </br>
**04B-preprocess-embedding-pca-v3-5thed.R** </br>
These script loads LTE-3 embeddings for CT, WHO, NCIT  Tumors and then performs PCA.WHO database 5th editions was used in this script. </br>

**05A-cluster-on-ADA2-embedding-Kmeans.R**</br>
This script computes Kmeans cluster using ADA-002 embeddings  and also computes silhouette index.WHO database all editions was used in this script. </br>
**05A-ADA2-embedding-Kmeans-5thed.R** </br>
This script computes Kmeans cluster using ADA-002 embeddings  and also computes silhouette index.WHO database 5th editions was used in this script. </br>

**05B-v3-embedding-Kmeans.R** </br>
This script computes Kmeans cluster using LTE-3 embeddings  and also computes silhouette index.WHO database all editions was used in this script. </br>
**05B-v3-embedding-Kmeans-5thed.R** </br>
This script computes Kmeans cluster using LTE-3 embeddings  and also computes silhouette index.WHO database 5th editions was used in this script. </br>

**06A-cluster-on-ADA-embedding-affinity.R** </br>
This script computes affinity propagation clustering using ADA-002 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database all editions was used in this script. </br>
**06A-cluster-on-ADA-embedding-affinity-5thed.R** </br>
This script computes affinity propagation clustering using ADA-002 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database 5th editions was used in this script.</br>

**06B-cluster-on-V3-embedding-affinity.R** </br>
This script computes affinity propagation cluster using LTE-3 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database all editions was used in this script.</br>
**06B-cluster-on-V3-embedding-affinity-5thed.R** </br>
This script computes affinity propagation cluster using LTE-3 embeddings. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.WHO database 5th editions was used in this script. </br>


**07A-annotate-cluster-result-NCIT-WHO.R**</br>
This script annotates Affinity propagation cluster results of ADA-002 embeddings. WHO database all editions was used in this script. </br>
**07A-annotate-cluster-result-NCIT-WHO-5thed.R**</br>
This script annotates Affinity propagation cluster results of ADA-002 embeddings. WHO database 5th editions was used in this script. </br>

**07B-annotate-cluster-result-V3-NCIT-WHO.R** </br>
This script annotates Affinity propagation cluster results of LTE-3 embeddings. WHO database all editions was used in this script. </br>
**07B-annotate-cluster-result-V3-NCIT-WHO-5thed.R** </br>
This script annotates Affinity propagation cluster results of LTE-3 embeddings. WHO database 5th editions was used in this script. </br>
**08-outlier-detection-embeddings.R** </br>
This script is used to detect if Affinity propagation cluster members are outliers using LOF and Isolation Forest.We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database all editions was used in this script </br>
**08-outlier-detection-embeddings-5thed.R** </br>
This script is used to detect if Affinity propagation cluster members are outliers using LOF and Isolation Forest.We perform this for clusters formed using both ADA-002 and LTE-3 embeddings. WHO database 5th editions was used in this script </br>

**09-cluster-reassignment-outlier.R** </br>
This script performs reannotates Affinity cluster  after outlier detection. We perform this for clusters formed using both ADA002 and V3 embeddings. WHO database all editions was used in this script. </br>
**09-cluster-reassignment-outlier-5thed.R** </br>
This script performs reannotates Affinity cluster  after outlier detection. We perform this for clusters formed using both ADA-002 and LTE-3 embeddings. WHO database 5th editions was used in this script. </br>

**10-assign-who-ncit-outlier-kmeans-editdistance-clustering.R**</br>
This script to detect outliers for embedding-based-Kmeans and edit distance based standardization. WHO database all editions was used in this script. </br>
**10-assign-who-ncit-outlier-kmeans-editdistance-clustering-5thed.R** </br>
This script to detect outliers for embedding-based-Kmeans and edit distance based standardization. WHO database 5th editions was used in this script. </br>

**11-os-embedding-euclidean-dist.R**</br>
This script is used to compute Euclidean distance matrices between tumor names in CTR, WHO and, NCIt  using embeddings obtained from non-Open AI models. For each embedding type, we identify the closest match WHO 5th edition, WHO all edition, and NCIt terms for ever CTR tumor name. </br>

**12-sample-CT-tumors-validation.R**</br>
Randomly Sample 1600 tumors from the CTR Tumor Names. </br>

**13-summarize-results.R** </br>
For each standardization method compute their standardization accuracy against the sampled tumor names for which at least one ground truth was found from the WHO system . </br>

**14-compute-euclidean-v3-threshold.R** </br>
Script is used for generating box-plots for comparing the standardization results of the methods LTE-3+Euclidean Dist and all-MiniLM-L12-v2+Euclidean Dist. </br>

**15-compute-euclidean-v3-threshold.R** </br>
Script is used to generate average silhouette score vs number of cluster plot in figure 5.</br>

**16-generate-heatmaps-mi-analysis.R** </br>
Script is used to perform mutual information analysis to identify three high accuracy standardization methods for performing majority vote to select standardization for a given CTR tumor name. The script also generates heatmap (figure 6 and 7) of pairwise mutual information among the high-accuracy methods.</br>

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
24. pdist </br>
