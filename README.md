# CT-Embedding-Paper

## Description
This repository contains the code, tables, and plots associated with the CT Embedding paper. The pipeline built in this repository does the following task: <br/>

1. Extract tumor names from the CT database if they are associated with an NCT ID and have an associated drug belonging to the categories  of Drug, Biological,Combination Product,Genetic. A total 50,410 diseases are extracted.<br/>

2. These 50410 diseases are flagged as tumors and non tumors by the pipeline, which are then further manually annotated pediatric and adult tumors. A total of 13,329 tumors are identified from the 50,410 diseases.<br/> 

3. We compute the distance of each  13,329 clinical trials tumors, 4720 WHO tumors, and 1395 NCIT tumors.  Distance metrics used are Levenshtein, Cosine , and Jarro-Winkler. </br>

4. We find the closest matching WHO term for each tumor for each distance metric and then also group the top 0.05% closest matching group of tumors. Each of the closest match terms are standardized to their closest matching WHO Term. </br>

5. We also use the distance matrices computed to perform 3 levels of nested affinity clustering and group the tumors. After grouping the tumors they are standardized to their closest matching WHO Term. </br>

6. We generate embeddings for each tumor terms (CT, WHO, NCIT) using Open AI's ADA2.0 and V-3 Large text-embedding models. We then compute the closest matching (Euclidean Distance) WHO terms for each tumor.  </br>

7. We also perform PCA on each of the embedding types and then run K-means and Affinity Clustering to group the tumors together. We refine the clusters by filtering outliers using isolation forest and local outlier factor.  </br>

8. After cluster refinement, each cluster is standardized to the WHO term that matches a majority of the members of that cluster.  
   
   
## Run Instructions
Following are the steps for running the pipeline: </br>
1. Clone this Github repository to your local machine </br> 
1. Navigate to the following website: and download the zip file labeled '' </br>
2. Unzip the file and store the files in each of the following directory of the cloned GitHub repository: </br>
	
CT_Embeddings_ADA2.csv			
CT_Embeddings_V3.csv			
NCIT_Embeddings_V3.csv		
WHO_Aggregate_ADA2.csv			
WHO_Terms_All_V3.csv		
combined_embedding_ada2_df.RData	
combined_embedding_ada2_df_5thed.RData	
combined_embedding_v3_df.RData		
combined_embedding_v3_df_5thed.RData
disease_transform_pca_ada2.csv 
dissimilarity_matrix_cosine.RData
dissimilarity_matrix_jw.RData
dissimilarity_matrix_lv.RData
disease_transform_pca_v3.csv
script-5B.RData
script-5b_5thed.Rdata
script5a_5thed.RData
script5a.RData

3. 
   


The results of this study is used to standardize the tumor names in CT database, so they can be used in a meaningful way for further downstream analysis of relating drugs and targets to tumors.

## Scripts

00-generate-ct-disease-file.R:  </br> 
This script loads data from clinical trials and select only the diseases with NCT ID , and associated with Intervention types of Drug, Biological,Combination Product,Genetic. Totally 50410 diseases are extracted. </br> 


01-generate-disease-annotation-for-manual-review.R </br>
This script annotates the 50K diseases automatically as cancer or not. </br>

02-calculate-edit-distance.R </br>
02-calculate-edit-distance-5thed.R</br>
This script loads the manually annotated disease file with pediatric and adult cancer annotation. <br/>

03-edit-distance-clustering.R </br>
03-edit-distance-clustering-5thed.R </br>
This Script computes edit distances between clinical trials tumors, WHO tumors, and NCIT tumors.Performs Affinity Cluster with 3 levels of nesting. </br>

04-preprocess-embedding-pca.R </br>
04A-preprocess-embedding-pca-ADA2-5thed.R
This script loads embeddings V3 and ADA2 for CT, WHO, NCIT  Tumors and then performs PCA. </br>

04B-preprocess-embedding-pca-v3.R </br>
04B-preprocess-embedding-pca-v3-5thed.R </br>


05A-cluster-on-ADA2-embedding-Kmeans.R </br>
05A-ADA2-embedding-Kmeans-5thed.R </br>
This script computes Kmeans cluster of ADA2 data and also computes silhouette index. </br>

05B-v3-embedding-Kmeans.R </br>
05B-v3-embedding-Kmeans-5thed.R </br>


05B-cluster-on-v3-embedding-Kmeans.R </br>
This script computes Kmeans cluster of V3 data and also computes silhouette index. </br>

06A-cluster-on-ADA-embedding-affinity.R </br>
06A-cluster-on-ADA-embedding-affinity-5thed.R </br>
This script computes affinity cluster of ADA2 data. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership. <br/>


06B-cluster-on-V3-embedding-affinity.R </br>
06B-cluster-on-V3-embedding-affinity-5thed.R </br>
This script computes affinity cluster of V3 data. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.


07A-annotate-cluster-result-NCIT-WHO.R </br>
07A-annotate-cluster-result-NCIT-WHO-5thed.R
This script annotates Affinity cluster results of ADA2 embeddings. </br>

07B-annotate-cluster-result-V3-NCIT-WHO.R </br>
07B-annotate-cluster-result-V3-NCIT-WHO-5thed.R </br>
This script annotates Affinity cluster results of V3 embeddings. </br>


08-outlier-detection-embeddings.R </br>
08-outlier-detection-embeddings-5thed.R </br>
This script is used to detect if Affinity cluster members are outliers using LOF and Isolation Forest on ADA2 and V3 embedding data </br>

09-cluster-reassignment-outlier.R </br>
09-cluster-reassignment-outlier-5thed.R </br>
This script performs Affinity cluster reassignment after outlier detection </br>

10-assign-who-ncit-outlier-kmeans-editdistance-clustering.R </br>
10-assign-who-ncit-outlier-kmeans-editdistance-clustering-5thed.R </br>
This script to detect outliers for Kmeans and editdistance based standardization </br>

11-generate-records-annotation.R </br>
This script is used to annotate the types of ground truth found for each of the 1600 tumors sampled. </br>

12--summarize-results.R </br>
Prints the accuracy of each standardization method. </br>

13-plot-script.R </br>
Script is used for generating silhouette plots. </br>
