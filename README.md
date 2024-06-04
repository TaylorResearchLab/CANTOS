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
   
   
   

The results of this study is used to standardize the tumor names in CT database, so they can be used in a meaningful way for further downstream analysis of relating drugs and targets to tumors.


