# CT-Embedding-Paper

## Description
This repository contains the code, tables, and plots associated with the CT Embedding paper. The pipeline built in this repository does the following task: <br/>

1. Extract tumor names from the CT database if they are associated with an NCT ID and have an associated drug belonging to the categories  of Drug, Biological,Combination Product,Genetic. A total 50,410 diseases are extracted.<br/>

2. These 50410 diseases are flagged as tumors and non tumors by the pipeline, which are then further manually annotated pediatric and adult tumors. A total of 13,329 tumors are identified from the 50,410 diseases.<br/> 

3. We compute the distance of each  13,329 clinical trials tumors, 4720 WHO tumors, and 1395 NCIT tumors.  Distance metrics used are Levenshtein, Cosine , and Jarro-Winkler. </br>

4. We find the closest matching WHO term for each tumor for each distance metric and then also group the closest 0.05% matching group of tumors.
   
   
   
2. Cluster tumor names using distance metrics and show it is not an effective method <br />
3. Cluster tumor names using methods such as KNN and show it is not an effective method <br />
4. Finally, use embeddings generated from OpenAI's ADA 2.0 and cluster using affinity propagation <br/>
5. For each cluster find the closest matching NCIT or WHO Term. <br />

The results of this study is used to standardize the tumor names in CT database, so they can be used in a meaningful way for further downstream analysis of relating drugs and targets to tumors.


