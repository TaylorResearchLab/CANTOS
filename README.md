# CT-Embedding-Paper

## Description
This repository contains the code, tables, and plots associated with the CT Embedding paper. The pipeline built in this repository does the following task: <br/>
1. Extract tumor names from the CT database <br/>.
2. Cluster tumor names using distance metrics and show it is not an effective method <br />
3. Cluster tumor names using methods such as KNN and show it is not an effective method <br />
4. Finally, use embeddings generated from OpenAI's ADA 2.0 and cluster using affinity propagation <br/>
5. For each cluster find the closest matching NCIT or WHO Term. <br />

The results of this study is used to standardize the tumor names in CT database, so they can be used in a meaningful way for further downstream analysis of relating drugs and targets to tumors.


