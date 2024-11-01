#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

Rscript 00-generate-ct-disease-file.R

Rscript 01-generate-disease-annotation-for-manual-review.R

Rscript 02-calculate-edit-distance.R
Rscript 02-calculate-edit-distance-5thed.R

Rscript 03-edit-distance-clustering.R
Rscript 03-edit-distance-clustering-5thed.R

Rscript 04A-preprocess-embedding-pca.R
Rscript 04A-preprocess-embedding-pca-ADA2-5thed.R

Rscript 04B-preprocess-embedding-pca-v3.R
Rscript 04B-preprocess-embedding-pca-v3-5thed.R

Rscript 05A-cluster-on-ADA2-embedding-Kmeans.R
Rscript 05A-ADA2-embedding-Kmeans-5thed.R

Rscript 05B-v3-embedding-Kmeans.R
Rscript 05B-v3-embedding-Kmeans-5thed.R

Rscript 06A-cluster-on-ADA-embedding-affinity.R
Rscript 06A-cluster-on-ADA-embedding-affinity-5thed.R

Rscript 06B-cluster-on-V3-embedding-affinity.R
Rscript 06B-cluster-on-V3-embedding-affinity-5thed.R

Rscript 07A-annotate-cluster-result-NCIT-WHO.R 
Rscript 07A-annotate-cluster-result-NCIT-WHO-5thed.R

Rscript 07B-annotate-cluster-result-V3-NCIT-WHO.R
Rscript 07B-annotate-cluster-result-V3-NCIT-WHO-5thed.R

Rscript 08-outlier-detection-embeddings.R
Rscript 08-outlier-detection-embeddings-5thed.R

Rscript 09-cluster-reassignment-outlier.R
Rscript 09-cluster-reassignment-outlier-5thed.R

Rscript 10-assign-who-ncit-outlier-kmeans-editdistance-clustering.R
Rscript 10-assign-who-ncit-outlier-kmeans-editdistance-clustering-5thed.R

Rscript 11-generate-records-annotation.R 
Rscript 12-summarize-results.R 
Rscript 13-plot-script.R
