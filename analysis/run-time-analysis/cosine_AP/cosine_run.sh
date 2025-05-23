#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

RScript cosine_AP_5th.R
RScript cosine_AP_all.R