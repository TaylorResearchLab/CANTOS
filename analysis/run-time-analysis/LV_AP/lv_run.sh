#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

RScript LV_AP_5th.R
RScript LV_AP_all.R