#!/usr/bin/env bash
set -e

# 1. Imputation 15MA by sampling 16MA and 14MA
awk -f- <<EOF <(gzip -dc ) <(gzip -dc)
EOF
# 2. Imputation 14L by sampling 15L and 13L
# 3. Imputation 14R by sampling 15L and 15R
# 4. Downsampling of duplciated samples
# 5. Discarding of low quality samples
