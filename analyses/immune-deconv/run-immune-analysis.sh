#!/bin/bash


set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo "analyses deconvoluted data..."
Rscript analyse-dists.R \
--fracs_file 'results/quantiseq_output.rds' \
--output_dir 'results'