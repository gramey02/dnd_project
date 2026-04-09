#!/bin/bash
#$ -N setup_scripts
#$ -cwd
#$ -o logs/out/setup_scripts.out
#$ -e logs/err/setup_scripts.err

SEARCH="Grace.Ramey@ucsf.edu"
REPLACE="<your email here>"
DIR=""
echo $DIR

# find "$DIR" -type f -exec sed -i "s|$SEARCH|$REPLACE|g" {} +

# echo "Done replacing '$SEARCH' with '$REPLACE' in $DIR"