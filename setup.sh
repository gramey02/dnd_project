#!/bin/bash

[SEARCH1="](mailto:SEARCH1=%22new.email@example.com)Grace.Ramey@ucsf.edu"
[REPLACE1="](mailto:REPLACE1=%22new.email@example.com)$2"

SEARCH2="/wynton/home/capra/gramey02/dnd_project”
REPLACE2="$1"

DIR="$1"

if [ -z "$DIR" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# Replace recursively in all files
find "$DIR" -type f -exec sed -i "s|$SEARCH1|$REPLACE1|g" {} +
find "$DIR" -type f -exec sed -i "s|$SEARCH2|$REPLACE2|g" {} +

echo "Done replacements in $DIR"
