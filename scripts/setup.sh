#!/bin/bash

SEARCH1="Grace.Ramey@ucsf.edu"
REPLACE1="<your_email_here>"
DIR="<dnd_project_filepath_here>"

find "$DIR" -type f -exec sed -i "s|$SEARCH1|$REPLACE1|g" {} +
echo "Done replacing '$SEARCH1' with '$REPLACE1' in $DIR"

SEARCH2="/wynton/home/capra/gramey02/dnd_project"
REPLACE2="<absolute_dnd_project_filepath"

find "$DIR" -type f -exec sed -i "s|$SEARCH2|$REPLACE2|g" {} +
echo "Done replacing '$SEARCH2' with '$REPLACE2' in $DIR"