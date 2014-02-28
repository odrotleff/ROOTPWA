#!/bin/bash
# start 4 fitting jobs simultaneusly!

IDS=( 1 2 3 4 5 6 7 8 9 )
#IDS=( 5 6 7 8 )
#IDS=( 9 10 11 12 )

for i in {1..6}; do
#for i in 0 ; do
export SGE_TASK_ID=${IDS[$i]}
echo "Starting Job $SGE_TASK_ID"
./rootfitBatch.sge > batch$SGE_TASK_ID.log &
done