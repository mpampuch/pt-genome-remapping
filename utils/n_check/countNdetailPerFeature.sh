#!/usr/bin/bash

featureList=$1

arr=()
while IFS= read -r line; do
   arr+=("$line")
done <$featureList

for i in "${arr[@]}"
do
    inFileName="${i}.fa"
    outFileName="${i}.ncounts"
    result="$(awk -F "\n" 'BEGIN{RS=">"} {subs=gsub(/N/,"",$2)} subs >=1 {print $1,subs}' $inFileName)"
    echo -e "$result" > $outFileName

done