#!/bin/bash
for file in RanGenDataBig/data/*
do
	if [ -f "$file" ]
	then
		l=$(echo "$file" | awk -F '_' '{print $5}')
		d=$(echo "$file" | awk -F '_' '{print $7}')
		q=$(echo "$file" | awk -F '_' '{print $9}')
		echo "$file"
		./my_fm $file -l $l -d $d -q $q >> result_after_AP_big
	fi
done