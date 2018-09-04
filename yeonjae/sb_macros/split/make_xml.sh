#!/bin/bash


for f in `ls *root`
do
	#echo $f
	arrIN=(${f//./ })
	tmp=${arrIN[0]} 
	echo $tmp

	cat template.xml >> full.xml
	sed -i "s/XXXX/${tmp}/g" full.xml
	sed -i "s/YYYY/${f}/g" full.xml

done
