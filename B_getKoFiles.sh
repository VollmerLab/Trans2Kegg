#!/bin/bash
mkdir -p ko
cut -f1 koList.tsv > koList2.tsv
while read ko
do
  wget -P ko http://rest.kegg.jp/get/$ko
done <koList2.tsv
