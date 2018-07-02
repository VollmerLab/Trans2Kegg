#!/bin/bash
wget -O orgs.tsv http://rest.kegg.jp/list/organism
cut -f2 orgs.tsv > orgs2.tsv
mkdir -p orgs
while read org
do
  wget -P orgs -O uniprot.$org.tsv http://rest.kegg.jp/conv/$org/uniprot
done <orgs2.tsv

