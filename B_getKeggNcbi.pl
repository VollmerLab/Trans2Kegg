#!/usr/bin/perl
use warnings;
use strict;

`wget -O orgs.tsv http://rest.kegg.jp/list/organism`;
open(ORGS,'<', 'orgs.tsv') or die $!;

while (<ORGS>){
	chomp();
	my @fields = split("\t", $_);
	my $org = $fields[1];
	`wget -O ncbi.$org.tsv http://rest.kegg.jp/conv/$org/ncbi-geneid`;
}

