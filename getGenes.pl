#!/usr/bin/perl
use warnings;
use strict;

my $org = $ARGV[0];
open(ORG, '<', "uniprot.$org.tsv") or die $!;

while(<ORG>){
	chomp();
	my $line = $_;
	my($sp, $kegg) = split("\t", $line);
	`wget http://rest.kegg.jp/get/$kegg`;
}
