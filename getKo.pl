#!/usr/bin/perl
use warnings;
use strict;

open(KO, '<', "listKo.txt") or die $!;

while(<KO>){
	chomp();
	my $ko = $_;
	`wget http://rest.kegg.jp/get/$ko`;
}
