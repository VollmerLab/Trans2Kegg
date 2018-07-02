#!/usr/bin/perl
use warnings;
use strict;

open( UNI,         '<', 'uniprot.tsv' )             or die $!;
open( SP,          '<', '../kTransSwissProt2.tsv' ) or die $!;
open( UNI_TO_KEGG, '>', 'uniToKegg.tsv' )           or die $!;

my %uniToKegg;
while (<UNI>) {
	chomp();
	my $line = $_;
	if ( $line =~ /up:(.*?)\s+(.*?)$/ ) {
		my $uniProt = $1;
		my $kegg    = $2;
		$uniToKegg{$uniProt}{$kegg}++;
	}
}

while (<SP>) {
	chomp();
	my $line     = $_;
	my @fields   = split( /\t|\|/, $line );
	my $kID      = $fields[0];
	my $sp       = $fields[2];
	my $slen     = $fields[4];
	my $alen     = $fields[5];
	my $ident    = $fields[6];
	my $evalue   = $fields[7];
	my $bitscore = $fields[8];
	if ( $uniToKegg{$sp} ) {

		foreach my $kegg ( keys %{ $uniToKegg{$sp} } ) {
			print UNI_TO_KEGG join(
				"\t",
				$kID, $sp, $slen, $alen, $ident, $evalue, $bitscore, $kegg
			  ),
			  "\n";
		}
	}
}
