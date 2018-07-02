#!/usr/bin/perl
use warnings;
use strict;

opendir( KO, 'ko' ) or die $!;
open( KEGG_KO, '>', 'keggToKO.tsv' ) or die $!;

while ( readdir KO ) {
	my $file = $_;
	open( KO_GENE, '<', "ko/$file" ) or die $!;
	my $def = 'NA';
	while (<KO_GENE>) {
		chomp();
		my $line = $_;
		if ( $line =~ /DEFINITION\s+(.*?)$/ ) {
			$def = $1;
		}
		if ( $line =~ /([A-Z]{3,3}):\s+(.*?)$/ ) {
			my $org        = lc($1);
			my $geneString = $2;
			#print $geneString, "\n";
			my @genes      = split(" ",$geneString);
			foreach my $gene (@genes) {
				if($gene =~ /(.*?)\(/){
					$gene = $1;
				}
				print KEGG_KO join( "\t", "$org:$gene", $file, $def ), "\n";
			}
		}
	}
	close(KO_GENE);
}
