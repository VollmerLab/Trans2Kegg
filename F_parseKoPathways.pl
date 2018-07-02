#!/usr/bin/perl
use warnings;
use strict;

opendir( KO, 'ko' ) or die $!;

open( KEGG_PATH, '>', '../koToPathway.tsv' ) or die $!;
print KEGG_PATH join( "\t", "Gene", "Path", "PathDesc" ), "\n";

while ( readdir KO ) {
	my $file = $_;
	open( KO_GENE, '<', "ko/$file" ) or die $!;
	my $def = 'NA';
	while (<KO_GENE>) {
		chomp();
		my $line = $_;
		if ( $line =~ /\s+(ko\d+)\s+(.*?)$/ ) {
			my $path = $1;
			my $desc = $2;
			$file =~ s/ko://;
			print KEGG_PATH join( "\t", $file, $path, $desc ), "\n";
		}
	}
	close(KO_GENE);
}
