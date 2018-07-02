#!/usr/bin/perl
use warnings;
use strict;

open( UNI_TO_KEGG, '<', 'uniToKegg.tsv' )  or die $!;
open( KEGG_TO_KO,  '<', 'keggToKO.tsv' )   or die $!;
open( UNI_TO_KO,   '>', 'uniToKO.tsv' )    or die $!;
open( TX2GENE,     '>', '../tx2gene.csv' ) or die $!;
open( KODESC,      '>', '../koDesc.tsv' )  or die $!;
open( COV,      '>', 'kCoverage.csv' )    or die $!;

my %keggToKO;
my %koToDesc;
my %tx2gene;

while (<KEGG_TO_KO>) {
	chomp();
	my $line = $_;
	my ( $kegg, $ko, $desc ) = split( "\t", $line );
	$ko =~ s/ko://;
	$keggToKO{$kegg}{$ko}++;
	$koToDesc{$ko} = $desc;
}

while (<UNI_TO_KEGG>) {
	chomp();
	my $line = $_;
	my ( $kID, $sp, $slen, $alen, $ident, $evalue, $bitscore, $kegg ) =
	  split( "\t", $line );
	if ( $keggToKO{$kegg} ) {
		foreach my $ko ( %{ $keggToKO{$kegg} } ) {
			my $desc = $koToDesc{$ko};
			my $cov  = $alen / $slen;
			$cov = sprintf( "%.2f", $cov );
            if ( $desc and ( ( $evalue <= 1e-5 ) and ( $cov > .5 ) ) ) {
                #if ( $desc and ( $evalue <= 1e-5 ) ) {
				print UNI_TO_KO join( "\t",
					$kID,  $sp,    $slen,   $alen,
					$cov,  $ident, $evalue, $bitscore,
					$kegg, $ko,    $desc ),
				  "\n";
				$tx2gene{$kID}{$ko}++;
				print COV join(",", $kID, $ko, $cov), "\n";
				
			}
			else {
				#print join("\t", $kegg, $ko),"\n";
			}
		}
	}
	else {
		#print join( "\t", $sp, $kegg ), "\n";
	}
}

my %printed;
foreach my $kID ( keys %tx2gene ) {
	my %kos = %{ $tx2gene{$kID} };
	foreach my $ko (
		sort { $kos{$b} <=> $kos{$a} }
		keys %kos
	  )
	{
		if ( not $printed{$kID} ) {
			print TX2GENE join( ",", $kID, $ko ), "\n";
			$printed{$kID}++;
		}
	}
}
print KODESC join( "\t", "Gene", "GeneDesc" ), "\n";
foreach my $ko ( keys %koToDesc ) {
	print KODESC join( "\t", $ko, $koToDesc{$ko} ), "\n";
	#print KODESC join( "\t", $ko, $koToDesc{$ko} ), "\n";
}
