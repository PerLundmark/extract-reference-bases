#!/usr/bin/perl -w
#
# A script to extract reference bases from fasta files for VCF export from Chiasma.
# The script use a csv-manifest file for the infinium chip as well as a reference fasa file 
# and creates a file with the list of reference bases and contig names for VCF generation.
#

use lib '/usr/lib/bioperl/BioPerl-1.6.1';

use Bio::SeqIO;

my $manifest = shift;
my $reference = shift;
my $outfile = shift;

my $seq_io = Bio::SeqIO->new(-file => $reference);
my @header_lines = (); #the header for the chiasma VCF file

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

#general part of VCF header
push @header_lines, '##fileformat=VCFv4.1';
push @header_lines, '##fileDate=' . sprintf("%04d:%02d:%02d", $year + 1900, $mon + 1, $mday);
push @header_lines, '##source=extract_reference_bases.pl_and_chiasma';
push @header_lines, '##reference=file://' . $reference;
push @header_lines, '##phasing=none';

# Make a hash of references to sequence objects for all sequences in the reference file
# and build the contig part of the header for the Chiasma VCF export.
#
my %sequences = ();
while (my $seq = $seq_io->next_seq){
	my $seq_name = $seq->id();
	print "Reading reference sequences\n";
	print $seq_name . "\n";
	print $seq->desc . "\n";
	print $seq->length . "\n";
	
	if (exists($sequences{$seq_name})){
		print "!!!  Warning: $seq_name was replaced by multiple occurences of the same sequence name !!!\n";
	}
	$sequences{$seq_name} = $seq;

	#header
	push @header_lines, '##contig=<ID=' . $seq_name . ',length=' . $seq->length . ',assembly=b37>';
}

push @header_lines, '##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\'Genotype\'>';

#Open the infinium manifest file
open MANIFEST, $manifest or die "Unable to open manifest!\n";

#skip header and set index for fields from the manifest to account for variable order in manifests.
my $idx_snp_name = -1;
my $idx_snp_chr = -1;
my $idx_snp_pos = -1;
my $idx_snp_var = -1;
my $idx_snp_ref_strand = -1;

while (my $row = <MANIFEST>){
	if ($row =~ /IlmnID/){
		my @header_fields = split /,/, $row, -1;
		for (0..$#header_fields){
			if ($header_fields[$_] eq 'Name'){
				$idx_snp_name = $_;
			}elsif($header_fields[$_] eq 'Chr'){
				$idx_snp_chr = $_;
			}elsif($header_fields[$_] eq 'MapInfo'){
				$idx_snp_pos = $_;
			}elsif($header_fields[$_] eq 'SNP'){
				$idx_snp_var = $_;
			}elsif($header_fields[$_] eq 'RefStrand'){
				$idx_snp_ref_strand = $_;
			}
		}
		last;
	}
}

#process output and print 
my $counter = 0;
open OUTFILE, ">", $outfile or die "Unable to open outfile!\n";

#print header info for VCF
print OUTFILE "#vcf_header\n";
foreach my $line (@header_lines){
	print OUTFILE $line . "\n";
}

#print ref base data
print OUTFILE "#refbase_data\n";
print OUTFILE join("\t", ('snp_name', 'snp_chr', 'snp_pos', 'snp_var', 'snp_ref_strand', 'ref_name', 'ref_comment', 'ref_base', 'ref_length')), "\n";
while ($row = <MANIFEST>){
	if($row =~ /Controls/){last;}
	
	chomp($row);
	$row =~ s/\r//;
	$row =~ s/\n//;
	my @fields = split(',', $row, -1);
	
	my $snp_name = $fields[$idx_snp_name];
	my $snp_chr  = $fields[$idx_snp_chr];
	my $snp_pos  = $fields[$idx_snp_pos];
	my $snp_var  = $fields[$idx_snp_var];
	my $snp_ref_strand = $fields[$idx_snp_ref_strand];
	
	my $current_ref = "";
	my $current_ref_desc = "";
	my $current_ref_length = "";
	my $current_ref_id = "";
	my $current_ref_base = '';
	
	#deal with markers in the infinium manifest with pos 0 and chr 0
	if ($snp_pos == 0){
		if (($snp_chr eq "0") || ($snp_chr eq "Unknown")){
			$current_ref_id = "Unknown";
			$current_ref_desc = "Position was given as 0 in Infinium manifest";
			$current_ref_length = 1;
			$current_ref_base = 'N';
		}else{
                        $current_ref = $sequences{$snp_chr eq "XY"?"X":$snp_chr};
			$current_ref_desc = $current_ref->desc;
			$current_ref_length = $current_ref->length;
			$current_ref_id = $current_ref->id();
			$current_ref_base = 'N';
		}
	}else{
		$current_ref = $sequences{$snp_chr eq "XY"?"X":$snp_chr};
		$current_ref_desc = $current_ref->desc;
		$current_ref_length = $current_ref->length;
		$current_ref_id = $current_ref->id();

		if ($snp_pos > $current_ref_length){
			$current_ref_base = 'N';
			print $snp_name . " has position " . $snp_pos . " outside the reference sequence " . $current_ref_id . ", length: " . $current_ref_length . "\n";
		}else{
			$current_ref_base = $current_ref->subseq($snp_pos, $snp_pos);
		}
	}

	print $counter++ , "\n";

	print OUTFILE join("\t", ($snp_name, $snp_chr, $snp_pos, $snp_var, $snp_ref_strand, $current_ref_id, $current_ref_desc, $current_ref_base, $current_ref_length));
	print OUTFILE "\n";
	
}
