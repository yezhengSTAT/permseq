###################################################################
#
#   Count number of reads aligned to each nucleotide, no extension and only reports position with at least on tags aligned to.
#
#   Command arguments: 
#   #   
#   - bed: Aligned uni-reads in BED format
#   - outfile: Name of output file
###################################################################
#!/usr/bin/perl
    use strict;
    use warnings;
    use Data::Dumper;

    my $bed = $ARGV[0];
    
    my $outfile=$ARGV[1];
    my %h; 
    my $c=0;
       open IN, "$bed" or die "Unable to open $bed";
    while (<IN>){

    my @data = split /\s+/, $_;
    #my $read_length=$data[2]-$data[1]+1;
    #$data[1] = $data[1] -$read_length +1 if ( $data[5] eq "-" );

    my $t=join(" ",$data[0],$data[1]);
    if(exists $h{$t}){
    $h{$t}=$h{$t}+$data[4]/1000;}else{
    $h{$t}=$data[4]/1000;
    }
 
    }

    close( IN );
open OUT, ">$outfile" or die "Cannot open $outfile\n";
foreach my $content (keys %h) {
my @temp=split /\s+/, $content;
print OUT "$temp[0]\t$temp[1]\t$h{$content}\n";

}
close(OUT);
