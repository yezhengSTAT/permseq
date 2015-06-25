###################################################################
#
#   Convert SAM to uni-reads BED format
#
#   Command arguments: 
#   #   
#   - sam: Aligned reads in SAM format
#   - outfile: Name of output file
###################################################################
#!/usr/bin/perl

    use strict;
    use warnings;
    use Data::Dumper;

    my $sam = $ARGV[0];
    my $outfile=$ARGV[1];
    my %h; 
    my $c=0;
       open IN, "$sam" or die "Unable to open $sam";
#read in SAM file
    while (my $line = <IN>){
      chomp($line);
      if ( $line =~ /^[^@].+/ ) { 
         my @data = split /\s+/, $line;
         my $t=$data[0];
         $h{$t}++ ; 
    
 
    }
   }

    close( IN );

open IN, "$sam" or die "Unable to open $sam";
open OUT, ">$outfile" or die "Cannot open $outfile\n";
# find uni-reads from SAM file and output in BED format
 my $str;
my $end;
while (my $line = <IN>) {
	chomp($line);

   
    if ( $line =~ /^[^@].+/ ) { 
    
    my ($t1, $bwflag, $chrt, $pos, $t2, $t3, $t4, $t5, $t6, $seq, $t7 ) =  split( /\s+/, $line );
    if ( $h{$t1}==1) {
       $pos = int($pos)-1;  #bed start from 0
       my $read_length = length $seq;
       $end=$pos+$read_length-1;
        
         if ( $bwflag & 4 or $bwflag & 512 or $bwflag & 1024 ) {
            # exclude invalid lines
                
         
        } elsif ( $bwflag & 16 ) {
            # reverse strand
            
            $str = "-";
            print  OUT "$chrt\t$pos\t$end\t$t1\t1000\t$str\n";
        } else {
            $str = "+";
            print  OUT "$chrt\t$pos\t$end\t$t1\t1000\t$str\n";
        }   
   
     
     
    }
   
    }

}


    

close IN;
close OUT
