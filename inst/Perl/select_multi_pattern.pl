###################################################################
#
#   Repartition genomes using a subset of data sets based on partitioned genomes using all data sets.
#
#   Command arguments: 
#   - patternnamefile: a vector of the unique pattern ID of selected data sets. 
#   - outfile: Name of output files
#   - prior_infile: Input file (directory + file name), the original genome clustering file based on all data sets.
#   - selectindex: the index of data sets used in the final model.

###################################################################


my ($patternnamefile,$outfile,$prior_infile,@selectindex)=@ARGV;




my @temp;
my @temp1;
my @index;
my %keys_index;
my $old_alpha;
my $new_alpha;
open IN, "$patternnamefile" or die "Unable to open $patternnamefile";

# the final group ID should be integers, convert pattern id to integers.
    my $line = <IN>; 

    chomp($line);

    my(@index) =  split( /\t/, $line );

for(my $i=0;$i<scalar(@index);$i++){

  $keys_index{$index[$i]}=$i;

}

close IN;
my @group_indicator;

open OUT, ">$outfile" or die "Cannot open $outfile\n";
 my @total_location;
  my @sum_chip;
 my $initial_pos=0;
 open IN, "$prior_infile" or die "Unable to open $prior_infile";
# read in original data and repartition the genome.
 while (my $line = <IN>) {


   chomp($line);
          my($data_chrt, @other) =  split( /\t/, $line );
            @group_indicator=();
           @group_indicator=$data_chrt;
          my $i=0;
          my $got_big=$other[2*$i+1];
          @temp=split('_',$got_big);
          for(my $j=0;$j<scalar(@selectindex);$j++){

            $temp1[$j]=$temp[$selectindex[$j]];}
         
          my $new_id= join('_',@temp1);
             $old_alpha=$new_id;
          #print $new_id."\n";
          my $pos=$other[2*$i];
          push @group_indicator, $pos;
          push @group_indicator,$keys_index{$new_id};
           for (my $i=1;$i<scalar(@other)/2;$i++){
             my $got_big=$other[2*$i+1];
             my $pos=$other[2*$i];

              @temp=split('_',$got_big);
              for(my $j=0;$j<scalar(@selectindex);$j++){

                 $temp1[$j]=$temp[$selectindex[$j]];}
         
             my $new_id= join('_',@temp1);
             $new_alpha=$new_id;
             #print $new_id."\n";
                if($new_alpha ne $old_alpha){
                push @group_indicator, $pos;## start from 0
                push @group_indicator ,$keys_index{$new_id};
                $old_alpha=$new_alpha;
           }else{
                $group_indicator[scalar(@group_indicator-2)]=$group_indicator[scalar(@group_indicator-2)]+$pos;
                
             }
           
          }


print OUT join("\t", @group_indicator);
#print join("\t", @group_indicator);
print OUT "\n";
}
close OUT;
close IN;
     
