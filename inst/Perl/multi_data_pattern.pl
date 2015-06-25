###################################################################
#
#   Partition genome based on multiple data sets.
#
#   Command arguments: 
#   #
#   - outfile: Name of output file.   
#   - chrlength: length of chromosome.
#   - genome_cluster_file: Names of Files, each encodes which segments of the genome are in which group based on data set.
###################################################################
use List::Util qw(first);
use List::Util qw(sum);
    

my ($outfile,$chrlength,@genome_cluster_file)=@ARGV;


my $chr;
my @dna;
my @temp;
my @temp1;
my %total_position;
my %chip_position;
my %sum_chip;
my %notback;
for(my $i=0;$i<scalar(@genome_cluster_file);$i++){
 $temp[$i]=0;
 $temp1[$i]=0;}
my $back=join("_", @temp);
# read in dnase data

 my @total_location;
  my @sum_chip;
 my $initial_pos=0;
# read in first file
 open IN, "$genome_cluster_file[0]" or die "Unable to open $genome_cluster_file[0]";
    my $line = <IN>; 
  $initial_pos=0;
  
	  chomp($line);
          my($data_chrt, @other) =  split( /\t/, $line );
          $chr=$data_chrt; 
          
                    for (my $i=0;$i<scalar(@other)/2;$i++){
             
            
             my $got_big=$other[2*$i+1];
             my $pos=$initial_pos+$other[2*$i]-1;
             if($got_big!=0){
               $temp[0]=$got_big;
               for (my $j = $initial_pos; $j <= $pos; $j++) {
                $notback{$j}++;
                
                 $dna[$j] = join("_", @temp) ;
               
                                             }
              }
              $initial_pos=$pos+1;
              }
close IN;


# read other files and partition genome.
  
 for(my $file_i=1;$file_i<scalar(@genome_cluster_file);$file_i++){
     open IN, "$genome_cluster_file[$file_i]" or die "Unable to open $genome_cluster_file[$file_i]";
     $initial_pos=0;   
      my $line = <IN> ;
      chomp($line);
      
    
        	 
          my($data_chrt, @other) =  split( /\t/, $line );
                      for (my $i=0;$i<scalar(@other)/2;$i++){
             
             
             my $got_big=$other[2*$i+1];
             my $pos=$initial_pos+$other[2*$i]-1;
             if($got_big!=0){
               for (my $j = $initial_pos; $j <= $pos; $j++) {
                 $notback{$j}++;
                 if(exists $dna[$j]){
                 
                 @temp= split( "_", $dna[$j] );
                 $temp[$file_i]=$got_big;
                 $dna[$j] =join("_", @temp) ;
                 #$total_position {$dna[$j]}++;  
                 }else{
                 @temp=@temp1;
                 $temp[$file_i]=$got_big;
                 $dna[$j] =join("_", @temp) ;
                  
                 } 
                 #print "$j\t$dna[$j]\n";
                            }
               }
              $initial_pos=$pos+1;
                 
         }
close  IN;
}

####write group id#######

my @group_indicator=$data_chrt;
my @selectposions=keys %notback;


my @sort_selectposions=sort { $a <=> $b }@selectposions;


open OUT, ">$outfile" or die "Cannot open $outfile\n";
     my  $old_alpha;
     my  $new_alpha;
     my $old_position=-1;###start from 0
     my $new_position;
     my $pos=$sort_selectposions[0];
     my $a=$dna[$sort_selectposions[0]];
           $new_alpha=$a;
      $old_alpha=$new_alpha;
       
   
      if ($pos!=0){
              
              if($new_alpha ne $back){
                push @group_indicator, $pos;## start from 0
                push @group_indicator ,$back;
                $old_position=$pos-1;
               
                }
           
                  }
       $new_position=$pos;

for(my $i=1;$i<scalar(@sort_selectposions);$i++){
	 $pos=$sort_selectposions[$i];
        
        my $a=$dna[ $pos];
        
        $new_alpha=$a;
        #print "$pos\t$new_position\t$old_position\t$new_alpha\t$old_alpha\n";
          if($pos!=$new_position+1&$old_alpha ne $back){
          push @group_indicator, $new_position-$old_position;
          push @group_indicator, $old_alpha;
          $old_position=$new_position;
          $new_position=$pos-1;
          $old_alpha=$back;
                                          }
         if($pos!=$new_position+1&$old_alpha eq $back){
          $new_position=$pos-1;         
                                         }
         
         if($new_alpha ne $old_alpha){
           
          push @group_indicator, $new_position-$old_position;
          push @group_indicator, $old_alpha;
          $old_position=$pos-1;
          $new_position=$pos;
                          }

          else{
           $new_position=$pos;
              }
          $old_alpha=$new_alpha;
}
close IN;

 $pos=$new_position;

 if($pos<$chrlength-1){
        if($new_alpha eq $back){
         push @group_indicator, $chrlength-$old_position-1; ## -1
         push @group_indicator, $back;} 
        else{ 
         push @group_indicator, $pos-$old_position;
         push @group_indicator, $new_alpha;
         push @group_indicator, $chrlength-$pos-1;##  -1
         push @group_indicator, $back;
             }

  if($pos==$chrlength-1){
 push @group_indicator, $chrlength-$old_position-1;
 push @group_indicator,$new_alpha;

}

     
                                     
}


print OUT join("\t", @group_indicator);
#print join("\t", @group_indicator);
print OUT "\n";

close OUT;
     

