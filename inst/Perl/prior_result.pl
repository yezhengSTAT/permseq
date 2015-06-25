###################################################################
#
#   Partition genome based on one data set
#
#   Command arguments: 
#   - ref_file: File name for chromosome info (chr size)
#   - dnase_dir: Directory of nucleotide-level counts files
#   - dnase_infile: nucleotide-level counts files
#   - outfile: File name of output
#   - data_chr: Included chromosomes
#   - q: dnaseThres cutoff. Genome partition is according to q
###################################################################
# from csem output
#uni read $5=1000
# 3 groups tag=0,1,2, did not print out

my ($ref_file,$dnase_dir, $dnase_infile,$outfile,$data_chr,@q)=@ARGV;

my @group_indicator=$data_chr;

# read in chromsome sizes
 open IN, "$ref_file" or die "Cannot open $ref_file\n";
    <IN>; # skip first line
   my $chr_len=<IN>;	
   chomp($chr_len);
   my @chr_list_len= split( /\s+/, $chr_len );
 
   my $chr_vec = <IN>;
   chomp($chr_vec);
   my @chr_list = split( /\s+/, $chr_vec );
   my %chr_ref;
   for( my $i=0;$i<scalar(@chr_list);$i++){
     $chr_ref{$chr_list[$i]}=$chr_list_len[$i];}

# read in nucleotide-level counts file of given chromosome
# partition genomes based on q values
   open OUT, ">$outfile" or die "Cannot open $outfile\n";



     my  $old_alpha;
     my  $new_alpha;
     my $old_position=0;###should be 0 before -1
     my $new_position;
     my @dna;
     my $dnase_infile_new=$dnase_dir.$data_chr.$dnase_infile;
     open IN, "$dnase_infile_new" or die "Unable to open $dnase_infile_new";
     my $line = <IN>;
     chomp($line);
      # procee read file, based on "format" option
      my ($chrt, $pos, $tag) =  split( /\t/, $line );
      my $a= find_location($tag,\@q);
      $new_alpha=$a;
      $old_alpha=$new_alpha;
       
      
      if ($pos!=1){
              
              if($new_alpha!=0){
                push @group_indicator, $pos-1;## start from 1
                push @group_indicator ,0;
                $old_position=$pos-1;
               
                }
           
                  }
       $new_position=$pos;

while (my $line = <IN>) {
	  chomp($line);

        my ($chrt, $pos, $tag) =  split( /\t/, $line );
        if($pos<=$chr_ref{$data_chr}){
        my $a= find_location($tag,\@q);
        $new_alpha=$a;
        
          if($pos!=$new_position+1&$old_alpha!=0){
          push @group_indicator, $new_position-$old_position;
          push @group_indicator, $old_alpha;
          $old_position=$new_position;
          $new_position=$pos-1;
          $old_alpha=0;
                                          }
         if($pos!=$new_position+1&$old_alpha==0){
          $new_position=$pos-1;         
                                         }
         
         if($new_alpha!=$old_alpha){
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
}
close IN;
 #open OUT1, ">diag.txt";
 $pos=$new_position;
 #print OUT1 $chr_ref{$data_chr};
 #print OUT1 $pos;
# if current position does not reach the chromosomsize
 if($pos<$chr_ref{$data_chr}){
        if($new_alpha==0){
         push @group_indicator, $chr_ref{$data_chr}-$old_position; ## no -1
         push @group_indicator, 0;} 
        else{ 
         push @group_indicator, $pos-$old_position;
         push @group_indicator, $new_alpha;
         push @group_indicator, $chr_ref{$data_chr}-$pos;## no -1
         push @group_indicator, 0;
             } 

     
                                     
}
if($pos==$chr_ref{$data_chr}){
         push @group_indicator, $pos-$old_position;
         push @group_indicator, $new_alpha;
}
  

print OUT join("\t", @group_indicator);
#print join("\t", @group_indicator);
print OUT "\n";
      close OUT;       

sub find_location {
    my ($tag,$q) =@_ ; 

my $got_big=scalar(@q);
if($tag<=$q[0]){
$got_big=0;}
else{
 for (my $i=1;$i<scalar(@q);$i++){
  if($tag<=$q[$i]){
       $got_big=$i; 
        last;}
     }
    }
return $got_big;

}

          
