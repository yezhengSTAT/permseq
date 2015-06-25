###################################################################
#
#   Calculate averaged ChIP counts at different read counts (or group defined by single data) defined in dnase_infile
#
#   Command arguments: 
#   #   
#   - chip: Uni-reads in BED format
#   - dnase_infile: File encodes which segments of the genome are in which group based on one data set.
#   - outfile: Name of output file
###################################################################
use List::Util qw(first);
use List::Util qw(sum);
    

my ($chip,$dnase_infile,$outfile)=@ARGV;



my @dna;


# read in dnase data

 my @total_location;
  my @sum_chip;
 my $initial_pos=0;
 open IN, "$dnase_infile" or die "Unable to open $dnase_infile";
    while (my $line = <IN>) {
        
         $initial_pos=0;
         @dna=();
	  chomp($line);
          my($data_chrt, @other) =  split( /\t/, $line );
                    for (my $i=0;$i<scalar(@other)/2;$i++){
             
             
             my $got_big=$other[2*$i+1];
             my $pos=$initial_pos+$other[2*$i]-1;
             if($got_big!=0){
               for (my $j = $initial_pos; $j <= $pos; $j++) {
                 $dna[$j] = $got_big;}
                            }
              $initial_pos=$pos+1;
              if ( exists $total_location[$got_big]){
                  $total_location[$got_big]=$total_location[$got_big]+$other[2*$i];}else{
                     $total_location[$got_big]=$other[2*$i];
                    

                                                     }
           # print "totalafter".$total_location[$got_big]."\n";

         }
      





 # read in chip

 open IN1, "$chip" or die "Cannot open $chip\n";
 my %data;
 my %data_chr_list;
 while (my $line = <IN1>) {
	chomp($line);
       
       my($chrt, $pos,$end,$id, $score,$str) =  split( /\t/, $line );
         
         if($chrt eq $data_chrt){
           
            if( exists $dna[$pos]){
                  
          
                  $sum_chip[$dna[$pos]]++; }else{
          
                 $sum_chip[0]++;
         }
                      }        
 
}
 close IN1;



}
close IN;
open OUT, ">$outfile" or die "Cannot open $outfile\n";

for ( my $loops = 0; $loops < scalar(@total_location); $loops++ )
{
      if ( !exists $sum_chip[$loops]){        
               $sum_chip[$loops]=0;
       }
}

 print OUT join("\t", @total_location)."\n";
 print OUT join("\t", @sum_chip)."\n";
close OUT;
