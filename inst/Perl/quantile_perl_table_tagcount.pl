###################################################################
#
#   Generate summary of counts distribution. Return the number of positions with the same counts.
#
#   Command arguments: 
#   - infile: Nucleotide-level counts data.
#   - ref_file: File name for chromosome info (chr size)
#   - outfile: Name of output file
#   - data_chr: included chromosomes
#
###################################################################
# from csem output
#uni read $5=1000
use POSIX;
my ($infile,$ref_file,$outfile,@data_chr)=@ARGV;
my %count_UR= ();
use Statistics::Descriptive;

# read in chromosome sizes
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

my $data_background=0;
my %data_tag=();

my $total_location=0;
for(my $i=0;$i<scalar(@data_chr);$i++){
  my $chr=$data_chr[$i];
  my $new_infile=$chr.$infile;
 my $data_temp=0;
  open IN, "$new_infile" or die "Cannot open $new_infile\n";
while (my $line = <IN>) {
	chomp($line);

    # procee read file count the number of positions with the same count
    # background records positions with counts<=2.
     my ($chrt, $pos, $tag) =  split( /\t/, $line );
    $tag=ceil($tag);
    $data_tag{$tag}++;
     $data_temp++;
                        }
    $data_background=$data_background+$chr_ref{$chr}-$data_temp;
    $total_location=$total_location+$chr_ref{$chr};
   }
 
  my @tag_dist;
  open OUT, ">$outfile" or die "Cannot open $outfile\n";
   foreach my $content(keys %data_tag){
     print OUT "$content\t$data_tag{$content}\n";}
  print OUT "$data_background\t$total_location\n";
   
