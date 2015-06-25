#!/usr/bin/perl

die "Usage : genRef.pl outF\n" unless scalar(@ARGV) == 1;

for ($i = 0; $i < 3; $i++) {
    $line = <STDIN>;
}

$size = 0;
(@names, @lens) = ();
while ($line = <STDIN>) {
    ++$size;
    chomp($line);
    ($seqn, $name, $len) = split(/[ \t]+/, $line);
    push(@names, $name);
    push(@lens, $len);
}

open(OUTPUT, ">$ARGV[0]");
print OUTPUT "$size\n";
print OUTPUT "@lens\n";
print OUTPUT "@names\n";
close(OUTPUT);


