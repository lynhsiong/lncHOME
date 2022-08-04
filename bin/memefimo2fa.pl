#!/usr/bin/perl -w

$meme=shift;
$fimo=shift;
$dir=shift;

$num=2000;

open FIMO, $fimo;
while(<FIMO>){
     chomp;
     next if (/^#/);
     @Ln=split(/\t/, $_);
     push @{$fimohash{$Ln[0]}}, $Ln[-1];
}

open MEME, $meme;
while(<MEME>){
     chomp;
     #if(/^MOTIF\s+(\w+)\s+MEME/){
     if(/^MOTIF\s+(\d+)/){
         $motif=$1;
         open $motif, ">$dir/$motif.txt";
         %hash=();
     }elsif (/^letter-probability matrix:.+ w=\s+(\d+)/){
           $width=$1;
           for ($i=1;$i<=$width;$i++){ 
             $newLine=<MEME>; chomp $newLine; $newLine=~s/^\s+//; $newLine=~s/\s+$//;
             ($a, $c, $g, $t)=split(/\s+/, $newLine);
             for (1..$num){
                 $hash{$_}.=&randomBASE($a, $c, $g, $t);
             }
           }
           foreach $key (keys %hash){
             print $motif $hash{$key},"\n";
           }

           if (defined $fimohash{$motif}){
                 @temparr=@{$fimohash{$motif}};
                 $temparr=scalar @temparr; $prop=1.5; $ratio=int($prop*$num/$temparr);    # $prop=0.5 means fimo seqs weight about 1/3 in final out motif matrix
																						  ### Add20170321, Motif Matrix 2000 seqs, How many fimo lines? when prop=x, final contribution of fimos is 1-1/(1+x)=x/(x+1)! 
                 for (1..$ratio){
                     print $motif join("\n", @temparr),"\n";
                 }
           }
     close $motif;
     }
}

sub randomBASE{
     $A=shift; $C=shift; $G=shift; $T=shift;
     $rand=rand($A+$C+$G+$T);
     if (($rand>=0) && ($rand<$A)){
        return "A";
     }elsif (($rand>=$A) && ($rand<($A+$C))){
        return "C";
     }elsif (($rand>=($A+$C)) && ($rand<($A+$C+$G))){
        return "G";
     }elsif (($rand>=($A+$C+$G)) && ($rand<=$A+$C+$G+$T)){
        return "T";
     }
}


#b       SEC     3490    3507    +       22.0991 2.98e-08                TAAGGATAATGGCCTCCA
#b       SEC     3580    3597    -       22.6937 1.68e-08                TAAAGAAAACGCGGTACA
#
#
# ACGT
#  0.000000        0.200000        0.100000        0.700000
#    0.800000        0.100000        0.000000        0.100000
