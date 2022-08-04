#!/usr/bin/perl -w 

$beforeM = shift;
$afterM = shift;

%hasha=(); %hashb=();

open BM, $beforeM;
while(<BM>){
     chomp;
     if(/^MOTIF\s+(\d+)/){
         $_motif=$1;
     }elsif (/^letter-probability matrix:.+ w=\s+(\d+)/){
           $_width=$1;
           for ($j=1;$j<=$_width;$j++){ 
             $_newLine=<BM>; chomp $_newLine; $_newLine=~s/^\s+//; $_newLine=~s/\s+$//;
             ($_a, $_c, $_g, $_t)=split(/\s+/, $_newLine);
             $hashb{$_motif}{$j}=&INF($_a, $_c, $_g, $_t);      
           }                  
     }
}
close BM;

open AM, $afterM;
while(<AM>){
     chomp;
     if(/^MOTIF\s+(\d+)/){
         $motif=$1;
     }elsif (/^letter-probability matrix:.+ w=\s+(\d+)/){
           $width=$1;
           for ($i=1;$i<=$width;$i++){ 
             $newLine=<AM>; chomp $newLine; $newLine=~s/^\s+//; $newLine=~s/\s+$//;
             ($a, $c, $g, $t)=split(/\s+/, $newLine);
             $hasha{$motif}{$i}=&INF($a, $c, $g, $t);      
           }                  
     }
}
close AM;

foreach $MotifID (sort keys %hasha){
     if (defined $hashb{$MotifID}){
        %hash_t_A=%{$hasha{$MotifID}};
        %hash_t_B=%{$hashb{$MotifID}};
        $sum_A=$sum_B=0;
        @info_A=values %hash_t_A; for (@info_A){ $sum_A+=$_; }
        @info_B=values %hash_t_B; for (@info_B){ $sum_B+=$_; }
        if ( ($sum_A >= $sum_B) || ($sum_A/scalar(@info_A) >=1) ){
           print " ",$MotifID;
        }
     }
}

#print "\n";


sub INF{
     $A=shift; $C=shift; $G=shift; $T=shift;
     $Ainf=($A==0)?0:$A*log($A)/log(2); $Cinf=($C==0)?0:$C*log($C)/log(2);
     $Ginf=($G==0)?0:$G*log($G)/log(2); $Tinf=($T==0)?0:$T*log($T)/log(2);
     $H= 2 + ($Ainf+ $Cinf + $Ginf + $Tinf);
     #$H= 2 + ( $A*log($A)/log(2) + $C*log($C)/log(2) + $G*log($G)/log(2)+ $T*log($T)/log(2) );
     ### info carried by a motif is compared to a background motif without info. (0.25, 0.25, 0.25, 0.25)
     return $H;
}
