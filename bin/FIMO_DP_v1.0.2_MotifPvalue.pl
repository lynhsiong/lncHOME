#!/usr/bin/perl -w

BEGIN {
  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts";
}

use Lncrnaevo::Myownperlmodule qw(fimolines2blocks motif_include motifpos getArrIdx);
use experimental qw(smartmatch);

use strict;
use List::Util qw(min max);
use Math::Trig;

my $infa  =shift;
my $infimo=shift;

my (%fa,$id);
open FA, $infa;
while(<FA>){
     chomp;
     if(/^>(\S+)/){
         $id=$1;
     }else{
         $fa{$id}.=$_;
     }
}

my %redundancy;
open FIMOrep, $infimo;
while(<FIMOrep>){
     next if (/^#/);
     chomp;
     my @F=split(/\t/,$_);
     $redundancy{$F[0]."\t".$F[4]."\t".$F[1]}=1;
     $redundancy{$F[1]}=0; 
}
close FIMOrep;

my @deRedundancy;
foreach my $Fk (keys %redundancy){if ($redundancy{$Fk} == 0){push @deRedundancy, $Fk};}
$redundancy{$deRedundancy[0]} = $deRedundancy[1]; $redundancy{$deRedundancy[1]} = $deRedundancy[0];

my %hash;
open FIMO, $infimo;
while(<FIMO>){
     next if (/^#/);
     chomp;
     my @Ln=split(/\t/,$_);
     if (defined $redundancy{$Ln[0]."\t".$Ln[4]."\t".$redundancy{$Ln[1]}}){
         push @{$hash{$Ln[1]}}, $_;
     }
}

my ($_hash, $__hash)=fimolines2blocks(\%hash, 0);
my %_hash=%{$_hash};
my %__hash=%{$__hash};


=head
my (%_hash,%__hash);
foreach my $key (keys %hash){
     my @temp=@{$hash{$key}}; 
     my @first=split(/\t/,$temp[0]);
     my ($sp,$st,$ed)=($first[1],$first[2],$first[3]);
     my $i="000";
     push @{$_hash{$key.$i}}, $temp[0];
     $__hash{$key.$i}=$sp."_".$st."_".$ed;
     foreach (1..$#temp){
     my @tempArray=split(/\t/,$temp[$_]);
       if (($tempArray[1] eq $sp) && ($tempArray[2] <= $ed)){
          push @{$_hash{$key.$i}}, $temp[$_];
          ($sp,$st,$ed)=($tempArray[1],$st,max($ed,$tempArray[3]));
          $__hash{$key.$i}=$sp."_".$st."_".$ed;
       }else{
          $i++;
          push @{$_hash{$key.$i}}, $temp[$_];
          $__hash{$key.$i}=$tempArray[1]."_".$tempArray[2]."_".$tempArray[3];
          ($sp,$st,$ed)=($tempArray[1],$tempArray[2],$tempArray[3]);
       }
     }
}
=cut

=head
foreach my $key (sort keys %__hash){
    print $key,"; ",$__hash{$key},"; ",join(" ",@{$_hash{$key}}),"\n"
}
=cut

my (%motif_fa,%Macro_seq);
my $motif_fa=\%motif_fa;
foreach my $key (sort keys %__hash){
     my ($keyid,$st,$ed)=split(/_/,$__hash{$key});
     if (defined($fa{$keyid})){
          $motif_fa->{$keyid}{$key}=substr($fa{$keyid},$st-1,$ed-$st+1);
     }
     $Macro_seq{$keyid}.=$key;
}

=head
foreach my $key (keys %motif_fa){
     my %thash=%{$motif_fa{$key}};
     foreach my $_key (keys %thash){
           print $_key,"\t",$thash{$_key},"\n";
     }
}
=cut

my (@sp,@motif_id_sp1,@motif_seq_sp1,@motif_id_sp2,@motif_seq_sp2, @FIR_motif_str, @SEC_motif_str, $share_motif);
@sp=(sort keys %motif_fa);
@motif_id_sp1=(sort {$a cmp $b} keys %{$motif_fa{$sp[0]}});
for (@motif_id_sp1){push @motif_seq_sp1, ${$motif_fa{$sp[0]}}{$_};}
@motif_id_sp2=(sort {$a cmp $b} keys %{$motif_fa{$sp[1]}});
for (@motif_id_sp2){push @motif_seq_sp2, ${$motif_fa{$sp[1]}}{$_};}
#print join(" ", @motif_id_sp1), "\n",join(" ", @motif_seq_sp1), "\n",join(" ", @motif_id_sp2), "\n",join(" ", @motif_seq_sp2), "\n";

#shared motif number between two arrays
for (@motif_id_sp1){ for(@{$_hash{$_}}){ push @FIR_motif_str, ${[split("\t", $_)]}[0]."_".${[split("\t", $_)]}[4]; } }
for (@motif_id_sp2){ for(@{$_hash{$_}}){ push @SEC_motif_str, ${[split("\t", $_)]}[0]."_".${[split("\t", $_)]}[4]; } }
$share_motif=intersect_two_redun_array(\@FIR_motif_str,\@SEC_motif_str);

my %Macro_matrix;
my $Macro_matrix=\%Macro_matrix;

if ($share_motif >=3){

#Macro_DP
my %_macro_matrix;
my $_macro_matrix=\%_macro_matrix;
for my $i (0..$#motif_id_sp1){
  for my $j (0..$#motif_id_sp2){
##  $_macro_matrix->{$motif_id_sp1[$i]}{$motif_id_sp2[$j]}=dp($motif_seq_sp1[$i],$motif_seq_sp2[$j]);
  $_macro_matrix->{$motif_id_sp1[$i]}{$motif_id_sp2[$j]}=dp($motif_id_sp1[$i],$motif_id_sp2[$j]);
#  print $motif_id_sp1[$i],"_",$motif_id_sp2[$j],"\n";
#  print $motif_seq_sp1[$i],"_",$motif_seq_sp2[$j],"\n";
  }
}

=head
my @axis;
foreach my $key1 (sort {$a cmp $b} keys %_macro_matrix){
  my %temp=%{$_macro_matrix{$key1}};
  print $key1,"\t";
  foreach my $key2 (sort {$a cmp $b} keys %temp){
    push @axis, $key2 if (!($key2~~@axis)); 
    print $temp{$key2},"\t";
  }
  print "\n";
}
print "\t",join("\t",@axis),"\n\n";
=cut

#Normalization
my ($_M, $_m)=(6, -6);   
my $slop=1; ######### bigger $slop, more step-like; when $slop=1, sigmoid min/max key -> value: -6 -> 0 and 6->1  ###########
my @t_array;
for (values %_macro_matrix){ push @t_array, (values %{$_});}    #print join(" ", @{[(values %{$_})]})
my ($_max, $_min)=(max(@t_array), min(@t_array));
($_max, $_min)=(10, -10);    #renew to constant values
##if ($_max!=$_min){
foreach my $k1 (sort {$a cmp $b} keys %_macro_matrix){
  my %temphash=%{$_macro_matrix{$k1}};
  foreach my $k2 (sort {$a cmp $b} keys %temphash){
    #$Macro_matrix->{$k1}{$k2}= sprintf "%.3f", tanh(($temphash{$k2}-$_min)/($_max-$_min)*($_M-$_m)+$_m);        ##### tanh transformation
    #$Macro_matrix->{$k1}{$k2}= sprintf "%.3f", 1/(1+exp(-$slop*(($temphash{$k2}-$_min)/($_max-$_min)*($_M-$_m)+$_m)))+$_m/($_M-$_m);   ##### sigmoid transformation
    $Macro_matrix->{$k1}{$k2}= sprintf "%.3f", 2*1/(1+exp(-$slop*((($temphash{$k2}-4)-$_min)/($_max-$_min)*($_M-$_m)+$_m)))-1;       ##### sigmoid transformation
  }
}
##}else{
##foreach my $k1 (sort {$a cmp $b} keys %_macro_matrix){
##  my %temphash=%{$_macro_matrix{$k1}};
##  foreach my $k2 (sort {$a cmp $b} keys %temphash){
##    $Macro_matrix->{$k1}{$k2}= sprintf "%.3f", $temphash{$k2};
##  }
##}
##}

=head
my @_axis;
foreach my $key1 (sort {$a cmp $b} keys %Macro_matrix){
  my %temp=%{$Macro_matrix{$key1}};
  print $key1,"\t";
  foreach my $key2 (sort {$a cmp $b} keys %temp){
    push @_axis, $key2 if (!($key2~~@_axis));
    print $temp{$key2},"\t";
  }
  print "\n";
}
print "\t",join("\t",@_axis),"\n\n";
=cut

#Micro_DP
#print $Macro_seq{$sp[0]},"\n",$Macro_seq{$sp[1]},"\n";
my ($out_result, @out_D2array) = DP($Macro_seq{$sp[0]},$Macro_seq{$sp[1]});

=head
print "\n\n";
for my $i ( 0 .. $#out_D2array ) { 
    for my $j ( 0 .. $#{$out_D2array[$i]} ) { 
        my $temp_score = sprintf "%.3f", $out_D2array[$i][$j]{score};
        print $temp_score,"_",$out_D2array[$i][$j]{pointer},"\t"; 
    } 
    print "\n";
} 
=cut

if ($out_result){
print $out_result,"\n";
}

}

print "\n\nMicroH_Finished";





sub intersect_two_redun_array{
    my $A=shift; my $B=shift; my $_common=0;
    my (%counta, %countb);
    for my $_a_ele (@{$A}){$counta{$_a_ele}++;}
    for my $_b_ele (@{$B}){$countb{$_b_ele}++;}
    foreach my $_ele_ (keys %counta){
        if(defined($countb{$_ele_})){ $_common+=min($counta{$_ele_},$countb{$_ele_}); }
    }
    return $_common;
}


sub dp{

###################   5	FIR	41	47	+	13.1143	4.39e-05		GCTTGGC     ###############
my ($match,$mismatch)=(1,-1);
my (@tempA1,@tempA2,@singleMotifA1,@singleMotifA2);
my $seq1id=shift; @tempA1=@{$_hash{$seq1id}};
my $seq2id=shift; @tempA2=@{$_hash{$seq2id}};

=head
my (@A1Motif,@A2Motif,%counta,%countb);
for (0..$#tempA1){ @singleMotifA1=split("\t", $tempA1[$_]); push @A1Motif, $singleMotifA1[0]."_".$singleMotifA1[4]; }
for (0..$#tempA2){ @singleMotifA2=split("\t", $tempA2[$_]); push @A2Motif, $singleMotifA2[0]."_".$singleMotifA2[4]; }
my $common=0;
        for my $a_ele (@A1Motif){$counta{$a_ele}++;}
        for my $b_ele (@A2Motif){$countb{$b_ele}++;}
        foreach my $ele (keys %counta){
            if(defined($countb{$ele})){$common+=min($counta{$ele},$countb{$ele});}
        }
my $score= ($common>0)?(($common)*$match):$mismatch;
return $score;
undef %counta; undef %countb;
=cut

#=head
my (%A1Motif,%A2Motif,%A1Prob,%A2Prob,%count,%A1count,%A2count,);
for (0..$#tempA1){ @singleMotifA1=split("\t", $tempA1[$_]); $A1Motif{$_}=$singleMotifA1[0]."_".$singleMotifA1[4]; $A1Prob{$_}=log($singleMotifA1[6])/log(10); }
for (0..$#tempA2){ @singleMotifA2=split("\t", $tempA2[$_]); $A2Motif{$_}=$singleMotifA2[0]."_".$singleMotifA2[4]; $A2Prob{$_}=log($singleMotifA2[6])/log(10); }

my ($SizeScale,$GapScale)=(1,0.5);   ### how to combine the Prob size (e.g. e-5) and Prob gap (or distance of two Prob values)  ###

foreach my $A1ranknum (sort keys %A1Motif){
     foreach my $A2ranknum (sort keys %A2Motif){
        if ($A1Motif{$A1ranknum} eq $A2Motif{$A2ranknum}){
            my $ProbGap=abs($A1Prob{$A1ranknum}-$A2Prob{$A2ranknum}); 
            my $ProbSize=($A1Prob{$A1ranknum}>$A2Prob{$A2ranknum})?abs($A2Prob{$A2ranknum}):abs($A1Prob{$A1ranknum});
            my $SizeMinusGapScore=$SizeScale*$ProbSize-$GapScale*$ProbGap;
            if (defined($count{$A1Motif{$A1ranknum}})){ if ($SizeMinusGapScore < $count{$A1Motif{$A1ranknum}}){
                                                           $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
            }else{ $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
        } else {
            my $A1ProbSize=abs($A1Prob{$A1ranknum});
            if (defined($A1count{"A1MaxMotif"})){ if ($A1ProbSize > $A1count{"A1MaxMotif"}){ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            }else{ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            my $A2ProbSize=abs($A2Prob{$A2ranknum});
            if (defined($A2count{"A2MaxMotif"})){ if ($A2ProbSize > $A2count{"A2MaxMotif"}){ $A2count{"A2MaxMotif"}=$A2ProbSize; }
            }else{ $A2count{"A2MaxMotif"}=$A2ProbSize; }
        }
     }
}

my $score;
if (!%count){
    if ((%A1count) && (%A2count)){
       #$score=$SizeScale*min($A1count{"A1MaxMotif"},$A2count{"A2MaxMotif"})-$GapScale*abs($A1count{"A1MaxMotif"}-$A2count{"A2MaxMotif"});
       $score=$SizeScale*max($A1count{"A1MaxMotif"},$A2count{"A2MaxMotif"});
       $score=$mismatch*$score;
    }
}else {
    foreach my $key (keys %count){
       $score+=$count{$key};
    }
    $score=$match*$score;
}

return sprintf "%.3f",$score;
undef %count; undef %A1count; undef %A2count; 
#=cut

}

=head
sub dp{

my $seq1=shift;
my $seq2_org=shift;
my $seq2_rev=reverse $seq2_org;
   $seq2_rev=~tr/ACGTacgtUu/TGCAtgcaAa/;
my ($MATCH,$MISMATCH,$GAP) = (2,-3,-2);

my @max_score;
for ($seq2_org,$seq2_rev){
my $seq2=$_;
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = 0;
    $matrix[0][$j]{pointer} = "none";
}
for(my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
}

for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);

        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);       
        if (uc($letter1) eq uc($letter2)) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }

        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;
        
        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next;
        }

        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        } else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        }

        push @max_score, $matrix[$i][$j]{score};
    }
}
}
 return max(@max_score);
}
=cut


sub DP{

my $seq1=shift;
my $seq2=shift;
my $GAP = 0;
#my $GAP = 2*1/(1+exp(-$slop*$_m))-1;
my $MATCH;

my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1)/6; $j++) {
    $matrix[0][$j]{score}   = 0;
    $matrix[0][$j]{pointer} = "none";
}
for(my $i = 1; $i <= length($seq2)/6; $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
}

my %cmp_hash;
for(my $i = 1; $i <= length($seq2)/6; $i++) {
    for(my $j = 1; $j <= length($seq1)/6; $j++) {
        my ($diagonal_score, $left_score, $up_score);

        my $letter1 = substr($seq1, ($j-1)*6, 6);
        my $letter2 = substr($seq2, ($i-1)*6, 6);
           $MATCH=$Macro_matrix->{$letter1}{$letter2};
           $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;

        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;

        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next;
        }

        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                if ($diagonal_score >0){
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        } else {
            if ($up_score >= $left_score) {
                if ($up_score >0){
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        }

        $cmp_hash{$i."_".$j} = $matrix[$i][$j]{score};
=head        
        if ($matrix[$i][$j]{score} > $max_score) {
            $max_i     = $i;
            $max_j     = $j;
            $max_score = $matrix[$i][$j]{score};
        }
=cut

    }
}

# trace-back
my $max_i     = 0;
my $max_j     = 0;
#my $max_score = 0;


my @pool;
my @align1;
my @align2;
my $result;
my $x;
my $y;

foreach my $_key_ (sort { $cmp_hash{$b} <=> $cmp_hash{$a} } keys %cmp_hash){

unless ($_key_~~ @pool){
($max_i, $max_j)=split("_", $_key_);
if($matrix[$max_i][$max_j]{pointer} eq "diagonal"){
my $j = $max_j;
my $i = $max_i;

while (1) {
    last if $matrix[$i][$j]{pointer} eq "none";
    push @pool, $i."_".$j;
    if ($matrix[$i][$j]{pointer} eq "diagonal") {
        push @align1, substr($seq1, ($j-1)*6, 6);
        push @align2, substr($seq2, ($i-1)*6, 6);
        $i--; $j--;
    }elsif ($matrix[$i][$j]{pointer} eq "left") {
        push @align1, substr($seq1, ($j-1)*6, 6);
        push @align2, "------";
        $j--;
    }elsif ($matrix[$i][$j]{pointer} eq "up") {
        push @align1, "------";
        push @align2, substr($seq2, ($i-1)*6, 6);
        $i--;
    }   
}

@align1 = reverse @align1;
@align2 = reverse @align2;

$result .= ">>>Score".$cmp_hash{$_key_}."\n".join("",@align1)."\n".join("",@align2)."\n\n";

my @motifID;
my $sys;
foreach (0..$#align1){
    $sys=$_;

    $result.=$align1[$sys].":";
    if (defined ($_hash{$align1[$sys]})){
    foreach (@{$_hash{$align1[$sys]}}){
        @motifID=split("\t",$_);
        $result.=$motifID[0]."_".$motifID[4].",";
        $x++;
    }
    }else{
        $result.="------";
    }

    $result.="\t\t\t";

    $result.=$align2[$sys].":";
    if (defined ($_hash{$align2[$sys]})){
    foreach (@{$_hash{$align2[$sys]}}){
        @motifID=split("\t",$_);
        $result.=$motifID[0]."_".$motifID[4].",";
        $y++;
    }
    }else{
        $result.="------";
    }

$result.="\n";

}

$result.="\n";

undef @align1;
undef @align2;
}
}

}

$result=(($x>=3)&&($y>=3))?$result:"";
return ($result, @matrix);

}


=head
sub fimolines2blocks{
  my (%t_hash,%t__hash);
  my $inhash=shift; my $thresholdratio=shift;
  my %inhash=%{$inhash};

  foreach my $key (keys %inhash){
     my @temp=@{$inhash{$key}};
     my @first=split(/\t/,$temp[0]);
     my ($sp,$st,$ed)=($first[1],$first[2],$first[3]);
     my $i="000";
     push @{$t_hash{$key.$i}}, $temp[0];
     $t__hash{$key.$i}=$sp."_".$st."_".$ed;
     foreach (1..$#temp){
     my @tempArray=split(/\t/,$temp[$_]);
       if (($tempArray[1] eq $sp) && ($tempArray[2] <= $ed)){
          if (motif_include($thresholdratio, \@tempArray, \@{$t_hash{$tempArray[1].$i}})){
             push @{$t_hash{$key.$i}}, $temp[$_];
             ($sp,$st,$ed)=($tempArray[1],min($st,$tempArray[2]),max($ed,$tempArray[3]));
             $t__hash{$key.$i}=$sp."_".$st."_".$ed;
          }else{
             $i++;
             push @{$t_hash{$key.$i}}, $temp[$_];
             $t__hash{$key.$i}=$tempArray[1]."_".$tempArray[2]."_".$tempArray[3];
             ($sp,$st,$ed)=($tempArray[1],$tempArray[2],$tempArray[3]);
          }
       }else{
          $i++;
          push @{$t_hash{$key.$i}}, $temp[$_];
          $t__hash{$key.$i}=$tempArray[1]."_".$tempArray[2]."_".$tempArray[3];
          ($sp,$st,$ed)=($tempArray[1],$tempArray[2],$tempArray[3]);
       }
     }
   }

  return (\%t_hash, \%t__hash);
}



sub motif_include{
  my $_ratio=shift;
  my $cand_arr=shift; my @cand_arr=@{$cand_arr};
  my ($cand_s, $cand_e)=($cand_arr[2], $cand_arr[3]);
  my $cand_len=$cand_e-$cand_s+1;
  my $target_arrlib=shift; my @target_arrlib=@{$target_arrlib};
  my ($overlaplen,$resultSignal); 

  for (@target_arrlib){
    my @temp=split(/\t/, $_);
    my ($target_s, $target_e)=($temp[2], $temp[3]);
    my $target_len=$target_e-$target_s+1;

#    $overlaplen=$target_e-$cand_s+1;
    if ($target_s<=$cand_s){
        $overlaplen=($target_e<=$cand_e)?$target_e-$cand_s+1:$cand_e-$cand_s+1;
    }else{
        $overlaplen=($target_e<=$cand_e)?$target_e-$target_s+1:$cand_e-$target_s+1;
    }

    if ( ($overlaplen/$cand_len>=$_ratio)||($overlaplen/$target_len>=$_ratio) ){
       $resultSignal.="1";
    }else{
       $resultSignal.="0";
    }
  }

  return ($resultSignal=~/0/)?0:1;
}
=cut

