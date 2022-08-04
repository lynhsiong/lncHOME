#!/usr/bin/perl -w

#################################################
# two kinds of usages:
#     a.  perl parsing_DP_FIMO.pl fimo_rox.sort DP_rox.sigmoid.new.out 0001 selected.fimo >./bbbbb.out
#     b.  perl parsing_DP_FIMO.pl fimo_rox.sort DP_rox.sigmoid.new.out >./bbbbb.out
#################################################

BEGIN {
  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts";
}

use Lncrnaevo::Myownperlmodule qw(eucdist intersect_two_redun_array fimolines2blocks motif_include);

use List::Compare;
use List::Util qw(min max);
use experimental qw(smartmatch);

$inFimo = shift;
$DP_results=shift;
if ($#ARGV==1){ $jobs =shift; $outfimo =shift; }     ###number 0001~9999

open FIMOrep, $inFimo;
while(<FIMOrep>){
     next if (/^#/);
     chomp;
     @F=split(/\t/,$_);
     $redundancy{$F[0]."\t".$F[4]."\t".$F[1]}=1;
     $redundancy{$F[1]}=0; 
}
close FIMOrep;

foreach $Fk (keys %redundancy){if ($redundancy{$Fk} == 0){push @deRedundancy, $Fk};}
$redundancy{$deRedundancy[0]} = $deRedundancy[1]; $redundancy{$deRedundancy[1]} = $deRedundancy[0];

open Fimo, "<$inFimo";
%_hashnew=();
while ($color = <Fimo>){
        next if ($color=~/^#/);
        chomp $color;
        @_t = split /\t/, $color;
        if (defined $redundancy{$_t[0]."\t".$_t[4]."\t".$redundancy{$_t[1]}}){
        push @{$_hashnew{$_t[1]}}, $color;
        }
}

open DP, "<$DP_results";
$num="0000";
while ($dp=<DP>){
    chomp $dp;
    next if (($dp=~/^\d/)||($dp=~/^$/));
    if ($dp=~/^>>>Score(\S+)/){
        $num++;
        $hashscore{$num}=$1;
        $align_up=<DP>; chomp $align_up; $align_down=<DP>; chomp $align_down;
        @align_up=$align_up=~/(.{1,6})/g; $up{$num}=\@{[@align_up]}; 
        #$up{$num}=[$align_up=~/(.{1,6})/g];
        @align_down=$align_down=~/(.{1,6})/g; $down{$num}=\@{[@align_down]};
    }elsif (($dp=~/:/)&&($dp!~/-:-/)){
      while ($dp=~/(\S+):(\S+)/g){
        $hashalign{$num."_".$1}=\@{[grep(!/^$/, @{[split(",", $2)]})]};
      }
    }
}


($hashFullLine, $hashb, $hasha, $hasha_prob)=fimolines2blocks(\%_hashnew, 2);
%hashFullLine=%{$hashFullLine};   ### full length of fimo lines (%t_hash)
%hashb=%{$hashb};		  ### coordinates of fimo blocks (%t__hash)
%hasha=%{$hasha};		  ### motifs contained in fimo blocks (%motif_hash)
%hasha_prob=%{$hasha_prob};       ### motif probabilities in fimo blocks (%motif_hash_prob)


foreach $key (sort {$a<=>$b} keys %hashscore){
    $Totalscore=$hashscore{$key};
    @a=@{$up{$key}}; @b=@{$down{$key}};
#print "fir=>",scalar(@a), "\tsec=>",scalar(@b),"\n";
    for (0..$#a){
      $firblock=$a[$_]; $secblock=$b[$_]; $temp_common_motif=""; $temp_common_motif_num=0;
      if ((defined $hashalign{$key."_".$firblock})&&(defined $hashalign{$key."_".$secblock})){
        ($temp_common_motif, $temp_common_motif_num)=intersect_two_redun_array($hashalign{$key."_".$firblock}, $hashalign{$key."_".$secblock}); 

        @AtempArray=@BtempArray=();
#print ">>>common motif\t=>\t", join(" ", @{$temp_common_motif}),"\n";
        for (@{$hashFullLine{$firblock}}){ 
		    $AtmpMotofId=${[split("\t", $_)]}[0]."_".${[split("\t", $_)]}[4];
#print ">>fir\t$AtmpMotofId\n";
		    if ($AtmpMotofId~~@{$temp_common_motif}){ 
#print $firblock, "=>", $AtmpMotofId, "\t=>OK!\n";
			    push @AtempArray, $_; 
				push @{$posA_s{$key."_".$Totalscore}}, ${[split("\t", $_)]}[2];
                push @{$posA_e{$key."_".$Totalscore}}, ${[split("\t", $_)]}[3]; 
			} 
		}
        for (@{$hashFullLine{$secblock}}){ 
		    $BtmpMotofId=${[split("\t", $_)]}[0]."_".${[split("\t", $_)]}[4];
#print ">>sec\t$BtmpMotofId\n";
		    if ($BtmpMotofId~~@{$temp_common_motif}){ 
#print $secblock, "=>", $BtmpMotofId, "\t=>OK!\n";
			    push @BtempArray, $_;  
                push @{$posB_s{$key."_".$Totalscore}}, ${[split("\t", $_)]}[2];
                push @{$posB_e{$key."_".$Totalscore}}, ${[split("\t", $_)]}[3]; 
			} 
		}

        push @{$hashSelectedFimoLine{$key."_".$Totalscore}}, join("\n",@AtempArray)."\n".join("\n",@BtempArray);

        #####################################
        ### IMPORTANT NOTE: the distant used in the following calculation is from $hashb (only merged according to the given "overlapping" threshold, ie. 0.8), but after filtering fimo lines in a block by common motif set, the previous coordinates shoud be renew!!!
        ####################################

		#print $firblock," <=> ",$secblock,"\n" if ($key == "0001");
        @distanta=split("_", $hashb{$firblock});$disA=$distanta[1];
        @distantb=split("_", $hashb{$secblock});$disB=$distantb[1];
        if    ($_>=1){ push @{$hashfir{$key."_".$Totalscore}}, $firblock."_".($disA-$normnalizeA);
		       push @{$hashsec{$key."_".$Totalscore}}, $secblock."_".($disB-$normnalizeB);  }
        elsif ($_==0){ $normnalizeA=$disA; push @{$hashfir{$key."_".$Totalscore}}, $firblock."_".0;
		       $normnalizeB=$disB; push @{$hashsec{$key."_".$Totalscore}}, $secblock."_".0; }

        ##### @distanta=split("_", $hashb{$a[$_]});$disA=$distanta[1];
        ##### @distantb=split("_", $hashb{$b[$_]});$disB=$distantb[1];
        ##### push @{$hashfir{$key."_".$Totalscore}}, $disA-$normnalizeA if ($_>=1);
        ##### push @{$hashsec{$key."_".$Totalscore}}, $disB-$normnalizeB if ($_>=1);

        for (@{$temp_common_motif}){
          $ele=$_; 
          if (defined($hasha_prob{$firblock."_".$ele}) && defined($hasha_prob{$secblock."_".$ele})){
            push @{$hashMotifProbA{$key."_".$Totalscore}}, map {$firblock."_".$ele."_".$_} @{$hasha_prob{$firblock."_".$ele}};
            push @{$hashMotifProbB{$key."_".$Totalscore}}, map {$secblock."_".$ele."_".$_} @{$hasha_prob{$secblock."_".$ele}};
          }
        }

#        $hashLen{$key."_".$Totalscore}++;
        push @{$hash{$key."_".$Totalscore}}, $temp_common_motif_num; 

        ##### if ($_==0){ $normnalizeA=$disA; push @{$hashfir{$key."_".$Totalscore}}, 0; } 
        ##### if ($_==0){ $normnalizeB=$disB; push @{$hashsec{$key."_".$Totalscore}}, 0; }

      }
    }
}

foreach $_key (sort keys %hash){

    $score=${[split("_",$_key)]}[1];
    @temp=@{$hash{$_key}}; $total=0; for (@temp){ $total+=$_; }
    @selected=grep{$_>=1}@temp; $max=max(@temp);
    ($A_s, $A_e, $B_s, $B_e)=( min(@{$posA_s{$_key}}), max(@{$posA_e{$_key}}), min(@{$posB_s{$_key}}), max(@{$posB_e{$_key}}) );
    @FIR_pos=@{$hashfir{$_key}};  @SEC_pos=@{$hashsec{$_key}};
    $euc_dist=eucdist(\@{[map {s/\w+_//; $_} @FIR_pos]}, \@{[map {s/\w+_//; $_} @SEC_pos]});
    ##if (($score>=2.0) || ($max>=3) || (scalar(@selected)>=3)){ if (!($hashfir{$_key}) || !($hashsec{$_key})){ die "ERR: $DP_results\t$_key\n"; } }
#print scalar(@selected)," ?= ",$hashLen{$_key},"\n";
if ((defined $hashMotifProbA{$_key})&&(defined $hashMotifProbB{$_key})&&(defined $hashfir{$_key})&&(defined $hashsec{$_key})){
    if (defined $jobs) {
      if (${[split("_",$_key)]}[0] eq $jobs){
         print $DP_results,"\t",join("\t", @{[split("_",$_key)]}),"\t",$euc_dist,"\t",$A_s,"_",$A_e,"\t",$B_s,"_",$B_e,"\t",$total,"\t",scalar(@selected),"\t",$max,"\n",join(" ",@{$hashMotifProbA{$_key}}),"\n",join(" ", @{$hashMotifProbB{$_key}}),"\n",join(" ",@{$hashfir{$_key}}),"\n",join(" ", @{$hashsec{$_key}}),"\n" if (($score>=2.0) || ($max>=3) || ($total>=3));
      }
    }else{
      print $DP_results,"\t",join("\t", @{[split("_",$_key)]}),"\t",$euc_dist,"\t",$A_s,"_",$A_e,"\t",$B_s,"_",$B_e,"\t",$total,"\t",scalar(@selected),"\t",$max,"\n",join(" ",@{$hashMotifProbA{$_key}}),"\n",join(" ", @{$hashMotifProbB{$_key}}),"\n",join(" ",@{$hashfir{$_key}}),"\n",join(" ", @{$hashsec{$_key}}),"\n" if (($score>=2.0) || ($max>=3) || ($total>=3));  ###@selected: maximun of syntenic block numbers   ; $max: maximun of merged motif number for individual blocks
    }

    if ((defined $jobs) && ($jobs eq ${[split("_",$_key)]}[0]) && (defined $hashSelectedFimoLine{$jobs."_".$score})){
       open OUT, ">$outfimo";  print OUT join("\n", @{$hashSelectedFimoLine{$jobs."_".$score}});     
       `cat $outfimo|sort -k2,2 -k3,3n -k4,4n >$outfimo.s`; `rm -rf $outfimo`;
    }
}

}




