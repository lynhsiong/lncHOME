#!/usr/bin/perl -w
use strict;

#my ($ori_fas) = $ARGV[0];
#my ($con_anc) = $ARGV[0];
my ($con_fas) = $ARGV[0];
#my ($out_inf) = $ARGV[1];
my $usage = "This script is to get the information of anchor motif.
usage: $0 <re_fas> <con_str> <out_inf>
";
die $usage if $#ARGV<1;

#open(ORI, $ori_fas)||die("open $ori_fas error!\n");
#open(ANC, $con_anc)||die("open $con_anc error!\n");
open(FAS, $con_fas)||die("open $con_fas error!\n");
open(OUT, ">C1_outmotif");
my @ancfile = @ARGV[1..$#ARGV];

my $i = 0;
my $j = 0;
my $k = 0;
my $r = 0;
my $t1 = 0;
my $t2 = 0;
my $mlen = 60;
my $weight = 0.0;
my $seqnum = 0;
my $flag = 0;
my $sen = "";
my $fas_file = "";
my $out_file1 = "";
my $out_file2 = "";
my $con_anc = "";
my $comm = "";
my $sta = "";
my $ke = "";
my $seq = "";
my @sent = "";
my @elemen = ();
my @info = ();
my @seq = ();
my @rest = ();
my @anchor1 = ();
my @anchor2 = ();
my @anchor = ();
my @lim = ();
my @match = ();
my @turn = ();
my @tem = ();
my @tinfo = ();
my @cord = ();
my %inf = ();
my @RAWcolor = ("red","yellow","blue","lightgreen","green","lightskyblue","brown","lightpink","purple","orange","black","navajowhite","orchid","rosybrown","springgreen","tomato","wheat","sienna","violet","darkgray","cyan","darkred","darkseagreen","darkblue","aquamarine","lime","PeachPuff","Peru","Linen","olive","gold","DarkKhaki","Ivory","greenyellow","beige","lime","honeydew","teal","azure","indigo","plum","navy","pink","crimson","Turquoise","MintCream","LimeGreen","OliveDrab","LemonChiffon","PapayaWhip","tan","DarkOrange","Chocolate","Coral","FireBrick","DimGray","Gainsboro","MistyRose","Tomato","SaddleBrown","Moccasin","Cornislk");


my $motif1 = "";
my $pval1 = 0;
my $motif2 = "";
my $pval2 = 0;
my $item;
my $sna = "";
my $sna1 = "";
my $sna2 = "";
my @lastre = (0)x6;
my @reg1 = ();
my @reg2 = ();
my @reg3 = ();

while($sen=<FAS>){
	chomp($sen);
	if(($sen =~m/^>/) && ($seq ne "")){
		@sent = split(//,$seq);
		@anchor = (-1)x($#sent+1);
		@anchor1 = (0)x($#sent+1);
		@anchor2 = ([])x($#sent+1);
		$j = 0;
		for($i=0; $i<=$#sent; $i++){
			if($sent[$i] ne "-"){
				$anchor[$j] = $i;
				$j++;
			}
		}
		push(@turn, [@anchor]);
		push(@match, [@anchor1]);
		push(@tinfo, [@anchor2]);
		$seq = "";
		$k++;
	}elsif($sen =~m/^>/ || $sen eq " "){
	}else{
		$seq = $seq.$sen;
	}
}
@sent = split(//,$seq);
@anchor = (-1)x($#sent+1);
@anchor1 = (0)x($#sent+1);
@anchor2 = ([])x($#sent+1);
$j = 0;
for($i=0; $i<=$#sent; $i++){
	if($sent[$i] ne "-"){
		$anchor[$j] = $i;
		$j++;
	}
}
push(@turn, [@anchor]);
push(@match, [@anchor1]);
push(@tinfo, [@anchor2]);

#1	2	180	202	15	10.04	1671_181
#1009	FIR	1396	1415	+	9.9292	5.79e-05		TTATTGTTTCATTCTTGTGA
#C1_1396_2789_C2_2166_3104.txt
$con_anc = $ancfile[0];
@cord = split(/\//,$con_anc);
@cord = split(/\./,$cord[$#cord]);
@cord = split(/_/,$cord[0]);
open(ANC, $con_anc)||die("open $con_anc error!\n");
while($sen=<ANC>){
	chomp($sen);
	@lim = split(/\t/,$sen);
	#@tem = split(/_/,$lim[6]);
	#$lim[6] = $tem[0];
	if($lim[4] eq "-"){
		next;
	}
	if($lim[1] eq "FIR"){
		$j = ${$turn[0]}[$lim[2]-$cord[1]];
		$r = $lim[2]-$cord[1];
		$sna = "0_".$lim[0]."_".$j;
		$inf{$sna} = [$lim[0],$lim[3]-$lim[2]+1,$lim[2]-$cord[1],$lim[3]-$cord[1],$lim[6],$lim[8]];
		push(@{${$tinfo[0]}[$j]},$lim[0]);
		if(${$match[0]}[$j] eq 0){
			${$match[0]}[$j] = $lim[0];
		}
		else{
			${$match[0]}[$j] = ${$match[0]}[$j]."_".$lim[0];
		}
	}else{
		$j = ${$turn[1]}[$lim[2]-$cord[4]];
		$r = $lim[2]-$cord[4];
		$sna = "1_".$lim[0]."_".$j;
		$inf{$sna} = [$lim[0],$lim[3]-$lim[2]+1,$lim[2]-$cord[4],$lim[3]-$cord[4],$lim[6],$lim[8]];
		push(@{${$tinfo[1]}[$j]},$lim[0]);
		if(${$match[1]}[$j] eq 0){
			${$match[1]}[$j] = $lim[0];
		}
		else{
			${$match[1]}[$j] = ${$match[1]}[$j]."_".$lim[0];
		}
	}

}
close ANC;
for ($i=1; $i<=$#ancfile; $i++){
	$con_anc = $ancfile[$i];
	@cord = split(/\//,$con_anc);
	@cord = split(/\./,$cord[$#cord]);
	@cord = split(/_/,$cord[0]);
	open(ANC, $con_anc)||die("open $con_anc error!\n");
	while($sen=<ANC>){
		chomp($sen);
		@lim = split(/\t/,$sen);
		#@tem = split(/_/,$lim[6]);
		#$lim[6] = $tem[0];
		if(($lim[4] eq "-") || ($lim[1] eq "FIR")){
			next;
		}else{
			$j = ${$turn[$i+1]}[$lim[2]-$cord[4]];
			$r = $lim[2]-$cord[4];
			$k = $i+1;
			$sna = $k."_".$lim[0]."_".$j;
			$inf{$sna} = [$lim[0],$lim[3]-$lim[2]+1,$lim[2]-$cord[4],$lim[3]-$cord[4],$lim[6],$lim[8]];
			push(@{${$tinfo[$i+1]}[$j]},$lim[0]);
			if(${$match[$i+1]}[$j] eq 0){
				${$match[$i+1]}[$j] = $lim[0];
			}
			else{
				${$match[$i+1]}[$j] = ${$match[$i+1]}[$j]."_".$lim[0];
			}
		}
		#$j = ${$turn[$lim[0]-1]}[$lim[2]-1];
		#${$match[$lim[0]-1]}[$j] = ${$match[$lim[0]-1]}[$j]."_".$lim[6];
		#$j = ${$turn[$lim[1]-1]}[$lim[3]-1];
		#${$match[$lim[1]-1]}[$j] = ${$match[$lim[1]-1]}[$j]."_".$lim[6];
	}
	close ANC;
}




for($i=0; $i<=$#{$match[0]}; $i++){
	if((${$match[0]}[$i] eq 0) && (${$match[1]}[$i] eq 0) &&(${$match[2]}[$i] eq 0)){
		next;
	}
	print OUT $i+1,"\t";
	for($j=0; $j<=$#match; $j++){
		print OUT ${$match[$j]}[$i],"\t";
	}
	print OUT "\n";
}
close FAS;
close OUT;

open(OUT1, ">C1_pmotif.txt");
my @lastre1 = (-1,-1,-1,-1,-1,-1);
my @lastre2 = (-1,-1,-1,-1,-1,-1);
@reg1 = ();
@reg2 = ();
@reg3 = ();
$flag = 0;
my $col = 0;

for($i=0; $i<=$#{$match[0]}; $i++){
	$flag = 0;
	if((${$match[0]}[$i] eq 0) && (${$match[1]}[$i] eq 0) &&(${$match[2]}[$i] eq 0)){
		next;
	}
	if(${$match[0]}[$i] ne 0){
		@tem = split(/_/,${$match[0]}[$i]);
		$motif1 = $tem[0];
		$sna = "0_".$motif1."_".$i;
		$pval1 = ${$inf{$sna}}[4];
		for($j=1; $j<=$#tem; $j++){
			$motif2 = $tem[$j];
			$sna = "0_".$motif2."_".$i;
			if(${$inf{$sna}}[4] < $pval1){
				$motif1 = $motif2;
				$pval1 = ${$inf{$sna}}[4];
			}
		}
		$sna1 = "0_".$motif1."_".$i;
		#if(grep {$_ eq $motif1} @{${$tinfo[1]}[$i]}){
		#	$sna2 = "1_".$motif1."_".$i;
		#	print OUT1 "1\t2\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
		#}
		for($r=0; $r<=$mlen; $r++){
			$k = $i-$r;
			@tem = split(/_/,${$match[1]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "1_".$motif1."_".$k;
				print OUT1 "1\t2\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				#$flag = 1;
				#@lastre2 = (${$inf{$sna1}}[2]+1,${$inf{$sna1}}[3]+1,${$inf{$sna2}}[2]+1,${$inf{$sna2}}[3]+1,0,0);
				if(($lastre1[1] <= ${$inf{$sna1}}[2]+1) || ($lastre1[3] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = ${$inf{$sna1}}[2]+1;
					$lastre1[1] = ${$inf{$sna1}}[3]+1;
					$lastre1[2] = ${$inf{$sna2}}[2]+1;
					$lastre1[3] = ${$inf{$sna2}}[3]+1;
				}else{
					$flag = 1;
					@lastre2 = (${$inf{$sna1}}[2]+1,${$inf{$sna1}}[3]+1,${$inf{$sna2}}[2]+1,${$inf{$sna2}}[3]+1,0,0);
				}
				last;
			}
			$k = $i+$r;
			if($k > $#{$match[0]}){
				last;
			}
			@tem = split(/_/,${$match[1]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "1_".$motif1."_".$k;
				print OUT1 "1\t2\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				#$flag = 1;
				#@lastre2 = (${$inf{$sna1}}[2]+1,${$inf{$sna1}}[3]+1,${$inf{$sna2}}[2]+1,${$inf{$sna2}}[3]+1,0,0);
				if(($lastre1[1] <= ${$inf{$sna1}}[2]+1) || ($lastre1[3] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					#@reg1 = (@reg1, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg2 = (@reg2, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = ${$inf{$sna1}}[2]+1;
					$lastre1[1] = ${$inf{$sna1}}[3]+1;
					$lastre1[2] = ${$inf{$sna2}}[2]+1;
					$lastre1[3] = ${$inf{$sna2}}[3]+1;
				}else{
					$flag = 1;
					@lastre2 = (${$inf{$sna1}}[2]+1,${$inf{$sna1}}[3]+1,${$inf{$sna2}}[2]+1,${$inf{$sna2}}[3]+1,0,0);
				}
				last;
			}
		}
		for($r=0; $r<=$mlen; $r++){
			$k = $i-$r;
			@tem = split(/_/,${$match[2]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "2_".$motif1."_".$k;
				print OUT1 "1\t3\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				if($flag == 2){
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$flag = 3;
				}elsif($flag == 1){
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$t1 = $lastre2[0];
					$t2 = $lastre2[1];
					#@reg1 = (@reg1, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = $lastre2[0];
					$lastre1[1] = $lastre2[1];
					$t1 = $lastre2[2];
					$t2 = $lastre2[3];
					#@reg2 = (@reg2, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$lastre1[2] = $lastre2[2];
					$lastre1[3] = $lastre2[3];
					$flag = 3;
				}elsif(($lastre1[1] <= ${$inf{$sna1}}[2]+1) || ($lastre1[5] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					#@reg1 = (@reg1, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = ${$inf{$sna1}}[2]+1;
					$lastre1[1] = ${$inf{$sna1}}[3]+1;
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$flag = 3;
				}
				#if($flag >= 2){
				#	if($col + 1 <= $#RAWcolor){
				#		$col = $col + 1;
				#	}
				#}
				last;
			}
			$k = $i+$r;
			if($k > $#{$match[0]}){
				last;
			}
			@tem = split(/_/,${$match[2]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "2_".$motif1."_".$k;
				print OUT1 "1\t3\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				if($flag == 2){
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$flag = 3;
				}elsif($flag == 1){
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$t1 = $lastre2[0];
					$t2 = $lastre2[1];
					#@reg1 = (@reg1, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = $lastre2[0];
					$lastre1[1] = $lastre2[1];
					$t1 = $lastre2[2];
					$t2 = $lastre2[3];
					#@reg2 = (@reg2, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$lastre1[2] = $lastre2[2];
					$lastre1[3] = $lastre2[3];
					$flag = 3;
				}elsif(($lastre1[1] <= ${$inf{$sna1}}[2]+1) || ($lastre1[5] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					#@reg1 = (@reg1, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg1, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[0] = ${$inf{$sna1}}[2]+1;
					$lastre1[1] = ${$inf{$sna1}}[3]+1;
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					$flag = 3;
				}
				#if($flag >= 2){
				#	if($col + 1 <= $#RAWcolor){
				#		$col = $col + 1;
				#	}
				#}
				last;
			}
		}
		if($flag >= 2){
			if($col + 1 <= $#RAWcolor){
				$col = $col + 1;
			}
		}
	}
	if(${$match[1]}[$i] ne 0){
		@tem = split(/_/,${$match[1]}[$i]);
		$motif1 = $tem[0];
		$sna = "1_".$motif1."_".$i;
		$pval1 = ${$inf{$sna}}[4];
		for($j=1; $j<=$#tem; $j++){
			$motif2 = $tem[$j];
			$sna = "1_".$motif2."_".$i;
			if(${$inf{$sna}}[4] < $pval1){
				$motif1 = $motif2;
				$pval1 = ${$inf{$sna}}[4];
			}
		}
		$sna1 = "1_".$motif1."_".$i;
		#if(grep {$_ eq $motif1} @{${$tinfo[1]}[$i]}){
		#	$sna2 = "1_".$motif1."_".$i;
		#	print OUT1 "1\t2\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
		#}
		for($r=0; $r<=$mlen; $r++){
			#@tem = split(/_/,${$match[1]}[$i+$r]);
			#@tem = split(/_/,${$match[2]}[$i]);
			$k = $i-$r;
			@tem = split(/_/,${$match[2]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "2_".$motif1."_".$k;
				print OUT1 "2\t3\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				if(($lastre1[3] <= ${$inf{$sna1}}[2]+1) || ($lastre1[5] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					#@reg2 = (@reg2, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[2] = ${$inf{$sna1}}[2]+1;
					$lastre1[3] = ${$inf{$sna1}}[3]+1;
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					#$flag = 3;
					if($col + 1 <= $#RAWcolor){
						$col = $col + 1;
					}
				}
				last;
			}
			$k = $i+$r;
			if($k > $#{$match[0]}){
				last;
			}
			@tem = split(/_/,${$match[2]}[$k]);
			if(grep {$_ eq $motif1} @tem){
				$sna2 = "2_".$motif1."_".$k;
				print OUT1 "2\t3\t",$sna1,"\t",$sna2,"\t",${$inf{$sna1}}[2],"\t",${$inf{$sna1}}[3],"\t",${$inf{$sna1}}[4],"\t",${$inf{$sna1}}[5],"\t",${$inf{$sna2}}[2],"\t",${$inf{$sna2}}[3],"\t",${$inf{$sna2}}[4],"\t",${$inf{$sna2}}[5],"\n";
				if(($lastre1[3] <= ${$inf{$sna1}}[2]+1) || ($lastre1[5] <= ${$inf{$sna2}}[2]+1)){
					$t1 = ${$inf{$sna1}}[2]+1;
					$t2 = ${$inf{$sna1}}[3]+1;
					$flag = 2;
					#@reg2 = (@reg2, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg2, [$t1."-".$t2.":",$motif1]);
					$t1 = ${$inf{$sna2}}[2]+1;
					$t2 = ${$inf{$sna2}}[3]+1;
					#@reg3 = (@reg3, $t1."-".$t2.":".$RAWcolor[$col]);
					push(@reg3, [$t1."-".$t2.":",$motif1]);
					$lastre1[2] = ${$inf{$sna1}}[2]+1;
					$lastre1[3] = ${$inf{$sna1}}[3]+1;
					$lastre1[4] = ${$inf{$sna2}}[2]+1;
					$lastre1[5] = ${$inf{$sna2}}[3]+1;
					#$flag = 3;
					if($col + 1 <= $#RAWcolor){
						$col = $col + 1;
					}
				}
				last;
			}
		}
	}
}
close OUT1;
my @reg = (@reg1,@reg2,@reg3);
my @uni = ();
my %unio = ();
for($i=0; $i<=$#reg; $i++){
	push(@uni, ${$reg[$i]}[1]);
}
$i = 0;
$j = 0;
for($i=0; $i<=$#uni; $i++){
	$sna = $uni[$i];
	if(exists $unio{$sna}){
		
	}else{
		$unio{$sna}=$RAWcolor[$j];
		if($j + 1 <= $#RAWcolor){
			$j = $j + 1;
		}
	}
}


open(OUT2, ">C1_colmotif.txt");
print OUT2 ">C1_1\n";
for($i=0; $i<=$#reg1; $i++){
	#print OUT2 $reg1[$i]." ";
	print OUT2 ${$reg1[$i]}[0].$unio{${$reg1[$i]}[1]}." ";
}
print OUT2 "\n";
print OUT2 ">C1_2\n";
for($i=0; $i<=$#reg2; $i++){
	#print OUT2 $reg2[$i]." ";
	print OUT2 ${$reg2[$i]}[0].$unio{${$reg2[$i]}[1]}." ";
}
print OUT2 "\n";
print OUT2 ">C1_3\n";
for($i=0; $i<=$#reg3; $i++){
	#print OUT2 $reg3[$i]." ";
	print OUT2 ${$reg3[$i]}[0].$unio{${$reg3[$i]}[1]}." ";
}
print OUT2 "\n";
print OUT2 ">C1_1c\n";
for($i=0; $i<=$#reg1; $i++){
	#print OUT2 $reg1[$i]." ";
	print OUT2 ${$reg1[$i]}[0].$unio{${$reg1[$i]}[1]}." ";
}
print OUT2 "\n";
print OUT2 ">C1_2c\n";
for($i=0; $i<=$#reg2; $i++){
	#print OUT2 $reg2[$i]." ";
	print OUT2 ${$reg2[$i]}[0].$unio{${$reg2[$i]}[1]}." ";
}
print OUT2 "\n";
print OUT2 ">C1_3c\n";
for($i=0; $i<=$#reg3; $i++){
	#print OUT2 $reg3[$i]." ";
	print OUT2 ${$reg3[$i]}[0].$unio{${$reg3[$i]}[1]}." ";
}
print OUT2 "\n";
close OUT2;

exit;
