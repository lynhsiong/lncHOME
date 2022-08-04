#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;
use experimental qw(smartmatch);
use Parallel::ForkManager;

#/Share/home/zhangqf2/xiongtl/misc/disk-ioz/TErelated/rmblast/rmblast/bin/dustmasker -in test1.fa -out test1.out -window 30 -level 10 -outfmt fasta
my ($human_motif_db) = $ARGV[0];
my ($mouse_motif_db) = $ARGV[1];
my ($proc) = $ARGV[2];

my $usage = "This script is to build the motif database for eight species.
usage: $0 <human_motif_db> <mouse_motif_db> <proc>
";
die $usage if $#ARGV<2;


my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
my ($hum_motif, $mou_motif) = &getmotif($human_motif_db, $mouse_motif_db);
my %hum_mot = %{$hum_motif}; my %mou_mot = %{$mou_motif}; my %com_mot = ();

foreach my $key (keys %hum_mot){
	$com_mot{$key} = 1;
}
foreach my $key (keys %mou_mot){
	if(exists $com_mot{$key}){
		$com_mot{$key} = 2;
	}else{
		$com_mot{$key} = 0;
	}
}

my $hspe = &hum2spe(); my $mspe = &mou2spe(); my $hmou = &mou2hum("/Share/home/zhangqf2/GroupShare/motif/Gencode/hum2mou.genepair.syn"); my $gene_seq = &hash_fa();
my %hum_spe = %{$hspe}; my %mou_spe = %{$mspe}; my %hum_mou = %{$hmou}; my %gene_fa = %{$gene_seq};

my $human_gene_fa = "/Share/home/zhangqf2/GroupShare/motif/Gencode/hum.lncGeneMerge.fa"; my $mouse_gene_fa = "/Share/home/zhangqf2/GroupShare/motif/Gencode/mou.lncGeneMerge.fa"; 


print "ALL input file ready!\n";
my $mp = Parallel::ForkManager->new($proc);

foreach my $key (sort {$a cmp $b} keys %com_mot){
	$mp->start and next;
	my $hmot_gfa; my $mmot_gfa;
	my $hmot_file; my $mmot_file; my $reinf; my $num;
	my $hgene; my $mgene;
	my %hmotif = (); my %mmotif = ();
	my $HHdir = substr($key, 0, 4);
	`mkdir -p ${HHdir}`;
	if ($com_mot{$key} eq 2){
		%hmotif = %{$hum_mot{$key}};
		%mmotif = %{$mou_mot{$key}};
		$hmot_file = $key."_hum";
		$mmot_file = $key."_mou";
		&pri_motif(\%hmotif, $key, "hum");
		&pri_motif(\%mmotif, $key, "mou");
		($reinf, $num) = &motif_define2(\%hmotif, 0.5, $hmot_file, $human_gene_fa);
		%hmotif = %{$reinf};
		&pri_motif(\%hmotif, $key, "hum");
		($reinf, $num) = &motif_define2(\%mmotif, 0.5, $mmot_file, $mouse_gene_fa);
		%mmotif = %{$reinf};
		&pri_motif(\%mmotif, $key, "mou");
		$hgene = &motif_hgene($hmot_file, $human_gene_fa);
		$mgene = &motif_hgene($mmot_file, $mouse_gene_fa);
		&motif_gfa($hgene, $mgene, \%hum_spe, \%mou_spe, $key, $gene_seq);
		&six_motif(\%mmotif, $key);
	}elsif($com_mot{$key} eq 1){
		%hmotif = %{$hum_mot{$key}};
		$hmot_file = $key."_hum";
		&pri_motif(\%hmotif, $key, "hum");
		($reinf, $num) = &motif_define2(\%hmotif, 0.5, $hmot_file, $human_gene_fa);
		%hmotif = %{$reinf};
		&pri_motif(\%hmotif, $key, "hum");
		$hgene = &motif_hgene($hmot_file, $human_gene_fa);
		&hum_mou_gfa($hgene, \%hum_mou, "mou", $key, $gene_seq);
		$mmot_gfa = $key."_mou.fa";
		&pri_motif(\%hmotif, $key, "mou");
		$mmot_file = $key."_mou";
		($reinf, $num) = &motif_define2(\%hmotif, 0.5, $mmot_file, $mmot_gfa);
		if($num eq 0){
			#next;
			`rm $mmot_file`;
			`rm $hmot_file`;
			`rm $mmot_gfa`;
		}else{
			%mmotif = %{$reinf};
			&pri_motif(\%mmotif, $key, "mou");
			($reinf, $num) = &motif_define2(\%mmotif, 0.5, $mmot_file, $mmot_gfa);
			%mmotif = %{$reinf};
			&pri_motif(\%mmotif, $key, "mou");
			$mgene = &motif_hgene($mmot_file, $mmot_gfa);
			&motif_gfa($hgene, $mgene, \%hum_spe, \%mou_spe, $key, $gene_seq);
			&six_motif(\%mmotif, $key);
			`rm $mmot_gfa`;
		}
	}else{
		%mmotif = %{$mou_mot{$key}};
		$mmot_file = $key."_mou";
		&pri_motif(\%mmotif, $key, "mou");
		($reinf, $num) = &motif_define2(\%mmotif, 0.5, $mmot_file, $mouse_gene_fa);
		%mmotif = %{$reinf};
		&pri_motif(\%mmotif, $key, "mou");
		$mgene = &motif_hgene($mmot_file, $mouse_gene_fa);
		&hum_mou_gfa($mgene, \%hum_mou, "hum", $key, $gene_seq);
		$hmot_gfa = $key."_hum.fa";
		&pri_motif(\%mmotif, $key, "hum");
		$hmot_file = $key."_hum";
		($reinf, $num) = &motif_define2(\%mmotif, 0.5, $hmot_file, $hmot_gfa);
		if($num eq 0){
			#next;
			`rm $mmot_file`;
			`rm $hmot_file`;
			`rm $hmot_gfa`;
		}else{
			%hmotif = %{$reinf};
			&pri_motif(\%hmotif, $key, "hum");
			($reinf, $num) = &motif_define2(\%hmotif, 0.5, $hmot_file, $hmot_gfa);
			%hmotif = %{$reinf};
			&pri_motif(\%hmotif, $key, "hum");
			$hgene = &motif_hgene($hmot_file, $hmot_gfa);
			&motif_gfa($hgene, $mgene, \%hum_spe, \%mou_spe, $key, $gene_seq);
			&six_motif(\%mmotif, $key);
			`rm $hmot_gfa`;
		}
	}
	`mv ${key}_* ${HHdir}`;
	print $key."finished!\n" if ( not $? );
	$mp->finish;
}
$mp->wait_all_children;


############################################
#
#The following contents is the sub function
#
#############################################

sub pri_motif{
	my $mot = shift; my $motif_name = shift; my $spe_name = shift;
	my %mot_inf = %{$mot};
	my $mot_file = $motif_name."_".$spe_name;
	open(OUT,">",$mot_file);
	my $len = ${$mot_inf{0}}[2]; my $i; my $new_var = sprintf("%.3f", ${$mot_inf{0}}[6]);
	print OUT ${$mot_inf{0}}[0], ${$mot_inf{0}}[1], " renew= ",${$mot_inf{0}}[3],"/",${$mot_inf{0}}[4],"/",${$mot_inf{0}}[5]," inf= ",$new_var;
	for($i=1; $i<=$len; $i++){
		$new_var = sprintf("#%.3f", &INF(${$mot_inf{$i}}[0], ${$mot_inf{$i}}[1], ${$mot_inf{$i}}[2], ${$mot_inf{$i}}[3]));
		print OUT $new_var;
	}
	print OUT "\n";
	for($i=1; $i<=$len; $i++){
		$new_var = sprintf(" %.6f  %.6f  %.6f  %.6f \n", ${$mot_inf{$i}}[0], ${$mot_inf{$i}}[1], ${$mot_inf{$i}}[2], ${$mot_inf{$i}}[3]);
		print OUT $new_var;
	}
	#my $i = length(keys %mot_inf);
	#print "MEME version 4.10.1 (Release date: Wed Mar 25 11:40:43 2015 +1000)\nstrands: +\n\nMOTIF\t$mot\nletter-probability matrix: alength= 4 w= $i\n";
	#foreach my $keys (sort {$a<=>$b} keys %mot_inf){
		#print ${$mot_inf{$keys}}[0],"\t",${$mot_inf{$keys}}[1],"\t",${$mot_inf{$keys}}[2],"\t",${$mot_inf{$keys}}[3],"\n";
		#my $new_var = sprintf("\s%.6f\s\s%.6f\s\s%.6f\s\s%.6f\s\n", ${$mot_inf{$keys}}[0],${$mot_inf{$keys}}[1],${$mot_inf{$keys}}[2],${$mot_inf{$keys}}[3]);
		#print OUT $new_var;}
	close OUT;
}

sub getmotif{
	my $hum_file = shift; my $mou_file = shift; 
	my $sen = ""; my $seq = ""; my $sna = "";
	my %hum_mot = (); my %mou_mot = ();
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $motif = "--"; my $width = 0; my $newLine; my $mver; my $hmot;
	my $i = 0; my $j = 0;
	open(FILE1, $hum_file)||die("open $hum_file error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		if($sen =~m/^MEME version/){
			$mver = $sen."\n";
		}elsif($sen =~m/^strands:/){
			$mver = $mver.$sen."\n\n";
		}
		elsif($sen =~m/^MOTIF\s+(\w+)/){
			$hmot = $sen."\n";;
			$motif = $1;
		}elsif($sen =~m/^letter-probability matrix:.+ w=\s+(\d+)/){
			$width = $1;
			$hmot = $hmot.$sen;
			${$hum_mot{$motif}}{0} = [$mver, $hmot, $width, 0, 0, 0, 0];
			for ($i=1; $i<=$width; $i++){ 
				$newLine=<FILE1>;
				chomp($newLine);
				@sen = split(/\s+/,$newLine);
				${$hum_mot{$motif}}{$i} = [$sen[1],$sen[2],$sen[3],$sen[4]]; 
			}
		}
	}
	close FILE1;
	open(FILE1, $mou_file)||die("open $mou_file error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		if($sen =~m/^MEME version/){
			$mver = $sen."\n";
		}elsif($sen =~m/^strands:/){
			$mver = $mver.$sen."\n\n";
		}elsif($sen =~m/^MOTIF\s+(\w+)/){
			$hmot = $sen."\n";;
			$motif = $1;
		}elsif($sen =~m/^letter-probability matrix:.+ w=\s+(\d+)/){
			$width = $1;
			$hmot = $hmot.$sen;
			${$mou_mot{$motif}}{0} = [$mver, $hmot, $width, 0, 0, 0, 0];
			for ($i=1; $i<=$width; $i++){ 
				$newLine=<FILE1>;
				chomp($newLine);
				@sen = split(/\s+/,$newLine);
				${$mou_mot{$motif}}{$i} = [$sen[1],$sen[2],$sen[3],$sen[4]]; 
			}
		}
	}
	close FILE1;
	return (\%hum_mot, \%mou_mot);
}

sub six_motif{
	my $inf = shift; my $motif_name = shift;
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my %inf = %{$inf};	
	my $mot_file; my $fa_file; my $reinf;
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $i; my $j; my $r; my $num;
	for($i=2; $i<=$#spe; $i++){
		$mot_file = $motif_name."_".$spe[$i];
		&pri_motif(\%inf, $motif_name, $spe[$i]);
		$fa_file = $motif_name."_".$spe[$i].".fa";
		($reinf, $num) = &motif_define2(\%inf, 0.5, $mot_file, $fa_file);
		%inf = %{$reinf};
		if($num eq 0){
			next;
			`rm $fa_file`;
			`rm $mot_file`;
		}else{
			&pri_motif(\%inf, $motif_name, $spe[$i]);
			($reinf, $num) = &motif_define2(\%inf, 0.5, $mot_file, $fa_file);
			%inf = %{$reinf};
			&pri_motif(\%inf, $motif_name, $spe[$i]);
			`rm $fa_file`;
		}		
	}
}

sub motif_define2{
	my $mot_mat = shift; my $p = shift; my $mo_file = shift; my $fa_file = shift; my $limit = 0.1; my $limit2 = 0.05;
	my %mat = %{$mot_mat}; my %remat = (); my %mot_seq = (); my %ori_mat = %mat; my %max_mat = %mat;
	my $fimo_out = $mo_file.".fimo";
	my $fimo_filter_fa = $fimo_out.".fa";
	my $fimo_filter = $fimo_out.".filter";
	my $sen = ""; my $seq = ""; my $sna = "";
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $num; my $fi_num; my $keys; my $w1; my $w2; my $sta; my $rem; my $max; my $ke; my $len = ${$mat{0}}[2];
	my $i = 0; my $j = 0; my $k = 0; 
	srand(time());
	foreach $keys (sort {$a<=>$b} keys %mat){
		$remat{$keys} = $mat{$keys};
	}
	`/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text $mo_file $fa_file > $fimo_out`;
	print "fimo finished!\n" if ( not $? );
	open(FILE1, $fimo_out)||die("open $fimo_out error!\n");
	open(OUT, ">", $fimo_filter_fa);
	$sen = <FILE1>;
	#pattern name	sequence name	start	stop	strand	score	p-value	matched sequence
	#1_MHYTD_10bp	NONMMUG000004	2594	2603	-	10.1212	4.85e-05		GTGTGTGTGG
	while($sen = <FILE1>){
		$i = $i + 1;
		chomp($sen);
		@sen = split(/\t/,$sen);
		$mot_seq{$i} = [$sen[5], $sen[6], $sen[8], $sen];
		print OUT ">",$i,"\n",$sen[8],"\n";
	}
	$fi_num = $i;
	close FILE1;
	close OUT;
	`/Share/home/zhangqf2/xiongtl/misc/disk-ioz/TErelated/rmblast/rmblast/bin/dustmasker -in $fimo_filter_fa -out $fimo_filter -window 30 -level 10 -outfmt fasta`;
	print "dustmasker finished!\n" if ( not $? );
	$num = 0;
	open(FILE1, $fimo_filter)||die("open $fimo_filter error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		if($sen =~m/>([0-9]+)$/){
			$i = $1;
			$sen = <FILE1>;
			chomp($sen);
			my $count = ($sen =~ s/([ATCG])/$1/g);
			my $lsen = length($sen);
			if($count/$lsen >= 0.6){
				$num = $num + 1;
			}else{
				delete $mot_seq{$i};
			}
		}
	}
	close FILE1;
	`rm $fimo_filter_fa`;
	`rm $fimo_filter`;
	
	%ori_mat = %mat;
	%max_mat = %mat;
	%remat = %mat;
	if( $num eq 0){
		${$remat{0}}[3] = $num;
		${$remat{0}}[4] = $num;
		${$remat{0}}[5] = $fi_num;
		#${$remat{0}}[6] = $rem;
		#return (\%remat, $num);
	}else{
		$k = 0;
		$w1 = 1 - $p;
		$w2 = $p/$num;
		$sta = &motif_inf(\%mat);
		
		$max = $sta;
		#print $sta,"\n";
		#&pri_motif(\%mat);
		foreach $keys (sort {${$mot_seq{$b}}[0] <=> ${$mot_seq{$a}}[0]} keys %mot_seq){
			$k = $k + 1;
			#%mat = %remat;
			@sen1 = split(//, uc(${$mot_seq{$keys}}[2]));
			#print $#sen1,"\n";
			#foreach $j (sort {$a<=>$b} keys %remat){
			for ($j=1; $j<=$len; $j++){
				$remat{$j} = [0, 0, 0, 0];
				#$k = rand();
	
				if($sen1[$j-1] eq "A"){
					${$remat{$j}}[0] = ${$remat{$j}}[0] + 1;
				}elsif($sen1[$j-1] eq "C"){
					${$remat{$j}}[1] = ${$remat{$j}}[1] + 1;
				}elsif($sen1[$j-1] eq "G"){
					${$remat{$j}}[2] = ${$remat{$j}}[2] + 1;
				}else{
					${$remat{$j}}[3] = ${$remat{$j}}[3] + 1;
				}
			}
			for ($j=1; $j<=$len; $j++){
				for ($i=0; $i<=3; $i++){
					${$remat{$j}}[$i] = $w2/($w1 + $w2) * ${$remat{$j}}[$i] + $w1/($w1 + $w2) * ${$mat{$j}}[$i];
				}
			}		
			$rem = &motif_inf(\%remat);
			#print $rem,"\n";
			#&pri_motif(\%remat);
			if (&motif_n50_inf(\%mat, \%remat) > 1){
				%max_mat = %remat;
				$max = $rem;
				%mat = %remat;
			}elsif(($rem >= $sta*(1-$limit)) && ($rem >= $max*(1-$limit)) && ($rem >= 0.5) && (&motif_n50_inf(\%ori_mat, \%remat) > (1-$limit2)) && (&motif_n50_inf(\%max_mat, \%remat) > (1-$limit2)) && (&motif_n50_per(\%ori_mat, \%remat) > (1-$limit2)) && (&motif_n50_per(\%max_mat, \%remat) > (1-$limit2))){
				%mat = %remat;
			}else{
				last;
			}
			#print OUT ${$mot_seq{$keys}}[3],"\n";
			$w1 = $w1 + $p/$num;
		}
		${$remat{0}}[3] = $k;
		${$remat{0}}[4] = $num;
		${$remat{0}}[5] = $fi_num;
		${$remat{0}}[6] = $rem;
	}
	`rm $fimo_out`;
	
	return (\%remat, $num);
}

sub hum2spe{
	#my @rel = ("hum2mou","hum2chc","hum2cow","hum2liz","hum2ops","hum2xtr","hum2zbf");
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my $filename;
	my %hum_spe = ();
	my $sen = "";
	my @sen = ();
	my $i;
	for ($i = 2; $i<=$#spe; $i++){
		$filename = "/Share/home/zhangqf2/GroupShare/motif/Gencode/hum2".$spe[$i].".genepair.syn";
		open(FILE1, $filename)||die("open $filename error!\n");
		while($sen = <FILE1>){
			chomp($sen);
			@sen = split(/\t/,$sen);
			$hum_spe{$spe[$i]."_".$sen[0]} = $sen[1];
		}
		close FILE1;
	}
	return \%hum_spe;
}

sub mou2spe{
	#my @rel = ("hum2mou","hum2chc","hum2cow","hum2liz","hum2ops","hum2xtr","hum2zbf");
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my $filename;
	my %mou_spe = ();
	my $sen = "";
	my @sen = ();
	my $i;
	for ($i = 2; $i<=$#spe; $i++){
		$filename = "/Share/home/zhangqf2/GroupShare/motif/Gencode/mou2".$spe[$i].".genepair.syn";
		open(FILE1, $filename)||die("open $filename error!\n");
		while($sen = <FILE1>){
			chomp($sen);
			@sen = split(/\t/,$sen);
			$mou_spe{$spe[$i]."_".$sen[0]} = $sen[1];
		}
		close FILE1;
	}
	return \%mou_spe;
}

sub mou2hum{
	my $filename = shift;
	my %mou_hum = ();
	my $sen = "";
	my @sen = ();
	open(FILE1, $filename)||die("open $filename error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		@sen = split(/\t/,$sen);
		$mou_hum{$sen[1]} = $sen[0];
		$mou_hum{$sen[0]} = $sen[1];
	}
	return \%mou_hum;
}

sub hash_fa{
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my $sen; my $sna; my $fa_file; 
	my $i; my $j; my $k;
	my %gene_fa = ();	
	for($i=0; $i<=$#spe; $i++){
		$fa_file = "/Share/home/zhangqf2/GroupShare/motif/Gencode/".$spe[$i].".lncGeneMerge.fa";
		open(FILE1, $fa_file)||die("open $fa_file error!\n");
		while($sen = <FILE1>){
			chomp($sen);
			if($sen =~m/>(.+)$/){
				$sna = $1;
				$sen = <FILE1>;
				#chomp($sen);
				${$gene_fa{$spe[$i]}}{$sna} = $sen;
			}
		}
		close FILE1;
	}
	return \%gene_fa;
}

sub motif_gfa{
	my $hgene = shift; my $mgene = shift; my $hum_spe = shift; my $mou_spe = shift; my $motif_name = shift; my $gene_data = shift; #my $fa_file = shift; my $out_file = shift;
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my @hg = @{$hgene}; my @mg = ($mgene); my %hs = %{$hum_spe}; my %ms = %{$mou_spe}; my %gene_fa = %{$gene_data};
	my $sen; my $sna; my $fa_file; my $out_file;
	my $i; my $j; my $k;
	my @sen = (); 
	my %sgene = ();
	for($j=2; $j<=$#spe; $j++){
		%sgene = ();
		for($i=0; $i<=$#hg; $i++){
			if(exists $hs{$spe[$j]."_".$hg[$i]}){
				$sna = $hs{$spe[$j]."_".$hg[$i]};
				$sgene{$sna} = 1;
			}
		}
		for($i=0; $i<=$#mg; $i++){
			if(exists $ms{$spe[$j]."_".$mg[$i]}){
				$sna = $ms{$spe[$j]."_".$mg[$i]};
				$sgene{$sna} = 1;
			}
		}
		#$fa_file = $spe[$j].".mergeGene.fa"
		$out_file = $motif_name."_".$spe[$j].".fa";
		#open(FILE1, $fa_file)||die("open $fa_file error!\n");
		open(OUT, ">", $out_file);
		foreach $k (keys %sgene){
			if(exists ${$gene_fa{$spe[$j]}}{$k}){
				print OUT ">",$k,"\n";
				print OUT ${$gene_fa{$spe[$j]}}{$k};
			}else{
				print "match error\n";
			}
		}
		close OUT;
		#close FILE1;
	}
}

sub hum_mou_gfa{
	my $gene = shift; my $hum_mou = shift; my $spename = shift; my $motif_name = shift; my $gene_data = shift; #my $fa_file = shift; my $out_file = shift;
	my @spe = ("hum", "mou", "cow", "ops", "chc", "liz", "xtr", "zbf");
	my @genel = @{$gene}; my %mou_hum = %{$hum_mou}; my %gene_fa = %{$gene_data};
	my $sen; my $sna; my $fa_file; my $out_file;
	my $i; my $j; my $k;
	my @sen = (); 
	my %sgene = ();
	for($i=0; $i<=$#genel; $i++){
		if(exists $mou_hum{$genel[$i]}){
			$sna = $mou_hum{$genel[$i]};
			$sgene{$sna} = 1;
		}
	}
	#$fa_file = $spe[$j].".mergeGene.fa"
	$out_file = $motif_name."_".$spename.".fa";
	#open(FILE1, $fa_file)||die("open $fa_file error!\n");
	open(OUT, ">", $out_file);
	foreach $k (keys %sgene){
		if(exists ${$gene_fa{$spename}}{$k}){
			print OUT ">",$k,"\n";
			print OUT ${$gene_fa{$spename}}{$k};
		}else{
			print "match error\n";
		}
	}
	close OUT;
	#close FILE1;
}


sub motif_hgene{
	my $mo_file = shift; my $fa_file = shift;
	my %hg = ();
	my $sen = ""; my $seq = ""; my $sna = ""; my $file = $mo_file.".fimo";
	my @sen = (); my @sen1 = (); my @sen2 = (); my @hgene = ();
	my $i = 0;
	`/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text $mo_file $fa_file > $file`;
	open(FILE1, $file)||die("open $file error!\n");
	$sen = <FILE1>;
	#pattern name	sequence name	start	stop	strand	score	p-value	matched sequence
	#1_MHYTD_10bp	NONMMUG000004	2594	2603	-	10.1212	4.85e-05		GTGTGTGTGG
	while($sen = <FILE1>){
		#$i = $i + 1;
		chomp($sen);
		@sen = split(/\t/,$sen);
		$hg{$sen[1]} = 1;
	}
	close FILE1;
	`rm $file`;
	@hgene = keys %hg;
	return \@hgene;
}

sub motif_inf{
	my $mot_inf = shift;
	my %mot_inf = %{$mot_inf};
	my $i = 0; my $j; my $key; my $sum = 0.0; my $len = ${$mot_inf{0}}[2];
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		#print &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]),"\t";
		$sum = $sum + &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]);
	}
	$sum = $sum/$i;
	#print $sum,"\n";
	return $sum;
}

sub INF{
    my $A=shift; my $C=shift; my $G=shift; my $T=shift;
    my $Ainf=($A==0)?0:$A*log($A)/log(2); my $Cinf=($C==0)?0:$C*log($C)/log(2);
    my $Ginf=($G==0)?0:$G*log($G)/log(2); my $Tinf=($T==0)?0:$T*log($T)/log(2);
    my $H= 2 + ($Ainf+ $Cinf + $Ginf + $Tinf);
    #$H= 2 + ( $A*log($A)/log(2) + $C*log($C)/log(2) + $G*log($G)/log(2)+ $T*log($T)/log(2) );
    ### info carried by a motif is compared to a background motif without info. (0.25, 0.25, 0.25, 0.25)
    return $H;
}


sub motif_n50_inf{
	my $inf1 = shift; my $inf2 = shift;
	my %mot_inf1 = %{$inf1}; my %mot_inf2 = %{$inf2}; my %pinf = ();
	my $i = 0; my $j; my $key; my $sum = 0.0; my $len = ${$mot_inf1{0}}[2]; my $s1 = 0.0; my $s2 = 0.0;
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		#print &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]),"\t";
		$pinf{$key}= [&INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]), &INF(${$mot_inf2{$key}}[0], ${$mot_inf2{$key}}[1], ${$mot_inf2{$key}}[2], ${$mot_inf2{$key}}[3])];
		$sum = $sum + &INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]);
	}
	$sum = $sum/2;
	foreach $key ( sort {${$pinf{$b}}[0] <=> ${$pinf{$a}}[0]} keys %pinf){
		$s1 = $s1 + ${$pinf{$key}}[0];
		$s2 = $s2 + ${$pinf{$key}}[1];
		if($s1 >= $sum){
			last;
		}
	}
	return $s2/$s1;
}

sub motif_n50_per{
	my $inf1 = shift; my $inf2 = shift;
	my %mot_inf1 = %{$inf1}; my %mot_inf2 = %{$inf2}; my %pinf = ();
	my $i = 0; my $j; my $key; my $sum1 = 0.0; my $sum2 = 0.0; my $len = ${$mot_inf1{0}}[2]; my $s1 = 0.0; my $s2 = 0.0;
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		#print &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]),"\t";
		$pinf{$key}= [&INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]), &INF(${$mot_inf2{$key}}[0], ${$mot_inf2{$key}}[1], ${$mot_inf2{$key}}[2], ${$mot_inf2{$key}}[3])];
		$sum1 = $sum1 + ${$pinf{$key}}[0];
		$sum2 = $sum2 + ${$pinf{$key}}[1];
	}
	#$sum = $sum/2;
	foreach $key ( sort {${$pinf{$b}}[0] <=> ${$pinf{$a}}[0]} keys %pinf){
		$s1 = $s1 + ${$pinf{$key}}[0];
		$s2 = $s2 + ${$pinf{$key}}[1];
		if($s1 >= $sum1/2){
			last;
		}
	}
	return ($s2/$sum2)/($s1/$sum1);
}


sub motif_ou_inf{
	my $inf1 = shift; my $inf2 = shift;
	my %mot_inf1 = %{$inf1}; my %mot_inf2 = %{$inf2};
	my $i = 0; my $j; my $key; my $sum = 0.0; my $len = ${$mot_inf1{0}}[2];
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		#print &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]),"\t";
		$sum = $sum + &INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]);
	}
	$sum = $sum/$i;
	#print $sum,"\n";
	return $sum;
}

################################################



