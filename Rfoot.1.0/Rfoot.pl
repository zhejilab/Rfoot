if ($#ARGV < 1) {
	print "##Rfoot analyses"."\n";
	print "usage: perl Rfoot.pl -i readFile -t transcriptFile -o outputFile [-w windownSize] [-r readnumCutoff] [-f translatedORFFile]"."\n";
	print "-i readFile: input read mapping file in SAM format;"."\n";
	print "-t transcriptFile: transcript annotation file in genePred format;"."\n";
	print "-o outputFile: output file reporting candidate non-ribosomal sites;"."\n";
	print "-w windownSize [optional]: scanning window size [optional], default: 60 nt;"."\n";
	print "-r readnumCutoff [optional]: cutoff of supported read number, default: 10;"."\n";
	print "-f translatedORFFile [optional]: file with ORFs and associated translated probabilities."."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("itowrf",\%args);
my $samfile=$args{i}; 
if (! $samfile) {
	print "No read mapping file"."\n";
    exit;
}
my $transfile=$args{t}; 
if (! $transfile) {
	print "No transcript annotation file"."\n";
    exit;
}
my $outputFile=$args{o}; 
if (! $outputFile) {
	print "No output file"."\n";
    exit;
}
my $win=$args{w};
my $readcoff=$args{r};
my $orffile=$args{f};

($win=60) if (!$win);
($readcoff=10) if (!$readcoff);

############
my %dist;
for (my $i=15; $i<=50; $i++) {
	$dist{$i}=int($i/2);
}

my %read;
open (IN, "$samfile");
while (<IN>) {
	chomp;
	if ($_ !~ /^@/) {
		my @s1=split /\t/, $_;
		my @s2=split /\D+/, $s1[5];
		my @s3=split /\d+/, $s1[5];
		my $len=0;
		for (my $i=0; $i<$#s3; $i++) {
			if ($s3[$i+1] eq 'M') {
				$len+=$s2[$i];
			}
		}
		if (exists ($dist{$len})) {
			my $asite=$dist{$len};
			my $loc=$s1[3];
			my $ra=0;
			my $ind=1;
			if ($s1[1]==0) {
				for (my $i=0; $i<$#s3 && $ind==1; $i++) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc+=$asite;
							$ind=0;
						} else {
							$loc+=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc+=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc+=$s2[$i];
					}
				}
				my $k=$s1[2].":"."+".":".($loc-1);
				$read{$k}++;
			} else {
				for (my $i=0; $i<$#s3; $i++) {
					$loc+=$s2[$i];
				}
				$loc--;
				for (my $i=($#s3-1); $i>=0 && $ind==1; $i--) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc-=$asite;
							$ind=0;
						} else {
							$loc-=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc-=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc-=$s2[$i];
					}
				}
				my $k=$s1[2].":"."-".":".($loc-1);
				$read{$k}++;
			}
		}
	}
}
close IN;

##################
my %tra;
if ($orffile) {
open (IN, "$orffile");
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			if ($j >= ($s1[5]-6) && $j < ($s1[6]+6)) {
				my $k=$s1[1].":".$s1[2].":".$j;
				$tra{$k}=1;
			}
		}
	}
}
}
close IN;

#################
my $ds=3;
my %filter;
open (AN, "$transfile");
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my $pa=0;
	my @tog;
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			my $k=$s1[1].":".$s1[2].":".$j; 
			if (exists ($read{$k}) && $read{$k}>=$readcoff) {
				$pa=1;
			}
			push @tog, $j;
		}
	}
	if ($pa==1) {
	my $ind=1;
	my %sel;
	for (my $sam=0; $sam<=($#tog-$win+1); $sam++) {
		my @post;
		for (my $j=0; $j<$win; $j++) {
			push @post, $tog[$sam+$j];
		}
		my $len1=0;
		my $tot=0;
		my @val;
		if ($s1[3] eq '-') {
			@post=reverse(@post);
		}
		for (my $m=0; $m<=$#post; $m++) {
			$len1++;
			my $k=$s1[1].":".$s1[2].":".$post[$m]; 
			if (exists ($read{$k})) {
				$tot+=$read{$k};
				$val[int(($len1-1)/$ds)]+=$read{$k};
			}
		}
		if ($tot >= $readcoff) {
			my $ent=0;
			my $ten=0;
			my $a=int(($len1+2)/$ds);
			my $b=$tot;
			my $t1=int(($a+$b-1)/$b);
			my @val2;
			for ($i=0; $i<=$#val; $i++) {
				if ($val[$i] > 0) {
					$val2[int($i/$t1)]+=$val[$i];
				}
			}
			for ($i=0; $i<=$#val2; $i++) {
				if ($val2[$i] > 0) {
					my $p=$val2[$i]/($tot);
					$ent+=$p*log(1/$p);
				}
			}
			my $t2=int($a/$t1);
			my $d1=int($b/$t2);
			my $d2=$b%$t2;
			my @va;
			for (my $i=0; $i<$t2; $i++) {
				$va[$i]=$d1;
			}
			for (my $j=0; $j<$d2; $j++) {
				$va[$j]++;
			}
			for (my $i=0; $i<=$#va; $i++) {
				if ($va[$i] > 0) {
					$p=$va[$i]/($tot);
					$ten+=$p*log(1/$p);
				}
			}
			my $per;
			if ($ten == 0) {
				$per=1;
			} else {
				$per=$ent/$ten;
			}
			if ($per > 0.5) {
				for (my $m=0; $m<=$#post; $m++) {
					my $k=$s1[1].":".$s1[2].":".$post[$m]; 
					if (exists ($read{$k})) {
						$sel{$k}=1;
					}
				}
			} else {
				$ind=0;
			}
		}
	}
	if ($ind==0) {
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			my $k=$s1[1].":".$s1[2].":".$j; 
			if (exists ($read{$k}) && ! exists ($sel{$k}) && ! exists ($tra{$k})) {
				my $k2=$s1[0].":".$s1[1].":".$s1[2].":".$j; 
				$filter{$k2}=$read{$k};
			}
		}
	}
	}
	}
}
close AN;

open (OUT, ">$outputFile");
print OUT "transcriptID"."\t"."chrom"."\t"."strand"."\t"."start"."\t"."end"."\t"."length"."\t"."read.num"."\t"."max.pos"."\t"."max.num"."\t"."positions"."\n";

my %cluster;
my $dis=5;
foreach my $key (sort keys %filter) {
	if (! exists ($cluster{$key})) {
		my @sep=split /:/, $key;
		my $start=$sep[3];
		my $end=$sep[3];
		my $si;
		my $ind=1;
		while ($ind==1) {
			$ind=0;
			for (my $i=($start-$dis); $i<=($end+$dis); $i++) {
				my $k1=$sep[0].":".$sep[1].":".$sep[2].":".$i;
				if (! exists ($cluster{$k1}) && exists ($filter{$k1})) {
					if ($i<$start) {
						$start=$i;
					} 
					if ($i>$end) {
						$end=$i;
					}
					$si.=$i.":".$filter{$k1}."|";
					$cluster{$k1}=1;
					$ind=1;
				}
			}
		}
		if ($si =~ /:/) {
			my $d=$end-$start+1;
			my @s1=split /\|/, $si;
			my $ma=0;
			my $mc;
			my $tot=0;
			for (my $i=0; $i<=$#s1; $i++) {
				my @s2=split /:/, $s1[$i];
				if ($s2[1] > $ma) {
					$ma=$s2[1];
					$mc=$s2[0];
				}
				$tot+=$s2[1];
			}
			if ($ma >= $readcoff) {
				print OUT $sep[0]."\t".$sep[1]."\t".$sep[2]."\t".$start."\t".$end."\t".$d."\t".$tot."\t".$mc."\t".$ma."\t".$si."\n";
			}
		}
	}
}
close OUT;

