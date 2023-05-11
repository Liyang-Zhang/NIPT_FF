#!/usr/bin/perl
#2207 ADD ERROR RATE
use strict;
use warnings;
use HTML::TableExtract;
use Getopt::Std;
use File::Basename;
use File::Find;

my $usage="Usage:\n\tperl $0 outfile\n
	Options:\n
    -d	input directory. default '.' .
    -s  If Scan is required.True for YES,False for NO. default 1 (YES) . 
    -t	Sequencing Type. default 'p' .
        'g' for WGS,'e' for WES,'p' for Panel sequencing.\n";    

my ($in_dir,$type,$scan);
my $GFF=0;
my %type=(p=>'Panel',e=>'Exon',g=>'Genome');
my %opts = (d => '.', t => 'p', s=>1);
my %opts1;
getopts(':d:t:s:', \%opts1);

unless (@ARGV){
	print $usage;
	exit;
}

for my $key (keys %opts1){
	my $val=$opts1{$key};
	$opts{$key}=$val;
}

for my $key (keys %opts){
	my $val=$opts{$key};
	if($key eq 'd'){
		$in_dir=$val;
		print "Input directory: $val .\n";
		}
	if($key eq 't' ){
		$val=lc $val;
		$type= $val;
		print "Sequencing Type: ", $type{$val},".\n";
		}
    if($key eq 's'){
        $scan=$val;
        print "Directory Scan is required.\n" if $scan;
        print "Sinple glob $in_dir for *_stats.\n" unless $scan;
        }
}

my @ladder=(1,20,30,50,100,300);
if($type eq 'g'){
	@ladder=(1,5,10,15,20,30);
}elsif($type eq 'e'){
	@ladder=(1,10,20,30,50,100,150);
}
my $OUT=shift;
my (@samples,%data);

my @dirs;
if($scan){
    find(sub {
		push @dirs,$File::Find::name if  -d and /_stats$/;
		}, $in_dir);
 #   print "sanning...\n";
    }else{
    @dirs=($in_dir);
#    print "No sanning...\n";
}

if(@dirs){
	print scalar @dirs," directory at  $in_dir .\n";
}else{
	die "No qualimap result directory ending with *_stats at $in_dir . \n";
}

for my $dir(@dirs){
    #my $sam=$dir;
    #$sam=~s/$in_dir\/?//;
	my $sam=basename $dir;
	$sam=~s/_stats$//;
	push @samples,$sam;
	my $html=$dir.'/qualimapReport.html';
	die "NO file:$html. May Qualimap run error." unless -e $html;
	my $te =HTML::TableExtract->new( );
	$te->parse_file($html);
	my @table=();
	for my $ts($te -> tables){
		my @cords= $ts->coords;
        if($cords[1]==0){
            for my $row($ts->rows){
                my $cmd=$$row[0];
                @table= $cmd =~/-gff/ ? (3,6,8,9,4):(2,4,6,7);
                $GFF=1 if $cmd =~/-gff/;
                $data{$sam}{SD}= $cmd =~/-sd/ ?1 :0 ;
            }
        }
        if($cords[1]==$GFF+2){
            my $flag=0;
            for my $row($ts->rows){
                $flag=1 if $$row[0] =~/^Reference size$/;
            }
            if($flag==0){
                @table=map ++$_, @table;
            }
        }
		for my $row ($ts->rows){			
			if($cords[1]==$table[0]){
                if($$row[0] eq 'Number of reads'){
                    my $str=$$row[1];
                    $str=~s/,//g;
                    $data{$sam}{Total}=sprintf("%.1f",$str/1000000);
                }
				if( $$row[0] eq 'Mapped reads'){
					my $str=$$row[1];
					$str=~s/%$//;
                    my ($num,$pct)=split /\//,$str;
					$data{$sam}{Map}=$pct;
				}
				if($$row[0] eq 'Duplicated reads (flagged)' or $$row[0] eq 'Duplicated reads (estimated)'){
					my $str=$$row[1];
					$str=~s/%$//;
                    my ($num,$pct)=split /\//,$str;
					$data{$sam}{Dup}=$pct;
				}
			}
			
			if($#table ==4 and $cords[1]==$table[4]){
				if($$row[0] eq 'Regions size/percentage of reference'){
					my $str=$$row[1];
					$str=~s/\s*\/.+$//;
					$data{$sam}{T_size}=$str;
					}
				if( $$row[0] eq 'Mapped reads'){
					my $str=$$row[1];
					$str=~s/%$//;
                    my ($num,$pct)=split /\//,$str;
					$data{$sam}{T_Map}=$pct;
				}
				if($$row[0] eq 'Duplicated reads (flagged)' or $$row[0] eq 'Duplicated reads (estimated)'){
					my $str=$$row[1];
					$str=~s/%$//;
                    my ($num,$pct)=split /\//,$str;
					$data{$sam}{T_Dup}=$pct;
				}
				
			}	
			if($cords[1]==$table[1]){
				$data{$sam}{T_Mean}=$$row[1] if $$row[0] eq 'Mean';
				if($$row[0] eq 'Mean'){
						my $str= $$row[1];
                        $str=~s/,//g;
						$str=sprintf "%.0f",$str;
						$data{$sam}{T_Mean}=$str;
					}
				}
			if($cords[1]==$table[2]){
				$data{$sam}{Insert}=$$row[1] if $$row[0] eq 'P25/Median/P75';
			}						
			if($cords[1]==$table[3]){
                if( $$row[0] eq 'General error rate'){
                    $data{$sam}{ERROR}=$$row[1];
                    $data{$sam}{ERROR}=~s/%$//;
                }
			}						
		}
        $data{$sam}{Dup}=exists $data{$sam}{Dup} ?$data{$sam}{Dup} :0 ;
        $data{$sam}{T_Dup}=exists $data{$sam}{T_Dup} ?$data{$sam}{T_Dup} :0 ;
	}
	
	my $fcov=$dir.'/raw_data_qualimapReport/coverage_histogram.txt';
	open my $fh,'<' ,$fcov or die "Cannot open file:$!";
	my (%covs,$total,%coverage);
    my $site0=0;
    my ($n20, $base0);
    my $Btarget=$data{$sam}{T_size};
    $Btarget=~s/,//g;
	while(<$fh>){
		chomp;
		next if /^#/;
		my ($read,$site)=(split /\t/)[0,1];
		$covs{$read}=$site;
        if($GFF){
            $Btarget=$Btarget-$site if $read==0 or $read==0.0;
            $n20=$read if $site0/$Btarget<0.2 and ($site0+$site)/$Btarget>=0.2;
            $site0+=$site if $read>0;
            $base0+=$site*$read;
        }
		$total+=$site;
	}
	close($fh);
	for my $read(keys %covs){
		for my $lad(@ladder){
			$coverage{$lad} += $covs{$read} if $lad < $read;
			}
	}
	$data{$sam}{Cover}=join ("\t", map {my $count=exists $coverage{$_} ? $coverage{$_} :0;
                                        sprintf "%.3f", $count/$total} @ladder);	
    my $cov30=$coverage{1} ? sprintf("%.3f",$coverage{30}/$coverage{1}) : 0;
    $data{$sam}{Cover} .= "\t".$cov30 ;
    $data{$sam}{Fold80}=sprintf("%.3f",$base0/($n20*$site0) ) if $GFF;
}

my (@keys, @outkeys);
if($GFF){
	@keys=qw(Total Map T_size T_Map T_Dup T_Mean Insert ERROR SD Fold80 Cover);
	@outkeys=qw(Total_Read(M) Map(%) T_size On_Target(%) T_Dup(%) T_Mean Insert_Size ERROR(%) SD Fold80);
}else{
	@keys=qw(Total Map Dup T_Mean Insert ERROR SD Cover);
	@outkeys=qw(Total_Read(M) Map(%) Dup(%) T_Mean Insert_Size ERROR(%) SD);
}
open my $wfh,'>' ,$OUT or die "Cannot write file:$!";
print $wfh join("\t",'Sample',@outkeys,map { '>'.$_.'X'} @ladder),"\tAdjust_30X\n";
for my $sam (sort {
                my ($da,$db);
                ($da=$data{$a}{T_size})=~s/,//g;
                ($db=$data{$b}{T_size})=~s/,//g;
                $da <=> $db or $a cmp $b;
                } @samples){
	print $wfh join("\t",$sam, map { $data{$sam}{$_} } @keys),"\n";
	
}
