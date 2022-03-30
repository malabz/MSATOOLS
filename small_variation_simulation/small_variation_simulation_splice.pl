#!/usr/bin/env perl

package GenerateAlignment;

=head1 NAME

Generate simulated multi-FASTA alignment data

=head1 DESCRIPTION

This script will allow you to generate simulated multi-FASTA alignment files.

=head1 LICENSE

Wellcome Trust Sanger Institute
Copyright (C) 2016  Wellcome Trust Sanger Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


1 随机生成0-1的3个小数，A T G 剩下的是C的比例
2 序列长1000bp，从0-999中随机抽取、A%的数字变成A，$DNA[0,2,4]="	

=head1 AUTHOR

Wellcome Trust Sanger Institute

=cut

# use Moose;
use Cwd qw(abs_path); 
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(shuffle);
use List::MoreUtils qw(uniq);
use List::Maker; 




my($output_filename, $num_sub, $num_ins, $num_del, $num_samples, $genome_length, $help);
GetOptions(
       'o|output=s'        => \$output_filename,
       'sub|num_sub=i'      => \$num_sub,
	   'ins|num_ins=i'      => \$num_ins,
	   'del|num_del=i'      => \$num_del,
       's|num_samples=i'   => \$num_samples,
       'l|genome_length=i' => \$genome_length,
	   'a|percentA=f'	=> \$pA,
	   'a|percentA=f'	=> \$pA,
	   'a|percentA=f'	=> \$pA,
		'h|help'            => \$help
   );
	 
	 if(!defined($output_filename) || !defined($num_sub) || !defined($num_samples) || !defined($genome_length) || defined($help))
	 {
     print <<USAGE;
Usage:   snp_sites_create_simulated_data [options]
Create a simulated multi-FASTA alignment file

Options: 
         -s INT    Number of samples
         -l INT    Number of bases in each genome
		 -sub INT    Number of substitution sites in the alignment
		 -ins INT    Number of insertion sites in the alignment
		 -del INT    Number of deletion sites in the alignment
         -o STR    Output file name
         -h        this help message

Example: Create an alignment with 20 samples, 1000 bases, 100 substitution sites, 10 insertion sites and 10 deletion sites in each genome

     small_variation_simulation -o output.aln -s 20 -l 1000 -sub 10 - ins 10 -del  10
USAGE
   exit();
	 }
   


my %number;
$number{"A"} = int( $genome_length * 0.25 );
$number{"T"} = int( $genome_length * 0.25 );
$number{"C"} = int( $genome_length * 0.25 );
$number{"G"} = $genome_length - $number{"A"} - $number{"T"} - $number{"C"};

my @a= <"A" xx $number{"A"}>;  #Maker xx
my @t= <"T" xx $number{"T"}>; 
my @c= <"C" xx $number{"C"}>; 
my @g= <"G" xx $number{"G"}>; 



my @dna=(@a, @t, @c, @g);
@dna = shuffle(@dna); ## 随机洗牌产生一条随机DNA序列

print "template DNA seq:\n", @dna,"\n";

my $out_seq0_io = Bio::SeqIO->new( -file => ">" . "template.fasta", -format => 'Fasta' );
$out_seq0_io->write_seq( Bio::Seq->new( -display_id => "sample_0", -seq => join("",@dna) ) );



## 复制模板序列生成等长度alignment
my @alignment; 
# my @tmp = @dna;
# for $i(0..$num_samples-1){
# 	$alignment[$i]=\@dna;
# }


## sub and del
for my $i (1..$num_samples-1){  

	my @delsub_sites = uniq_randompositions($num_del+$num_sub,$genome_length);
	my @delsites = @delsub_sites[0..$num_del-1];
	my @subsites = @delsub_sites[$num_del..$num_del+$num_sub-1];
	print "delsites: ", join("-",@delsites),"\t";
	print "subsites: ", join("-",@subsites),"\n";
	my @tmp = @dna;
	$alignment[$i]=\@tmp; #拷贝一份防止来回修改
	for my $subsite ( @subsites){  
		$alignment[$i]->[$subsite] = substitution($dna[$subsite]);
	}  
	for my $delsite ( @delsites){  
		$alignment[$i]->[$delsite] = '-';
	}
	print @{$alignment[$i]},"\n";
}  




## insertion
# 统一插入"-"，稀疏插入位点插入具体的碱基（预存hash）,后重新计算插入位点
my %hash_insites;
my @all_insites;

for my $i (1..$num_samples-1){

	my @ins_sites =  uniq_randompositions($num_ins,$genome_length);
	
	push(@all_insites, @ins_sites);
	
	$hash_insites{$i} = \@ins_sites;
	$hash_insbases{$i}= [insertions($num_ins)];
}

# print join("-",@all_insites),"\n";

# print @all_insite

# 插入位点全部都插入"-" perl -e '$dna = "ATCGATCG"; substr($dna,1,0,"RR");;print $dna'
my @uniq_insites = sort {$a <=> $b} uniq @all_insites;
print join("-",@uniq_insites),"\n\n";

# print @uniq_insites;
my @template_align=@dna;
my $num_uniq_insites=0;
for my $ins (@uniq_insites){
	my @template_align =  [splice (@template_align,$ins+$num_uniq_insites,0,"-")];
	$num_uniq_insites++;
}
print @template_align,"\n";

my $out_seq_io = Bio::SeqIO->new( -file => ">" . $output_filename, -format => 'Fasta' );
$out_seq_io->write_seq( Bio::Seq->new( -display_id => "sample_0", -seq => join("",@template_align )) );


for my $i (1..$num_samples-1){
	my $num_uniq_insites=0;
	my $index_insites=0;
	my @alignment_i = @{$alignment[$i]} ;
	print @{$hash_insites{$i}},"\n";
	for my $ins (@uniq_insites){
		if (grep $ins == $_,  @{$hash_insites{$i}}){

			splice(@alignment_i,$ins+$num_uniq_insites,0,$hash_insbases{$i}->[$index_insites]);
			$index_insites++;  #插入一个index增加一位 max num_ins
		}else{
			splice(@alignment_i,$ins+$num_uniq_insites,0,"-");
		}
		
		$num_uniq_insites++;	
	}

	print join("",@alignment_i ),"\n";
	$out_seq_io->write_seq( Bio::Seq->new( -display_id => "sample_$i", -seq => join("",@alignment_i) ));
}

#子程序：生成突变的位点index
sub uniq_randompositions
{
    my ($number,$length) = @_; 
    my %h; my @h;
    while ($#h+1 < $number) {
        my $r = int rand($length);
        if (!$h{$r}){
			push @h, $r;
		};
        $h{$r} = 1;
		
    }
	# my $h =  join("\t",@h),"\n";
    return @h;
}


# 随机替换成3中不同碱基中的一种
sub substitution
{
	my($base) = @_;
	my @subs = $base eq "insert" ?  qw(A T G C) : grep( !($base eq $_), qw(A T G C));
	my $newbase = $subs[rand @subs];
	return $newbase
}

# 插入的碱基list
sub insertions
{
	my($num_ins) = @_;
	my @insertions;
	for my $i (0..$num_ins-1){
		push @insertions, substitution('insert')
	}
	return @insertions
}
