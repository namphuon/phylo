#!/usr/bin/perl

package Phylo;

#
#This perl script performs lots of modular tasks.  It is used to generate statistical scores, eventually mask alignments, estimate phylogenetic trees, and score them.
########

#cat green_genes/SUBPROBLEMconservativetree | sed -e 's/[0-9][0-9]*/\!/g' -e 's/[^\!]//g'|tr -d '\n'|wc -c
use strict;
use warnings;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::stat;
use Time::localtime;
#use Statistics::Distributions;
#use Bio::Phylo::IO;
#use Bio::TreeIO;
use POSIX qw(ceil floor);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
no strict 'refs';

##Global variables, programs and where they live
my $raxml = "/projects/sate3/namphuon/bin/raxmlHPC-PTHREADS-git-June9-gcc";
my $raxmlV = "/projects/sate3/namphuon/bin/raxmlHPC-PTHREADS-git-June9-gcc";
my $raxmlPC = "/projects/sate3/namphuon/bin/raxmlHPC";
my $msa_score="/u/namphuon/classes/cs394C/project/COS_LRM_v2.01/msa_set_score";
my $msa_scoreV = "/u/namphuon/classes/cs394C/project/scoring/msa_set_score";
my $cos="/u/namphuon/classes/cs394C/project/COS_LRM_v2.01/cos.pl";
my $sateV = "/u/namphuon/programs/SATE2/temp/msaml/sate/run_sate.py";
my $sate = "/u/namphuon/programs/SATE2/msaml/sate/run_sate.py";
my $config = "/u/namphuon/programs/SATE2/msaml/sate/db_large_custom.cfg";
my $readseq = "/projects/sate3/namphuon/scripts/readseq.jar";
my $gblock = "/projects/sate3/namphuon/bin/Gblocks_0.91b/Gblocks";
my $wd = "/projects/sate3/namphuon";
my $fpfn = "/u/namphuon/bin/getFpFn.py";
my $sumfpfn = "/u/namphuon/bin/getSumOfFpFnRf.py";
my $perl = "/u/namphuon/localperl/bin/perl";
my $hmmer = $ENV{'HMMER'};

my $res_pair = "_res_pair.scr";
my $sequence_name_file = "seq_names.txt";
my $HoT_file = "hot_H.fasta";
my $script_directory = "/projects/sate3/namphuon/scripts";
my $bin = "/projects/sate3/namphuon/bin";
if (defined $ENV{'BIN'}) {
  $bin = $ENV{'BIN'};
}
my $rose = "/projects/sate3/namphuon/rose/bin/rose";
my $python_bin = "/u/namphuon/bin";

my $python = "/lusr/bin/python";

my $testing = 0;
my $machine = 1;
#my $temp_dir = "/tmp";
#my $temp_dir = "/projects/sate4/namphuon/tmp";
my $temp_dir = $ENV{'TEMP_DIR'}; 
my $phylo_lab = $ENV{'PHYLO_LAB'};
my $tnt = "$bin/tnt/tnt";
my $tnt_file = "$bin/tnt/tnt_properties.tnt";
my $superfine = "$wd/superfine/reup-1.0/reup/scripts/runReup.py";
my $mrp_paup = "$wd/superfine/reup-1.0/reup/scripts/runMRP.py";

my $fastmrp = "/projects/sate9/namphuon/programs/mrpmatrix/mrp.jar";
my $fastsp = $ENV{'FAST_SP'};
my $fasttree = "$bin/FastTree";
my $blast64 = "/projects/sate7/tools/ncbi-blast-2.2.25+/bin/";
my $blast32 = "/projects/sate7/tools/ncbi-blast-2.2.25+-32bit/bin/";
my $megan = "/projects/sate7/tools/megan/MEGAN";

sub get_final_alignment_tree_pasta {
  my $log_file = $_[0];
  my $tree_line = Phylo::trim("".`grep "Writing resulting tree" $log_file`);
  my $alignment_line = Phylo::trim("".`grep "Writing resulting alignment" $log_file`);
  
  if ($tree_line eq "" or $alignment_line eq "") {
    return (undef,undef);
  }
  
  my @tree_results = split(/\s+/, $tree_line);
  my @alignment_results = split(/\s+/, $alignment_line);
  
  if (not -e $tree_results[-1] or not -e $alignment_results[-1] or -s $tree_results[-1] == 0 or -s $alignment_results[-1] == 0) {
    return (undef,undef);
  }
  return ($alignment_results[-1], $tree_results[-1]);
}

sub separate_blast_fragments {
  my $input_file = $_[0];
  my $blast_file = $_[1];
  my $output_file = $_[2];
  
  
}

sub split_fasta {
  my $file = $_[0];
  my $output = $_[1];
  my $size = $_[2];
  
  my %fragments = %{read_fasta_file($file, 0)};
  my $counter = 0;
  my $idx = 0;
  my %temp = ();
  foreach my $key (keys %fragments) {
    $temp{$key}=$fragments{$key};
    $counter++;
    if ($counter % ($size+1) == $size) {
      write_alignment(\%temp, "$output.$idx");
      $idx++;
      %temp = ();
    }    
  }
  write_alignment(\%temp, "$output.$idx");
}

sub generate_sequence_stats {
  my %sequences = %{Phylo::read_fasta_file($_[0],0)};
  my $output_file = $_[1];
  my $model = $_[2];
  
  open(OUTPUT, ">$output_file");
  foreach my $key (keys %sequences) {
    my $str = $sequences{$key};
    $str =~ s/-//g;
    print OUTPUT sprintf("$model,$key,%d\n",length($str));
  }
  close(OUTPUT);
}

#Add fake branch lengths
sub add_branch_lengths {
  my $input = $_[0];
  my $output = $_[1];
  my $results = `cat $input`;
  $results = Phylo::trim($results);
  $results =~ s/(.),/$1:1.0,/g;
  $results =~ s/\)/:1.0\)/g;  
  `echo "$results" > $output`;
}

#finds more recently modified file
sub get_most_recent_file {
  my @files = @{$_[0]};
  my $mrf = undef;
  my $mrt = -1;
  foreach my $file (@files) {
    my $timestamp = stat($file)->mtime;
    if ($timestamp > $mrt) {
      $mrt = $timestamp;
      $mrf = $file;
    }
  }
  return $mrf;
}

#reads supports from a file
sub read_supports {
  my $file = $_[0];
  my @supports = ();
  open(INPUT, "<$file");
  while (my $line = <INPUT>) {
    my @results = $line =~ m/\)(\d+\.?\d*):/g;
    push(@supports, @results);
  }
  close(INPUT);
  #print "@supports\n";
  return \@supports;
}

#reads branch lengths from a file
sub read_branch_lengths {
  my $file = $_[0];
  my @supports = ();
  open(INPUT, "<$file");
  while (my $line = <INPUT>) {
    my @results = $line =~ m/:(\d+\.?\d*)/g;
    push(@supports, @results);
  }
  close(INPUT);
  #print "@supports\n";
  return \@supports;
}


#Read memory log file, returns running time in seconds
sub read_memory_time_log {
  my $log = $_[0];
  my $result = `grep -i Time $log`;
  $result =~ m/(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
  if (not defined $1) {
    print "Can't read $log\n";
    return undef;
  }
  my $total = $1+$2+$3+$4;
  return $total;
}

#Reads a condor file
sub read_condor_log_time_log {
  my $log = $_[0];
  open(INPUT, "<$log");
  my $total = 0;
  while (my $line = <INPUT>) {
    if ($line =~ m/\(1\) Normal termination \(return value 0\)/) {
      my ($day,$hr,$min,$sec) = (0,0,0,0);
      $line = <INPUT>;
      $line = <INPUT>;
      $line = <INPUT>;
      $line =~ m/Usr\s+(\d+)\s+(\d+):(\d+):(\d+).*Sys\s+(\d+)\s+(\d+):(\d+):(\d+)/;
      $day+=$1+$5; $hr+=$2+$6; $min+= $3+$7; $sec+= $4+$8;
      $line =~ m/Usr\s+(\d+)\s+(\d+):(\d+):(\d+).*Sys\s+(\d+)\s+(\d+):(\d+):(\d+)/;
      $day+=$1+$5; $hr+=$2+$6; $min+= $3+$7; $sec+= $4+$8;
      my $time = $sec+$min*60+$hr*60*60+$day*60*60*24;
      if ($time > $total) {
        $total = $time;
      }
    }
  }
  close(INPUT);
  return $total;
}


#Finds the taxa that differ in an alignment and tree
sub find_missing {
  my $input_tree = rel2abs($_[0]);
  my %input_alignment = %{read_fasta_file(rel2abs($_[1]), 0)};
  
  my $tree_str = `cat $input_tree`;
  my @taxa = @{get_taxa($tree_str)};
  my %taxa_map = ();
  foreach my $taxon (@taxa) {
    if (exists $input_alignment{$taxon}) {
      delete $input_alignment{$taxon};
    } else {
      $taxa_map{$taxon} = $taxon;        
    }
  }   
  if (scalar (keys %taxa_map) > 0) {
    my @extra = keys %taxa_map;
    print "Tree contains the following extra taxa: @extra\n";
  } 
  if (scalar (keys %input_alignment) > 0) {
    my @extra = keys %input_alignment;
    print "Alignment contains the following extra taxa: @extra\n";
  }
}

sub average_length {
  my $input_alignment = $_[0];
  my %sequences = %{read_fasta_file($input_alignment,0)};
  my $counts = 0;
  my $lengths = 0;
  foreach my $sequence (values %sequences) {
    $counts++;
    $sequence =~ s/-//;
    $lengths+=length($sequence);
  } 
  $lengths/=$counts;
  print "$counts $lengths\n";
}

sub average_support {
  my $input_tree = rel2abs($_[0]);
  open(INPUT, "<$input_tree");
  my $line = <INPUT>;
  my @results = $line =~ m/(\d+):/g;
  my $total = 0;
  foreach my $result (@results) {
    $total+=$result;
  }
  return $total/scalar @results;
}

#Fast prune tree
sub fast_prune_tree {
  my $input_tree = rel2abs($_[0]);
  my $output_tree = rel2abs($_[1]);
  my $names_list = rel2abs($_[2]);
  my $strip = $_[3];
  
  if (defined $strip) {
    $strip = "false";
  } else {
    $strip = "";
  }
  
  #names_list is list of names of sequences to keep, delimited by space
  my $currdir = `pwd`;
  $currdir = Phylo::trim($currdir);
  chdir("$phylo_lab");
  `java -Xmx2g phylolab/prune/FastPrune $input_tree $output_tree $names_list $strip`;
  chdir($currdir);
}


sub getFastFPFN {
  my $true_tree = $_[0];
  my $estimated_tree= $_[1];
  my $output_file = $_[2];
  
  if (not -e $true_tree or not -e $estimated_tree) {
    exit;
  }

  my $ref_temp = get_temp_file() . "_";
  my $est_tmp= "$ref_temp.estimated";

  `/lusr/bin/perl /projects/sate9/namphuon/programs/phylo/CompareTree.pl -tree $true_tree -versus $estimated_tree -output $ref_temp > $ref_temp.tmp`;
  `/lusr/bin/perl /projects/sate9/namphuon/programs/phylo/CompareTree.pl -versus $true_tree -tree $estimated_tree -output $est_tmp > $est_tmp.tmp`;

  if (not -e $ref_temp) {
    print "Can't find reference result\n";
    exit;
  }  
  #echo "3"
  my $shared=`awk '{print \$1}' $ref_temp`;
  my $total_ref=`awk '{print \$2}' $ref_temp`;
  my $total_est=`awk '{print \$2}' $est_tmp`;
  
  my $fp=($total_est-$shared)/$total_est;
  my $fn=1-($shared/$total_ref);
  my $rf=($total_ref+$total_est-2*$shared)/($total_ref+$total_est);
  open(OUTPUT, ">$output_file");
  print OUTPUT "($fp, $fn, $rf)\n";
  `rm $ref_temp*`;    
}

sub getFPFN_induced {
  my $true_tree = $_[0];
  my $estimated_tree = $_[1];  
  my $output_file = $_[2];
  
  my $file = get_temp_file();
  if (not -e $true_tree or not -e $estimated_tree or -s $estimated_tree == 0) {
		return;
  }
  
  open(INPUT, "<$true_tree");
  my $true = <INPUT>;
  close(INPUT);
  
  open(INPUT, "<$estimated_tree");
  my $estimated = <INPUT>;
  close(INPUT);
  
  my %trues = map {$_=>$_} @{Phylo::get_taxa($true)};
  my %estimates = map {$_=>$_} @{Phylo::get_taxa($estimated)};
  
  my %common = %{intersection(\%trues,\%estimates)};
  my @names = keys %common;
  local $" = " ";
  open(OUTPUT, ">$file.names");
  print OUTPUT "@names";
  close(OUTPUT);
  fast_prune_tree($true_tree, "$file.true","$file.names");
  fast_prune_tree($estimated_tree, "$file.estimated","$file.names");
  getFastFPFN("$file.true","$file.estimated",$output_file);
  `rm $file.*`;
}



#Prunes an intput tree based on a list of taxa to keep/remove
sub prune_tree {
  my $input_tree = rel2abs($_[0]);
  my $output_tree = rel2abs($_[1]);
  my @taxa_name = @{$_[2]};
  my $keep = $_[3];

  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }    

  if (defined $keep) {
    $keep = "Y";
  } else {
    $keep = "";
  } 
  open(TEMP, ">$temp_file");
  foreach my $name (@taxa_name) {
    print TEMP ">$name\n-\n";
  }
  close(TEMP);
  
  my $prune = "/projects/sate3/namphuon/taxon2/metagen/scripts/python/prune_tree.py";
  print `/lusr/bin/python $prune $temp_file $input_tree $output_tree $keep`;
  `rm $temp_file -rf`;
}

sub generate_random_tree {
  my $input_tree = rel2abs($_[0]);
  my $output_tree = rel2abs($_[1]);
  my $edges = $_[2];  
  
  my $curr_directory = `pwd`;
  $curr_directory = trim($curr_directory);
  chdir("/projects/sate2/kliu/repository/sate");
  #`setenv CLASSPATH .:/projects/sate2/kliu/repository/cli/commons-cli-1.0/commons-cli-1.0.jar`;

  `java gsp.util.PEdgeContract -i $input_tree -o $output_tree.unrefined -k $edges`;
  print "java gsp.util.PEdgeContract -i $input_tree -o $output_tree.unrefined -k $edges\n";
  randomly_refine("$output_tree.unrefined", "$output_tree");    
  #`rm $output_tree.unrefined`;
  
  chdir($curr_directory);  
}

sub write_mapping {
    my %hash = %{$_[0]};
    my $output_file = $_[1];
    
    open(OUTPUT, ">$output_file");
    foreach my $key (keys %hash) {
  print OUTPUT "$key,$hash{$key}\n";    
    }
    close(OUTPUT);
}

sub read_r_table {
  my $file = $_[0];
  my @data = ();
  my %keys = ();

  open(INPUT, "<$file");
  my $line = <INPUT>;
  $line =~ s/"//g;
  $line = Phylo::trim($line);  
  my @results = split(/\s+/, $line);
  my $counter = 0;
  foreach my $result (@results) {
    $keys{$result}=$counter;
    $counter++;
  }
  
  while ($line = <INPUT>) {
    $line = Phylo::trim($line);
    my @results = split(/\s+/, $line);
    shift(@results);
    push(@data, \@results);
  }
  close(INPUT);
  return (\%keys, \@data);
}

sub read_mapping {
    my $input_file = $_[0];
    my $delimiter = $_[1];
    if (not defined $delimiter) {
      $delimiter = ",";
    }
    my %hash = ();
    open(INPUT, "<$input_file");
    my $line = <INPUT>;    
    #print "$input_file\n";
    while (defined $line) {  
      my @results = split($delimiter, trim($line));
      if (not defined $results[0] or not @results) {
        print "$input_file problem\n";
        next;
      }
      if (defined $results[0] and (scalar @results == 2)) {
        $hash{trim($results[0])}=trim($results[1]);
      } else {
        $hash{trim($results[0])} = \@results;
      } 
      $line = <INPUT>;    
    }
    close(INPUT);
    return \%hash;
}


sub sum_hash {
   my @hashes = @{$_[0]};
   
   my %sum_hash = ();
   foreach my $hash_ref (@hashes) {
     foreach my $key (keys %{$hash_ref}) {
       if (not exists $sum_hash{$key}) {
   $sum_hash{$key} = 0;
       }
       if (looks_like_number($hash_ref->{$key})) {
   $sum_hash{$key}+= $hash_ref->{$key};
       } else {
   print STDERR "$hash_ref->{$key} does not look like a number!\n";
       } 
     }
   }
   return \%sum_hash;
}


sub get_temp_file {
  my $curr_temp_dir = $_[0];
  if (not defined $curr_temp_dir) {
     #my $result = `echo \$DISPLAY`;
     #if ($result =~ m/0\.0/) {
     if (not -e $temp_dir) {
       $curr_temp_dir = "/projects/sate9/namphuon/temp/";
     } else {
       $curr_temp_dir = $temp_dir;
     }     
  }

  my $temp_file = $curr_temp_dir . "/". Phylo::get_random_name($curr_temp_dir);
  while (-e $temp_file) {
    $temp_file = $curr_temp_dir . "/".Phylo::get_random_name($curr_temp_dir);
  }
  return $temp_file;  
}

############ START METAGENOMICS ###################
sub run_megan {
  my $input_file = rel2abs($_[0]);
  my $frag_file = rel2abs($_[1]);  
  my $output_file = rel2abs($_[2]);
  my $input_mapping = rel2abs($_[3]);
  
  my ($filename, $basename, $empty) = fileparse($output_file);  
  make_directory($basename);
  
  my $temp_file = get_temp_file() . ".MEGAN";
  open(MEGAN, ">$temp_file");
  if (defined $input_mapping) {
    print MEGAN "load synonymsfile='$input_mapping';\n";
  }  
    
  print MEGAN "import blastfile='$input_file' fastafile='$frag_file' meganfile='$temp_file.rma' maxmatches=100 minscore=35.0 toppercent=20.0 winscore=0.0 minsupport=5 mincomplexity=0.0 useseed=false usekegg=false blastformat=BLASTXML;\n";    
  print MEGAN "select nodes=all;\n\n\nuncollapse subtrees;\n\n\n\nselect nodes=all;\n\n\n\n";
  print MEGAN "select nodes=all;\nuncollapse subtrees;\nselect nodes=all;\n";
  print MEGAN "select nodes=all;\nuncollapse subtrees;\nselect nodes=all;\n";
  print MEGAN "select nodes=all;\nuncollapse subtrees;\nselect nodes=all;\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nselect nodes=all;\n\n";
  print MEGAN "select nodes=all;\nuncollapse subtrees;\nselect nodes=all;\n";
  print MEGAN "export what=CSV format=readname_taxonid separator=tab file='$temp_file.out';\n";  
  print MEGAN "quit;\n";
  close(MEGAN);    
  my $results = `/projects/sate7/tools/megan/MEGAN +g < $temp_file`;
  if (-e "$temp_file.out") {
    `mv $temp_file.out $output_file`;
  } else {
    print "No such file!\n";
    print "$results\n";
    exit;
  } 
  `rm $temp_file $temp_file.rma -rf`;
}

sub convert_megan_to_tab {
  my $input_file = rel2abs($_[0]);
  my $output_file = rel2abs($_[1]);
  my @levels = ("species", "genus", "family", "order", "class", "phylum");

  my $taxonomy_file = "/projects/sate8/metagen/metadata/all_taxon.taxonomy";
  my ($taxon_ref, $level_ref, $key_ref) = read_taxonomy_mapping($taxonomy_file, "lower");
  my %taxonomy = %{$taxon_ref};
  my %level_map = %{$level_ref};
  my %key_map = %{$key_ref};
  
  my %assignments = %{Phylo::read_mapping($input_file, "\\s+")};  
  my ($filename, $basename, $empty) = fileparse($output_file);  
  make_directory($basename);  
  
  my $temp_file = get_temp_file() . ".output";
  open(OUTPUT, ">$temp_file");
  print OUTPUT "#Gene\t#Species\t#Genus\t#Family\t#Order\t#Class\t#Phylum\n";
  my @keys = sort (keys %assignments);
  foreach my $key (@keys) {
    print OUTPUT $key;
    if (not exists $taxonomy{$assignments{$key}}) {      
      #print "WTF?! $key:$assignments{$key}?!\n";
      print OUTPUT "\tNA\tNA\tNA\tNA\tNA\tNA";      
    } else {
      my @results = @{$taxonomy{$assignments{$key}}};
      foreach my $level (@levels) {
  if ($results[$key_map{$level}] eq "") {
    print OUTPUT "\tNA";
  } else {
    print OUTPUT "\t$results[$key_map{$level}]";
  } 
      }
    } 
    print OUTPUT "\n";
  }
  close(OUTPUT);
  `mv $temp_file $output_file`;
}

sub read_taxonomy_mapping {
  my $taxonomy_file = $_[0];
  my $lower = $_[1];
  if (not defined ($lower)) {
    $lower = 0;
  }
  my @level_names = ("gene", "species", "genus", "family", "order", "class", "phylum");
  my %level_maps = ("gene", {}, "species", {}, "genus", {}, "family", {}, "order", {}, "class", {}, "phylum", {});  
  open(TAXONOMY, "<$taxonomy_file");
  
  #First line used to initialize the mapping  
  my $line = <TAXONOMY>;
  my @keys = $line =~ m/\"([^\"]*)\"/g;
  #split(",",Phylo::trim($line));
  my %key_map = ();
  for my $idx (0..(scalar(@keys)-1)) {
    $key_map{$keys[$idx]} = $idx;
  }
  
  #Now fill up taxonomy, level map
  my %taxonomy = ();
  my %level_map = ();
  $line = <TAXONOMY>;
  while (defined $line) {
    $line = Phylo::trim($line);
    my @results = $line =~ m/\"([^\"]*)\"/g;
    if ($lower) {
      $results[3] = lc $results[3];
    }
    $taxonomy{$results[0]} = \@results;    
    if (scalar(@results) != scalar(@keys)) {
      print "WTF!\n@results\n" . scalar(@results) . " " . scalar(@keys) . "\n";
      exit;
    }    
    for my $idx (0..(scalar(@results)-1)) {
      if ($results[$idx] eq "") {
  next;
      }
      if (not exists $level_map{$keys[$idx]}) {
  $level_map{$keys[$idx]} = {};
      } 
      if (not exists $level_map{$keys[$idx]}->{$results[$idx]}) {      
  $level_map{$keys[$idx]}->{$results[$idx]} = {};
      } 
      $level_map{$keys[$idx]}->{$results[$idx]}->{$results[0]} = $results[0];      
    }    
    $line = <TAXONOMY>;
  }

  return (\%taxonomy, \%level_map, \%key_map);
}

sub get_taxonomy_name_id_mapping {
  my %taxonomy = %{$_[0]};
  my %id_name_map = ();
  foreach my $id (keys %taxonomy) {
    $id_name_map{lc($taxonomy{$id}->[3])} = $id;
  }
  return \%id_name_map;
}


sub make_blast_database {
  my $input_file = rel2abs($_[0]);  
  my $output_file = rel2abs($_[1]);
  my $input_mapping = rel2abs($_[2]);
  
  my $blast = get_blast_version();
  
  if (not defined $input_mapping) {
    $input_mapping = "";
  } else {
    $input_mapping = "-taxid_map $input_mapping";
  } 
  
  my ($filename, $basename, $empty) = fileparse($output_file);    
  make_directory($basename);
  
  my $temp_file = get_temp_file() . "/";
  make_directory($temp_file);     
  
  `$blast/makeblastdb -in $input_file -dbtype nucl $input_mapping -out $temp_file/$filename`;
  make_directory($basename);
  `mv $temp_file/* $basename/`;
  `rm $temp_file -rf`;    
}

sub blast_search {
  my $input_file = rel2abs($_[0]);
  my $database = rel2abs($_[1]);  
  my $output_file = rel2abs($_[2]);
  
  my $blast = get_blast_version();
  
  my ($filename, $basename, $empty) = fileparse($output_file);  
  my $temp_file = get_temp_file();
  make_directory($temp_file);
    
  make_directory($basename);
  `$blast/blastn -outfmt 5 -query $input_file -out $temp_file/$filename -db $database`;
  `mv $temp_file/$filename $output_file`;
  `rm $temp_file -rf`;
}

sub get_chip_version {
  if (-d "/scratch/cluster/namphuon/") {
    return "64";
  } else {
    return "32";
  } 
}

sub get_blast_version {
 if (get_chip_version() eq "64") {
   return $blast64;
 } else {
   return $blast32;
 }
}

#Splits a json file into multiple json
sub split_json {
  my $input_json = $_[0];
  my $size = $_[1];
  my $output = $_[2];
  
  my $temp_file = get_temp_file();
  make_directory($temp_file);

  `python \$SEPP_HOME/scripts/python/split_jsons.py -i $input_json -o $temp_file/test.json -s $size`;

  my @json_files = glob("$temp_file/test.json.*");  
  foreach my $file (@json_files) {
    $file =~ /\.(\d+)/;
    `mv $file $output.$1.json`;
  }
  `rm $temp_file -rf`;
}

#Split out a guppy sing file into its original trees
sub split_out_tree_file {
  my $input_sing = $_[0];
  my $backbone_tree = $_[1];
  my $output_dir = $_[2];
  
  my $temp_file = get_temp_file();
  Phylo::make_directory($temp_file);
  Phylo::make_directory($output_dir);
    
  open(BASE, "<$backbone_tree");
  my $tree = <BASE>;
  my %original_taxa = map {$_ => $_} @{Phylo::get_taxa($tree)};
  close(BASE);
  
  open(TREES, "<$input_sing");
  while (my $line = <TREES>) {
    my %new_taxa = map {$_ => $_} @{Phylo::get_taxa($line)};
    my %difference = %{Phylo::difference(\%new_taxa, \%original_taxa)};
    my @names = keys %difference;
    if (scalar @names != 1 or $names[0] =~ m/_/) {
      my @results = split(/_/, $names[0]);
      foreach my $result (@results) {
        my $copy = $line;
        $copy =~ s/$names[0]/$result/g;
        open(OUT, ">$temp_file/$result.tree");
        print OUT "$copy";
        close(OUT);        
      }
    } else {
      open(OUT, ">$temp_file/$names[0].tree");
      print OUT "$line";
      close(OUT);
    }
  }
  close(TREES);
  `mv $temp_file/*.tree $output_dir`;
  `rm $temp_file -rf`;
}

#Merges two or more tab files, taking the LCA that they agree upon
sub merge_tabs {
  my @estimate_files = @{$_[0]};
  my $output_file = rel2abs($_[1]);  
  
  my %tabs = ();
  my %fragments = ();
  my $counter = 0;
  foreach my $estimate_file (@estimate_files) {
    my %estimate_mapping = %{Phylo::read_multiple_mapping($estimate_file, "\t", 1)};
    foreach my $key (keys %estimate_mapping) {
      $fragments{$key}++;
    }
    $tabs{$counter} = \%estimate_mapping;
    $counter++;
  }
  
  my %output_mapping = ();
  foreach my $fragment (sort keys %fragments) {
    my $defined = 0;    
    foreach my $idx (0..($counter)) {      
      if (not defined $tabs{$idx}->{$fragment}) {
        next;
      }      
      
      #If fragment hasn't been defined, then find first occurrance 
      if (not $defined) {
        $output_mapping{$fragment} = $tabs{$idx}->{$fragment};        
        $defined = 1;
        next;
      }
      
      #Now find LCA of the new estimated and 
      my @estimate = @{$tabs{$idx}->{$fragment}};
      my @merge = @{$output_mapping{$fragment}};      
      my $lca = 6;
      while ($lca > 0) {
        if (((lc $merge[$lca]) ne (lc $estimate[$lca]))) {
          foreach my $miss (1..$lca) {
            $output_mapping{$fragment}->[$miss] = "UNCLASSIFIED";
          }
          last;
        }
        $lca--;
      }
    }
  }
  
  write_multiple_mapping(\%output_mapping, $output_file, "\t", "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum");    
  return;
}

sub score_tab {
  my $estimate_tab = rel2abs($_[0]);
  my $tax_tab = rel2abs($_[1]);
  my $output_file = rel2abs($_[2]);
  my $frags_file = $_[3];
  
  
  if (not -e $estimate_tab or not -e $tax_tab) {
    print "No such input files\n";
    exit;
  }
  
  my %taxa_mapping = %{Phylo::read_multiple_mapping($tax_tab, "\t", 1)};
  my %estimate_mapping = %{Phylo::read_multiple_mapping($estimate_tab, "\t", 1)};
  my @keys = ();
  if (not defined $frags_file) {
    @keys = keys %taxa_mapping;
  } elsif ($frags_file eq "estimated") {      
    @keys = keys %estimate_mapping; 
  } else {    
  my %frags = %{Phylo::read_fasta_file(rel2abs($frags_file), 0)};
    @keys = keys %frags;
  }
  
  my %total = ("species",0.,"genus",0.,"family",0.,"order",0,"class",0.,"phylum",0.);
  my %classified = ("species",0.,"genus",0.,"family",0.,"order",0,"class",0.,"phylum",0.);
  my %hit = ("species",0.,"genus",0.,"family",0.,"order",0,"class",0.,"phylum",0.);
  my @levels = ("species", "genus", "family", "order", "class", "phylum");

  foreach my $key (@keys) {
    if ($key eq "fragment") {
      next;
    }
    if (not exists $taxa_mapping{$key} or not defined $taxa_mapping{$key}) {
      print STDERR "$key in $estimate_tab not found in true mapping file $tax_tab\n";
      return;
    }
    my @true_results = @{$taxa_mapping{$key}};
    my @predicted_results = ();
    if (not exists $estimate_mapping{$key}) {
      @predicted_results = ("$key", "NA", "NA", "NA", "NA", "NA", "NA");
    } else {
    @predicted_results = @{$estimate_mapping{$key}};    
  } 
  foreach my $idx (0..5) {
    if (not defined $true_results[$idx+1]) {
      print "Failed to find $idx in @true_results\n, $estimate_tab $tax_tab\n";
      return;
    }
    if ($true_results[$idx+1] ne "NA") {
      $total{$levels[$idx]}++;
    } else {
      next;
    } 
    if ((uc $predicted_results[$idx+1] ne "NA" and uc $predicted_results[$idx+1] ne "UNCLASSIFIED") and not ($predicted_results[$idx+1] =~ /#/)) {
      $classified{$levels[$idx]}++;
      if ((lc $predicted_results[$idx+1]) eq (lc $true_results[$idx+1])) {
        $hit{$levels[$idx]}++;
      }
    }
  }
  }
  
  my $temp_file = get_temp_file();
  
  open(OUTPUT, ">$temp_file");
  foreach my $stat ("precision", "classified_count", "classified_perc", "total") {  
      foreach my $idx (0..5) {      
        my $value = "";
        if ($stat eq "precision") {
          if ($classified{$levels[$idx]} == 0) {
          $value = 0;
        } else {
          $value = $hit{$levels[$idx]}/$classified{$levels[$idx]}*100;
        }
        } elsif ($stat eq "classified_count") {
          $value = $classified{$levels[$idx]};
        } elsif ($stat eq "classified_perc") {
          if ($total{$levels[$idx]} == 0) {
          $value = 0;
        } else {
            $value = $classified{$levels[$idx]}/$total{$levels[$idx]}*100;
        }
        } elsif ($stat eq "total") {
            $value = $total{$levels[$idx]};
        }
        print OUTPUT sprintf("$stat\t$levels[$idx]\t%1.5f\n", $value);
      }
  }
  close(OUTPUT);  
  my ($filename, $basename, $empty) = fileparse($output_file);  
  mkdir($basename);
  if (-e $temp_file) {
    `mv $temp_file $output_file`;
  }
  return (\%hit, \%classified, \%total);
}

sub convert_stats {
  my $stats_file = $_[0];
  my $output_file = $_[1];
  my $method = $_[2];
  my $gene = $_[3];
  my $length = $_[4];
  my $prefix = $_[5];
  my $threshold = $_[6];
  my $leaveout_level = $_[7];
  
  my ($hit, $classify, $total) = Phylo::read_stats_file("$stats_file");
  my %stats;
  my $temp_file = get_temp_file();  
  open(OUTPUT, ">$temp_file");  
  foreach my $class ("species", "genus", "family", "order", "class", "phylum") {
    foreach my $stat ("precision", "classified_count", "classified_perc", "total") {
      my $value = "";
      if ($stat eq "precision") {
        $value = $hit->{$class}/($classify->{$class}+0.0000000000001);
      } elsif ($stat eq "classified_count") {
        $value = $classify->{$class};
      } elsif ($stat eq "classified_perc") {
        $value = 100*$classify->{$class}/($total->{$class}+.0000000000000000000000001);
      } elsif ($stat eq "total") {
        $value = $total->{$class};
      }
      print OUTPUT "$method $gene $length $leaveout_level $prefix\_$leaveout_level $threshold $stat $class $value\n";
    }
  }
  close(OUTPUT);
  if (-e $temp_file) {
    `mv $temp_file $output_file`;
  }
}

sub aggregate_stats_file {
  my @files = @{$_[0]};
  my $output_file = $_[1];
  my %hit = ();
  my %classified = ();
  my %total = ();
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  foreach my $file (@files) {
    if (not -e "$file") {
      next;
    } 
    
    my ($h1, $c1, $t1) = read_stats_file("$file");  
    my @t = ($h1, \%hit);  
    %hit = %{sum_hash(\@t)};
    @t = ($c1, \%classified);  
    %classified = %{sum_hash(\@t)};
    @t = ($t1, \%total);  
    %total = %{sum_hash(\@t)};      
  }
  
  open(OUTPUT, ">$output_file");
  foreach my $level (@levels) {
    if (not exists $classified{$level} or $classified{$level} == 0) {
      $classified{$level} = .00000000001;
    }
    if (not exists $total{$level} or $total{$level} == 0) {
      $total{$level} = .00000000001;
    }  
    if (not exists $hit{$level} or $hit{$level} == 0) {
      $hit{$level} = .00000000000;
    }    
    print sprintf("$level\tPrec:%0.3f\tSens:%0.3f\n", $hit{$level}/$classified{$level}, $hit{$level}/$total{$level});  
    print OUTPUT sprintf("$level,precision,%f\n", $hit{$level}/$classified{$level});
    print OUTPUT sprintf("$level,sensitivity,%f\n", $hit{$level}/$total{$level});
  }
  close(OUTPUT);  
}

# sub read_stats_file_2 {
#   my $file = $_[0];
#   open(INPUT, "<$file");
#   my $line = <INPUT>;
#   my %stats = ();
#   while (defined $line) {    
#     $line = Phylo::trim($line);
#     my @results = split(/,/, $line);
#     if (not defined $stats{$results[0]}) {
#       $stats{$results[0]} = {};
#     }
#     $stats{$results[0]}->{$results[1]}=$results[2];
#     $line = <INPUT>;
#   }
#   close(INPUT);
#   return \%stats;
# }

sub read_stats_file {
  my $stats_file = $_[0];
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  my %hit = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %total = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %classify = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  foreach my $level (@levels) {
    my $p = `grep "precision" $stats_file | grep "$level" | awk {'print \$3'}`;    
    my $c = `grep "classified_count" $stats_file | grep "$level\t" | awk {'print \$3'}`;    
    my $t = `grep "total" $stats_file | grep "$level" | awk {'print \$3'}`;
    
    my $pc = $p*$c;
    $hit{$level}+=$pc;
    $total{$level}+=$t;
    $classify{$level}+=$c;    
  }
  return (\%hit, \%classify, \%total);
}

sub make_db {
  my $json_file = rel2abs($_[0]);
  my $refpkg = rel2abs($_[1]);
  my $outputdb = rel2abs($_[2]);
  my $options = $_[3];
  my $new = $_[4];
  
  if (not defined $options) {
    $options = "";
  }
  
  my ($filename, $basename, $empty) = fileparse($outputdb);  
  if (not -e $basename and not -d $basename) {
    make_directory($basename);
  }    
  if (not defined $new or $new eq "") {
    `/projects/sate3/namphuon/taxon2/metagen/scripts//bash/make-taxonomic-db.sh $outputdb $refpkg $json_file $options`;      
  } else {
    `/projects/sate3/namphuon/taxon2/metagen/scripts//bash/make-taxonomic-db.sh.latest $outputdb $refpkg $json_file $options`;  
  } 
}

sub generate_classification_new_tipp {
  my $input = rel2abs($_[0]);
  my $output = $_[1];
  my $threshold = $_[2];
  my $taxonomy_file = $_[3];
  my $optional = $_[4];
  
  my %level_map = ("root", 6, "species",0, "genus", 1, "family", 2, "order", 3, "class", 4, "phylum", 5);
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  my ($taxon_ref, $level_ref, $key_ref) = read_taxonomy_mapping($taxonomy_file, "lower");
  my %names = ();

  open(INPUT, "<$input");
  open(OUTPUT, ">$output");
  print OUTPUT "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n";
  my $old_name = "";
  my $old_probability = "";
  my $old_id = ""; 
  my $old_rank = ""; 
  local $" = "\t";
  while (my $line = trim("".<INPUT>)) {    
    my @results = split(/,/, $line);
    my ($name, $id, $rank, $probability) = ($results[0], $results[1], $results[3], $results[4]);
    $names{$name} = $name;
    if ($name ne $old_name) {
      if ($old_name ne "") {
        my @lineage = @{$taxon_ref->{$old_id}};
        my @output_line = ($old_name);
        foreach my $level (@levels) {
          my $clade = $lineage[$key_ref->{$level}];
          if ($clade eq "") {
            $clade = "NA";
          }
          push(@output_line, $clade);
        }
        print OUTPUT "@output_line\n";
      }
      $old_name = $name;
      $old_rank = "root";
      $old_probability = 1;
      $old_id = 1;
    }
    
    if (defined $level_map{$rank} and ($level_map{$old_rank} > $level_map{$rank}) and ($probability > $threshold)) {
      $old_rank = $rank;
      $old_probability = $probability;
      $old_id = $id;
    } elsif (defined $level_map{$rank} and ($level_map{$old_rank} == $level_map{$rank}) and ($probability > $old_probability)) {
      $old_rank = $rank;
      $old_probability = $probability;
      $old_id = $id;    
    }
  }
  my @lineage = @{$taxon_ref->{$old_id}};
  my @output_line = ($old_name);
  foreach my $level (@levels) {
    my $clade = $lineage[$key_ref->{$level}];
    if ($clade eq "") {
      $clade = "NA";
    }
    push(@output_line, $clade);
  }
  print OUTPUT "@output_line\n";  
  close(OUTPUT);
  close(INPUT);
  
  if (defined $optional) {
    my %seq = %{Phylo::read_fasta_file($optional, 0)};    
    open(OUTPUT, ">>$output");
    foreach my $key (sort(keys %seq)) {
      if (not defined $names{$key}) {
        print OUTPUT "$key\tNA\tNA\tNA\tNA\tNA\tNA\n";
      }
    }
    close(OUTPUT);
  }
}


sub generate_classification {
  my $database = rel2abs($_[0]);
  my $output = $_[1];
  my $threshold = $_[2];
  my $fragments = $_[3];
  
  if (not defined $threshold or $threshold eq "") {
    $threshold = "";
  } else {
    $threshold = "-t $threshold";
  } 
  
  if (not defined $fragments or $fragments eq "") {
    $fragments = "";
  } else {
    $fragments = "-f $fragments";
  } 
  
  `/lusr/bin/python /projects/sate3/namphuon/taxon2/metagen/sepp/classification/classify.py $threshold $fragments -d $database -o $output`;  
}

sub read_meta {
  my $metafile = rel2abs($_[0]);
  open(INPUT, "<$metafile");
  my $line = <INPUT>;
  $line = <INPUT>;
  my $read_sets = 1;
  my %sequence_sets = ();
  my %fragment_assignments = ();
  
  while (defined $line) {
    if ($line =~  m/Fragment assignments/) {
      $read_sets = 0;
      $line = <INPUT>;
    }
    $line =~ m/(\d+):\s+(.*)/;
    my $idx = Phylo::trim($1);
    my $remaining = Phylo::trim($2);
    if ($read_sets) {
      my @align_names = split(/\s+/, $remaining);
      $sequence_sets{$idx} = \@align_names;            
    } else {
      my @frag_names = split(/\s+/, $remaining);
      my @temp_array = ();
      foreach my $frag_name (@frag_names) {
        $fragment_assignments{$frag_name} = $idx;
      }      
    }    
    $line = <INPUT>;
  }
  close(INPUT);
  return (\%sequence_sets, \%fragment_assignments);
}

sub read_multiple_mapping {
    my $input_file = $_[0];
    my $delimiter = $_[1];
    my $skip = $_[2];
    if (not defined $delimiter) {
      $delimiter = ",";
    }
    my %hash = ();
    if (not -e $input_file) {
      print "$input_file does not exist!\n";
      exit;
    }
    open(INPUT, "<$input_file");
    my $line = <INPUT>;    
    my $first = "";
    #If we say skip, skip the first line
    if (defined $skip and $skip) {
      $first = $line;
      $line = <INPUT>;
    }
    while (defined $line) {  
      my @results = split($delimiter, trim($line));
      if (defined $results[0] and defined $results[1]) {
        $hash{trim($results[0])}=\@results;
      }
      $line = <INPUT>;    
    }
    close(INPUT);    
    return \%hash;
}

sub write_multiple_mapping {
  my %mapping = %{$_[0]};
  my $output_file = $_[1];
  my $delimiter = $_[2];
  my $header = $_[3];
  
  if (not defined $delimiter) {
    $delimiter = " ";
  }
  local $" = $delimiter;
  open(OUTPUT, ">$output_file");
  if (defined $header) {
    print OUTPUT "$header\n";    
  }
  my @keys = sort keys %mapping;
  foreach my $key (@keys) {
    print OUTPUT "@{$mapping{$key}}\n";
  }
  close(OUTPUT);
}

sub mafft_add_fragments {
  my $backbone_alignment = $_[0];
  my $fragments = $_[1];
  my $output = $_[2];
  my $options = $_[3];
  
  if (not defined $options) {
    $options = "";
  }
  my $temp_file = get_temp_file();
  `/projects/sate9/namphuon/programs/mafft-6.956-with-extensions/bin/bin/mafft --add $fragments $options $backbone_alignment > $temp_file`;
  my %alignment = %{read_fasta_file($temp_file)};
  if ((scalar keys %alignment) > 0) {
    `mv $temp_file $output`;
  }
}

sub hmmr_global_profile {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $options = $_[2];
  
  if (not defined $options) {
    $options = "";
  }
  
  my_cmd("$hmmer/hmmbuild $options $output_file $input_file");
}

sub hmmr_profile {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $format = $_[2];
  my $options = $_[3];
  
  if (not defined $format) {
    $format = "--informat afa";
  } elsif ($format eq "stockholm") {
    $format = "";
  } else {
    $format = "--informat afa";
  } 
  
  if (not defined $options) {
    $options = "";
  }
  
  my_cmd("$hmmer/hmmbuild $options --symfrac 0.0 $format $output_file $input_file");
}

sub hmmr_align {
  my $input_file = $_[0];
  my $hmmr_file = $_[1];  
  my $output_file = $_[2];
  my $options = $_[3];
  
  my_cmd("$hmmer/hmmalign $options -o $output_file $hmmr_file $input_file");
}

sub convert_fasta_to_sto {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
   `/lusr/bin/python /projects/sate3/namphuon/taxon2/metagen/scripts/python/convert_to_sto.py $input_file $output_file`;
}

sub convert_fasta_to_fastq {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $mapping_output = $_[2];
  
  if (not -e $input_file) {
    print "Failed to find $input_file\n";
    exit();
  }
  
  my %reads = %{Phylo::read_fasta_file("$input_file",0)};  
  open(OUTPUT, ">$output_file");
  foreach my $read (keys %reads) {
    my $quality = "a" x length($reads{$read});
    print OUTPUT "\@S$read A812M3ABXX:1:1:1:1#0\n$reads{$read}\n+\n$quality\n";
  }
  close(OUTPUT);
  
}

sub filter_reads_by_length {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $length = $_[2];

  my %fragments = %{Phylo::read_fasta_file($input_file,0)};
  foreach my $key (keys %fragments) {
    if (length($fragments{$key}) < $length) {
      delete $fragments{$key};
    }
  }
  Phylo::write_alignment(\%fragments, $output_file);
}

sub hmmr_filter {
  my $input_align_file = $_[0];
  my $hmm_file = $_[1];
  my $output_file = $_[2];
  
  #Read in alignment file
  my %aln = %{read_fasta_file($input_align_file, 0)};
  
  #Read hmm search file results
  open(HMM_FILE, "<$hmm_file");
  my $line = <HMM_FILE>;
  my %good_seq = ();
  while (defined $line) {
    if ($line =~ m/^>>\s+([^\s]+)/) {
      $good_seq{$1} = $aln{$1};
    }
    $line = <HMM_FILE>;
  }  
  write_alignment(\%good_seq, $output_file);
  close(HMM_FILE);
}

sub hmmr_search {
  my $input_file = $_[0];
  my $hmm_file = $_[1];  
  my $output_file = $_[2];
  my $elim = $_[3];
  my $options = $_[4];
  
  if (not defined $elim) {
    $elim = 0.00001    
  }
  if (not defined $options) {
    $options = "";
  }
  
  my_cmd("$hmmer/hmmsearch --noali $options -E $elim -o $output_file $hmm_file $input_file");
}

sub get_hmmr_query_scores {
  my $input_file = $_[0];
  my $fragment_file = $_[1];

  my $fragments = undef;
  if (defined $fragment_file) {
    $fragments = Phylo::read_fasta_file($fragment_file, 0);
  }
  
  my %output_scores = ();
  open(INPUT, "<$input_file");
  my $line = <INPUT>;
  my $start_match = 0;
  while (defined $line) {
    $line = Phylo::trim($line);
    if ($line =~ m/inclusion threshold/) {
      $line = <INPUT>;
      next;
    }    
    
    if ($line =~ m/Scores for complete sequences/ and not $start_match) {
      $line = <INPUT>;
      $line = <INPUT>;
      $line = <INPUT>;
      $start_match = 1;
    } elsif ($start_match and $line ne "") {
      my @results = split(/\s+/, $line);
      my $name = $results[8];
      if (not defined $name) {
        print "Reading $input_file\n";                
        exit;
      }
      if ((defined $fragments and exists $fragments->{$name}) or not defined $fragments) {
        $output_scores{"$name\tbit_score"} = $results[1];
        $output_scores{"$name\te_value"} = $results[0];
      }
    } elsif ($start_match) {
      last;
    }
    $line = <INPUT>;
  }
  close(INPUT);
  return \%output_scores;
}

############ END METAGENOMICS ###################

#Returns the sp score, input is the output of lobster
sub get_sp_lobster {
    my $score_file = $_[0];
    if (not -e $score_file) {
  return -1;
    }
    
    open(FILE, "<$score_file");
    my $line = <FILE>;
    while (defined($line)) {
  $line =~ m/SP=([^;]+);/;
  if (defined($1)) {
      return $1;
  }
  $line = <FILE>;
    }
    close(FILE);
    return -1;
}

sub get_alignment_statistics {
    my $stats_file = $_[0];
    my $directory = $_[1];
    
    if (defined $directory) {
      my %stats = ();
      foreach my $stats_file (@{$stats_file}) {
  if (not -e $stats_file) {
      next;
  }
  open(FILE,"<$stats_file");
  my $line =  <FILE>;
  while (defined($line)) {
      $line =~ m/([\S]+)\s+\|\s+([^\s]+)/;
      if (defined($1) && defined($2)) {
    if (defined $stats{$1}) {
      push(@{$stats{$1}},$2);
    } else {
      $stats{$1}=[$2];
    }
      }
      $line = <FILE>;
  }
  #print %stats;
  close(FILE);  
      }      
      return \%stats;      
    } else {
      my %stats = ();
      if (not -e $stats_file) {
    return \%stats;
      }
      
      open(FILE,"<$stats_file");
      my $line =  <FILE>;
      while (defined($line)) {
    $line =~ m/([\S]+)\s+\|\s+([^\s]+)/;
    if (defined($1) && defined($2)) {
        $stats{$1}=$2;
    }
    $line = <FILE>;
      }
      #print %stats;
      close(FILE);
      return \%stats;
    }
}

#Returns an upper triangular sparse matrix containing distances
sub calc_distance_matrix {
  my $input_file = $_[0];  
  my $base = $_[1];  
  
  if (not defined $base) {
    $base = 2;
  }
  
  #Read the FASTA file
  my %name_map = %{read_fasta_file($input_file)};  
  my @names = sort(keys %name_map);
  my $taxa = scalar (@names);
  my $columns = length($name_map{$names[0]});
  
  #TODO: Can do this much more efficiently
  my @distanceMatrix = ();
  for (my $i_idx=0; $i_idx < $taxa; $i_idx++) {
    my @row = ();
    my $curr_taxa = $name_map{$names[$i_idx]};
    for (my $j_idx=$i_idx+1; $j_idx < $taxa; $j_idx++) {
      my $next_taxa = $name_map{$names[$j_idx]};
      my $distance = 0;
      for (my $k_idx = 0; $k_idx < $columns; $k_idx++) {
        $a = substr( $curr_taxa, $k_idx , 1 );
        $b = substr( $next_taxa, $k_idx , 1 );        
        if ($a eq "?" or $b eq "?" and $a ne $b) {
          $distance+=0.5;
        }
        elsif ($a ne $b) {
          $distance+=1.0;
        }
      }
      #Compute the corrected distance for CF model
      if ($distance/$columns * $base/($base-1.0) >= 1.0) {        
        $distance = -($base-1.0)/$base * log((1.0)/($columns*$base));
      }
      else {
        $distance = -($base-1.0)/$base * log(1.0-$base/($base-1.0)*($distance/$columns));
      }
      push(@row, $distance);
    }
    push(@distanceMatrix, \@row);
  }
  return (\@names, \@distanceMatrix);
}

sub create_stat_file {
  my $column_file = $_[0];
  my $keep_file = $_[1];
  my $output_file = $_[2];
  my $meta = $_[3];
  
  if (not -e "$column_file" or not -e "$keep_file") {
    print "$column_file or $keep_file is missing\n";
    exit;
  }
  open(INPUT, "<$keep_file");
  my $line = Phylo::trim("".<INPUT>);
  my @results = split(/\s+/, $line);  
  close(INPUT);

  my %keep_sites = map {$_=>$_} @results;

  my %original = ();  
  open(INPUT, "<$column_file");
  
  my @totals = (0,0,0);
  my @totals_masked = (0,0,0);
  #site,correct,incorrect,total
  $line = <INPUT>;
  while ($line = <INPUT>) {
    $line = Phylo::trim($line);
    my @results = split(/,/, $line);
    $original{$results[0]} = \@results;
    $totals[0]+= $results[1];
    $totals[1]+= $results[2];
    $totals[2]+= $results[3];

    if (defined $keep_sites{$results[0]}) {
      $totals_masked[0]+= $results[1];
      $totals_masked[1]+= $results[2];
      $totals_masked[2]+= $results[3];
    }
  }
  close(INPUT);
  
  my $temp_file = Phylo::get_temp_file();
  open(OUTPUT, ">$temp_file");
  #dataset,name,method,alignment,rep,pre_total_correct_homo,pre_total_incorrect_homo,post_total_correct_homo,post_total_incorrect_homo
  print OUTPUT "$meta,$totals[0],$totals[1],$totals_masked[0],$totals_masked[1]";
  close(OUTPUT);
  `mv $temp_file $output_file`;
}

sub find_masked_columns {
  my $original_alignment_file = $_[0];
  my $masked_alignment_file = $_[1];
  my $output_file = $_[2];
  
  my %original_alignment = %{Phylo::read_fasta_file($original_alignment_file,0)};
  my %masked_alignment = %{Phylo::read_fasta_file($masked_alignment_file,0)};

  my @keys = sort keys %original_alignment;

  if (length($masked_alignment{$keys[0]}) == 0) {
    open(OUTPUT, ">$output_file");
    print OUTPUT "";
    close(OUTPUT);
  }

  my %col_original_alignment = ();
  my %col_masked_alignment = ();
  my $length_original = 0;
  my $length_masked = 0;
  foreach my $key (@keys) {
    my @temp = split(//, uc $original_alignment{$key});
    $original_alignment{$key} = \@temp;
    $length_original = scalar @temp;
    foreach my $idx (0..(scalar @temp - 1)) {
      $col_original_alignment{$idx}.=$temp[$idx];
    }
    if (not defined $masked_alignment{$key}) {
      print "Error! $masked_alignment_file\n";
      exit;
    }

    @temp = split(//, uc $masked_alignment{$key});
    $masked_alignment{$key} = \@temp;
    $length_masked = scalar @temp;
    foreach my $idx (0..(scalar @temp - 1)) {
      $col_masked_alignment{$idx}.=$temp[$idx];
    }
  }
  
  my $original_idx = 0;
  my @sites_keep = ();
  foreach my $idx (0..($length_masked-1)) {
    my $found = 0;
    while (not $found and $original_idx < $length_original) {
      if ($col_original_alignment{$original_idx} eq $col_masked_alignment{$idx}) {
        push(@sites_keep,$original_idx);        
        $found = 1;
      }
      $original_idx++;
    }
  }
  open(OUTPUT, ">$output_file");
  print OUTPUT "@sites_keep";
  close(OUTPUT);
  
  return \@sites_keep;
}

#Given an alignment map and output file, write to file
sub write_alignment {
  my %aln = %{$_[0]};
  my $output_file = $_[1];
    
  open(OUTPUT, ">$output_file");
  foreach my $key (sort keys %aln) {
    print OUTPUT ">$key\n$aln{$key}\n";
  }
  close(OUTPUT);
}

sub run_fastsp {
  my $reference_alignment = $_[0];
  my $estimated_alignment = $_[1];  
  my $output_file = $_[2];
  my $column_file = $_[3];
  
  my $file = get_temp_file();
  my $result = undef;
  if (defined $column_file) {
    chdir("/projects/sate9/namphuon/programs/FastSP/src");
    $result = `java phylolab/alg/sp/FastSP -r $reference_alignment -e $estimated_alignment -O $column_file > $file`;
  } else {
    $result = `java -Xmx5g -jar $fastsp -r $reference_alignment -e $estimated_alignment -ml > $file`;
  }
  if (defined $output_file) {
    #open(OUTPUT, ">$output_file");
    #print OUTPUT $result;
    #close(OUTPUT);
    `mv $file $output_file`;
  } else {
    rm $file;
  }
  
  my %stats = ();
  foreach my $stat ("SPFN", "SPFP", "SP-Score", "Compression", "TC") {
    $result =~ m/$stat\s+(\d+\.\d+)/;
    $stats{$stat} = $1;
  }
  return \%stats;
}

sub run_fastsp_induced {
  my $reference_alignment = $_[0];
  my $estimated_alignment = $_[1];  
  my $output_file = $_[2];
  my $column_file = $_[3];  
  
  my $file = get_temp_file();
  if (not -e $reference_alignment or not -e $estimated_alignment or -s $estimated_alignment == 0) {
		return;
  }
  
  
  my %ref = %{read_fasta_file("$reference_alignment")};
  my %est = %{read_fasta_file("$estimated_alignment")};

  my %common = %{intersection(\%ref,\%est)};  
  my %new_ref = map {$_=>$ref{$_}} keys %common;
  my %new_est = map {$_=>$est{$_}} keys %common;
  write_alignment(\%new_ref, "$file.ref");
  write_alignment(\%new_est, "$file.est");
  run_fastsp("$file.ref","$file.est",$output_file,$column_file);
  `rm $file.*`;
}


sub run_fastsp_log {
  my $reference_alignment = $_[0];
  my $estimated_alignment = $_[1];  
  my $output_file = $_[2];
  my $optional = $_[3];
  
  if (not defined $optional) {
    $optional = "";
  }
    
  my $file = get_temp_file();
  my $results = `java -Xmx4g -jar $fastsp $optional -r $reference_alignment -e $estimated_alignment 2> $file`;
    
  $results = Phylo::trim(`grep "SPFN" $file`."");
  if ($results ne "") {
    `mv $file $output_file`;
  }
}

sub read_fastsp_log {
  my $log_file = $_[0];
  my $results = Phylo::trim(`grep Exception $log_file`."");
  my %stats = ();
  if ($results ne "") {
    return \%stats;
  }
  open(INPUT, "<$log_file");
  my $skip = 1;
  while (my $line = <INPUT>) {
    if ($line =~ m/MaxLenNoGap=/) {
      $skip = 0;
    } elsif($skip == 0) {
    } else {
      next;
    }
    if ($line =~ m/Number of shared homologies:\s+(\d+[^\s]*)/) {
      $stats{"shared_homologies"} = $1; 
    } elsif ($line =~ m/Number of homologies in the reference alignment:\s+(\d+[^\s]*)/) {
      $stats{"reference_homologies"} = $1; 
    } elsif ($line =~ m/Number of homologies in the estimated alignment:\s+(\d+[^\s]*)/) {
      $stats{"estimated_homologies"} = $1; 
    } elsif ($line =~ m/NumSeq=/) {
      $line =~ m/NumSeq=\s+(\d+)/;
      $stats{"total_sequences"} = $1;
      $line =~ m/LenRef=\s+(\d+)/;
      $stats{"length_ref"} = $1;
      $line =~ m/LenEst=\s+(\d+)/;
      $stats{"length_est"} = $1; 
    } elsif ($line =~ m/Number of correctly aligned columns:\s+(\d+[^\s]*)/) {
      $stats{"correctly_aligned_columns"} = $1; 
    } elsif ($line =~ m/Number of aligned columns in ref\. alignment:\s+(\d+[^\s]*)/) {
      $stats{"total_aligned_columns"} = $1; 
    } else {
      $line =~ m/(.*)\s+(\d+\.\d+[^\s]*)/;
      my ($key,$value) = ($1,$2);
      if (not defined $key) {
        next;
      }
      $stats{"$1"} = $2; 
    } 
  }
  close(INPUT);
  return \%stats;
}

sub read_fastsp_column_file {
  my $input_file = $_[0];
  open(INPUT, "<$input_file");
  my $line = <INPUT>;
  my %stats = ();
  while ($line = <INPUT>) {
    $line = trim($line);
    my @results = split(/,/, $line);
    $stats{$results[0]}=\@results;
    print "@results\n";
  }
  close(INPUT);
  return \%stats;
}

sub read_fastsp_result {
  my $input_file = $_[0];
  if (not -e $input_file or -s $input_file == 0) {
    return {};
  }
  my %stats = ();
  open(INPUT, $input_file);
  my $columns = 0;
  while (my $line = <INPUT>) {
    my @results = split(/\s+/, $line);
    if (scalar @results != 2 and $line !~ m/correctly aligned/) {
      next;
    }
    $stats{trim($results[0])} = trim($results[1]);
    if ($line =~ m/correctly aligned/) {
      $line =~ m/:\s+(\d+)/;
      $stats{"total_tc"} = $1;
    }
  }
  #my $results = `cat $input_file`;
  #my %stats = ();
  #foreach my $stat ("SPFN", "SPFP", "SP-Score", "Compression", "TC") {
  #  $results =~ m/$stat\s+([^\s]+)/;
  #  $stats{$stat} = $1;
  #}
  return \%stats;
}

sub write_distance_matrix {
  my @names = @{$_[0]};
  my @distanceMatrix = @{$_[1]};
  my $output_file = $_[2];
  my $format = $_[3];
  
  if (not defined $format) {
    $format = "nexus";
  }
  
  open(OUTPUT, ">$output_file");
  my $numTaxa = scalar @names;
  if ($format eq "nexus" ) {
    
    print OUTPUT ("#NEXUS\nbegin taxa;\n");
    print OUTPUT sprintf("\tdimensions ntax = %d;\n", $numTaxa);
    my $name_string = "";
    foreach my $name (@names) {
      $name_string = $name_string . " " . $name;
    } 
    print OUTPUT sprintf("\ttaxlabels %s;\n", $name_string);
    print OUTPUT ("end;\n\n");
    
    print OUTPUT ("begin distances;\n");
    print OUTPUT ("\tformat triangle=upper diag=n labels;\n");
    print OUTPUT ("\tmatrix\n");
      
    for (my $i = 0; $i < $numTaxa-1; $i++) {
      print OUTPUT sprintf("\t\t%s %s\n", $names[$i],  "@{$distanceMatrix[$i]}");
    }
    print OUTPUT sprintf("\t\t%s;\n", $names[$numTaxa-1]);
    print OUTPUT ("\nend;\n");    
  } elsif ($format eq "phylip") {
    #Build the matrix
    my @squareMatrix = ();
    for (my $i = 0; $i < $numTaxa; $i++) {
      my @distances = @{$distanceMatrix[$i]};
      my @row = ();
      for (my $j = 0; $j <= $i; $j++) {
        if ($j == $i) {
          push(@row, 0);
        } else {
          push(@row, $squareMatrix[$j]->[$i]);          
        } 
      }
      if ($i != ($numTaxa-1)) {
      foreach my $distance (@distances) {
        push(@row, $distance);
      }}
      push(@squareMatrix, [@row]);
    }
    print OUTPUT "$numTaxa\n";
    for (my $i = 0; $i < $numTaxa; $i++) {
      print OUTPUT "$names[$i]\n@{$squareMatrix[$i]}\n";
    }    
  }
  close(OUTPUT);
}

sub run_nj_fastme {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $base = $_[2];
  
  if (not defined $base) {
    $base = 2;
  }
  
   my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
   while (-e $temp_file) {
     $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
   }

  my ($a, $b) = Phylo::calc_distance_matrix($input_file, $base);
  Phylo::write_distance_matrix($a, $b, $temp_file, "phylip");
  my_cmd("$bin/fastme -i $temp_file -o $output_file");  
  Phylo::my_cmd("$perl -p -i -e \"s/\\)[^,;\)\(]*/\)/g\" $output_file");
  my_cmd("rm $temp_file");
}

sub run_nj_paup {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $base = $_[2];
  
  if (not defined $base) {
    $base = 2;
  }
  
   my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
   while (-e $temp_file) {
     $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
   }

  my ($a, $b) = Phylo::calc_distance_matrix($input_file, $base);
  Phylo::write_distance_matrix($a, $b, $temp_file, "nexus");
  open(OUTPUT, ">>$temp_file");
  
  print OUTPUT "begin paup;\n\tlog file = $temp_file.log replace;\n\tdset distance=user;\n\tnj;\n\tsavetrees file=$output_file format=phylip brlens=yes replace;\nend;";
  close(OUTPUT);  
  my_cmd("$bin/paup -n $temp_file");  
  #Phylo::convert_paup_newick("$temp_file.tree", $output_file);
  #Phylo::my_cmd("$perl -p -i -e \"s/\\[.*\\]//g\" $output_file");  
  my_cmd("rm $temp_file $temp_file.log");
}

sub combine_alignment_partitions_protein {
  my %genes = %{$_[0]};
  my %partitions = %{$_[1]};
  my %partition_model = %{$_[2]};
  my $alignment_file = $_[3];
  my $partition_file = $_[4];
  
  open(OUTPUT, ">$partition_file");
  my %concat = ();
  my $idx = 0;
  my $length = 1;
  foreach my $partition (keys %partitions) {    
    local $" = "_";
    my $name = "@{$partitions{$partition}}";
    print OUTPUT "$partition_model{$partition}, combined_$idx=$length-";
    $idx++;
    foreach my $gene (@{$partitions{$partition}}) {
      my $curr_length = 0;
      foreach my $seq (keys %{$genes{$gene}}) {
        $concat{$seq}.=$genes{$gene}->{$seq};
        $curr_length = length($genes{$gene}->{$seq});
      }
      $length+=$curr_length;      
    }
    print OUTPUT "$length\n";
    $length++;    
  }
  close(OUTPUT);
  Phylo::write_alignment(\%concat, $alignment_file);
}

sub read_raxml_info_rates {
  my $input_file = $_[0];
  my $line = Phylo::trim(`grep "alpha\\[0\\]" $input_file | head -n1`."");
  if ($line eq "") {
    return undef;
  }
  $line =~ m/alpha\[0\]: (\d+\.\d+) rates\[0\] ac ag at cg ct gt: (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+)/;
  my @rates = ($1, $2, $3, $4, $5, $6, $7);
  return \@rates;
}

sub split_out_by_codon {
  my $input_file = $_[0];
  my $output_prefix = $_[1];
  
  my %sequences = %{Phylo::read_fasta_file($input_file)};
  my @positions = ({},{},{});
  local $" = "";
  foreach my $key (keys %sequences) {
    my @results = split(//,$sequences{$key});
    my @codon_1 = map {$results[$_ * 3]} (0.. (scalar @results / 3)-1);
    my @codon_2 = map {$results[$_ * 3 + 1]} (0.. (scalar @results / 3)-1);
    my @codon_3 = map {$results[$_ * 3 + 2]} (0.. (scalar @results / 3)-1);
    $positions[0]->{$key} = "@codon_1";
    $positions[1]->{$key} = "@codon_2";
    $positions[2]->{$key} = "@codon_3";
  }
  foreach my $idx (1..3) {
    Phylo::write_alignment($positions[$idx-1], "$output_prefix.$idx.fasta");
  }
}

sub remove_third_codon_position {
  my $input_phylip = $_[0];
  my $output_phylip = $_[1];
  my $input_partition = $_[2];
  my $output_partition = $_[3];
  
  open(IALIGN, "<$input_phylip");
  open(OALIGN, ">$output_phylip");
  my $line = <IALIGN>;
  my @results = split(/\s+/, $line);
  print OALIGN sprintf("%d %d\n", $results[0], ($results[1]/3)*2);
  while (my $line = <IALIGN>) {
    my @first = split(/\s+/, $line);
    print OALIGN $first[0] . "\n";    
    $line = $first[1];
    @results = split(//,$line);
    my $new_align = "";
    for(my $idx = 0; $idx < scalar(@results); $idx++) {
      if ($idx % 3 == 2) {
        #$new_align = $new_align . "-";
        next;
      }
      $new_align = $new_align . $results[$idx];
    }
    print OALIGN $new_align . "\n";    
  }
  close(IALIGN);    
  close(OALIGN);  

  if (defined $input_partition) {
    open(IPART, "<$input_partition");
    open(OPART, ">$output_partition");
    while (my $line = <IPART>) {
      #DNA,4048=1-1659
      $line = Phylo::trim($line);
      $line =~ m/=\s*(\d+)-(\d+)/;
      my ($start, $end) = ($1, $2);
      my $replace = "$start-$end";
      $start = int($start/3*2)+1;
      $end = int($end/3*2);
      $line=~s/$replace/$start-$end/;
      print OPART $line . "\n";    
    }
    close(IPART);    
    close(OPART);  
  }
}

sub remove_all_gap_columns_raxml {
  my $input_phylip = rel2abs($_[0]);
  my $output_phylip = rel2abs($_[1]);
  my $model = $_[2];
  my $input_partition = $_[3];
  my $output_partition = $_[4];
  my $temp_dir = get_temp_file();
  make_directory($temp_dir);
  my $cmd = "-m $model -f c -n test -s $input_phylip -w $temp_dir/";  
  if (defined $input_partition) {
    $input_partition = rel2abs($input_partition);
    $output_partition = rel2abs($output_partition);
    $cmd = $cmd . " -q $input_partition ";
  }
  `raxmlHPC-git-June2-gcc $cmd`;
  if (-e "$input_phylip.reduced") {
    `mv $input_phylip.reduced $output_phylip`;
  }
  
  if (defined $input_partition and -e "$input_partition.reduced") {
    `mv $input_partition.reduced $output_partition`;
  }  
  `rm $temp_dir -rf`;
}

sub remove_all_gap_columns {
    my $alignment_file = $_[0];
    my $output_file = $_[1];
    my $format = $_[2];
    if (not defined $format) {
    $format = "fasta";
    }
    `/lusr/bin/python /projects/sate3/namphuon/taxon2/metagen/scripts/python/removeGaps.py $alignment_file $output_file $format`;
}

sub calc_average_pairwise_distances {
  my $alignment_file = $_[0];  
  my %average = ();
  
  my $temp_file = get_temp_file();
  calc_pairwise_distances($alignment_file, $temp_file);
  open(INPUT, "<$temp_file");
  while (my $line = trim("".<INPUT>)) {
    my @results = split(/,/, $line);
    if (not defined $average{$results[0]}) {
      $average{$results[0]} = [0,0];
    }
    $average{$results[0]}->[0]+=$results[2];
    $average{$results[0]}->[1]++;
    
    $average{$results[1]}->[0]+=$results[2];
    $average{$results[1]}->[1]++;
  }
  foreach my $key (keys %average) {
    $average{$key}->[0]/=$average{$key}->[1];
  }
  `rm $temp_file`;
  return \%average;
}

sub fast_calc_pairwise_distances {
  my $alignment_file = $_[0];
  my $output_file = $_[1];
  
  my %alignment = %{Phylo::read_fasta_file($alignment_file)};  
  my @names = keys %alignment;
  
  #Convert all characters to lowercase, dashes to character with ASCII value greater than 300
  #XOR this character with any other A-Za-z
  my $special_char = chr(300);
  foreach my $name (@names) {
    $alignment{$name} = lc $alignment{$name};
    $alignment{$name} =~ s/-/$special_char/g;
  }  
  
  my $num_sequences = (scalar @names)-1;
  open(OUTPUT,">$output_file");
  my $lines = 0;
  foreach my $idx_x (0..$num_sequences) {
    foreach my $idx_y (($idx_x+1)..$num_sequences) {
      my $xor = ($alignment{$names[$idx_x]} ^ $alignment{$names[$idx_y]});
      my $or = ($alignment{$names[$idx_x]} | $alignment{$names[$idx_y]});
      my $distance = $xor =~ tr/\001-\255//;
      my $cols = $or =~ tr/\001-\255//;
      print OUTPUT sprintf("$names[$idx_x],$names[$idx_y],%d,%d\n", $distance,$cols);
      if ($lines % 5000 == 0) {
        print STDERR "$lines\n";
      }
      $lines++;      
    }
  }
  close(OUTPUT);
}
 
sub calc_pairwise_distances {
  my $alignment_file = $_[0];
  my $output_file = $_[1];
  
  my %alignment = %{Phylo::read_fasta_file($alignment_file)};
  my @names = keys %alignment;
  my $num_sequences = (scalar @names)-1;
  open(OUTPUT,">$output_file");
  foreach my $idx_x (0..$num_sequences) {
    foreach my $idx_y (($idx_x+1)..$num_sequences) {
      my $hd = hamming_distance($alignment{$names[$idx_x]}, $alignment{$names[$idx_y]});
      if ($hd == -1) {
        next;
      }
      print OUTPUT "$names[$idx_x],$names[$idx_y],$hd\n";
    }
  }
  close(OUTPUT);
}


sub hamming_distance {
  my $str_1 = $_[0];
  my $str_2 = $_[1];
  
  my $len = length ($str_2);
  my $num_mismatch = 0;
  my $non_gap = 0;
  
  for (my $i=0; $i<$len; $i++) {
    if (substr($str_1, $i, 1) eq "-" or substr($str_2, $i, 1) eq "-") {
      next;
    }
    $non_gap++;
    ++$num_mismatch if lc substr($str_1, $i, 1) ne lc substr($str_2, $i, 1);
  }
  if ($non_gap eq 0) {
    return -1;
  }
  return $num_mismatch/$non_gap;
}

sub calc_stat {
    my $alignment_file = $_[0];
    my $out_file = $_[1];
    my $curr_dir = `pwd`;
    
    my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
    while (-e $temp_file) {
  $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
    }  
    
    remove_all_gap_columns($alignment_file, $temp_file);
    if ($machine == 0) {      
  my_cmd("java -cp \"$bin/distance\" AlignmentStatistics $temp_file $out_file");    
    } else {
  my_cmd("java -cp \"$bin/distanceV\" AlignmentStatistics $temp_file $out_file");  
    }
    `rm $temp_file`;
}

sub write_pairwise_distances {
  my $backbone = $_[0];
  my $output = $_[1];  
  
  my %stats = %{Phylo::calc_average_pairwise_distances($backbone)};
  my @averages = map {$stats{$_}->[0]} keys %stats;
  my ($mean, $std, $counts) = Phylo::mean_std(\@averages);
  
  open(OUTPUT, ">$output");
  foreach (sort { ($stats{$a}->[0] cmp $stats{$b}->[0])} keys %stats) {
    print OUTPUT "$_,". $stats{$_}->[0] .",". ($stats{$_}->[0]-$mean)/$std."\n";
  }
  close(OUTPUT);
}

sub read_alignment_stats {
  my $file = $_[0];
  my %map = ();
  
  open(INPUT, "<$file");
  while (my $line = <INPUT>) {
    my @results = split(/\|/, $line);
    $results[0] = trim($results[0]);
    $results[1] = trim($results[1]);
    $map{$results[0]} = $results[1];
  }
  close(INPUT);
  return \%map;
}

sub calculate_mrp_score {
   #Tree file
   my $input_file = $_[0];
   #MRP Matrix
   my $mrp_file = $_[1];  
  #Score result
   my $output_file = $_[2];
 
   my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
   while (-e $temp_file) {
     $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
   }  
   convert_newick_to_mrp_tnt($input_file, $temp_file, "$temp_file.map", "map");
   convert_mrp_to_alignment($temp_file, "$temp_file.fasta", "fasta");  
   estimate_ml_tree_fasta_fastree("$temp_file.fasta", "$output_file.temp", "-nt -cat 1");
   remap_newick("$output_file.temp", $output_file, "$temp_file.map", "remap");
  my_cmd("rm $output_file.temp $temp_file.map $temp_file $temp_file.fasta");
}

sub generate_mrp_tnt {
   my $source_file = $_[0];
   my $mrp_file = $_[1];  
  my $map_file = $_[2];  

   my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
   while (-e $temp_file) {
     $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
   }  
  convert_newick_to_mrp_tnt($source_file, $temp_file, "$temp_file.map", "map");
  Phylo::my_cmd("mv $temp_file $mrp_file");
  Phylo::my_cmd("mv $temp_file.map $map_file");
  
  #TODO:REMAP
}

sub get_alignment_stats {
  my $file_name = $_[0];
  my $result = trim(my_cmd("$python $script_directory/mask.py -f $file_name -l"));  
  return $result;
}

sub measure_memory {  
  my $file_name = $_[0];
  my $pid = fork();
  if ($pid == 0) {    
    if (-e $file_name) {
      my_cmd("rm $file_name");
    }
    my_cmd("echo \"$$\n\" > $file_name");    
    my $idx = 0;
    while ($idx < 100000) {
      #print "Working $idx\n";
      my_cmd("top -b -n 1 -u namphuon >> $file_name");
      sleep(20);
      $idx++;
    }
    exit(0);
  } else {
    return $pid;
  }
}

sub end_memory {
  my $pid = $_[0];
  my_cmd("kill -9 $pid");
}

sub call {
  shift;
  my $name = shift;
  &{$name}(@_);
}

sub foo {
  foreach my $i (@_) {
    print "Foo:". $i . "\n";
  }
}

sub count_nt_pairs {
  my $fasta_file = $_[0];
  my $string = my_cmd("$python $script_directory/mask.py -z -f $fasta_file");
  return trim($string);
}

sub get_sp_score {
  my $score_file = $_[0];
  my $score = 0;

  open(SCORE,"<$score_file");
  my $line = <SCORE>;
  while (defined($line)) {    
    if($line =~ m/\s*#MEAN_RES_PAIR_SCORE\s*(\d+\.\d+)\s*#MEAN_COL_SCORE\s*(\d+\.\d+)/) {
      $score = $1;
      last;
      
    }
    $line = <SCORE>;
  }
   close(SCORE);
  return $score;
}

sub get_bali_score {
  my $ref_file = $_[0];
  my $test_file = $_[1];

  my $test_name = "test" .  int(rand()*100000) . ".msf";
  my $ref_name = "ref" .  int(rand()*100000) . ".msf";
  convert_fasta_to_msf($ref_file, $ref_name);
  convert_fasta_to_msf($test_file, $test_name);
  my $result = my_cmd("$bin/bali_score $ref_name $test_name");
  my_cmd("rm $ref_name $test_name");
  return $result;
}

sub get_bali_score_to_lobster {
  my $ref_file = $_[0];
  my $test_file = $_[1];
  my $output_file = $_[2];

  my $test_name = "test" .  int(rand()*100000) . ".msf";
  my $ref_name = "ref" .  int(rand()*100000) . ".msf";
  my $temp_file = "temp" .  int(rand()*100000);
  convert_fasta_to_msf($ref_file, $ref_name);
  convert_fasta_to_msf($test_file, $test_name);

  my $result = my_cmd("$bin/bali_score $ref_name $test_name > $temp_file");
  my_cmd("rm $ref_name $test_name");

  #Now convert it to a lobster score since I fucked up
  #TODO: Fix this
  open(FILE, "<$temp_file");
  my $line = <FILE>;
  my $sp_score = undef;
  while (defined($line)) {
    $line =~ m/SP\s*score\s*=\s*([\d]+\.[\d]+)/;
    if (defined($1)) {
      $sp_score = $1;
      last;
    }
    $line = <FILE>;
  }
  if (defined($sp_score)) {
    my_cmd("echo \"SP=$sp_score;\" >> $temp_file");
    my_cmd("mv $temp_file $output_file");
  } else {
    my_cmd("rm $temp_file");
  }   
  close(FILE);
}

sub get_lobster_score {
  my $ref_file = $_[0];
  my $test_file = $_[1];
  my $output_file = $_[2];

  my $test_name = "test" .  int(rand()*100000);
  my $ref_name = "ref" .  int(rand()*100000);
  my $output_name = "output". int(rand()*10000);
  to_upper($ref_file, $ref_name);  
  to_upper($test_file, $test_name);
  #print `pwd`;
  #print ("$bin/lobster -score $test_name -ref $ref_name > $output_name");
  my_cmd("$bin/lobster -score $test_name -ref $ref_name > $output_name");
  my_cmd("mv $output_name $output_file");
  my_cmd("rm $ref_name $test_name");
}

sub run_gappy_remove {
  my $input_file = $_[0];
  my $threshold = $_[1]; #% of alignment to remove
  my $output = $_[2];
  
  my $temp_file = get_temp_file();
  convert_fasta_to_phylip($input_file, "$temp_file.phylip");
  `$perl /projects/sate3/namphuon/scripts/prunekpercent_5Apr.pl $temp_file.phylip $threshold $temp_file.output`;
  phylip_to_fasta("$temp_file.output", "$temp_file.fasta", 1);
  if (check_alignment("$temp_file.fasta") eq "") {
    `mv $temp_file.fasta $output`;
  }
  if (defined $temp_file and $temp_file ne "") {
    `rm $temp_file*`;
  }
}


sub remove_gap_msa {
  my $input_alignment = $_[0];
  my $threshold = $_[1];
  my $output_alignment = $_[2];
  my $output_character = $_[3];
  
  if (not defined $output_character) {
    $output_character = "";
  }
  
  my %alignment = %{Phylo::read_fasta_file($input_alignment)};
  my @gap_count = ();
  my $taxa = scalar keys %alignment;
  my @names = keys %alignment;
  my $length = length $alignment{$names[0]};
  my $min_threshold =  ceil($taxa*$threshold);
  my $idx = 0;
  while ($idx < $length) {
    push(@gap_count, 0);
    $idx++;
  }
  
  foreach my $key (keys %alignment) {
    my @results = split(//, $alignment{$key});
    $idx = 0;
    while ($idx < $length) {
      if ($results[$idx] eq '-') {
        $gap_count[$idx]++;
      }
      $idx++;
    }    
  }
  
  foreach my $key (keys %alignment) {
    my @results = split(//, $alignment{$key});    
    $idx = 0;
    while ($idx < $length) {
      if ($gap_count[$idx] >= $min_threshold) {
        $results[$idx] = $output_character;
      }
      $idx++;
    }
    local $" = "";
    $alignment{$key} = "@results";
  }
  write_alignment(\%alignment, $output_alignment);  
}

#Removes all insertions into backbone alignment
sub mask_extended_alignment {
  my $backbone_alignment_file = $_[0];
  my $extended_alignment_file = $_[1];
  my $output_alignment_file = $_[2];
  
  my %backbone_alignment = %{Phylo::read_fasta_file($backbone_alignment_file)};
  my %extended_alignment = %{Phylo::read_fasta_file($extended_alignment_file)};
  
  my @full_length = keys %backbone_alignment;
  foreach my $key (keys %extended_alignment) {
    my @results = split(//,$extended_alignment{$key});
    $extended_alignment{$key} = \@results;
  }
  
  my @full_length_seq = @{$extended_alignment{$full_length[0]}};
  my @to_mask = ();
  for (my $idx = 0; $idx < scalar @full_length_seq; $idx++) {
    if ($full_length_seq[$idx] eq "-") {
      my $remove = 1;
      foreach my $key (@full_length) {
        if ($extended_alignment{$key}->[$idx] ne "-") {
          $remove = 0;
          last;
        }
      }
      if ($remove == 1) {
        push(@to_mask, $idx);
      }
    }
  }
  
  %extended_alignment = %{Phylo::read_fasta_file($extended_alignment_file)};
  my %masked = %{mask_sites(\%extended_alignment, \@to_mask)};
  Phylo::write_alignment(\%masked, $output_alignment_file);
}

sub mask_sites_gappy {
  my $input = $_[0];
  my $output = $_[1];
  my $threshold = $_[2];
  
  my %input_alignment = %{Phylo::read_fasta_file($input)};   
  my @column_count = ();
  foreach my $key (keys %input_alignment) {
    my @results = split(//,$input_alignment{$key});
    $input_alignment{$key} = \@results;
    foreach my $idx (0..((scalar @results) - 1)) {
      if ($results[$idx] ne "-") {
        $column_count[$idx]++;
      }
    }
  }
  
  foreach my $key (keys %input_alignment) {
    my @results = @{$input_alignment{$key}};
    my $seq = "";
    my $counter = 0;
    foreach my $idx (0..(@results-1)) {
      if ($column_count[$idx] > $threshold) {
        $seq = $seq . $results[$idx];
      }
    }
    $input_alignment{$key} = $seq;
  }
  Phylo::write_alignment(\%input_alignment, $output);
}

sub mask_sites {
  my %alignment = %{$_[0]};
  my @sites_array = @{$_[1]}; #List of sites to remove
  
  my %sites = map {$_ => $_} @sites_array;  
  my %masked_alignment = ();
  foreach my $key (keys %alignment) {
    my @results = split(//, $alignment{$key});
    my $seq = "";
    my $counter = 0;
    foreach my $idx (0..(@results-1)) {
      if (not defined $sites{$idx}) {
        $seq = $seq . $results[$idx];
      }
    }
    $masked_alignment{$key} = $seq;
  }
  return \%masked_alignment;
}

sub run_aliscore {
  my $input_file = rel2abs($_[0]);
  my $output_file = $_[1];
  my $options = $_[2];
  
  my $temp_file = get_temp_file();
  my ($filename, $basename, $empty) = fileparse($input_file); 
  make_directory($temp_file);
  if (not defined $options) {
    $options = "";
  }
  `ln -s $input_file $temp_file/input.fasta`;
  `$perl /projects/sate3/namphuon/bin/aliscore/Aliscore.02.2.pl $options -i $temp_file/input.fasta`;
  
  #Read the sites to be masked
  my $mask_file = "$temp_file/input.fasta_List_l_all.txt";
  if (not -e $mask_file) {
    $mask_file = "$temp_file/input.fasta_List_random.txt";
  }
  open(INPUT, "<$mask_file");
  my $line = <INPUT>;
  $line = trim($line);
  my @results = split(/\s+/, $line);
  foreach my $idx (0..(@results-1)) {
    $results[$idx]--;
  }
  close(INPUT);
  
  my $alignment = read_fasta_file($input_file, 1);
  my $masked_alignment = mask_sites($alignment, \@results);
  write_alignment($masked_alignment, $output_file);
  print STDERR "sites: @results\n";
  `rm $temp_file -rf`;
  return 1;
}

sub score_zorro {
  my $input_file = $_[0];
  my $output_file = $_[1];
  `/projects/sate9/namphuon/programs/zorro/zorro_linux_x86_64 $input_file > $output_file`;
}


sub run_zorro {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $cut_off = $_[2];
  if (not defined $cut_off) {
    $cut_off = 4;
  }
  my $temp_file = get_temp_file();
  `/projects/sate9/namphuon/programs/zorro/zorro_linux_x86_64 -sample $input_file > $temp_file`;
  
  #Read the sites to be masked
  open(INPUT, "<$temp_file");
  my $counter = 0;
  my @to_mask = ();
  while (my $line = <INPUT>) {
    $line = trim($line);
    if ($line <= $cut_off) {
      push(@to_mask, $counter);
    }
    $counter++;
  }
  close(INPUT);
  
  my $alignment = read_fasta_file($input_file, 1);
  my $masked_alignment = mask_sites($alignment, \@to_mask);
  write_alignment($masked_alignment, $output_file);
  print STDERR "sites: @to_mask\n";
  `mv $temp_file $output_file.score`;
}

sub gblock_msa {
  my $msa_file = $_[0];
  my $output_file = $_[1];
  my $options = $_[2];
  my $taxa = $_[3];
  
  my $temp_file = get_random_name();
  my_cmd("cp $msa_file $temp_file");    
  

  if (not defined($taxa)) {
    $taxa = trim(my_cmd("$python $script_directory/mask.py -x -f $msa_file"));
  }
  
  if (defined($options) and ($options eq "relaxed")) {
    $options = "-t=d -b2=" . (int($taxa)*.5625) . " -b3=10 -b4=5 -b5=H";
  } elsif (defined($options) and ($options eq "strict")) {
    $options = "-t=d";
  } elsif (not defined($options)) {
    $options = "-t=d";
  }
  
  my_cmd("$gblock $temp_file $options");
  my_cmd("mv $temp_file". "-gb $output_file");
  my_cmd("rm $temp_file". "-gb.htm");
  my_cmd("rm $temp_file");
}


sub run_noisy {
    my $ref_file = $_[0];
    my $output_file = $_[1];
    my $gap_option = $_[2];
    my $gap_option_str = "";
    if (defined($gap_option)) {
      $gap_option_str = " $gap_option ";
    }
    my $temp_file = get_random_name();
    my_cmd("cp $ref_file $temp_file.fasta");    
    if ($machine) {
  my_cmd("$bin/noisyV $gap_option_str $temp_file.fasta");
    } else {
  my_cmd("$bin/noisy $gap_option_str $temp_file.fasta");
    } 
    
    #my_cmd("rm temp.file temp_sta.gr temp_typ.eps");
    fix_noisy_fasta("$temp_file" . "_out.fas", $output_file);    
    my_cmd("rm $temp_file" . "_*");
    my_cmd("rm $temp_file.fasta");
}

sub run_bmge {
    my $ref_file = $_[0];
    my $output_file = $_[1];
    my $type = $_[2];
    my $type_str = "DNA";
    if (defined($type)) {
      $type_str = $type;
    }

    my $temp_file = get_random_name();
    my_cmd("java -jar $bin/BMGE.jar -i $ref_file -t $type_str -of $temp_file");
    my_cmd("mv $temp_file $output_file");
}

sub run_guidance {
  print "$perl www/Guidance/guidance.pl --seqFile /projects/sate3/namphuon/experiments/length_500/R0/RAXML_CLUSTALW/initial.alignment --msaProgram CLUSTALW --seqType nuc --outDir temp";
}

sub run_trimal {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $command = $_[2];
  if (not defined $command) {
    $command = "-automated1 ";
  }
  my $temp_file = get_random_name();
  if ($machine == 0) {
    my_cmd("$bin/trimal -in $input_file -out $temp_file $command");
  } else {
    my_cmd("$bin/trimalV -in $input_file -out $temp_file $command");
  } 
  my_cmd("$perl -p -i -e \"s/\\s+\\d+\\s+\\w+//g\" $temp_file");
  my_cmd("mv $temp_file $output_file");
}

sub verify_alignment {
  my $fasta_file = $_[0];
  my $result = my_cmd("$python $script_directory/mask.py -v -f $fasta_file");
  if (trim($result) eq trim("Invalid")) {
    return 0;
  }
  return 1;
}

sub verify_tree {
  my $tree_file = $_[0];
  my $result = my_cmd("$python $script_directory/phylo.py -t -i $tree_file");
  if (trim($result) eq trim("Invalid")) {    
    return 0;
  }
  if (-s $tree_file == 0) {
    return 0;
  }
  return 1;
}

sub unalign_sequence {
  my $input = $_[0];
  my $output = $_[1];
  my_cmd("cp $input $output");
  my_cmd("$perl -p -i -e \"s/-/\"\"/g\" $_[1]");
}


sub align_clustalw {
  my $input = $_[0];
  my $output = $_[1];
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }
  
  if ($machine == 0) {
    my_cmd("$bin/clustalw2 -align -infile=$input -outfile=$temp_file -output=fasta");
  } else {
    my_cmd("$bin/clustalw2 -align -infile=$input -outfile=$temp_file -output=fasta");
  } 
  if (-e $temp_file) {
    `mv $temp_file $output`;
  }  
}

sub align_probcons {
  my $input = $_[0];
  my $output = $_[1];

  my_cmd("$bin/nwaffine -d $input -m pw.mat -0 0 -1 1 -2 1 -n 1 -u 1 -f upgma.tree");
  my_cmd("$bin/paup -u -n pw.mat");
  my_cmd("python $script_directory/phylo.py -c -i upgma.tree -o upgma.tree.newick");
  my_cmd("$perl -p -i -e \"s/\\[.*\\]//g\" upgma.tree.newick");
  my_cmd("$bin/probcons -g upgma.tree.newick $input > $output");
  my_cmd("rm upgma.tree pw.mat upgma.tree.newick");  
}

sub fix_trees {
  my $input = $_[0];
  my $output = $_[1];
  
  open(TREE_FILE, "<$input");
  my $line = <TREE_FILE>;
  my $tree = "";
  while (defined $line) {
    $tree = $tree.$line;
    my $line = <TREE_FILE>;
  }  
  close(TREE_FILE);
}

sub fix_source_trees {
    my $input = $_[0];
    my $output = $_[1];
    
    `cp $input $output`;
    `$perl -p -i -e "s/\\s+//g" $output`;
    `$perl -p -i -e "s/'//g" $output`;
    `$perl -p -i -e "s/;/\n;/g" $output`;
}
######################################
# R scripts
######################################
#Runs scatterplot
sub scatterplot_test {
  my @sample_one = @{$_[0]};
  my @sample_two = @{$_[1]};
  my $file = $_[2];
  my $extra_ref = $_[3];
  my %extra;
  if (not defined $extra_ref) {
    %extra = ();
  } else {
    %extra = %{$extra_ref};
  } 
  my $extra_string = "";
  foreach my $key (keys %extra) {
    $extra_string = $extra_string . ", $key=\"$extra{$key}\"";
  }
  
  #Check same length
  my $length = scalar(@sample_one);
  if ($length != scalar(@sample_two)) {
    return -1;
  }
  my @valid_one = ();
  my @valid_two = ();
  
  #Keep only those that are valid
  foreach my $idx (0..($length-1)) {
    if (looks_like_number($sample_one[$idx]) and looks_like_number($sample_two[$idx])) {
      push(@valid_one, $sample_one[$idx]);
      push(@valid_two, $sample_two[$idx]);
    } 
  }
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }
  $temp_file = "temp";
  
  open(RFILE, ">$temp_file.r");
  print RFILE r_string_array("x", \@valid_one) . "\n";
  print RFILE r_string_array("y", \@valid_two) . "\n";
  print RFILE "postscript(file=\"$file\", paper=\"special\", width=10, height=10, horizontal=FALSE)\n";
  print RFILE "plot(x, y $extra_string)\n";
  close(RFILE);
  my $results = `cat $temp_file.r | R --vanilla`;  
  #`rm $temp_file.r`;
  
  return 1;
}

#Runs paired_t_test using R, returns p value
sub spearman_test {
  my @sample_one = @{$_[0]};
  my @sample_two = @{$_[1]};
  my $method = "pearson";
  my $name = "cor";
  if (defined $_[2] and ($_[2] eq "spearman")) {
    $name = "rho";
    $method = "spearman";
  }
  
  
  #Check same length
  my $length = scalar(@sample_one);
  if ($length != scalar(@sample_two)) {
    return -1;
  }
  my @valid_one = ();
  my @valid_two = ();
  
  #Keep only those that are valid
  foreach my $idx (0..($length-1)) {
    if (looks_like_number($sample_one[$idx]) and looks_like_number($sample_two[$idx])) {
      push(@valid_one, $sample_one[$idx]);
      push(@valid_two, $sample_two[$idx]);
    } 
  }
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }  
  #$temp_file = "temp";
  open(RFILE, ">$temp_file.r");
  print RFILE r_string_array("x", \@valid_one) . "\n";
  print RFILE r_string_array("y", \@valid_two) . "\n";
  print RFILE "cor.test(x,y, method=\"$method\")\n";
  #print `cat $temp_file.r`;  
  close(RFILE);
  my $results = `cat $temp_file.r | R --vanilla`;  
  `rm $temp_file.r`;
  $results =~ m/p-value\s+([^\s])\s+([^\s]+)/;
  my $pvalue = $2;
  $results =~ m/sample estimates:\s+$name\s+([^\s]+)\s+/;
  my $rho = $1;
  #print $results . "\n";  
  return ($rho, $pvalue);
}

#Runs paired_t_test using R, returns p value
sub paired_t_test {
  my @sample_one = @{$_[0]};
  my @sample_two = @{$_[1]};
  
  #Check same length
  my $length = scalar(@sample_one);
  if ($length != scalar(@sample_two)) {
    return -1;
  }
  my @valid_one = ();
  my @valid_two = ();
  
  #Keep only those that are valid
  foreach my $idx (0..($length-1)) {
    if (looks_like_number($sample_one[$idx]) and looks_like_number($sample_two[$idx])) {
      push(@valid_one, $sample_one[$idx]);
      push(@valid_two, $sample_two[$idx]);
    } 
  }
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }  
  open(RFILE, ">$temp_file.r");
  print RFILE r_string_array("x", \@valid_one) . "\n";
  print RFILE r_string_array("y", \@valid_two) . "\n";
  print RFILE "t.test(x,y=y, paired=TRUE)\n";
  #print `cat $temp_file.r`;
  close(RFILE);
  my $results = `cat $temp_file.r | R --vanilla`;  
  `rm $temp_file.r`;
  $results =~ m/p-value\s+([^\s])\s+([^\s]+)/;
  my $pvalue = $2;
  $results =~ m/mean of the differences\s+([^\s]+)\s+/;
  my $mean = $1;
  print $results . "\n";
  return ($mean, $pvalue);
}

sub r_string_array {
  my $name = $_[0];
  my @array = @{$_[1]};
  local $" = ',';
  my $str = "$name = c(@array)";
  return $str;
}

sub randomly_refine {
  my $input = $_[0];
  my $output = $_[1];
  
   my $temp_file = get_temp_file();
  `cp $input $temp_file`;  
  `perl -p -i -e "s/\\)\\d+\.\\d+/\\)/g" $temp_file`;
  
  my_cmd("$bin/random_refine_malloc -t $temp_file > $temp_file.done");  
  `mv $temp_file.done $output`;
  `rm $temp_file`;
  #my_cmd("cp $input $temp_file.temp");
  #my_cmd("$perl -p -i -e \"s/\\s+//g\" $temp_file.temp");
  #remove_weights("$temp_file.temp", "$temp_file.temp_weightless");
  #remap_newick("$temp_file.temp_weightless", "$temp_file.temp_weightless.mapped", "$temp_file.map", "map");  
  #my_cmd("$bin/random_refine_malloc -t $temp_file.temp_weightless.mapped > $temp_file.done");
  #remap_newick("$temp_file.done", "$output", "$temp_file.map", "remap");  
  #my_cmd("rm $temp_file.done $temp_file.temp $temp_file.temp_weightless $temp_file.map $temp_file.temp_weightless.mapped");  
}

sub align_pynast {
  my $backbone_alignment_file = $_[0];
  my $query_sequences = $_[1];
  my $output_file = $_[2];
  my $options = $_[3];
  #-m clustal mafft
  my $temp_name = get_temp_file();
  if (not defined $options ) {
    $options = "-m mafft -l 100 -p 10";
  }
  
  #Options   
  `/lusr/bin/python //projects/sate7/tools/pynast/PyNAST-1.1/scripts/pynast -i $query_sequences -t $backbone_alignment_file -a $temp_name.fasta -g $temp_name.log -f $temp_name.failure $options`;
  
  if (-e "$temp_name.fasta") {
   my %sequences = %{Phylo::read_fasta_file("$temp_name.fasta", 0)};
   my %renamed_alignment = ();
   foreach my $key (keys %sequences) {
     my @results = split(/\s+/,$key);
     $renamed_alignment{$results[0]} = $sequences{$key};
     
   }
   Phylo::write_alignment(\%renamed_alignment, "$temp_name.fasta");
   if (-e "$temp_name.fasta") {
    `mv $temp_name.fasta $output_file`;
   } 
   if (-e "$temp_name.log") {
    `mv $temp_name.log $output_file.log`;
   }
   if (-e "$temp_name.failure") {
    `mv $temp_name.failure $output_file.failure`;
   }
  }
  my $results = `ls $temp_name*`;
  if (trim($results) ne "") {
    `rm $temp_name*`; 
  }  
}

sub align_prank {
  my $input = $_[0];
  my $mafft_tree = $_[1];
  my $output = $_[2];

  my $temp_file = get_temp_file();
  my_cmd("$bin/prank -o=$temp_file -d=$input -t=$mafft_tree +F -matinitsize=5 -uselogs");
  my_cmd("mv $temp_file.best.fas $output");
}

sub generate_alignment {
  my $input = $_[0];
  my $output = $_[1];
  my $method = $_[2];
  my $options = $_[3];
  
  my $temp_file = get_temp_file();
  if ($method eq "mafft") {    
    `/projects/sate8/namphuon/tools/mafft/mafft-6.902-with-extensions/core/mafft --quiet --ep 0.123 $options $input > $temp_file`;
  } elsif ($method eq "opal") {
    `java -Xmx4G -jar /projects/sate3/namphuon/bin/opal.jar --in $input --out $temp_file`;
  } elsif ($method eq "clustalw") {
    `/projects/sate3/namphuon/bin/clustalw-64bit -align -infile=$input -outfile=$temp_file -output=fasta $options`;
  } elsif ($method eq "clustalo") {
    `/projects/sate9/namphuon/programs/clustalomega/bin/bin/clustalo $options -i $input -o $output`;
  } elsif ($method eq "muscle") {
    my $total = Phylo::trim(`grep ">" $input |wc -l`."");
    my $extra = "";
    if ($total > 3000) {
      $extra = "-maxiters 2";
    } else {
      $extra = "-maxiters 1000"
    }
    `/projects/sate7/tools/bin/muscle3.8.31_i86linux64 $extra -in $input -out $temp_file`;
  }
  if (-e $temp_file and -s $temp_file != 0) {
    `mv $temp_file $output`;
  } else {
    `rm $temp_file`;
  }
}

sub align_mafft {
  my $input = $_[0];
  my $output = $_[1];
  my $parttree = $_[2];

  my $temp_file = get_temp_file();
  
  if (defined $parttree) {
    my_cmd("/projects/sate7/tools/mafft-6.864-with-extensions/bin/mafft --ep 0.123 --parttree --retree 2 --quiet --partsize 1000 $input > $temp_file");
  } else {
    #my_cmd("$bin/mafft --localpair --maxiterate 1000 --quiet $input > $output");
    my_cmd("/projects/sate7/tools/mafft-6.864-with-extensions/bin/mafft --ep 0.123 --localpair --maxiterate 1000 --quiet $input > $temp_file");
    #my_cmd("/projects/sate7/tools/mafft-6.864-with-extensions/bin/mafft --ep 0.123 --retree 2 --quiet $input > $temp_file");
  }
  if (check_alignment($temp_file) eq "") {
    `mv $temp_file $output`;
  }
  #my_cmd("/u/namphuon/bin/mafft/mafft --localpair --maxiterate 1000 --quiet $input > $output");
}

sub align_opal {
  my $input = $_[0];
  my $output = $_[1];
  
  my_cmd("java -Xmx4g -jar $bin/opal.jar --in $input --out $output");
}

sub align_opal_mem {
  my $input = $_[0];
  my $output = $_[1];
  
  my_cmd("java -Xmx8g -jar $bin/opal.jar --in $input --out $output");
}

sub align_muscle {
  my $input = $_[0];
  my $output = $_[1];
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }  

  my_cmd("$bin/muscle -in $input -out $temp_file -quiet");
  if (check_alignment($temp_file) eq "") {
    `mv $temp_file $output`;
  }
}

sub convert_paup_newick {
  my $input_file = $_[0];
  my $output_file = $_[1];  
  my_cmd("$python $script_directory/phylo.py -c -i $input_file -o $output_file");
  my_cmd("$perl -p -i -e \"s/\\[.*\\]//g\" $output_file");
}

sub read_rnasim_nexus_tree_alignment {
  my $input_file = $_[0];
  my $output_tree = $_[1];
  my $output_alignment = $_[2];
  
  open(INPUT, "<$input_file");
  my $read_alignment = 0;
  
  while (my $line = <INPUT>) {
    if ($line =~ m/begin trees/) {
      $line = <INPUT>;
      $line =~ m/=\s+(.*)\s+/;
      my $tree = $1;
      $tree =~ s/_I\d+//;
      open(OUTPUT, ">$output_tree");
      print OUTPUT $tree;
      close(OUTPUT);
    }
    if ($line =~ m/begin characters/) {
      $read_alignment = 1;
      $line = <INPUT>;
      $line = <INPUT>;
      $line = <INPUT>;
      open(OUTPUT, ">$output_alignment");
      next;
    }    
    if ($read_alignment == 1 and $line =~ m/;/) {
      close(OUTPUT);
      close(INPUT);
      last;
    } elsif ($read_alignment == 1) {
      if ($line =~ m/^_/) {
        next;
      }
      $line = Phylo::trim($line);
      my @results = split(/\s+/, $line);
      $results[0] = Phylo::trim($results[0]);
      $results[1] = Phylo::trim($results[1]);
      print OUTPUT ">$results[0]\n$results[1]\n";
    } 
  }  
}

sub convert_nexus_newick {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  my $line = "";
  my $tree_block = 0;
  while (defined $line) {
    $line = <INPUT>;
    if (not defined $line) {
      next;
    }
    if (not $tree_block and $line =~ m/begin trees/i) {
      $tree_block = 1;
      next;
    }
    if ($tree_block and $line =~ m/end\s*;/) {
      $line = undef;
      next;
    }
    if ($tree_block and $line =~ m/^\s*tree\s+[^=]*\s*=/) {
      $line =~ m/(\([^;]*;)/;
      if (defined $1) {
        print OUTPUT $1 . "\n";
      } else {
        my $tree = <INPUT>;
        print OUTPUT $tree;
      }
    }
  }
  close(INPUT);
  close(OUTPUT);
  return "";
}


sub remove_seq {
  my $input_file = $_[0];
  my $output_file = $_[1];

  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  my $sequence = "";
  my $line = <INPUT>;
  while (defined($line)) {
    if($line =~ />/) {      
      $line  =~ m/>([\S]*)[\s]*[\s]*([\S]*)/;
      if (($1 =~ /SEQ[\d]+/) or ($1 =~ /TheTree/)) {
        $line = <INPUT>;
        while (not ($line =~ />/)) {
          if (not defined($line)) {
            last;
          }
          $line = <INPUT>;
        }
        next;
      } else {
        print OUTPUT ">$1\n$2";
        $line = <INPUT>;
      }      
    } else {
      print OUTPUT $line;
      $line = <INPUT>;
    } 
    
    
  }
  close(INPUT);
  close(OUTPUT);
}

sub generate_alignment_files {
  my $input_file = $_[0];
  my $output_dir = $_[1];
  my $dict_ref = get_dictionary_from_rose($input_file);
  my $current = trim(my_cmd("pwd"));  
  chdir($output_dir);  
  print my_cmd("$rose $current/$input_file");
  write_tree_file($dict_ref, "rose.tt");
  remove_seq("test.fas", "rose.fas");
  remove_seq("test.fa", "rose.aln.true.fasta");
  my_cmd("rm test.fa test.fas test.tree");
  write_rose_file($dict_ref, "rose.internal.model");
  chdir($current);
}

sub generate_rose_files {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $sequence_length = $_[2];
  
  my %options = ("OutputFilebase","\"test\"", "AlignmentFormat", "\"FASTA\"", "StdOut", "False");
  
  my $dict_ref = change_sequence_length($input_file, $sequence_length);
  write_rose_file($dict_ref, $output_file, \%options);  
}

#Changes the sequence length of a rose internal model file and 
#outputs the new version out.  
#Args input file, sequence length
sub change_sequence_length {
  my $input_file = $_[0];  
  my $sequence_length = $_[1];
  
  my $new_sites = draw_sites($sequence_length);
  my $ref_dict = get_dictionary_from_rose($input_file);
  @$ref_dict{"TheMutationProbability"}=$new_sites;
  @$ref_dict{"SequenceLen"}=$sequence_length;
  return $ref_dict;
}

#N draws to a gamma distribution, returns a string in rose ready format
sub draw_sites {
  my $draws = $_[0];
  my $sites = "[" . rgamma(1.0);
  
  my $counter = 1;
  while ($counter < $draws) {
    $sites = $sites . "," . sprintf("%2.20f", rgamma(1.0));
    $counter = $counter+1;
  }  
  return $sites . "]";
}

#Draw from random beta, code from Kevin
sub rbeta {
  my $a=shift;
  my $b=shift;
  $b=$b-1;
  my $i=1;
  while(1) {
      my $ri=rand();
      my $x=$ri**(1/$a);
      my $ri1=rand();
      my $y=$ri1**(1/$b);
      if ($x+$y<=1) {return $x;}
      $i=$i+2;
  }
}

#Draw from random gamma, code from Kevin
sub rgamma {
  my $a=shift;
  my $z=0;
  my $i;
    if (int($a)==$a) {    # integral
      for ($i=0;$i<$a;$i++) {
          $z=$z+rexp(1);
      }
      return $z;
  }
      for ($i=1;$i<=$a;$i++) {$z+=rand();}
  my $x=rbeta($a-int($a),2-$a+int($a));
  my $y=-log(rand()*rand());
  return ($z+$x*$y);
}

# from numerical recipes
sub rexp {
  my $avg=shift;
  my $dum=0.0;
  do { $dum=rand(1); } while($dum==0);
  return -log($dum)*$avg;
}

# from numerical recipes
sub gaussian_rand {
    my $mean = $_[0];
    my $std = $_[1];
    
    if (not defined $mean) {
      $mean = 0;
    } if (not defined $std) {
      $std = 1;
    } 
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return ($std*$g1)+$mean;
}

#Writes a dictionary back into rose ready format
sub write_rose_file {
  my $ref_dict = $_[0];
  my $output = $_[1];
  my $optional_ref = $_[2];

  my %dict = %$ref_dict;
  if (defined($optional_ref)) {
    my %optional = %$optional_ref;
    foreach my $key (keys %optional) {
      $dict{$key} = $optional{$key};      
    }
  }
    
  open (OUTPUT, ">$output");
  foreach my $keys (keys %dict) {
    print OUTPUT "$keys = $dict{$keys}\n";
  }
  close(OUTPUT);
}

sub write_tree_file {
  my $ref_dict = $_[0];
  my $output = $_[1];
  
  my %dict = %$ref_dict;
  open (OUTPUT, ">$output");
  print OUTPUT $dict{"TheTree"};
  close(OUTPUT);
}

sub get_dictionary_from_rose {
  my $input = $_[0];
  my %dictionary = ();

  open(ROSE_FILE,"<$input");  
  my $field_value;
  my $field_key;
  while (my $line = <ROSE_FILE>) {
    if($line =~ /#/) {      
      next;
    }    
    elsif($line =~ /=/) {      
      $line =~ m/[\s]*([\S]+)[\s]*=[\s]*(.*)/;    
      my ($e1, $e2) = ($1, $2);
      if (!defined($field_key)) {
        $field_key = $e1;
        $field_value = $e2;
      } else {
        #print "$field_key:$field_value\n";
        $dictionary{$field_key} = $field_value;
        $field_key = $e1;
        $field_value = $e2;
      }
    } else {
      $field_value = $field_value . $line;
    }
  }
  $dictionary{$field_key} = $field_value;
  close(ROSE_FILE);
  return \%dictionary;
}

sub run_sate {
  my $alignment_file = $_[0];
  my $config_file = $_[1];
  my $machine_t = $_[2];
  if ((defined($machine_t) and $machine_t == 1) || $machine == 1) {
    return my_cmd("$python $sateV -i $alignment_file $config_file");
  }
  return my_cmd("$python $sate -i $alignment_file $config_file");
}

sub run_sate2 {
  my $alignment_file = $_[0];
  my $config_file = $_[1];
  my $output_alignment = $_[2];
  my $output_tree = $_[3];  
  my $directory = $_[4];  
  
  my $temp_file = get_temp_file();
  print "Temp directory: $temp_file\n";
    
  print my_cmd("python2.7 /projects/sate3/namphuon/bin/satesrc-v2.2.3-2012May22/sate-core/run_sate.py -i $alignment_file -o $temp_file -j satejob $config_file");
  
  my ($filename, $basename, $empty) = fileparse($alignment_file); 
  if (check_alignment("$temp_file/satejob.marker001.$filename") ne "") {
    #`rm $temp_file* -rf`;
    return;
  }
  
  if (-e "$temp_file/satejob.tre") {
    my_cmd("cp $temp_file/satejob.tre $output_tree");    
  } else {
    #my_cmd("rm $temp_file -rf");
    return;
  } 
  if (-e "$temp_file/satejob.marker001.$filename") {
    my_cmd("cp $temp_file/satejob.marker001.$filename $output_alignment ");
  } else {
    my_cmd("cp $temp_file/satejob.marker001*.aln $output_alignment");
  } 
  `mv $temp_file $directory`;  
}


sub run_sate_2 {
  my $alignment_file = rel2abs($_[0]);
  my $config_file = rel2abs($_[1]);
  my $output_tree = rel2abs($_[2]);
  my $output_alignment = rel2abs($_[3]);
  my $directory = rel2abs($_[4]);

  if (defined $directory) {      
      my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
      while (-e $temp_file) {
        $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
      }
      my $name = Phylo::get_random_name($temp_dir);
      find_replace("$config_file", "$temp_file", {"<output_directory>", "$directory"});
      my_cmd("/projects/sate7/tools/bin/sate.2.1.0  -i $alignment_file -j $name $temp_file");
      my ($filename, $basename, $empty) = fileparse($alignment_file);  
      my_cmd("cp $basename/$name.tre $output_tree");
      if (-e "$basename/$name.marker001.$filename") {
        my_cmd("cp $basename/$name.marker001.$filename $output_alignment ");
      } elsif (-e "$basename/$name.marker001.$filename.aln") {      
        my_cmd("cp $basename/$name.marker001.$filename.aln $output_alignment ");      
      } else {
        my $align = trim(my_cmd("ls $basename/$name*.aln"));
        my_cmd("cp $align $output_alignment");      
      }
      my_cmd("rm $temp_file $basename/$name*");      
  }
   else {
    return my_cmd("/lusr/bin/python /projects/sate3/namphuon/sate2/sate-core/run_sate.py -i $alignment_file -j sate $config_file");
   }
   
}

sub run_sate_3 {
  my $alignment_file = $_[0];
  my $config_file = $_[1];
  my $commands = $_[2];  
  if (not defined $commands) {
    $commands = "";
  }
  print "/lusr/bin/python /projects/sate3/namphuon/sate2/sate-core/run_sate.py $commands -i $alignment_file -j sate $config_file\n";  
  return my_cmd("/lusr/bin/python /projects/sate3/namphuon/sate2/sate-core/run_sate.py $commands -i $alignment_file -j sate $config_file");
}
      

#Finds intersection of two maps
sub intersection {
  my %set1 = %{$_[0]};
  my %set2 = %{$_[1]};
  
  my %intersection = ();
  foreach my $key (keys %set1) {
    if (exists $set2{$key}) {
      $intersection{$key}=$key;
    }
  }
  return \%intersection;
}

#Finds difference of two maps
sub difference {
  my %set1 = %{$_[0]};
  my %set2 = %{$_[1]};
  
  my %difference = %set1;
  foreach my $key (keys %set2) {
    if (exists $difference{$key}) {
      delete $difference{$key};
    }
  }
  return \%difference;
}


#Finds union of two maps
sub union {
  my %set1 = %{$_[0]};
  my %set2 = %{$_[1]};
  
  my %map = map { $_ => 1 } (keys %set1, keys %set2);
#   foreach my $key (%map) {
#     print "$key\t";
#   }
  return \%map;
}


#Get all taxa in newick tree
sub get_taxa {
  my $str = $_[0];
  my @matches = ($str =~ m/[(,]([^(,:)]+)/g);
  my @taxa = ();
  foreach my $taxon (@matches) {
    if (defined $taxon) {
  push(@taxa, trim($taxon));
    }
  }
  return \@taxa;
}

#Check to see if the tree and alignment is on the same set of taxa
sub tree_alignment_on_same_taxa {
  my $tree_file = $_[0];
  my %alignment = %{read_fasta_file($_[1])};
  open(INPUT, "<$tree_file");
  my $tree = <INPUT>;
  close(INPUT);
  
  my @taxa = @{get_taxa($tree)};
  my $tree_total = scalar @taxa;
  my $align_total = scalar keys %alignment;
  if ($tree_total ne $align_total) {
    return 0;
  }
  foreach my $taxon (@taxa) {
    if (not defined $alignment{$taxon}) {
      return 0;
    }
  }
  return 1;
}

#Check to see if the trees are on the same set of taxa
sub trees_on_same_taxa {
  my $tree_1 = $_[0];
  my $tree_2 = $_[1];
  
  my $tree_str_1 = `cat $tree_1`;
  my $tree_str_2 = `cat $tree_2`;
  
  my @taxa_1 = @{get_taxa($tree_str_1)};
  my @taxa_2 = @{get_taxa($tree_str_2)};
  my %taxa_map_1 = map {$_=>$_} @taxa_1;
  my %taxa_map_2 = map {$_=>$_} @taxa_2;
  foreach my $key (keys %taxa_map_1) {
    if (not defined $taxa_map_2{$key}) {
      return 0;
    } else {
      delete $taxa_map_2{$key};
    } 
  }
  my @keys_2 = keys %taxa_map_2;
  if (scalar @keys_2 != 0) {
    return 0;
  }
  return 1;
}



# Fisher yates shuffle to generate in place random permutation
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

sub getEdgeDistances {
  my $true_tree = $_[0];
  my $estimated_tree = $_[1];
  my $output_file = $_[2];

  my $results = `java -jar \$PHYLO_DIR/phylonet_v2_4/phylonet_v2_4.jar rf -m $true_tree -e $estimated_tree`;
  $results =~ m/([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)/;
  my ($fn, $fp, $te, $tm) = ($1, $2, $3, $4);
  
  if (defined $fn and defined $output_file) {
    my $temp_file = get_temp_file();
    open(OUTPUT, ">$temp_file");
    print OUTPUT $results;
    close(OUTPUT);
    `mv $temp_file $output_file`;
  }
  return ($fn, $fp, $te, $tm);
}

sub getFPFN {
  my $true_tree = $_[0];
  my $estimated_tree = $_[1];
  my $output_file = $_[2];
  if (not -e $true_tree) {
      print "Missing true tree $true_tree\n";
      return;
  } elsif (not -e $estimated_tree) {
      print "Missing estimated tree $estimated_tree\n";
      return;
  }
  
  my $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
      $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }
  
  to_upper($true_tree, $temp_file . "true_upper"); 
  my $temp_true = $temp_file . "true_upper";
  to_upper($estimated_tree, $temp_file . "estimated_upper");
  my $temp_estimated = $temp_file . "estimated_upper";
  
  `$perl -p -i -e "s/'//g" $temp_true`;
  `$perl -p -i -e "s/'//g" $temp_estimated`;  
  
  remove_weights($temp_true, "$temp_true.1");
  remove_weights($temp_estimated, "$temp_estimated.1");
  my $results = `java -Xmx4g -jar /projects/sate3/namphuon/bin/phylonet_v2_4/phylonet_v2_4.jar rf -m $temp_true.1 -e $temp_estimated.1`;
  $results =~ m/([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)/;
  my ($fn, $fp, $te, $tm) = ($1, $2, $3, $4);
  my_cmd("rm $temp_true* $temp_estimated*");
  if ((not defined $fn) or (not defined $fp) or (not defined $tm) or (not defined $te)) {
    return undef;
  }
  #print "$results\n$fp $fn $te $tm \n";
  my $result = sprintf("(%0.6f, %0.6f, %0.6f)", $fp/$te, $fn/$tm, ($fp+$fn)/($tm+$te));
  if (defined $output_file) {
      open(OUTPUT, ">$output_file");
      print OUTPUT $result;
      close(OUTPUT);
  }

  return $result;
}



sub getSumFPFN {
  my $estimated_tree = $_[0];
  my $source_trees = $_[1];
  my $output_file = $_[2];
  my $averaged = $_[3];
  if (defined $averaged) {
    $averaged = "-a True";    
  } else {
    $averaged = "";
  } 
  if (not -e $estimated_tree) {
    print "Missing $estimated_tree\n";
    return -1;
  } elsif (not -e $source_trees) {
    print "Missing $source_trees\n";
    return -1;
  }   
  
  my $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }  
  
  to_upper($estimated_tree, $temp_file . "estimated_upper"); 
  to_upper($source_trees, $temp_file . "source_upper");  
   my $result = my_cmd("$python $sumfpfn -t $temp_file"."estimated" . "_upper -s $temp_file"."source" . "_upper $averaged");
  my_cmd("rm $temp_file"."estimated" . "_upper $temp_file"."source" . "_upper ");
  $result =~ s/None/0.0/g;
  if (defined $output_file and $output_file ne "") {
    open(OUTPUT, ">$output_file");
    print OUTPUT $result;
    close(OUTPUT);
  }    
  return $result;
}

sub get_tree_stats {
  my $stats_file = $_[0];
  if (not -e $stats_file) {
    return undef;
  }
  my $results = `cat $stats_file`;
  Phylo::trim($results);  
  $results =~ m/\(([^,]+),\s*([^,]+),\s*([^,]+)\)/;
  return ($1, $2, $3);
}

sub remove_3rd_codon_position {
  my $fasta_file = $_[0];
  my $output_file = $_[1];
  my $partition_file = $_[2];
  my $output_partition = $_[3];
  
  my %alignment = %{Phylo::read_fasta_file($fasta_file)};
  my $length = 0;
  foreach my $key (keys %alignment) {
    my $line = $alignment{$key};
    my @results = split(//,$line);
    my $new_align = "";
    for(my $idx = 0; $idx < scalar(@results); $idx++) {
      if ($idx % 3 == 2) {
        next;
      }
      $new_align = $new_align . $results[$idx];
    }
    $alignment{$key}=$new_align;
  }
  Phylo::write_alignment(\%alignment, $output_file);
  if (defined $partition_file) {
    open(IPART, "<$partition_file");
    open(OPART, ">$output_partition");
    while (my $line = <IPART>) {      
      $line = Phylo::trim($line);
      $line =~ m/=\s*(\d+)-(\d+)/;
      my ($start, $end) = ($1, $2);
      my $replace = "$start-$end";
      $start = int($start/3*2)+1;
      $end = int($end/3*2);
      $line=~s/$replace/$start-$end/;
      print OPART $line . "\n";    
    }
    close(IPART);    
    close(OPART);  
  }
}


sub convert_fasta_to_phylip {
  my $fasta_file = $_[0];
  my $output_file = $_[1];

  my %sequences = %{read_fasta_file($fasta_file, 0)};
  my $temp_file = get_temp_file();
  open(OUTPUT, ">$temp_file");
  my @names = sort keys %sequences;
  print OUTPUT scalar @names . " " . length($sequences{$names[0]}) . "\n";
  foreach my $name (@names) {
    print OUTPUT "$name $sequences{$name}\n"
  }
  close(OUTPUT);
  `mv $temp_file $output_file`;
  #my_cmd("java -cp $readseq run -o $output_file"." -f12 $fasta_file");
}

sub phylip_to_fasta {
  my $phylip_file = $_[0];
  my $fasta_file = $_[1];
  my $split_it = $_[2];
  my $temp_file = get_temp_file();

  open(INPUT, "<$phylip_file");
  open(OUTPUT,">$temp_file");
  my $line = <INPUT>;
  while (defined ($line = <INPUT>)) {
    if (trim($line) eq "") {
      next;
    }
    if (not defined $split_it) {
        print OUTPUT ">$line";
        $line = <INPUT>;
        print OUTPUT $line;
    } else {
      my @results = split(/\s+/,$line);
      print OUTPUT ">$results[0]\n$results[1]\n";      
    }
  }  
  close(INPUT);
  close(OUTPUT);
  `mv $temp_file $fasta_file`;
}

sub phylip_to_fasta_interleave {
  my $phylip_file = $_[0];
  my $fasta_file = $_[1];
  
  open(INPUT, "<$phylip_file");  
  my $line = <INPUT>;
  $line =~ m/(\d+)\s+(\d+)/;
  my ($taxa, $length) = ($1,$2);
  my @taxon_index = ();
  my $counter = 0;
  my %alignment = ();
  local $" = "";
  while (defined ($line = <INPUT>)) {
    if ($counter < $taxa) {
      my @results = split(/\s+/,$line);
      my $name = shift(@results);      
      $line = "@results";
      $line =~ s/\s+//g;
      push(@taxon_index, $name);
      $counter++;
      $alignment{$name}=$line;
    } else {
      $line =~ s/\s+//g;
      if ($line eq "") {
        next;
      }
      $alignment{$taxon_index[$counter % $taxa]}.=$line;
      $counter++;
    }
  }
  foreach my $key (keys %alignment) {
    if (length($alignment{$key}) != $length) {
      print "$key is not right length!\n";
      exit;
    }
  }
  Phylo::write_alignment(\%alignment, $fasta_file);
}

sub convert_to_fasta {
  my $fasta_file = $_[0];
  my $output_file = $_[1];

  my_cmd("java -cp $readseq run -o $output_file"." -f8 $fasta_file");
}

sub convert_fasta_to_paup {
  my $fasta_file = $_[0];
  my $output_file = $_[1];

  my_cmd("java -cp $readseq run -o $output_file"." -f17 $fasta_file");
}

sub read_stockholm_file {
  my $stockholm_file = $_[0];
  open(INPUT, "<$stockholm_file");
  my %seq;
  while (<INPUT>) {
      next unless /\S/;
      next if /^\s*\#/;
      if (/^\s*\/\//) {}
      else {
        chomp;
        my ($name, $seq) = split;
        $seq=~s/\./-/g;
        $seq{$name} .= $seq;
      }
  }
  return \%seq;
}

sub convert_stockholm_to_fasta {
  my %alignment = %{read_stockholm_file($_[0])};
  write_alignment(\%alignment, $_[1]);
}

sub convert_fasta_to_msf {
  my $fasta_file = $_[0];
  my $output_file = $_[1];
    
  my_cmd("java -cp $readseq run -o $output_file"." -f15 $fasta_file");
}

sub compare_trees {
  my $input1 = $_[0];
  my $input2 = $_[1];
  my $output = $_[2];
  `java -jar /projects/sate3/namphuon/bin/vtd.jar $input1 $input2 > $output`;
}


sub compare_files {
  my $fasta_file = $_[0];
  my $output_file = $_[1];
  return my_cmd("$python $script_directory/mask.py  -c -f $fasta_file -o $output_file");
}

sub mask_msa {
  my $msa_file = $_[0];
  my $score_file = $_[1];
  my $threshold = $_[2];
  my $output_file = $_[3];
  
  return my_cmd("$python $script_directory/mask.py -m -f $msa_file -s $score_file -t $threshold -o $output_file");
}

sub mask_nt_msa {
  my $msa_file = $_[0];
  my $score_file = $_[1];
  my $threshold = $_[2];
  my $mask = $_[3];
  my $output_file = $_[4];

  return my_cmd("$python $script_directory/mask.py -m -f $msa_file -s $score_file -t $threshold -p $mask -o $output_file");
}

sub condor_generate_scores {
  my $fasta_file = $_[0];
  my $alt_files = $_[1];
  my $name = $_[2];
  my $machine_t = $_[3];
  my $dir = $_[4];
  my $rc = "";
    
  #generate_file_list("temp_file_list.txt", $alt_files);
  
  mkdir("$temp_dir/$dir");
  mkdir("temp_$dir");

  to_upper($fasta_file, "temp_$dir.fasta");
  my $counter = 0;
  foreach my $file (@{$alt_files}) {
    #print my_cmd("pwd");
    my_cmd("cp $file temp_$dir/file_$counter");
    to_upper("temp_$dir/file_$counter", "temp_$dir/xfile_$counter");    
    $counter = $counter + 1;
  }
  my_cmd("rm temp_$dir/file_*");
  
  #Need to fix this, figure out why using uvanimor doesn't work
  if ( $machine_t == 1 || $machine == 1) {
    $rc=my_cmd($msa_scoreV . " " . "temp_$dir.fasta" . " " . "$temp_dir/$dir/temp" . " -d " . "temp_$dir");  
  } 
  else {    
    $rc=my_cmd($msa_score . " " . "temp_$dir.fasta" . " " . "$temp_dir/$dir/temp" . " -d " . "temp_$dir");      
  }  

  #my_cmd("rm temp_file_list.txt");
  my_cmd("rm -rf temp_$dir temp_$dir.fasta $temp_dir/$dir/temp_msa.scr $temp_dir/$dir/temp_col_col.scr $temp_dir/$dir/temp_res_pair.scr $temp_dir/$dir/temp_res_pair_seq.scr $temp_dir/$dir/temp_res_pair_seq_pair.scr");  
  my_cmd("mv $temp_dir/$dir/temp_res_pair_col.scr $temp_dir/$dir/$name" . "_res_pair_col.scr");
  my_cmd("mv $temp_dir/$dir/temp_res_pair_res.scr $temp_dir/$dir/$name" . "_res_pair_res.scr");
}

sub condor_scores {
  my $fasta_file = $_[0];
  my $alt_files = $_[1];
  my $name = $_[2];
  my $machine_t = $_[3];
  my $dir = $_[4];
  my $rc = "";
    
  #generate_file_list("temp_file_list.txt", $alt_files);
  
  mkdir("temp_$dir");

  to_upper($fasta_file, "temp_$dir.fasta");
  my $counter = 0;
  foreach my $file (@{$alt_files}) {
    #print my_cmd("pwd");
    my_cmd("cp $file temp_$dir/file_$counter");
    to_upper("temp_$dir/file_$counter", "temp_$dir/xfile_$counter");    
    $counter = $counter + 1;
  }
  my_cmd("rm temp_$dir/file_*");
  
  #Need to fix this, figure out why using uvanimor doesn't work
  if ( $machine_t == 1 || $machine == 1) {
    $rc=my_cmd($msa_scoreV . " " . "temp_$dir.fasta" . " " . "temp" . " -d " . "temp_$dir");  
  } 
  else {    
    $rc=my_cmd($msa_score . " " . "temp_$dir.fasta" . " " . "temp" . " -d " . "temp_$dir");      
  }  

  #my_cmd("rm temp_file_list.txt");
  my_cmd("rm -rf temp_$dir temp_$dir.fasta temp_col_col.scr temp_res_pair.scr temp_res_pair_seq.scr temp_res_pair_seq_pair.scr");  
  my_cmd("mv temp_res_pair_col.scr $dir" . "/temp_res_pair_col.scr");
  my_cmd("mv temp_res_pair_res.scr $dir" . "/temp_res_pair_res.scr");
  my_cmd("mv temp_msa.scr $dir" . "/temp_msa.scr");
}

sub calc_average_sp {
  my $model = $_[0];
  my $file = $_[1];
  
  my @replicates = (<$model/R*>);
  my $average_sp = 0;
  my $count = 0;
  foreach my $replicate (@replicates) {    
    if (not -e ("$replicate/" . $file)) {
      next;
    } else {
      $average_sp = $average_sp + Phylo::get_sp_score(("$replicate/" . $file));      
      $count++;
    }
  }
  if ($count == 0) {
    return -1;
  }
  return $average_sp/$count;
}

sub get_directories {  
  my @array = split(/\//, $_[0]);
  return \@array;
}

sub generate_scores {
  my $fasta_file = $_[0];
  my $alt_files = $_[1];
  my $type = $_[2];
  my $machine_t = $_[3];
  my $rc = "";
    
  #generate_file_list("temp_file_list.txt", $alt_files);
  
  my $random_number = Phylo::get_random_name();
  mkdir("temp_$random_number");

  to_upper($fasta_file, "temp_$random_number.fasta");    
  my $counter = 0;
  foreach my $file (@{$alt_files}) {
    #print my_cmd("pwd");
    my_cmd("cp $file temp_$random_number/file_$counter");
    to_upper("temp_$random_number/file_$counter", "temp_$random_number/xfile_$counter");    
    $counter = $counter + 1;
  }
  #my_cmd("rm temp_$random_number/file_*");
  
  #Need to fix this, figure out why using uvanimor doesn't work
  if ( $machine_t == 0 || $machine == 0) {
    #$rc=my_cmd($msa_score . " " . $fasta_file . " " . $type . " -f " . "temp_file_list.txt");
    $rc=my_cmd($msa_score . " " . "temp_$random_number.fasta" . " " . $type . " -d " . "temp_$random_number");  
  } 
  else {
    #$rc=my_cmd($msa_scoreV . " " . $fasta_file . " " . $type . " -f " . "temp_file_list.txt");
    $rc=my_cmd($msa_scoreV . " " . "temp_$random_number.fasta" . " " . $type . " -d " . "temp_$random_number");
  }  

  #my_cmd("rm temp_file_list.txt");
     #my_cmd("rm -rf temp_$random_number temp_$random_number.fasta");
}

sub generate_stats {
  my $file_name = $_[0];
  my $true_name = $_[1];

  #Now open both files and find out which NT pairs are TP vs FP and bin
  #them as data
  my %positives = ();
  my %negatives = ();
  my %columns = ();

  #Now process the files
  open(ESTIMATE_FILE, "< $file_name" . "_res_pair.scr");
  open(TRUE_FILE, "< $true_name" . "_res_pair.scr");
  open(COL_FILE, "< $file_name" . "_res_pair_col.scr");
  
  #Throw away first line
  my $estimate_support_line = <ESTIMATE_FILE>;
  my $true_support_line = <TRUE_FILE>;
  my $col_support_line = <COL_FILE>;

  
  $col_support_line = <COL_FILE>;
  $col_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)/;
  my $col = $1;
  my $score = $2;
  my $correct = 0;
  my $total = 0;

  while($estimate_support_line = <ESTIMATE_FILE>) {
    $true_support_line = <TRUE_FILE>;

    $estimate_support_line = trim($estimate_support_line);    
    $estimate_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)/;
    my ($e1, $e2, $e3, $e4) = ($1, $2, $3, $4);

    $true_support_line = trim($true_support_line);
    $true_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)/;
    my ($t1, $t2, $t3, $t4) = ($1, $2, $3, $4);
    
    #For each pair, find the support index and whether it is
    #a false positive or true positive
    if (!defined($t4) || !defined($e4)) {
      #print "NOT DEFINED!\n$sate_support_line\n$msa_support_line\n$sopfn_line\n";
      last;
    }
  
    if (($t1 ne $e1) || ($t2 ne $e2) || ($t3 ne $e3)) {
      print "ERROR READING FILE!\n";
      last;
    }

    if ($t1 ne $col) {
      my @temp_list = ($correct,$total,$score);
      $columns{$col} = \@temp_list;
      $correct = 0;
      $total = 0;

      $col_support_line = <COL_FILE>;
      $col_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)/;
      $col = $1;
      $score = $2;
    }

    if ($t4 > 0) {
      if (defined($positives{substr($e4,0,5)})) {
        $positives{substr($e4,0,5)} = $positives{substr($e4,0,5)}+1;
      } else {
        $positives{substr($e4,0,5)} = 1;
      }
      $correct++;
      $total++;
    } else {
      if (defined($negatives{substr($e4,0,5)})) {
        $negatives{substr($e4,0,5)} = $negatives{substr($e4,0,5)}+1;
      } else {
        $negatives{substr($e4,0,5)} = 1;
      }
      $total++;
    }
  }
  close(ESTIMATE_FILE);
  close(TRUE_FILE);
  close(COL_FILE);

  return (\%positives, \%negatives, \%columns);
}

sub generate_sumfp {
  my $file_name = $_[0];
  my $true_name = $_[1];
  my $machine_t = $_[2];  
  
  my @true_name = ($true_name);
  generate_scores($file_name, \@true_name, "fp_temp1", $machine_t);

  #Now open both files and find out which NT pairs are TP vs FP and bin
  #them as data
  my $positives = 0.0;
  my $negatives = 0.1;

  #Now process the files
  open(TRUE_FILE, "<fp_temp1_res_pair.scr");        
  
  #Throw away first line  
  my $true_support_line = <TRUE_FILE>;
  while($true_support_line = <TRUE_FILE>) {
    $true_support_line = trim($true_support_line);
    $true_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)/;
    my ($t1, $t2, $t3, $t4) = ($1, $2, $3, $4);
    
    #For each pair, find the support index and whether it is
    #a false positive or true positive
    if (!defined($t4)) {
      print "NOT DEFINED!\n$true_support_line\n";
      print "$t1 $t2 $t3 $t4\n";
      last;
    }
  
    if ($t4 > 0) {
      $positives = $positives + 1.0;
    } else {
      $negatives = $negatives + 1.0;
    }
  }
  close(TRUE_FILE);
  my_cmd("rm fp_temp1*.scr");
  return ($positives, $negatives, ($positives/($negatives+$positives)));
}

sub generate_nt_pairing_stats {
  my $file_name = $_[0];
  my $true_name = $_[1];
  my $file_list = $_[2];
  my $machine_t = $_[3];  
  
  my @true_name = ($true_name);
  generate_scores($file_name, $file_list, "phy_temp", $machine_t);
  generate_scores($file_name, \@true_name, "fp_temp", $machine_t);

  #Now open both files and find out which NT pairs are TP vs FP and bin
  #them as data
  my %positives = ();
  my %negatives = ();

  #Now process the files
  open(ESTIMATE_FILE, "< phy_temp_res_pair.scr");
  open(TRUE_FILE, "<fp_temp_res_pair.scr");        
  
  #Throw away first line
  my $estimate_support_line = <ESTIMATE_FILE>;
  my $true_support_line = <TRUE_FILE>;
  while($estimate_support_line = <ESTIMATE_FILE>) {
    $true_support_line = <TRUE_FILE>;

    $estimate_support_line = trim($estimate_support_line);    
    $estimate_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)/;
    my ($e1, $e2, $e3, $e4) = ($1, $2, $3, $4);

    $true_support_line = trim($true_support_line);
    $true_support_line =~ m/[\s]*([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)/;
    my ($t1, $t2, $t3, $t4) = ($1, $2, $3, $4);
    
    #For each pair, find the support index and whether it is
    #a false positive or true positive
    if (!defined($t4) || !defined($e4)) {
      print "NOT DEFINED!\n$estimate_support_line\n$true_support_line\n";
      print "$t1 $t2 $t3 $t4 $e1 $e2 $e3 $e4\n";
      last;
    }
  
    if (($t1 ne $e1) || ($t2 ne $e2) || ($t3 ne $e3)) {
      print "ERROR READING FILE!\n";
      last;
    }
    
    if ($t4 > 0) {
      if (defined($positives{substr($e4,0,5)})) {
        $positives{substr($e4,0,5)} = $positives{substr($e4,0,5)}+1;
      } else {
        $positives{substr($e4,0,5)} = 1;
      } 
    } else {
      if (defined($negatives{substr($e4,0,5)})) {
        $negatives{substr($e4,0,5)} = $negatives{substr($e4,0,5)}+1;
      } else {
        $negatives{substr($e4,0,5)} = 1;
      }
    }
  }
  close(ESTIMATE_FILE);
  close(TRUE_FILE);

  #my_cmd("rm phy_temp*.scr fp_temp*.scr");
  print "Printing values!\n";
  foreach my $t (keys %positives) {
    print "pos{$t} = " . $positives{$t} . "\n";
  }
  foreach my $t (keys %negatives) {
    print "neg{$t} = " . $negatives{$t} . "\n";
  }

  return (\%positives, \%negatives);
}

sub generate_file_list {
  my $file_name = $_[0];
  my $file_list = $_[1];
  open FILES, ">$file_name";    
  foreach my $file (@{$file_list}) {
    print FILES $file ."\n";
  }
  close(FILES);
}

sub get_bootstrap_tree {
  my $raxml_file = rel2abs($_[0]);
  my $threshold = $_[1];
  my $output = rel2abs($_[2]);
  
  `$perl /projects/sate3/namphuon/bin/calculate_bootstrap_tree_from_bipartitions.pl -i $raxml_file -o $output -t $threshold`
}


#This function finds a sequence in an alignment and returns the starting
#index and length.  If it's not found, it returns -1,-1.
sub find_substring {
  my $substring = lc $_[0];
  my $alignment = lc $_[1];
  
  #First unalign alignment
  my $unaligned = $alignment;
  $unaligned =~ s/-//g;
  
  #Next find if sequence is in the unaligned sequence
  my $result = index($unaligned, $substring);
  
  #result is not found, return (-1,-1)
  if ($result == -1) {
      return (-1,-1);
  }
  
  #Now find the substring in the alignment, start by finding prefix, substring, postfix in unaligned
  my $prefix = substr($unaligned,0,$result);
  my $fix = substr($unaligned,$result,(length $substring));
  my $postfix = substr($unaligned,(length $substring)+$result);
  
  #Find index of last letter of prefix and first letter of postfix in alignment.  Do this by counting
  my $start = substr($fix,0,1);
  my $end = substr($fix,-1,1);
  my $temp = "$prefix$start";
  my $start_count = ( $temp =~ s/$start//g);

  $temp = "$prefix$fix";
  my $end_count = ($temp =~ s/$end//g);
  
  #Now loop through original alignment, counting each occurance of a till we hit the specified count
  my $idx = index($alignment,$start);
  my $counts = 1;
  while ($counts < $start_count) {
    $idx = index($alignment,$start,$idx+1);    
    $counts++;
  }
  my $start_idx = $idx;

  #Now loop through original alignment, counting each occurance of a till we hit the specified count
  $idx = index($alignment,$end);
  $counts = 1;
  while ($counts < $end_count) {
    $idx = index($alignment,$end,$idx+1);    
    $counts++;
  }
  my $end_idx = $idx;
  return ($start_idx,$end_idx-$start_idx+1);
}

sub super_tree_nj {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $method = $_[2];
  
  if (not defined $method) {
    $method = "fastme";
  }
  
  my $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }
    
  compute_fast_mrp($input_file, $temp_file, "FASTA");
  if ($method eq "paup") {
    run_nj_paup($temp_file, $output_file, 2);
  } elsif ($method eq "fastme") {
    run_nj_fastme($temp_file, $output_file, 2);
  } 
  
  my_cmd("rm $temp_file");
}

sub super_tree_superfine {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $options = $_[2];
  
  my $dir = "-d $temp_dir   ";
  if (not defined $options) {
    $options = "-r gmrp";
  }
  print("$python $superfine $dir $options -o $output_file $input_file\n");
  my_cmd("$python $superfine $dir $options -o $output_file $input_file");
}

sub super_tree_fasttree {
  my $input_file = $_[0];
  my $output_file = $_[1];

  my $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }
  convert_trees_to_alignment($input_file, "$temp_file", "fasta");
  estimate_ml_tree_fasta_fastree("$temp_file", "$output_file", "-nt -cat 1");
  my_cmd("rm $temp_file");
}

sub super_tree_mrl {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $tree_method = $_[2];
  my $log_file = $_[3];
  my $special_options = $_[4];
  my $reweight = $_[5];
  
  if (not defined $log_file or $log_file eq "") {
        $log_file = "";
  } else {
        $log_file = " -log $log_file";
  }
  if (not defined $special_options) {
      $special_options = "";
    }
  
  my $temp_file = get_temp_file();  
  if ($tree_method eq "fasttree") {
    compute_fast_mrp($input_file,"$temp_file.aln", "FASTA", "-dna");
  } else {
    compute_fast_mrp($input_file,"$temp_file.aln", "PHYLIP", "");
  }
  
  if ($tree_method eq "raxml") {
    my $temp = get_temp_file();
    if ($special_options eq "weights") {
          `java -jar \$FASTMRP/wmrp.jar $input_file $temp_file.nexus NEXUS`;
          my $weight_line = `grep weight $temp_file.nexus`;  
          my @edges = $weight_line =~ m/(\d+\.\d+):\d+/g;
          foreach my $idx (0..(scalar @edges)-1) {
            if (defined $reweight and $reweight eq "1") {
              $edges[$idx]*=100;
            }
            $edges[$idx] = int($edges[$idx]);
          }
          open(OUTPUT, ">$temp_file.weights");
          print OUTPUT "@edges";
          close(OUTPUT);
          $special_options = "-a $temp_file.weights";
        }

    my_cmd("mkdir $temp");    
    my $align_file = File::Spec->rel2abs("$temp_file.aln");
    my $cur_dir = trim(my_cmd("pwd"));        
    chdir($temp);
    my_cmd("$bin/raxmlHPC -m BINCAT -n tree -s $align_file $special_options");
    chdir($cur_dir);
    my_cmd("mv $temp/RAxML_result.tree $output_file");
    my_cmd("rm $temp -rf");
    my_cmd("rm $temp_file.*");
  } else {    
    estimate_ml_tree_fasta_fastree("$temp_file.aln", "$output_file", "-gtr -noboot -nosupport -nt $log_file");
    my_cmd("rm $temp_file.*");
  }
}

sub super_tree_tnt {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $template = $_[2];  
  
  if (not defined $template) {
    $template = "quick.template";
  }
  $temp_dir = "/scratch/cluster/namphuon/goloboff_analyses/temp/tnt/";
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }
  #$temp_file = "/scratch/cluster/namphuon/goloboff_analyses/temp/tnt/goloboff_raxml_nt_tnt_large";
  #convert_newick_to_mrp_tnt($input_file, $temp_file, "$temp_file.map", "map");  
  remap_newick($input_file, $temp_file, "$temp_file.map", "map");
  `java -jar /$fastmrp $temp_file $temp_file.nexus NEXUS`;
  #`perl -p -i -e "s/ROOT\\s+[^\\s]+\\s*//" $temp_file`;
    
  my ($filename, $basename, $empty) = fileparse($output_file);  
  find_replace("/projects/sate3/namphuon/bin/tnt/$template", "$temp_file.tnt", {"<source_tnt>", "$temp_file.nexus","<output>", $filename.".temp", "<dir>", $basename});  
  my_cmd("$tnt proc $temp_file.tnt\\; quit;");    
  convert_nexus_newick("$output_file.temp", "$output_file.temp2");
  convert_nexus_newick("$output_file.temp.trees", "$output_file.temp3");
  
  #Write a random tree with best score
  my $random_tree = read_trees("$output_file.temp3", "random");
  open(OUTPUT, ">$output_file.random_tree");
  print OUTPUT $random_tree;
  close(OUTPUT);
  
  #Get greedy consensus tree
  print my_cmd("$python /u/namphuon/bin/getGreedyConsensus.py -s $output_file.temp3 -t $output_file.temp3 -d $temp_dir> $output_file.temp4");  
  convert_nexus_newick("$output_file.temp4", "$output_file.temp5");
  
  #Remap all trees
  `$perl -p -i -e "s/'//g" $output_file.temp2`;
  `$perl -p -i -e "s/'//g" $output_file.random_tree`;
  `$perl -p -i -e "s/'//g" $output_file.temp5`;
  remap_newick("$output_file.temp2", "$output_file", "$temp_file.map", "remap");
  remap_newick("$output_file.random_tree", "$output_file.best_tree", "$temp_file.map", "remap");
  remap_newick("$output_file.temp5", "$output_file.majority", "$temp_file.map", "remap");  
  
  #Delete everything
  my_cmd("rm $temp_file $temp_file.nexus $temp_file.map $temp_file.tnt $output_file.temp $output_file.temp2 $output_file.temp3 $output_file.temp4 $output_file.temp5 $output_file.temp.trees $output_file.random_tree");
}

sub read_trees {
  my $tree_file = $_[0];
  my $tree_idx = $_[1];
  my $tree_string = `cat $tree_file`;
  
  my @trees = ($tree_string =~ m/([^;]*;)/g);
  
  if (defined $tree_idx and $tree_idx eq "random") {
    return $trees[int(rand(scalar(@trees)))];
  }
  return \@trees;
}

sub calc_mp_score {
  my $input_tree = $_[0];
  my $input_alignment = $_[1];
  my $output_file = $_[2];
  my $fix_mrp_file = $_[3];  
  
  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }  
  print $temp_file . "\n";
  
  if (defined $fix_mrp_file and $fix_mrp_file eq "alignment") {    
    Phylo::my_cmd("cp $input_alignment $temp_file.aln");
    Phylo::my_cmd("$perl -p -i -e \"s/'//g\" $temp_file.aln");
  } else {
    Phylo::my_cmd("cp $input_alignment $temp_file.aln");    
  }   
  
  if (defined $fix_mrp_file and $fix_mrp_file eq "tree") {
    remove_whitespace_names($input_tree, $temp_file.".fixed");
    $input_tree = $temp_file.".fixed";
  }
  `cp $input_tree $temp_file.tree.temp`;
  `$perl -p -i -e "s/'//g" $temp_file.tree.temp`;
  remap_newick("$temp_file.tree.temp", "$temp_file.tree", "$temp_file.map", "map");    
  remap_nexus("$temp_file.aln", "$temp_file.alignment", "$temp_file.map", "remap");
  
  
  `cat $temp_file.alignment > $temp_file`;
  `echo "begin trees;\ntree tnt = [&U]" >> $temp_file`;
  `cat $temp_file.tree >> $temp_file`;
  `echo "\nend;\n" >> $temp_file`;  
  #`rm $temp_file.tree $temp_file.alignment $temp_file.map $temp_file.aln`;
  `$perl -p -i -e "s/'//g" $temp_file`;

   find_replace("/projects/sate3/namphuon/bin/tnt/score.template", "$temp_file.tnt", {"<source_tnt>", "$temp_file","<output>", "$temp_file.score"});
   print "$tnt proc $temp_file.tnt\\; quit;\n";
   my_cmd("$tnt proc $temp_file.tnt\\; quit;");
   if (-e "$temp_file.score") {
      `mv $temp_file.score $output_file`;
    }
   my_cmd("rm $temp_file.* $temp_file");
}


sub super_tree_fasttree_tnt {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $options = $_[2];
  
  if (not defined $options) {
    $options = "-gtr -nt -cat 1";
  }

  my $temp_file = $temp_dir . "/". Phylo::get_random_name($temp_dir);
  while (-e $temp_file) {
    $temp_file = $temp_dir . "/".Phylo::get_random_name($temp_dir);
  }
  
  convert_newick_to_mrp_tnt($input_file, $temp_file, "$temp_file.map", "map");
  convert_mrp_to_alignment($temp_file, "$temp_file.fasta", "fasta","root");  
  estimate_ml_tree_fasta_fastree("$temp_file.fasta", "$output_file.temp", $options);
  remap_newick("$output_file.temp", $output_file, "$temp_file.map", "remap");
}

sub super_tree_paup {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $iter = $_[2];
  my $options = $_[3];
  if (not defined $options) {
    $options = "";
  }
  
  if (not defined ($iter)) {
    $iter = 100;
  }
  my_cmd("$python $mrp_paup -i $input_file -o $output_file -n $iter $options");
}

sub estimate_tree {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $temp = "";
  if ($machine == 0 ) {
    $temp = $raxml;
  } else {
    $temp = $raxmlV;
  } 
  
  my_cmd("$temp -m GTRCAT -n tree -s $input_file");
  my_cmd("mv RAxML_result.tree $output_file");
  my_cmd("rm RAxML_*");
}

sub estimate_tree_fasta {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $options = $_[2];
  my $temp = "";
  if (not defined $options) {
    $options = "";
  }
  if ($machine == 0 ) {
    $temp = $raxml;
  } else {
    $temp = $raxmlV;
  } 
  my $temp_name = "temp" . int(rand()*100000);
  my_cmd("mkdir $temp_name");
  chdir($temp_name);
  convert_fasta_to_phylip("../" . $input_file, "../" . $input_file . ".phylip");
  my_cmd("$temp $options -m GTRCAT -n tree -s ../$input_file" . ".phylip");
    
  my_cmd("mv RAxML_result.tree ../$output_file");
  my_cmd("rm RAxML_* ../$input_file".".phylip");
  chdir("..");
  my_cmd("rm $temp_name -rf");  
}

sub fasta_to_phylip {
  my $input = $_[0];
  my $output = $_[1];
  my $outf = " > ".$output;
  open (F, " < ".$input);
  my $taxa = -1;
  my @taxonNames = ();
  my @sequences = ();

  my $line = "";
  while($line = <F>)
  {
  if($line =~ />/)
  {
  $taxa++;
  my $name = $line;
  $name = trim($name);
  $name =~ s/[\s+|\(|\)|\,|;]//g;
  $name =~ s/,//g;
  $name =~ s/>//g;
  $taxonNames[$taxa] = $name;
  }
  else
  {
  my $seq = $line;
  $seq =~ s/\s+//g;
  if (not defined $sequences[$taxa]) {
      $sequences[$taxa] = "";
  }
  $sequences[$taxa] = $sequences[$taxa].$seq;
  }
  }

  close(F);

  for(my $i = 0; $i <= $taxa; $i++)
  {
  #print $taxonNames[$i]." ".(length($sequences[$i]))."\ n";
  }

  my $s = $taxa + 1;
  my $bp = length($sequences[0]);

  open (F, $outf);

  print F $s." ".$bp."\n";

  for(my $i = 0; $i <= $taxa; $i++)
  {
  print F $taxonNames[$i]." ".$sequences[$i]."\n";
  }
}

sub estimate_ml_tree_fasta {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $temp = "";
  if ($machine == 0 ) {
    $temp = $raxml;
  } else {
    $temp = $raxmlV;
  } 

  convert_fasta_to_phylip($input_file, $input_file . ".phylip");
  my_cmd("$temp -m GTRCAT -n tree -s $input_file" . ".phylip");
  
  
  my_cmd("mv RAxML_result.tree $output_file");
  my_cmd("rm RAxML_* $input_file".".phylip");
}
sub estimate_ml_tree_fasta_fastree {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $commands = $_[2];
  my $temp = "";
  if ($machine == 0 ) {
    $temp = "/projects/sate9/namphuon/programs/sate/satesrc-v2.2.7-2013Feb15/sate-tools-linux/fasttreeMP";
  } else {
    $temp = "/projects/sate9/namphuon/programs/sate/satesrc-v2.2.7-2013Feb15/sate-tools-linux/fasttreeMP";
  }
  if (not -e $input_file) {
    return;
  }
  
  my $temp_file = get_temp_file();
  if (defined($commands) and $commands ne "") {    
    my_cmd("$temp $commands $input_file > $temp_file");
  } else {
    my_cmd("$temp -gtr -nt $input_file > $temp_file");
  }     
  #if (-e $temp_file and verify_tree($temp_file)) {
  if (-e $temp_file and -s $temp_file != 0) {
    `mv $temp_file $output_file`;
  } else {
    `rm $temp_file`;
  }
  #} else {
  #  `rm $temp_file`;
  #}
}

sub get_random_name {
  my $target_dir = $_[0];
    
  if (not defined $target_dir) {
    $target_dir = "";
  }
  local $"="";
  my @name_str = map { ("a".."z", "A".."Z", 0..9)[rand 62] } 1..20;
  my $name = "temp@name_str";
  while (-e "$target_dir$name") {
    @name_str = map { ("a".."z", 0..9)[rand 36] } 1..20;
    $name = "temp@name_str";
  }
  return $name;
}

sub estimate_ml_tree_fasta_condor {
  my $input_file = rel2abs($_[0]);
  my $output_file = rel2abs($_[1]);
  my $keep = $_[2];
  my $options = $_[3];
  
  if (not -e $input_file) {
    return -1;
  }
  if (not defined $keep) {
        $keep = 0;
  }

  if (not defined $options) {
    $options = "";
  }
  my $model = "-m GTRCAT";
  if ($options =~ m/-m/) {
      $model = "";
    }
  my $temp = "";
  if ($machine == 0 ) {
    $temp = $raxml;
  } elsif ($options =~ m/-T/) {
    $temp = "/projects/sate3/namphuon/bin/raxmlHPC-PTHREADS-git-June9-gcc";
  } else {
    $temp = "/projects/sate3/namphuon/bin/raxmlHPC-git-June2-gcc";
  }
  
  if (not -e $input_file) {
    print "$input_file does not exist.\n";
    return;
  }
  
  

  my $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_name) {
    $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }

  my_cmd("mkdir $temp_name");

  convert_fasta_to_phylip($input_file, "$temp_name/temp.phylip");
  my $cur_dir = trim(my_cmd("pwd"));
  chdir($temp_name);
  print "$temp $model -n tree -s temp.phylip $options";
  print my_cmd("$temp $model -n tree -s temp.phylip $options");
  chdir($cur_dir);

  if (-e "$temp_name/RAxML_bestTree.tree") {
      my_cmd("mv $temp_name/RAxML_bestTree.tree $output_file");
  } elsif (-e "$temp_name/RAxML_result.tree") {
    my_cmd("mv $temp_name/RAxML_result.tree $output_file");
  }
  if ($keep) {
      my_cmd("mv $temp_name/RAxML_info.tree $output_file.RAxML_info");
      my_cmd("mv $temp_name/RAxML_log.tree $output_file.RAxML_log");
      my_cmd("mv $temp_name/RAxML_parsimonyTree.tree $output_file.RAxML_parsimonyTree");      
  } 
  my_cmd("rm $temp_name -rf");
}

sub estimate_ml_tree_raxml {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $temp = "";
  
  if (not -e $input_file) {
    print "$input_file does not exist.\n";
    return;
  }

  my $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_name) {
    $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }
  my_cmd("mkdir $temp_name");

  my $cur_dir = trim(my_cmd("pwd"));
  my ($filename, $basename, $empty) = fileparse($input_file);
  my_cmd("ln -s $input_file $temp_name/$filename");
  
  chdir($temp_name);
  print `pwd`;
  print `ls`;
  my_cmd("rm RAxML_*");
  print "/projects/sate3/namphuon/bin/raxmlHPC -m BINCAT -n tree -s $filename\n";
  my_cmd("/projects/sate3/namphuon/bin/raxmlHPC -m BINCAT -n tree -s $filename");
        chdir($cur_dir);
  my_cmd("mv $temp_name/RAxML_result.tree $output_file");
  #my_cmd("rm $temp_name -rf");
}

sub calc_ml_score {
  my $input_tree = $_[0];
  my $input_alignment = $_[1];
  my $output_file = $_[2];
  my $temp = "";
  
  if (not -e $input_tree or not -e $input_alignment) {
    print "$input_tree or $input_alignment does not exist.\n";
    return;
  }

  my $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  while (-e $temp_name) {
    $temp_name = $temp_dir . "/"  .Phylo::get_random_name($temp_dir);
  }
  #$temp_name = "fuck";
  my_cmd("mkdir $temp_name");
  randomly_refine($input_tree, "$temp_name/temp.tree");
  my_cmd("cp $input_alignment $temp_name/temp.alignment");
  Phylo::my_cmd("$perl -p -i -e \"s/'//g\" $temp_name/temp.alignment");

  my $cur_dir = trim(my_cmd("pwd"));
  
  chdir($temp_name);
  print `pwd`;
  my_cmd("rm RAxML_*");
  my_cmd("/projects/sate3/namphuon/bin/raxmlHPC -m BINGAMMA -n tree -s temp.alignment -t temp.tree -f e");
  print "/projects/sate3/namphuon/bin/raxmlHPC -m BINGAMMA -n tree -s temp.alignment -t temp.tree -f e";
        chdir($cur_dir);
  my_cmd("mv $temp_name/RAxML_log.tree $output_file");  
  my_cmd("rm $temp_name -rf");
}

sub remove_whitespace_names {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  open(OUTPUT, ">$output_file");
  my $line = `cat $input_file`;
  my @matches = ($line =~ /([^']*)('([^']*)')([^']*)/g);
   my $counter = 0;
   foreach my $match (@matches) {     
    if ($counter % 4 != 1 and $counter % 4 != 2) {
      print OUTPUT $match;
    } elsif ($counter % 4 != 2) {
      $match=~s/\s+/_/g;
      $match=~s/'//g;
      print OUTPUT $match;
    } 
     $counter++;
   }  
   close(OUTPUT);
}

sub convert_newick_to_nexus {
  my $input_file = $_[0];
  my $output_file = $_[1];
    
  my $temp_name = get_random_name($temp_dir);
  convert_newick_to_tnt($input_file, $temp_dir."/".$temp_name);
  convert_to_mrp($temp_dir."/".$temp_name, $output_file);
  my_cmd("rm $temp_dir/$temp_name");
}

sub compute_fast_mrp {
  my $source_trees = $_[0];
  my $output_file = $_[1];
  my $output_format = $_[2];
  my $options = $_[3];
  
  if (not defined $options) {
    $options = "";
  }
  
  print `java -jar $fastmrp $source_trees $output_file $output_format $options`;  
}

sub compute_fast_mrl {
  my $source_trees = $_[0];
  my $output_file = $_[1];
  
  compute_fast_mrp($source_trees, "$output_file.fasta", "FASTA", "-dna");
  print `$fasttree -gtr -nt -nosupport -nocat $output_file.fasta > $output_file`;
}

sub convert_newick_to_mrp_tnt {    
  my $source_trees = $_[0];      
  my $output_file = $_[1];
  my $remap_file = $_[2];
  my $ram = $_[3];
  
  if (not defined($ram)) {
    $ram = 16000;
  }
  
  my $temp = $temp_dir . "/" .get_random_name($temp_dir);
  remap_newick($source_trees, $temp, "$remap_file", "map");  
  convert_newick_to_tnt($temp, "$temp.data.tnt");
  my ($filename, $basename, $empty) = fileparse($output_file);
  find_replace("/projects/sate3/namphuon/bin/tnt/mrp.template", "$temp.tnt", {"<source_tnt>", "$temp.data.tnt","<ram>", "$ram", "<output>", $filename, "<dir>", $basename});  
  my_cmd("$tnt proc $temp.tnt\\; quit;");  
  my_cmd("rm $temp $temp.data.tnt $temp.tnt");
}

sub convert_to_mrp {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  my $runner = get_random_name();
  
  my_cmd("echo \"proc $tnt_file; proc $input_file; mrp; export $output_file; quit;\" > $runner");
  my_cmd("$tnt proc $runner");
  my_cmd("rm $runner");  
}

sub estimate_nj_tree_fasta {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  my $temp_name = "temp.nexus" . int(rand()*100000);
  my $temp_tree = $temp_name . ".tree";
  my $temp_tree_nexus = $temp_name . ".tree.nexus";
  Phylo::convert_fasta_to_paup("$input_file", $temp_name);
  #Now to generate the PAUP Block and append it to the file
  
  #Now use 90% of chars?
  #[Ratchet parameters: NumChar = 2  Inclusion Percentage = 0.250000 (0 characters)  Replicates = 100  Final search = no]
  open(NEXUS, ">>$temp_name");

  #Now generate PAUP BLOCK, starting with initializer
  print NEXUS "\nbegin paup;\n\tset autoclose = yes warntree = no warnreset = no notifybeep = no monitor = yes taxlabels = full;\n\tlog file = ratchet.log replace;\n\tnj;\n\tsavetrees file = $temp_tree replace = yes format = altnex;\n\n\tlog stop;\nend;\nquit warntsave = no;";
  close(NEXUS);

  my_cmd("$bin/paup -n $temp_name");  
  Phylo::convert_paup_newick($temp_tree, $output_file);
  Phylo::my_cmd("$perl -p -i -e \"s/\\[.*\\]//g\" $output_file");
  my_cmd("rm $temp_name* ratchet.log");
}

sub estimate_mp_tree_fasta {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  my $temp_name = "temp.nexus" . int(rand()*100000);
  my $temp_tree = $temp_name . ".tree";
  my $temp_tree_nexus = $temp_name . ".tree.nexus";
  Phylo::convert_fasta_to_paup("$input_file", $temp_name);
  #Now to generate the PAUP Block and append it to the file
  
  #First step is to get the number of characters
  open(NEXUS, "<$temp_name");
  my $num_chars = 0;
  while (my $line = <NEXUS>) {
    if ($line =~ m/NCHAR[\W]*=[\W]*([\d]+)/i) {
      $num_chars = int($1);
    }
  }
  close(NEXUS);
  if (int($num_chars) == 0) {
    my_cmd("rm $temp_name* ratchet.log");
    print "In: ". `pwd`;
    print "Failed to convert $input_file to $output_file\n";
    return;
    #die;
  }
  my $used_chars = int(.25*$num_chars);

  #Now use 90% of chars?
  #[Ratchet parameters: NumChar = 2  Inclusion Percentage = 0.250000 (0 characters)  Replicates = 100  Final search = no]
  open(NEXUS, ">>$temp_name");
  print NEXUS "[Ratchet parameters: NumChar = $num_chars  Inclusion Percentage = 0.250000 ($used_chars characters)  Replicates = 100  Final search = no]";

  #Now generate PAUP BLOCK, starting with initializer
  print NEXUS "\nbegin paup;\n\tset autoclose = yes warntree = no warnreset = no notifybeep = no monitor = yes taxlabels = full;\n\tlog file = ratchet.log replace;\n\tset criterion = parsimony;\n\tpset collapse = no;\n";

  print NEXUS "\n\t[!][!*** Replicate 0 (initial tree) ***]\n\thsearch addseq = random nreps = 1 rseed = " .int(rand()*10000) . " swap = TBR multrees = no dstatus = 60;\n\tsavetrees file = $temp_tree format = altnex replace;\n\tsavetrees file = $temp_tree_nexus format = nexus replace;\n";
  
  #Now generate replicates
  for (my $i=1; $i<=100; $i++) {
     print NEXUS "\n\t[!][!*** Replicate #$i ***]\n\tweights 2: ;\n\thsearch start = current swap = TBR multrees = no dstatus = 60;\n\tweights 1: all;\n\thsearch start = current swap = TBR multrees = no dstatus = 60;\n\tsavetrees file = $temp_tree format = altnex append;\n\tsavetrees file = $temp_tree_nexus format = nexus append;\n"
  } 

  #Finally, genereate consensus tree
  print NEXUS "\n\t[!][!*** Determining consensus trees ***]\n\tset MaxTrees = 201;\n\tgettrees file = $temp_tree allblocks = yes warntree = no;\n\tset criterion = parsimony;\n\tcondense collapse = no deldupes = yes;\n\tfilter best = yes;\n\tcontree all / strict = yes treefile = $temp_tree.smrp replace;\n\tcontree all / majrule = yes strict = no treefile = $temp_tree.mmrp replace;\n\tcontree all / majrule = yes strict = no le50 = yes treefile = $temp_tree.gmrp replace;\n\tsavetrees file = $temp_tree replace = yes format = altnex;\n\n\tlog stop;\nend;\nquit warntsave = no;";
  close(NEXUS);

  my_cmd("$bin/paup -n $temp_name");  
  Phylo::convert_paup_newick("$temp_tree.gmrp", $output_file);
  #Phylo::convert_paup_newick($temp_tree, $output_file);
  Phylo::my_cmd("$perl -p -i -e \"s/\\[.*\\]//g\" $output_file");
  my_cmd("rm $temp_name* ratchet.log");
}

sub find_replace {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $patterns_ref = $_[2];
  my %patterns = %$patterns_ref;

  open(INPUT,$input_file) or die "Cannot open file: $input_file\n";
  open(OUTPUT,">$output_file");

  while(my $line = <INPUT>){
    while ((my $key, my $value) = each(%patterns)){
      if (!defined($key) or (!defined($value))) {
        print "K: " . $key . " V: " . $value . "\n";
      }
           $line =~ s/$key/$value/g;
    }
    print OUTPUT $line;
  }
  
  close INPUT;
  close OUTPUT; 
}

sub check_alignment {
  my $fasta_file = $_[0];
  if (not -e $fasta_file) {
    print "File does not exist\n";
    return "File does not exist\n";
  }
  open(FILE, "<$fasta_file");  
  my $line = "";
  $line = <FILE>;
  my $length = -1;  
  while (defined $line and $line =~ />([^\n\r\f]*)/) {
    my $sequence = "";
    my $name = $1;
    $line = <FILE>;    
    while (defined $line and not $line =~ />([^\n\r\f]*)/) {
      $line =~ s/\s//g;
      $sequence = $sequence . trim($line);
      $line = <FILE>;
    }
    if ($length == -1) {
      $length = length($sequence);
    } else {
      if ($length != length($sequence)) {
        print "Error in sequence length: $name\n";
        return "Error in sequence length: $name\n";
      }
    }
    if (length(trim($sequence)) == 0) {
      print "Error empty alignment\n";
      return "Error empty alignment\n";
    
    }         
    if (not $sequence =~ /[^-?]/) {      
      print "Error all indel characters\n";
      return "Error all indel characters\n";      
    }
  }
  if ($length == -1) {
        print "Error in length, empty\n";
        return "Error in length, empty\n";
  }
  if (defined $line) {
    print "Error in format\n";
    return "Error in format\n";
  }
  close(FILE);
  return "";
}

sub subtract {
  my @array1 = @{$_[0]};
  my @array2 = @{$_[1]};
  
  my $len1 = scalar @array1;
  my $len2 = scalar @array2;
  
  if ($len1 != $len2) {
    return undef;
  }
  
  my @new_array = ();
  foreach my $i (0..($len1-1)) {
      push(@new_array, $array1[$i]-$array2[$i]);
  }
  return \@new_array;
}

sub read_phylip_file {
  my $input_file = $_[0];
  my $verify = $_[1];
  my $interleaved = $_[2];
  if (not defined $verify) {
    $verify = 1;
  }
  
  my $temp_file = get_temp_file();
  if (defined $interleaved and $interleaved == 1) {
    phylip_to_fasta_interleave("$input_file ", "$temp_file.fasta");
  } else {
    phylip_to_fasta("$input_file ", "$temp_file.fasta", 1);
  }
  my %result = %{read_fasta_file("$temp_file.fasta", $verify)};
  `rm $temp_file.fasta`;
  return \%result;
}

sub create_consensus_sequence {
  my %alignment = %{$_[0]};  
  my $gap_threshold = $_[1];  
  
  if (not defined $gap_threshold) {
    $gap_threshold = .50;
  }
  
  my @site_counts = ();
  my $taxa = scalar keys %alignment;
  my @names = keys %alignment;
  my $length = length $alignment{$names[0]};
  my $idx = 0;
  while ($idx < $length) {
    push(@site_counts, {}); 
    $idx++;    
  }
  my $counter = 0;
  foreach my $key (keys %alignment) {
    print "$counter\n";
    $counter++;
    my @results = split(//, $alignment{$key});
    $idx = 0;
    while ($idx < $length) {
      $site_counts[$idx]->{uc $results[$idx]}++;
      $idx++;
    }   
  }
  my $threshold = $gap_threshold*$taxa;
  my $consensus = "";
  foreach my $map (@site_counts) {
    my $largest = "";
    my $largest_count = -1;
    foreach my $key (keys %{$map}) {
      if ($largest_count == $map->{$key} and $key eq "-") {
        next;
      }
      if ($map->{$key} >= $largest_count) {
        $largest = $key;
        $largest_count = $map->{$key};
      }
    }
    if (not defined $map->{'-'} or $map->{'-'} <= $threshold) {
      $consensus.="$largest";
    } else {
      $consensus.="-";
    }
  }
  return $consensus;
}

sub read_fasta_file {
  my $input_file = $_[0];
  my $verify = $_[1];
  my %name_map = ();
  
  if (not defined $verify) {
    $verify = 1;
  }
  
  if (not -e $input_file) {
    print "$input_file does not exist!\n";
    exit;
  }
  open(FILE, "<$input_file");  
  my $line = "";
  $line = <FILE>;
  my $length = -1;  
  while (defined $line and $line !~ />([^\n\r\f]*)/) {
    $line = <FILE>;
  }
  while (defined $line and $line =~ />([^\n\r\f]*)/) {
    my $sequence = "";
    my $name = $1;
    $line = <FILE>;    
    while (defined $line and not $line =~ />([^\n\r\f]*)/) {
      $line =~ s/\s//g;
      $sequence = $sequence . trim($line);
      $line = <FILE>;
    }
    if ($length == -1 and $verify) {
      $length = length($sequence);
    } else {
      if ($length != length($sequence) and $verify) {
        print "Error in sequence length $name \n";
        return "Error in sequence length $name\n";
      }
    }
    if (length(trim($sequence)) == 0 and $verify) {
      print "Error empty alignment\n";
      return "Error empty alignment\n";
    
    }         
    if (not $sequence =~ /[^-?]/ and $verify) {      
      print "Error all indel characters $name\n";
      return "Error all indel characters $name\n";      
    }
    $name_map{$name} = $sequence;
  }
  close(FILE);
  return \%name_map;
}


sub find_replace_string {
  my $input = $_[0];
  my $patterns_ref = $_[1];
  my %patterns = %$patterns_ref;

  while ((my $key, my $value) = each(%patterns)){
    if (!defined($key) or (!defined($value))) {
      print "K: " . $key . " V: " . $value . "\n";
    }
    $input =~ s/$key/$value/g;
  }
  return $input;
}


sub to_upper {
  my $input = $_[0];
  my $output = $_[1];
  my_cmd("tr '[:lower:]' '[:upper:]' < $input > $output");
}

sub to_lower {
  my $input = $_[0];
  my $output = $_[1];
  my_cmd("tr '[:upper:]' '[:lower:]' < $input > $output");
}

#Converts a string array to upper case or lower case (default)
sub change_case_array {
  my @input_array = @{$_[0]};
  my $case = $_[1];
  
  if (not defined $case) {
    $case = "lc";
  }
  
  foreach my $idx (0..(scalar(@input_array) - 1)) {
    if ($case eq "lc") {
      $input_array[$idx] = lc $input_array[$idx];
    } else {
      $input_array[$idx] = uc $input_array[$idx];
    }     
  }
}

sub make_directory {
  my $directories = $_[0];
  my $dir_ref = get_directories($directories);
  my @dirs = @{$dir_ref};
  my $cur_dir = my_cmd('pwd');
  my $counter = 0;
  my $undef = 0;
  my $first = 1;
  
  umask 0007;
  foreach my $dir (@dirs) {
    if ($dir eq "" and $first == 1) {
        $undef = 1;
        $first = 0;
        next
    }
    $first = 0;
    
    if ($undef == 1 and substr($directories, 0, 1) =~ m/\//) {
        $dir = "/" . $dir;
        $undef = 0;
    }
    if ($dir eq "") {
        next;
    }
    if ((! -d $dir) && (! -e $dir)) {      
      mkdir($dir);
      chdir($dir);
      $counter++;      
    } elsif (-e $dir){      
      chdir($dir);
      $counter++;
    } 
    else {
      last;
    } 
  }  
  chdir(trim($cur_dir));  
}

sub replace_undef {
  my (@numbers) = @{$_[0]};
  my $replacement = $_[1];
  for my $idx (0..(scalar @numbers - 1)) {
    if (not defined $numbers[$idx]) {
      $numbers[$idx] = $replacement;
    }
  }
  return \@numbers;
}

sub mean_std {
    my (@numbers) = @{$_[0]};
    my $ignore_neg = $_[1];
    
    #Prevent division by 0 error in case you get junk data
    #print "Numbers  " . @numbers . "\n";
    return (undef, undef) unless(scalar(@numbers));
     
    if (defined($ignore_neg)) {
      my @new_number = ();
      foreach my $idx (0..(@numbers-1)) {
      if ($numbers[$idx] eq "nan") {
        next;
      }
            if ($numbers[$idx] >= 0) {
        push(@new_number, $numbers[$idx]);  
            }
      }
  @numbers = @new_number;
    }
    
    #Remove all undefined numbers
    my @new_number = ();
    foreach my $idx (0..(@numbers-1)) {
  if (not defined($numbers[$idx]) or $numbers[$idx] eq "nan" or not looks_like_number($numbers[$idx])) {
    next;
  }

      if (defined $numbers[$idx] ) {
    push(@new_number, $numbers[$idx]);  
  }
    }
    @numbers = @new_number;
    
    return (undef, undef) if (scalar(@numbers) == 0);

    # Step 1, find the mean of the numbers
    my $total1 = 0;
    foreach my $num (@numbers) {
  $total1 += $num;
    }
    my $mean1 = $total1 / (scalar @numbers);    

    # Step 2, find the mean of the squares of the differences
    # between each number and the mean
    my $total2 = 0;
    foreach my $num (@numbers) {
  $total2 += ($mean1-$num)**2;
    }
    my $mean2 = $total2 / (scalar @numbers);    

    # Step 3, standard deviation is the square root of the
    # above mean
    my $std_dev = sqrt($mean2);    

    #print $mean1 . " " . $std_dev . "\n";    
    return ($mean1, $std_dev, scalar @numbers);
}

sub binomial {
  my @results = @{$_[0]};
  my $positive = 0;
  my $negative = 0;
  foreach my $result (@results) {
    if (defined $result and $result > 0) {
      $positive++;
    } elsif (defined $result and $result < 0) {
      $negative++;
    } 
  }
  return $positive;
}

sub has_branch_length {
  my $input_tree = $_[0];
  my $tree_string = "";
  if (not -e $input_tree) {
    return -1;
  }
  open(TREE,"<$input_tree");
  my $line = <TREE>;
  while (defined $line) {
    $tree_string = $tree_string . $line;
    $line = <TREE>;  
  }     
   my $forest = Bio::Phylo::IO->parse(-format => 'newick', -string => $tree_string);
  my $tree = $forest->first;
  close(TREE);
  if ($tree->is_cladogram()) {
    return 0;
  } else {    
    return 1;
  }  
}

sub has_branch_length_fast {
      my $input_tree = $_[0];
  my $tree_string = "";
  if (not -e $input_tree) {
    return -1;
  }
  open(TREE,"<$input_tree");
  my $line = <TREE>;
  while (defined $line) {
    $tree_string = $tree_string . $line;
    $line = <TREE>;  
  }     
  close(TREE);
  if ($tree_string !~ /:/) {
    return 0;
  } else {    
    return 1;
  }
}

sub remove_internal_names {
  my $tree_file = $_[0];
  my $output_file = $_[1];
  
  open(INPUT, "<$tree_file");  
  open(OUTPUT, ">$output_file");  
  my $idx = 0;
  while (my $line = <INPUT>) {
    $line =~ s/\[\d+\]//g;
    print OUTPUT $line;
  }
  close(INPUT);  
  close(OUTPUT);
  
}

sub remove_supports {
  my $tree_file = $_[0];
  my $output_file = $_[1];
  
  open(INPUT, "<$tree_file");  
  open(OUTPUT, ">$output_file");  
  my $idx = 0;
  while (my $line = <INPUT>) {
    $line =~ s/\)[^\:]+:/):/g;
    print OUTPUT $line;
  }
  close(INPUT);  
  close(OUTPUT);
}


sub map_names {
  my $input_file = rel2abs($_[0]);
  my $output_file = rel2abs($_[1]);
  my $map_file = rel2abs($_[2]);
  my $directions = $_[3];
  my $type = $_[4];
  
  my %map = %{Phylo::read_mapping($map_file)};
      
  if ($directions eq "remap") {
    %map = reverse(%map);    
  } 
  if ($type eq "fasta") {
    my %alignment = %{Phylo::read_fasta_file($input_file)};
    my %new_alignment = ();
    foreach my $key (keys %alignment) {
      if (not exists $map{$key}) {
        print "Failed to find $key in map when remapping alignment\n";
        exit;
      } else {
        $new_alignment{$map{$key}}=$alignment{$key}
      }
    }
    Phylo::write_alignment(\%new_alignment, $output_file);
  } else { 
    my $tree = `cat $input_file`;
    $tree = Phylo::trim($tree);
    my $new_tree = "";
    my $reading_name = 0;
    my $name = "";
    for (my $key = 0; $key < length($tree); $key++) {
      my $letter = substr ($tree, $key, 1);
      if ($reading_name and $letter =~ m/[\(]/){
        $new_tree = $new_tree . $letter;
        next;
      } elsif ($letter =~ m/[\(,]/ and not $reading_name) {
        $reading_name = 1;
        $new_tree = $new_tree . $letter;
        next
      } elsif ($reading_name and $letter =~ m/[,:\)\[\]]/) {
        $reading_name = 0;
        if (not exists $map{$name}) {
          print "Failed to find $key in map when remapping tree\n";
          exit;
        }
        $new_tree = $new_tree . $map{$name} . $letter;
        $name = "";
        next;
      }
      if (not $reading_name) {
        $new_tree = $new_tree . $letter;
        next;
      } elsif ($reading_name and $letter =~ m/[a-zA-Z0-9]/) {
        $name = $name . $letter;
        next;
      } else {
        print "WTF!\n";
        exit;
      }
    }
    open(OUTPUT, ">$output_file");
    print OUTPUT "$new_tree";
    close(OUTPUT);
  }
}

sub rename_fasta {
  my $input_fasta = $_[0];
  my $prefix = $_[1];
  my $output_file = $_[2];
  my $output_map = $_[3];

  my %fasta = %{Phylo::read_fasta_file($input_fasta,0)};
  my $counter = 0;
  my $temp_file = Phylo::get_temp_file();
  open(MAP, ">$temp_file.map");
  open(OUTPUT, ">$temp_file.output");
  foreach my $key (keys %fasta) {
    print OUTPUT ">$prefix$counter\n$fasta{$key}\n";
    print MAP "$prefix$counter\t$key\n";
    $counter++;
  }
  close(OUTPUT);
  close(MAP);
  `mv $temp_file.map $output_map`;
  `mv $temp_file.output $output_file`;
}

sub rename_tree {
  my $tree_file = $_[0];
  my $map_file = $_[1];
  my $output_file = $_[2];
  my $direction = $_[3];
    
  my %map = %{Phylo::read_mapping($map_file, "\\s+")};  
  my %reverse = reverse %map;
  if (defined $direction and $direction < 0) {
    %map = reverse(%map);
    %reverse = reverse %map;
  }
    
  foreach my $key (keys %map) {
    #$map{$key} = $map{$key};
  }
  
  my $temp = get_temp_file();  
  open(INPUT, "<$tree_file");
  open(OUTPUT, ">$temp");
  while (my $result = <INPUT>) {
    $result = Phylo::trim($result);
    
    #Check for names
    my @taxa = @{get_taxa($result)};
    my $bad = 0;
    foreach my $taxon (@taxa) {
      if (not defined $map{$taxon}) {
        print "$taxon is not found in mapping file!\n";
        #$bad = 1;
        $map{$taxon} = $taxon;
      }
    }    
    if ($bad == 1) {
      close(INPUT);
      close(OUTPUT);
      exit;
    }
    $result =~ s/([,\(])([^(,:\)]+)/$1$map{$2}/g;
    print OUTPUT "$result\n";
  }
  close(INPUT);
  close(OUTPUT);
      
  `mv $temp $output_file`;
}

sub remap_newick {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $rename_file = $_[2];
  my $direction = $_[3];
  
  my %name_map = ();
  my @matches = ();
  my $temp = "$temp_dir/".get_random_name($temp_dir);
  remove_weights($input_file, $temp);
  my $file = "";
    
  if ($direction eq "remap") {
    open(INPUT, "<$rename_file");    
    my $results = "";
    my $str = <INPUT>;
    while (defined $str) {
      #print "str :" .  $str;
      $results = $results . $str;
      $str = <INPUT>;
    }
    close(INPUT);
    my ($name_map_ref);
    eval $results;    
    %name_map = %{$name_map_ref};
    %name_map = reverse(%name_map);
  }
    
  #Read input file
  open(INPUT, "<$temp");
  open(OUTPUT, ">$output_file");
  my $line = "";
  my $output_line = "";
  my $counter = 0;
  while (defined $line) {  
    $output_line = $line;    
    if (@matches = $line =~ /[\(,]\s*([^,\);\(]+)/smg) {
      foreach my $match (@matches) {
        my $name = trim($match);
        $name =~ s/'//g;
        if (not defined $name_map{$name}) {
          $name_map{$name} = "tt" . $counter . "f";
          $counter++;                    
        }
        $output_line =~ s/([\(,])\s*$name\s*([,\);:])/$1$name_map{$name}$2/g  
      }
    }
    print OUTPUT $output_line;
    $line = <INPUT>;    
  }
  close(INPUT);
  close(OUTPUT);
  
  if ($direction eq "map") {
    open(OUTPUT,">$rename_file");
    my $stra = Data::Dumper->Dump([ \%name_map ], [ '$name_map_ref' ]);    
    print OUTPUT $stra;
    close(OUTPUT);
  } 
  `rm $temp`;
}

sub remap_nexus {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $rename_file = $_[2];
  my $direction = $_[3];
  
  my %name_map = ();
  my @matches = ();
  my $temp = "$temp_dir/".get_random_name($temp_dir);
  my $file = "";
    
  if ($direction eq "remap") {
    open(INPUT, "<$rename_file");    
    my $results = "";
    my $str = <INPUT>;
    while (defined $str) {
      #print "str :" .  $str;
      $results = $results . $str;
      $str = <INPUT>;
    }
    close(INPUT);
    my ($name_map_ref);
    eval $results;    
    %name_map = %{$name_map_ref};
    #%name_map = reverse(%name_map);
  }
#   foreach my $key (keys %name_map) {
#     print "$key $name_map{$key}\n";
#   }
#     
  #Read input file
  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  my $line = "";
  my $output_line = "";
  my $counter = 0;
  my $matrix = 0;  
  while (defined $line) {  
    $output_line = $line;
    if (not $matrix and $output_line =~ m/\s*matrix\s*\n/) {
      $matrix=1;
    } elsif ($matrix and not ($output_line =~ m/end\s*;\s*/)) {
      $output_line =~ m/s*([^\s]+)\s+([^\s]+)/;
      if (defined $1) {
        my $taxa_name = trim($1);        
        $taxa_name =~ s/'//g;
        if (defined $name_map{$taxa_name}) {
          $output_line =~ s/$taxa_name/$name_map{$taxa_name}/;
        }
      }
    } elsif ($matrix and $output_line =~ m/end\s*;\s*/) {
      print OUTPUT $output_line;
      $line = <INPUT>;
      while (defined $line) {
        print OUTPUT $output_line;
        $line = <INPUT>;      
      }
      last;
    }
    print OUTPUT $output_line;
    $line = <INPUT>;
  }
  close(INPUT);
  close(OUTPUT);
  
  if ($direction eq "map") {
    open(OUTPUT,">$rename_file");
    my $stra = Data::Dumper->Dump([ \%name_map ], [ '$name_map_ref' ]);    
    print OUTPUT $stra;
    close(OUTPUT);
  } 
}

    
sub remap_newick_back {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $rename_file = $_[2];
  my $direction = $_[3];
  
  my %name_map = ();
  my @matches = ();
  my $temp = "$temp_dir/".get_random_name($temp_dir);
  remove_weights($input_file, $temp);
  my $file = "";
  #Read input file
  open(INPUT, "<$temp");
  my $line = "";
  while (defined $line) {  
    $file = $file . $line;
    $line = <INPUT>;
  }
  close(INPUT);
  
  if ($direction eq "map") {
    my %names = ();
    my $counter = 0;
    if (@matches = $file =~ /[\(,]\s*([^,\);\(]+)/smg) {
      foreach my $match (@matches) {
        my $name = trim($match);
        
        if (not defined $name_map{$name}) {
          $name_map{$name} = "tt" . $counter . "f";
          $counter++;
        }
      }
    }
    open(OUTPUT,">$rename_file");
    my $stra = Data::Dumper->Dump([ \%name_map ], [ '$name_map_ref' ]);    
    print OUTPUT $stra;
    close(OUTPUT);
  } elsif ($direction eq "remap") {
    open(INPUT, "<$rename_file");    
    my $results = "";
    my $str = <INPUT>;
    while (defined $str) {
      #print "str :" .  $str;
      $results = $results . $str;
      $str = <INPUT>;
    }
    close(INPUT);
    my ($name_map_ref);
    eval $results;    
    %name_map = %{$name_map_ref};
    %name_map = reverse(%name_map);
  } 

  foreach my $key (keys %name_map) {
    $file =~ s/([\(,])\s*$key\s*([\(\),;:])/$1$name_map{$key}$2/g;
  }

  open(OUTPUT, ">$output_file");  
  print OUTPUT $file;
  close(OUTPUT);
}

sub remove_weights {
  my $input_file = $_[0];
  my $output_file = $_[1];
  open(INPUT, "<$input_file");
  my $tree = "";
  my $line = "";
  while (defined $line) {
    $tree = $tree . $line;
    $line = <INPUT>;
  }
  close(INPUT);
  
  open(OUTPUT,">$output_file");
  $tree =~ s/\s*:[^,\(\);]*//g;
  $tree =~ s/\s*\)\s*[\d]+\.?[\d]*\s*[^,\(\);]*/\)/g;  
  print OUTPUT $tree;
  close(OUTPUT);
}
    
sub convert_newick_to_tnt {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  my %taxa_names = ();
  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  
  my $tree = "";
  my $line = "";
  my $trees = "";
  my $counter = 0;
  while (defined $line) {
    $tree = "";
    $line = <INPUT>;
    while ((defined($line) and $line !~ /;/)) {
      $tree = $tree . $line;
      $line = <INPUT>;
    }
    if (defined($line)) {
      $tree = "tree t$counter = [&U]\n". $tree . $line;
    }    
    $tree =~ s/:[\d]+\.?[\d]*//g;
    $tree =~ s/\)[\d]+\.?[\d]*/\)/g;
    #$tree =~ s/[\d]+//;
    $trees = $trees . $tree;
    
    my @matches;
    if (@matches = $tree =~ /[\(,\)]([^\(,\);]+)/smg) {
      foreach my $match (@matches) {
        my $name = trim($match);
        if (length($name) > 0) {
          $taxa_names{$name} = $name;
        }
      }
    }
    $counter = $counter + 1;
  }
  my $num_taxa = scalar(keys %taxa_names);
  print OUTPUT "#NEXUS begin data;\ndimensions ntax=$num_taxa nchar=1;\nmatrix\n";
  foreach my $key (keys %taxa_names) {
    print OUTPUT $key . " 0\n";
  }
  print OUTPUT ";\nend;\nbegin trees;\n";
  print OUTPUT $trees;
  print OUTPUT "\nend;";
  close(INPUT);
  close(OUTPUT);  
}
    
sub convert_trees_to_mrp {
  my $input_file = $_[0];
  my $output_file = $_[1];
  
  if (not defined($input_file) or not defined($output_file)) {
    return;
  }
  my_cmd("$python $python_bin/printMRPMatrix.py -i $input_file -o $output_file");
  Phylo::my_cmd("$perl -p -i -e \"s/characters/data/\" $output_file");  
}

sub convert_trees_to_mrp_tnt {
  my $source_trees = $_[0];
  my $output_file = $_[1];  
  
  my $temp = $temp_dir . "/" .get_random_name($temp_dir) . "fuck";
  remap_newick($source_trees, $temp, "$temp.map", "map");  
  convert_newick_to_tnt($temp, "$temp.data.tnt");
  my ($filename, $basename, $empty) = fileparse("$temp.mrp");
  find_replace("/projects/sate3/namphuon/bin/tnt/mrp.template", "$temp.tnt", {"<source_tnt>", "$temp.data.tnt", "<output>", $filename, "<dir>", $basename});
  my_cmd("$tnt proc $temp.tnt\\; quit;");  
  my_cmd("rm $temp $temp.data.tnt $temp.tnt");
  map_mrp_helper("$temp.mrp", "$output_file", "$temp.map", 1, 1);  
  ($filename, $basename, $empty) = fileparse($temp);  
  my_cmd("rm $temp.map $temp.mrp $filename.mrp.log");  
  
  
}

sub map_mrp_helper {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $map_file = $_[2];
  my $reverse_map = $_[3];
  my $root = $_[4];
  my %name_map;
  
  #Load the map
  open(INPUT, "<$map_file");    
  my $results = "";
  my $str = <INPUT>;
  while (defined $str) {
    #print "str :" .  $str;
    $results = $results . $str;
    $str = <INPUT>;
  }
  close(INPUT);
  my ($name_map_ref);
  eval $results;    
  %name_map = %{$name_map_ref};
    
  if ($reverse_map) {
    %name_map = reverse(%name_map);
  }
  
  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  my $line = "";
  my $reading = 0;
  my $taxa = 0;
  my $character = 0;
  while (defined($line)) {
    $line = <INPUT>;
    if (not defined($line)) {
      next;
    }    
    if ($line =~ m/dimensions\s+ntax\s*=\s*([\d]+)\s+nchar\s*=\s*([\d]+)/) {
      $taxa = (int($1)-1);
      $character = $2;
      $line =~ s/ntax=([\d]+)/ntax=$taxa/g
    }
    if ($line =~ m/ROOT/ and defined $root) {
      next;
    }
    
    if ($line =~ m/([^\s]+)\s+([01\?]+)\s+/) {
      my $taxa = $1;
      my $alignment = $2;
      
      if (defined $name_map{$taxa}) {
        $taxa = $name_map{$taxa};
      }
      print OUTPUT "$taxa $alignment\n";
    } else {
      print OUTPUT $line;
    }
  }
  if ($taxa == 0) {
    my_cmd("rm $output_file");
  }
  close(INPUT);
  close(OUTPUT);
}


sub convert_trees_to_alignment {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $type = $_[2];  
  
  if (not defined($input_file) or not defined($output_file)) {
    return;
  }
    
  my_cmd("$python $python_bin/printMRPMatrix.py -i $input_file -o $output_file.temp");  
  convert_mrp_to_alignment("$output_file.temp", $output_file, $type);
  my_cmd("rm $output_file.temp");
}

sub convert_mrp_to_alignment {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $type = $_[2];
  my $root = $_[3];
  
  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  my $line = <INPUT>;
  my $reading = 0;
  my $taxa = 0;
  my $character = 0;
  while (defined($line)) {
    if ($line =~ m/dimensions/) {
      $line =~ m/ntax\s*=\s*([\d]+)\s*nchar\s*=\s*([\d]+)/;
      $taxa = $1;
      $character = $2;
      #If phylip format, need to print number character/taxa
      #However, TNT adds root node to MRP matrix, so need to be aware of that.
      if ($type eq "phylip") {      
        if (defined $root) {
          $taxa = $taxa -1;
        }
        print OUTPUT "$taxa $character\n";
      }
    }
    if ($line =~ m/matrix/ and not $reading) {
      $reading = 1;          
    } elsif ($line =~ m/;/ and $reading){      
      $reading = 0;        
      last;
    } elsif ($reading) {      
      $line =~ m/'?(.*)'?\s+([^\s]+)/;
      if (defined($1)) {
        my $name = trim($1);
        my $sequence = trim($2);
        if ($type eq "fasta") {
          $sequence =~ s/0/A/g;
          $sequence =~ s/1/T/g;
          $sequence =~ s/\?/-/g;
        }
        if (not defined $root or (defined $root and ($name !~ /root/i))) {
          if ($type eq "fasta") {
            print OUTPUT ">$name\n$sequence\n";
          } elsif ($type eq "phylip") {
            print OUTPUT "$name $sequence\n";
          } 
        } else {
        } 
      }
    } 
    $line = <INPUT>;
  }
  close(INPUT);
  close(OUTPUT);
}

sub fix_noisy_fasta {
  my $input_file = $_[0];
  my $output_file = $_[1];

  open(INPUT,$input_file) or die "Cannot open file: $input_file\n";
  open(OUTPUT,">$output_file");

  while(my $line = <INPUT>){
    $line =~ s/>\s+/>/g;    
    print OUTPUT $line;
  }
  
  close INPUT;
  close OUTPUT; 
}

sub create_condor_runs {
  my $bin = $_[0];
  my @lines = @{$_[1]};
  my $output_name = $_[2];
  my $max_size = $_[3];
  
  my $header = <<END;
+Group = "GRAD"
+Project = "COMPUTATIONAL_BIOLOGY"
+ProjectDescription = "sate shit"
Universe = vanilla
Requirements = Arch == "X86_64" && Memory >= 4000 && InMastodon
executable = $bin
getEnv=True
END
  my $str = $header;
  my $idx = 0;
  my $counter = 0;
  while ($counter < scalar @lines) {
    $str.= $lines[$counter];
    $counter++;
    if (($counter % $max_size) == ($max_size-1)) {
      open(OUTPUT, ">$output_name.$idx");
      print OUTPUT $str;
      close(OUTPUT);
      $str=$header;
      $idx++;
    }
  }
  if ($str ne $header) {
    open(OUTPUT, ">$output_name.$idx");
    print OUTPUT $str;
    close(OUTPUT);    
  }
}

sub make_condor_file {
  my $output_file = $_[0];
  my $executable = $_[1];  
  my @arguments = @{$_[2]};
  my @output = @{$_[3]};
  my @error = @{$_[4]};
  my @log = @{$_[5]};
  my $options = $_[6];
  my $split = $_[7];  

  my $counter = 0;
  my $fidx = 0;
  my $file = $output_file;
  if (defined $split) {        
        $file = $output_file . ".$fidx";
  }
  if (not defined $options) {
        $options = "";
  }
  open(OUTPUT, ">$file");
  print OUTPUT "+Group = \"GRAD\"\n+Project = \"COMPUTATIONAL_BIOLOGY\"\n+ProjectDescription = \"phylogeny runs\"\n";
  print OUTPUT "Universe = vanilla\nRequirements = Arch == \"X86_64\" && InMastodon && Memory >= 4000\ngetenv = True\n$options\n";
  print OUTPUT sprintf("executable = %s\n\n", $executable);
  my $total = scalar @arguments;
  foreach my $idx (0..($total-1)) {
        if ($split and $counter % 150 == 149) {
        $fidx++;
        close(OUTPUT);
        open(OUTPUT, ">$output_file.$fidx");
        print OUTPUT "+Group = \"GRAD\"\n+Project = \"COMPUTATIONAL_BIOLOGY\"\n+ProjectDescription = \"phylogeny runs\"\n";
        print OUTPUT "Universe = vanilla\nRequirements = Arch == \"X86_64\" && InMastodon && Memory >= 4000\ngetenv = True\n$options\n";
        print OUTPUT sprintf("executable = %s\n\n", $executable);                
        }
        print OUTPUT sprintf("Arguments = %s\n", $arguments[$idx]);
        print OUTPUT sprintf("Output = %s\n", $output[$idx]);
        print OUTPUT sprintf("Error = %s\n", $error[$idx]);
        print OUTPUT sprintf("Log = %s\n", $log[$idx]);
        print OUTPUT "Queue\n\n";
        $counter++;
  }
  close(OUTPUT);  
}

sub get_tree_from_rose {
  my $input = $_[0];
  my $output = $_[1];
}

################################################
#----------------------------------------------#
################################################
sub my_cmd {
  my $print_err = "";
  my $cmdstr=$_[0];
  
  if ($testing == 1) {
    print $cmdstr . "\n";
    $print_err = "2>&1";
  }
  
  my ($package, $filename, $line) = caller;  
  my $rc=`$cmdstr $print_err`;
   if ($testing == 1) {
    print $rc . "\n";
  }
  return $rc;
}#sub my_cmd
################################################

#Trims whitespace from strings
sub trim {
  my $string = $_[0];
  if (not defined($string)) {
      return undef;
  }
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}
1;
