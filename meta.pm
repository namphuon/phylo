package Meta;

#
# Helper functions to run all metagen scripts, useful stuff for leave out experiments
#
######

use phylo;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use Data::Dumper;
use File::Path;
use FileHandle;
use strict;

my @all_genes = ( 
"dnaG",  "nusA",  "pyrG",  "rplC",  "rplF",  "rplM",  "rplS",  "rpsB",  "rpsI",  "rpsM",
"frr",   "pgk",   "rplA",  "rplD",  "rplK",  "rplN",  "rplT",  "rpsC",  "rpsJ",  "rpsS",
"infC",  "pyrg",  "rplB",  "rplE",  "rplL",  "rplP",  "rpmA",  "rpsE",  "rpsK",  "smpB");
my @all_lengths = ("100bp", "300bp");
my $base_dir = $ENV{'METAGEN'};
my $taxonomy_file = "$base_dir/data/taxonomy/all_taxon.taxonomy";
my $mapping_file = "$base_dir/data/taxonomy/species.mapping";
my @levels = ("species", "genus", "family", "order", "class", "phylum");
my $perl = $ENV{'PERL'};
my $metaphyler = $ENV{'METAPHYLER'};
my $metaphlan = $ENV{'METAPHLAN'};
my $metaphlan_db = $ENV{'METAPHLAN_DB'};
my $blast_db = $ENV{'BLAST_DB'};
my $blast = $ENV{'BLAST'};

sub run_metaphyler {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];
  
  my $temp_name = Phylo::get_temp_file();
  
  my %sequences = %{Phylo::read_fasta_file($fragment_file,0)};
  my ($total,$lengths) = (0,0);
  foreach my $length (values %sequences) {
    $lengths+=length($length);
    $total++;
  }
  my $blast = "blastn";
  if ($lengths/$total > 150) {
    $blast = "blastx"
  }
  
  
  `mkdir -p $output_directory`;
  `ln -s $fragment_file $temp_name`;
  `$perl $metaphyler $temp_name $blast $output_directory/$prefix 16`;
  if (defined $temp_name and -e $temp_name) {
    `rm $temp_name`;
    `mv $temp_name.* $output_directory/`;
  }
}

sub run_metaphlan {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];
  
  `mkdir -p $output_directory`;
  `python $metaphlan --blastdb $metaphlan_db $fragment_file $output_directory/$prefix  --nproc 16`;
  `mv $fragment_file.outfmt6.txt $output_directory`;
}

sub run_motu {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];
  chdir($output_directory);
  `perl \$WORK/programs/motu/mOTUs.pl $fragment_file --processors=16 --output-directory=$output_directory`;
  
  if (-e "$output_directory/rownames") {
    if (not -e "$output_directory/species.abundance" and -e "$output_directory/NCBI.species.abundances.gz") {
      my $file = "$output_directory/NCBI.species.abundances.gz";
      my $true = Phylo::trim(`ls -l $file |awk '{print \$11}'`."");
      $true =~ s/scaled\.species\.fraction/scaled.taxaid.fraction/;
      `cp $true $output_directory/species.abundance.gz`;
      chdir("$output_directory/");
      `gzip -d species.abundance.gz`;
    }
    if (-e "$output_directory/species.abundance" and not -e "$output_directory/species.abundance.stats") {
      convert_motu_to_classification("$output_directory/species.abundance","$output_directory/species.abundance.stats", $prefix);
    }
  }  
}

sub get_metaphlan_name {
  my $name = lc $_[0];
  my $level = $_[1];
  my %mapping = ("genus", {"baumannia", "baumannia",'ensifer','sinorhizobium'}, "species", {'ensifer_meliloti','sinorhizobium meliloti', "pelodictyon_clathratiforme", "pelodictyon_phaeoclathratiforme", "marinobacter_aquaeolei", "marinobacter_hydrocarbonoclasticus", "alkalilimnicola_ehrlichei", "alkalilimnicola_ehrlichii", "baumannia_cicadellinicola", "candidatus_baumannia_cicadellinicola", "thermococcus_kodakaraensis", "thermococcus_kodakarensis", "clostridium_difficile", "[clostridium]_difficile", "chlorobium_phaeovibrioides", "chlorobium_phaeovibrioides", "chlorobium_vibrioformis",  "prosthecochloris_vibrioformis", "dehalococcoides_ethenogenes", "dehalococcoides_mccartyi"}, "class", {"chlamydiae", 'chlamydiia',"magnetococci", "magnetococci", "cyanophyceae","cyanophyceae","flavobacteria", "flavobacteriia", "sphingobacteria", "sphingobacteriia", "thermi", "deinococci", "spirochaetes","spirochaetia"}, "phylum", {"thermi", "deinococcus-thermus", "nanoarchaea", "nanoarchaea"}, "order", {"aciduliprofundales", "aciduliprofundales", "nanoarchales", "nanoarchales", "synechococcales","synechococcales","leptospirales", "spirochaetales"}, "family", {"aciduloprofundaceae", "aciduloprofundaceae", "nanoarchaceae", "nanoarchaceae", "methanocullaceae", "methanomicrobiaceae", "gloeobacteraceae", "gloeobacteraceae", "cyanobacteriaceae", "cyanobacteriaceae", "phormidiaceae", "phormidiaceae", "borreliaceae", "spirochaetaceae", "thermoanaerobacterales_family_iii_incertae_sedis", "thermoanaerobacterales_family_iii._incertae_sedis","clostridiales_family_xi_incertae_sedis", "clostridiales_family_xi._incertae_sedis", "treponemaceae", "spirochaetaceae","synechococcaceae", "synechococcaceae", "acaryochloridaceae", "acaryochloridaceae"});
  
  if (defined $mapping{$level}->{$name}) {
    return $mapping{$level}->{$name};
  }
  return $name;
}
  


sub bin_hmmer_fragments {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];

  my $temp_name = Phylo::get_temp_file();
  `mkdir -p $temp_name`;
  
  Phylo::make_directory($output_directory);
  my @alignments = ();
  my %gene_mapping = ();
  push(@all_genes, "16S_bacteria");
  push(@all_genes, "16S_archaea");
  foreach my $gene (@all_genes) {
    push(@alignments, "$base_dir/data/genes/$gene/results/RAXML_SATE/sate.fasta");
    if (not -e "$base_dir/data/genes/$gene/results/RAXML_SATE/sate.fasta") {
      exit;
    }
    $gene_mapping{"$base_dir/data/genes/$gene/results/RAXML_SATE/sate.fasta"} = $gene;
  }

  my %forward = %{Phylo::read_fasta_file($fragment_file, 0)};
  my $temp_directory = Phylo::get_temp_file();

  Phylo::make_directory($temp_directory);
  my %reverse = ();
  my $idx = 0;
  foreach my $sequence (keys %forward) {
    print "$idx\n";
    $idx++;
    $reverse{$sequence} = Meta::reverse_dna($forward{$sequence});
  }

  Phylo::write_alignment(\%reverse, "$temp_directory/reverse.fas");
  my %forward_results = %{Meta::get_best_hmm_score(\@alignments, "$fragment_file")};
  my %reverse_results = %{Meta::get_best_hmm_score(\@alignments, "$temp_directory/reverse.fas")};

  my %genes = ();
  foreach my $key (keys %forward) {
    my $forward_gene = $gene_mapping{$forward_results{$key}->[0]};
    my $forward_value = $forward_results{$key}->[1];
    
    my $reverse_gene = $gene_mapping{$reverse_results{$key}->[0]};
    my $reverse_value = $reverse_results{$key}->[1];
    
    my $gene = $forward_gene;
    if ($gene eq "None" and $reverse_gene ne "None") {
      $gene = $reverse_gene;
    }
    elsif ($gene ne "None" and $reverse_gene ne "None") {
      if ($reverse_value > $forward_value) {
        $gene = $reverse_gene;
        $forward{$key} = $reverse{$key};
      }
    }

    if (not defined $genes{$gene}) {
      $genes{$gene} = {};
    }
    $genes{$gene}->{$key} = $forward{$key};
  }

  foreach my $key (keys %genes) {
    Phylo::write_alignment($genes{$key}, "$output_directory/$prefix.$key.fas");
  }
  `rm $temp_directory -rf`;
}

sub bin_metaphyler_fragments {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];
  
  my $temp_name = Phylo::get_temp_file();
  `mkdir -p $temp_name`;
  `$perl $metaphyler $fragment_file blastn $temp_name/blast 1`;
  
  my %seq2marker = %{Phylo::read_mapping("$base_dir/data/simulatedReads/seq2marker.tab", "\\s+")};
  my %archaea = %{Phylo::read_fasta_file("$base_dir/data/genes/16S_archaea/16S_archaea.fas", 0)};
  my %bacteria = %{Phylo::read_fasta_file("$base_dir/data/genes/16S_bacteria/16S_bacteria.fas", 0)};  
  foreach my $key (keys %archaea) {
    $seq2marker{$key} = "16S_archaea";
  }
  foreach my $key (keys %bacteria) {
    $seq2marker{$key} = "16S_bacteria";
  }

  `mkdir $output_directory/`;
  my %genes = ();
  my %frags = %{Phylo::read_fasta_file($fragment_file,0)};
  open(INPUT, "<$temp_name/blast.blastn");
  while (my $line = <INPUT>) {
    my @results = split(/\s+/, $line);
    my $gene = $seq2marker{$results[1]};
    if (not defined $genes{$gene}) {
      $genes{$gene} = {};
    }
    $genes{$gene}->{$results[0]}=$frags{$results[0]};
  }
  close(INPUT);
  foreach my $key (keys %genes) {
    Phylo::write_alignment($genes{$key}, "$output_directory/$prefix.$key.fas");
    my @profiles = </scratch/cluster/namphuon/metagen/profiles/$key/profile.*>;
    fix_directions_gene("$output_directory/$prefix.$key.fas", "$output_directory/$prefix.$key.fas", \@profiles);
  }
  `rm $temp_name -rf`;
}

sub bin_blast_fragments {
  my $fragment_file = $_[0];
  my $output_directory = $_[1];
  my $prefix = $_[2];

  my $temp_name = Phylo::get_temp_file();
  `mkdir -p $temp_name`;
  `$blast -db  $blast_db -outfmt 6 -query $fragment_file -out $temp_name/blast.out -max_target_seqs 1`;

  my %seq2marker = %{Phylo::read_mapping("$base_dir/data/simulatedReads/seq2marker.tab", "\\s+")};
  my %archaea = %{Phylo::read_fasta_file("$base_dir/data/genes/16S_archaea/16S_archaea.fas", 0)};
  my %bacteria = %{Phylo::read_fasta_file("$base_dir/data/genes/16S_bacteria/16S_bacteria.fas", 0)};  
  foreach my $key (keys %archaea) {
    $seq2marker{$key} = "16S_archaea";
  }
  foreach my $key (keys %bacteria) {
    $seq2marker{$key} = "16S_bacteria";
  }

  `mkdir -p $output_directory/`;
  my %genes = ();
  my %frags = %{Phylo::read_fasta_file($fragment_file,0)};
  open(INPUT, "<$temp_name/blast.out");
  while (my $line = <INPUT>) {
    my @results = split(/\s+/, $line);
    my $gene = $seq2marker{$results[1]};
    if (not defined $genes{$gene}) {
      $genes{$gene} = {};
    }
    $genes{$gene}->{$results[0]}=$frags{$results[0]};
  }
  close(INPUT);
  foreach my $key (keys %genes) {
    Phylo::write_alignment($genes{$key}, "$output_directory/$prefix.$key.fas");
    my @profiles = <$base_dir/data/profiles/$key/profile.*>;
    fix_directions_gene("$output_directory/$prefix.$key.fas", "$output_directory/$prefix.$key.fas", \@profiles);
  }
  `rm -rf $temp_name`;
}

sub fix_directions_gene {
  my $fragment_file = $_[0];
  my $output_file = $_[1];
  my @hmmer_models = @{$_[2]};
  
  my %forward_results = %{Meta::get_best_hmm_score_from_profile(\@hmmer_models, $fragment_file)};
  my %forward = %{Phylo::read_fasta_file("$fragment_file", 0)};
  
  my %backward = ();
  foreach my $sequence (keys %forward) {
    $backward{$sequence} = Meta::reverse_dna($forward{$sequence});
  }
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::write_alignment(\%backward, "$temp_name");    
    
  my %backward_results = %{Meta::get_best_hmm_score_from_profile(\@hmmer_models, "$temp_name")};
  `rm $temp_name`;
    
  foreach my $sequence (keys %forward) {
    if ((not defined $forward_results{$sequence} or $forward_results{$sequence}->[1] eq "NA") and $backward_results{$sequence}->[1] ne "NA") {
      $forward{$sequence} = $backward{$sequence};
    } elsif ((not defined $backward_results{$sequence} or $backward_results{$sequence}->[1] eq "NA") and $forward_results{$sequence}->[1] ne "NA") {
    } elsif ($backward_results{$sequence}->[1] > $forward_results{$sequence}->[1]) {
      $forward{$sequence} = $backward{$sequence};
    } elsif ($backward_results{$sequence}->[1] <= $forward_results{$sequence}->[1]) {
    } else {
      delete $forward{$sequence};
    }
  }
  Phylo::write_alignment(\%forward, $output_file);
}

sub convert_16S_classification_to_nbci {
  my %classifications = %{$_[0]};
  my $input_taxonomy = $_[1];
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$taxonomy_file", "lower");
  my ($staxon_ref, $slevel_ref, $skey_ref) = Phylo::read_taxonomy_mapping("$input_taxonomy", "lower");
  
  my %name_mapping = ();
  foreach my $key (keys %{$taxon_ref}) {
    $name_mapping{$taxon_ref->{$key}->[$key_ref->{'tax_name'}]} = $key;
  }  
  
  foreach my $key (keys %classifications) {
    my @results = @{$classifications{$key}};    
    my $idx = 1;
    foreach my $level (@levels) {
      if ($results[$idx] eq "NA") {
        $idx++;
        next;
      }
      if ($level eq "species") {
        $results[1] = "NA";
        $idx++;
        next;
      }
      my $name = $staxon_ref->{$results[$idx]}->[$skey_ref->{'tax_name'}];
      $name =~ s/$level\_//;
      $name =~ s/_/ /g;
      if (defined $name_mapping{$name}) {
        $results[$idx] = $name_mapping{$name};
      } else {
        $results[$idx] = "NA";
      }
      $idx++;
    }
    $classifications{$key} = \@results;
  }
  return \%classifications;
}

sub combine_mapping {
  my @input_files = @{$_[0]};
  my $mapping = {};  
  foreach my $file (@input_files) {
    my %classifications = %{Phylo::read_multiple_mapping($file, "\t", 1)};
    foreach my $key (keys %classifications) {
      $mapping->{$key} = $classifications{$key};
    }
  }
  return $mapping;
}

sub generate_nbc_classification {
  my $fragment_file = $_[0];
  my @results = @{$_[1]};
  my $output = $_[2];  
  my $optional = $_[3];  
  my $taxonomy = $_[4];
  my $translation_file = $_[5];
  
  if (not defined $taxonomy) {
    $taxonomy = "$taxonomy_file";
  }

  if (not defined $translation_file) {
    $translation_file = "$base_dir/data/nbc/translate.txt";
  }
  
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$taxonomy", "lower");
  
  
  my %fragments = %{Phylo::read_fasta_file($fragment_file, 0)};
  my %best = ();
  my %translation = %{Phylo::read_mapping("$translation_file", ",")};
  my $idx = 0;
  foreach my $result (@results) {
    #my $temp_file = Phylo::get_temp_file();
    #`grep -v '\\-inf' $result > $temp_file`;
    my %mapping = %{Phylo::read_mapping($result, " ")};    
    $result =~ m/-15-(.*).txt/;
    my $genome = $1;
    foreach my $key (keys %mapping) {
      if (not defined $best{$key}) {
        $best{$key} = [$key, -9999999999999, "NA"];
      }    
      if ($mapping{$key} eq "-inf") {
        next;
      } else {
        if ($best{$key}->[1] < $mapping{$key}) {
          $best{$key} = [$key, $mapping{$key}, $translation{$genome}];
        }
      }
    }
    #`rm $temp_file`;
    $idx++;
    if ($idx % 5 == 0) {
      print "$idx " . scalar @results . "\n";
    }

  }
  open(OUTPUT, ">$output");
  print OUTPUT "#Name\t#Species\t#Genus\t#Family\t#Order\t#Class\t#Phylum\n";
  foreach my $key (keys %best) {
    if ($best{$key}->[2] eq "NA") {
      print OUTPUT "$key\tNA\tNA\tNA\tNA\tNA\tNA\n";
    } else {
      my $threshold = -23.7*length($fragments{$key})+490;
      if ($best{$key}->[1] >= $threshold or defined $optional) {
        print OUTPUT "$key";
        if (not defined $taxon_ref->{$best{$key}->[2]}) {
          next;
        }
        my @line = @{$taxon_ref->{$best{$key}->[2]}};
        foreach my $level (@levels) {
          my $clade = $line[$key_ref->{$level}];
          if ($clade eq "") {
            $clade = "NA";
          }
          print OUTPUT "\t$clade";          
        }
        print OUTPUT "\n";
      }
    }
  }
  close(OUTPUT);  
}

sub run_phylopythia {
  my $input_fragment = $_[0];
  my $output_result = $_[1];
  my $temp_dir = Phylo::get_temp_file();
  mkdir($temp_dir);
  `ln -s $input_fragment $temp_dir/temp.fas`;  
  `ruby /projects/sate8/namphuon/tools/phylopythia/svm_struct_phylo_distr/scripts/predict.rb $temp_dir/temp.fas /projects/sate9/namphuon/programs/phylopythia/project2/sample_config.txt`;
  if (-e "$temp_dir/temp.fas.nox.fna.PP.out") {
    my %classifications = %{read_phylopythia_classification("$temp_dir/temp.fas.nox.fna.PP.out", 1, $taxonomy_file)};
    Phylo::write_multiple_mapping(\%classifications, $output_result, "\t", "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum");    
  } else {
  }
  `rm $temp_dir -rf`;
}

sub read_phylopythia_classification {
  my $input_file = $_[0];
  my $convert = $_[1];
  my $taxonomy = $_[2];
  
  if (not defined $taxonomy) {
    $taxonomy = $taxonomy_file;
  }
  
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$taxonomy", "lower");
  my %name_mapping = ();
  if (not defined $convert or $convert == 0) {
    foreach my $key (keys %{$taxon_ref}) {
      $name_mapping{$key} = $key;
    }
  } else {
    foreach my $key (keys %{$taxon_ref}) {
      $name_mapping{$taxon_ref->{$key}->[$key_ref->{'tax_name'}]} = $key;
    }  
  }
  
  my %mapping = ();
  open(INPUT, "<$input_file");
  my $line = <INPUT>;
  my $begin = 0;
  while ($line = <INPUT>) {
    if (not $begin and $line =~ m/root/) {
	 $begin = 1;
    } elsif ($begin) {
	 my @classification = ("name", "NA", "NA", "NA", "NA", "NA", "NA");
	 $line = lc Phylo::trim($line);
	 my @results = split(/\t/, $line);
	 $classification[0] = $results[0];
	 foreach my $idx (3..7) {
	   if ($idx >= scalar @results) {
		last;
	   }
	   if ($results[$idx] ne "") {		
		if (defined $convert and $convert == 1) {
		  if (defined $name_mapping{$results[$idx]}) {
		    $results[$idx] = $name_mapping{$results[$idx]}
		  } else {
		    print "Unable to find $results[$idx] of clade $levels[8-$idx]\n";
		  } 
		}
		$classification[9-$idx] = $results[$idx];
	   }
	 }
	 $mapping{$results[0]} = \@classification;
    }
  }
  return \%mapping;
}

#Runs the RDP classifier
sub run_rdp_classifier {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $convert = $_[2];
  #my $taxonomy = $_[2];
  my $temp_file = Phylo::get_temp_file();
  `java -jar /projects/sate9/namphuon/programs/rdp_classifier_2.4/rdp_classifier-2.4.jar -f fixrank g 16srrna -o $temp_file -q $input_file`;
  convert_rdp_classification_file($temp_file, $output_file, 0.80, $convert);  
}

#Reads the RDP classifier output to form taxonomy
sub convert_rdp_classification_file {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $threshold = $_[2];
  my $convert = $_[3];
  
  my ($bacteria_taxon_ref, $bacteria_level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$base_dir/data/RDP_10_28/bacteria.taxonomy", "lower");
  my ($archaea_taxon_ref, $archaea_level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$base_dir/data/RDP_10_28/archaea.taxonomy", "lower");
  
  my %bacteria_name_mapping = ();
  foreach my $key (keys %{$bacteria_taxon_ref}) {
    $bacteria_name_mapping{$bacteria_taxon_ref->{$key}->[$key_ref->{'tax_name'}]} = $key;
  }
  my %archaea_name_mapping = ();
  foreach my $key (keys %{$archaea_taxon_ref}) {
    $archaea_name_mapping{$archaea_taxon_ref->{$key}->[$key_ref->{'tax_name'}]} = $key;
  }

  my %levels_map = ("species",0, "genus",1, "family",2, "class",3, "order",4, "phylum",5);
  open(RDP, "<$input_file");
  my $line = <RDP>;
  open(OUTPUT, ">$output_file");
  print OUTPUT "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n";
  local $" = "\t";
  while (defined $line) {    
    $line = Phylo::trim($line);      
    $line =~ s/"//g;    
    my @results = split(/\t+/, $line);
    my $ranks = ["NA", "NA", "NA", "NA", "NA", "NA"];
    my $idx = 1;
    while ($idx < scalar @results) {            
      my $level = $results[$idx+1];
      my $name = "$level\_". lc $results[$idx];
      $name =~ s/\s+/_/g;
      $name =~ s/\//_/g;
      my $support = $results[$idx+2];
      
      my $type = "bacteria";
      if (defined $levels_map{$level} and $support >= $threshold) {
        if (not defined $bacteria_name_mapping{"$name"} and not defined $archaea_name_mapping{"$name"}) {
          print "Unable to find $name\n$line\n";
          foreach my $i (0..5) {
            $ranks->[$i] = -1; 
          }
          last;
        }
        if (not defined $bacteria_name_mapping{"$name"}) {
          $ranks->[$levels_map{$level}] = $archaea_name_mapping{"$name"};
          if (defined $convert) {
            $ranks->[$levels_map{$level}] = $archaea_taxon_ref->{$archaea_name_mapping{"$name"}}->[$key_ref->{'tax_name'}];
            $ranks->[$levels_map{$level}] =~ s/^([^_]+_)//g
          }
        } else {
          $ranks->[$levels_map{$level}] = $bacteria_name_mapping{"$name"};
          if (defined $convert) {
            $ranks->[$levels_map{$level}] = $bacteria_taxon_ref->{$bacteria_name_mapping{"$name"}}->[$key_ref->{'tax_name'}];
            $ranks->[$levels_map{$level}] =~ s/^([^_]+_)//g
          }
        }
      }
      $idx+=3;
    }
    print OUTPUT "$results[0]\t@{$ranks}\n";
    $line = <RDP>;
  }
  close(OUTPUT);
  close(RDP);
}

sub get_ranks_from_node_dmp {
	my %names = %{$_[0]};
	my $node_dmp = "$base_dir/data/taxonomy/nodes.dmp";
	my %map = ();
	open(INPUT, "<$node_dmp");
	while (my $line = <INPUT>) {
		my @results = split(/\s*\|\s*/, $line);
		$map{$results[0]} = \@results;
	}
	close(INPUT);
	
	my $counter = 0;
	my %levels_map = map {$_=>$counter++} @levels;
	
	foreach my $key (keys %names) {				
		my $ranks = ["NA", "NA", "NA", "NA", "NA", "NA"];
		my $id = $key;
		while ($id != 1) {
			my $row = $map{$id};
			if (defined $levels_map{$row->[2]}) {
				$ranks->[$levels_map{$row->[2]}] = $id;
				if ($row->[2] eq "phylum") {
					last;
				}
			}
			$id = $row->[1];			
		}
		$names{$key} = $ranks;
	}
	return \%names;
}

sub generate_bowtie_mapping {  
  my $seq_to_marker = $_[0];
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping($taxonomy_file, "lower");
  my $mapping = Phylo::read_multiple_mapping($seq_to_marker, "\t");
  foreach my $key (keys %{$mapping}) {
    $key =~ m/(.+)_\d+_\d+/;
    my $name = $1;
    if (not defined $name) {
      print "Error parsing $name\n";
      exit;
    }
    my $results = $mapping->{$key};
    delete $mapping->{$key};      
    $mapping->{$name} = $results;
  }
  foreach my $fungus ("Ca21chr1", "Ca21chr2", "Ca21chr3", "Ca21chr4", "Ca21chr5", "Ca21chr6", "Ca21chr7", "Ca21chrR", "Ca19-mtDNA") {
    my @line = ($fungus);
    foreach my $level (@levels) {
      my $clade = $taxon_ref->{5476}->[$key_ref->{$level}];
      if ($clade eq "") {
        $clade = "NA";
      }
      push(@line, $clade)
    }
    $mapping->{$fungus} = \@line;
  }  
  return $mapping;
}

sub read_bowtie {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $seq_to_marker = $_[2];  
    
  my $mapping = undef;  
  if (defined $seq_to_marker) {
    $mapping = generate_bowtie_mapping($seq_to_marker);
  }

  open(INPUT, "<$input_file");
  open(OUTPUT, ">$output_file");
  local $" = "\t";
  my $current = undef;
  while (my $line = <INPUT>) {      
    if ($line =~ m/^@/) {
      next;
    }      
    my @results = split(/\t/, $line);    
    #Skip if fragment doesn't match anything
    if ($results[2] =~ m/\*/) {
      next;
    }
    if ($results[2] =~ m/^Ca/) {
      my @copy = @{$mapping->{$results[2]}};
      $copy[0] = $results[0];
      print OUTPUT "@copy\n";
      next;
    }
    my @source = split(/\|/, $results[2]);    
    my $frag = $source[3];
    $frag =~ s/\.\d+//g;
    if (not defined $mapping) {
      print OUTPUT "$results[0]\t$frag\n";
    } else {
      if (defined $mapping->{$frag}) {
        if (not defined $current) {
          my @copy = @{$mapping->{$frag}};
          $current = [$results[0], \@copy];
          $current->[1]->[0] = $results[0];
        } else {
          if ($current->[0] ne $results[0]) {
            print OUTPUT "@{$current->[1]}\n";
            $current->[0] = $results[0];
            my @copy = @{$mapping->{$frag}};
            $current->[1] = \@copy;
            $current->[1]->[0] = $results[0];
          } else {
            foreach my $idx (1..6) {
              if ($mapping->{$frag}->[$idx] ne $current->[1]->[$idx]) {
                $current->[1]->[$idx] = "NA";
              } else {
                last;
              }
            }
          }
        }
      } else {
        print STDERR "Unable to resolve $line\n";
      }
    }
  }
  print OUTPUT "@{$current->[1]}\n";
  close(INPUT);
  close(OUTPUT);
}

sub get_classification_counts {
  my $classification_file = $_[0];
  my $filter_na = $_[1];
  #my %classifications = %{Phylo::read_multiple_mapping($classification_file, "\t", 1)};
  my %counts = ();  
  open(INPUT, "<$classification_file");
  my $line = <INPUT>;
  while ($line = <INPUT>) {
    $line = Phylo::trim($line);
    my @results = split(/\t/, $line);
    if (defined $filter_na and $filter_na == 1) {
      if ($line =~ /NA\tNA\tNA\tNA\tNA\tNA/ or $line =~ /-1\t-1\t-1\t-1\t-1\t-1/) {
        next;
      }
    }
    foreach my $idx (0..(scalar @levels-1)) {
      if (not defined $counts{$levels[$idx]}) {
        $counts{$levels[$idx]} = {};
      }
      if ($results[$idx+1] eq "") {
        $results[$idx+1] = 'NA';
      }
      $counts{$levels[$idx]}->{$results[$idx+1]}++;      
    }    
  }
  # foreach my $frag (keys %classifications) {
    # my @results = @{$classifications{$frag}};
    # foreach my $idx (0..(scalar @levels-1)) {
      # if (not defined $counts{$levels[$idx]}) {
        # $counts{$levels[$idx]} = {};
      # }
      # $counts{$levels[$idx]}->{$results[$idx+1]}++;      
    # }
  # }
  return \%counts;
}

#Returns the best hmm score for each query sequence
sub get_best_hmm_score {
  my @alignment_files = @{$_[0]};
  my $query_file = $_[1];
  my %query_sequences = %{Phylo::read_fasta_file($query_file, 0)};
  my @query_names = sort keys %query_sequences;
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  my %best_scores = ();
  foreach my $idx (0..(scalar @alignment_files-1)) {    
    Phylo::hmmr_profile($alignment_files[$idx], "$temp_name/profile.$idx");
    #Phylo::hmmr_search($query_file, "$temp_name/profile.$idx", "$temp_name/search.$idx", "100000", "--max");
    Phylo::hmmr_search($query_file, "$temp_name/profile.$idx", "$temp_name/search.$idx", "1", "");
    my %scores = %{Phylo::get_hmmr_query_scores("$temp_name/search.$idx", "$query_file")};
    
    foreach my $query (@query_names) {
      if (not defined $best_scores{$query} and defined $scores{"$query\tbit_score"}) {
        $best_scores{$query} = [$alignment_files[$idx], $scores{"$query\tbit_score"}, $scores{"$query\te_value"}]; 
      } elsif (defined $scores{"$query\tbit_score"}) {
        if ($best_scores{$query}->[1] < $scores{"$query\tbit_score"}) {
          $best_scores{$query} = [$alignment_files[$idx], $scores{"$query\tbit_score"}, $scores{"$query\te_value"}]; 
        }
      }
    }
  }
  
  my @missing = keys %{Phylo::difference(\%query_sequences, \%best_scores)};
  foreach my $miss (@missing) {
    $best_scores{$miss} = ["None", "NA", "NA"];
  }
  `rm $temp_name -rf`;
  return \%best_scores;
}

#Returns the best hmm score for each query sequence
sub get_best_hmm_score_from_profile {
  my @profiles = @{$_[0]};
  my $query_file = $_[1];
  my %query_sequences = %{Phylo::read_fasta_file($query_file, 0)};
  my @query_names = sort keys %query_sequences;
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  my %best_scores = ();
  foreach my $idx (0..(scalar @profiles-1)) {    
    Phylo::hmmr_search($query_file, "$profiles[$idx]", "$temp_name/search.$idx", "1", "--cpu 8");
    my %scores = %{Phylo::get_hmmr_query_scores("$temp_name/search.$idx", "$query_file")};
    
    foreach my $query (@query_names) {
      if (not defined $best_scores{$query} and defined $scores{"$query\tbit_score"}) {
        $best_scores{$query} = [$profiles[$idx], $scores{"$query\tbit_score"}, $scores{"$query\te_value"}]; 
      } elsif (defined $scores{"$query\tbit_score"}) {
        if ($best_scores{$query}->[1] < $scores{"$query\tbit_score"}) {
          $best_scores{$query} = [$profiles[$idx], $scores{"$query\tbit_score"}, $scores{"$query\te_value"}]; 
        }
      }
    }
  }
  
  my @missing = keys %{Phylo::difference(\%query_sequences, \%best_scores)};
  foreach my $miss (@missing) {
    $best_scores{$miss} = ["None", "NA", "NA"];
  }
  `rm $temp_name -rf`;
  return \%best_scores;
}


sub get_best_hmm_score_from_stats {
  my $input_fragments = $_[0];
  my $stats_file = $_[1];
  my $output_directory = $_[2];
  
  Phylo::make_directory($output_directory);
    
  my %forward = %{Phylo::read_fasta_file($input_fragments, 0)};
  my %reverse = ();
  foreach my $key (keys %forward) {
    $reverse{$key} = Meta::reverse_dna($forward{$key});
  }
        
  #Load file
  my ($in) = new FileHandle "<$stats_file";
  local($/) = "";
  my($str) = <$in>;
  close $in;
  
  my($scores_ref);
  eval $str;
  my (%hash) = %$scores_ref;
  my %genes = ();
  my $forward = 0;
  my $backward = 0;
  my $total = 0;
  foreach my $key (keys %hash) {
    my $best_gene = "None";
    my $best_value = -100000000000000;
    my $direction = "None";
    foreach my $element (@{$hash{$key}}) {
      if ($best_value < $element->[1]) {
        $best_gene = $element->[0];
        $best_value = $element->[1];
        $direction = $element->[3];
      }
    }
    $best_gene =~ m/profiles\/([^\/]+)\//;
    my $gene = $1;
    if (not defined $genes{$gene}) {
      $genes{$gene} = {};
    }
    if ($direction eq "forward") {
      $genes{$gene}->{$key} = $forward{$key};
      $forward++;
    } elsif ($direction eq "reverse") {
      $genes{$gene}->{$key} = $reverse{$key};
      $backward++;
    } else {
      print "Unable to match gene!\n";        
      exit;
    }      
  }
  my $matched = $forward+$backward;
  print "Forward: $forward Backward: $backward Matched: $matched Total: $total\n";
  foreach my $gene (keys %genes) {
    Phylo::write_alignment($genes{$gene}, "$output_directory/$gene.fas");
  }  
}


sub get_taxon_name_mapping {
  my %taxonomy = %{$_[0]};
  my %name_mapping = ();
  foreach my $key (keys %taxonomy) {
    $name_mapping{lc $taxonomy{$key}->[3]}=$key;
  }
  return \%name_mapping;
}

sub refine_mappings {
  my $initial_file = $_[0];
  my $refined_file = $_[1];
  my $output_file = $_[2];
  my $option = $_[3];
  
  if (not defined $option) {
    $option = "restrict";
  }

  if (not -e $initial_file) {
    `cp $refined_file $output_file`;
    return;
  }
  my %initial_mapping = %{Phylo::read_multiple_mapping(rel2abs($initial_file), "\t", 1)};
  my %refined_mapping = %{Phylo::read_multiple_mapping(rel2abs($refined_file), "\t", 1)};
  my %merged_mapping = ();

  if ($option eq "restrict") {
    my %difference = %{Phylo::difference(\%initial_mapping, \%refined_mapping)};
    foreach my $key (keys %difference) {
      delete $initial_mapping{$key};
    }
  }
  
  my %keys = %{Phylo::union(\%initial_mapping, \%refined_mapping)};
  
  foreach my $key (keys %keys) {
    #If the key is in one file, but not the other, take the existing one.
    if (exists $initial_mapping{$key} and not exists $refined_mapping{$key} or not defined $refined_mapping{$key}) {
      $merged_mapping{$key} = $initial_mapping{$key};
      next;
    } elsif (not exists $initial_mapping{$key} and exists $refined_mapping{$key} or not defined $initial_mapping{$key}) {
      $merged_mapping{$key} = $refined_mapping{$key};
      next;
    }        

    my @initial = @{$initial_mapping{$key}};
    my @refined = @{$refined_mapping{$key}};
    #Now, find out the lowest level classified both placements
    my $initial_idx = 7;
    my $refined_idx = 7;
    
    foreach my $idx (1..6) {
      if (($initial[$idx] ne "NA" and $initial[$idx] ne "UNCLASSIFIED") and $idx <= $initial_idx) { 
        $initial_idx = $idx;
      }
      if (($refined[$idx] ne "NA" and $refined[$idx] ne "UNCLASSIFIED") and $idx <= $refined_idx) { 
        $refined_idx = $idx;
      }      
    }
    
    #No initial placement
    if ($initial_idx == 7) {
      $merged_mapping{$key} = $refined_mapping{$key};
      next;      
    }
    
    #Now check if initial placement is more specific than refined placement, if so, reject.    
    if ($initial_idx <= $refined_idx) {
      $merged_mapping{$key} = $initial_mapping{$key};
      next;
    }        
    
    #Now check if initial placement is more less specific than refined placement, if so accept, else reject
    if ((lc $initial[$initial_idx]) ne (lc $refined[$initial_idx])) {
      $merged_mapping{$key} = $initial_mapping{$key};
    } else {
      $merged_mapping{$key} = $refined_mapping{$key};
    }    
  }
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  open(OUTPUT, ">$temp_name/output.tab");
  print OUTPUT "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n";
  my @frags = sort(keys %merged_mapping);
  foreach my $frag (@frags) {
    my @results = @{$merged_mapping{$frag}};
    my $len = scalar @results;
    print OUTPUT "$results[0]";
    foreach my $result (1..($len-1)) {
      print OUTPUT "\t$results[$result]";
    }
    print OUTPUT "\n";
  }
  close(OUTPUT);
  `mv $temp_name/output.tab $output_file`;
  `rm $temp_name -rf`;
}

sub read_stats_file_helper {
  my $stats_file = $_[0];
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  my %precision = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %total = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %classified_perc = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %classified_count = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %stats = ("total", \%total, "classified_perc", \%classified_perc, "precision", \%precision, "classified_count", \%classified_count);
  open(INPUT, "<$stats_file");
  my $line = <INPUT>;
  while (defined $line) {
    $line = Phylo::trim($line);
    my @results = split(/\s+/, $line);
    $stats{$results[0]}->{$results[1]}=$results[2];
    $line = <INPUT>;
  }
  close(INPUT);
  return (\%stats);
}



sub read_stats_file {
  my $stats_file = $_[0];
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  my $stats = Meta::read_stats_file_helper($stats_file);
  my %hit = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %total = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %classify = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  foreach my $level (@levels) {
    my $p = $stats->{'precision'}->{$level};
    my $c = $stats->{'classified_count'}->{$level};
    my $t = $stats->{'total'}->{$level};
    
    my $pc = $p/100.0*$c;
    $hit{$level}+=$pc;
    $total{$level}+=$t;
    $classify{$level}+=$c;    
  }
  return (\%hit, \%classify, \%total);
}


sub rename_classifications {
  my $input_classification = rel2abs($_[0]);
  my $output_classification = rel2abs($_[1]);
  my $taxonomy = $_[2];
  my $direction = $_[3];
  if (not defined $taxonomy) {
    $taxonomy = $taxonomy_file;
  }
  
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping($taxonomy, "lower");
  my %taxon_map = %{$taxon_ref};
  my %level_map = %{$level_ref};
  my %key_map = %{$key_ref};
  my %name_map = %{Meta::get_taxon_name_mapping(\%taxon_map)};
  if ($direction eq "reverse") {
    %name_map = reverse(%name_map);
  }

  my %input_map = %{Phylo::read_mapping($input_classification, "\\t")};
  
  local $"="\t";
  my $temp_name = Phylo::get_temp_file();
  open(OUTPUT, ">$temp_name");
  
  if (defined $input_map{"fragment"}) {
    my @initial = @{$input_map{"fragment"}};
    print OUTPUT "@initial\n";
    delete $input_map{"fragment"};
  }
  
  my @frags = sort keys %input_map;    
  foreach my $frag (@frags) {
    my @results = @{$input_map{$frag}};
    foreach my $idx (1..6) {
      if (defined $name_map{lc $results[$idx]}) {
        $results[$idx] = $name_map{lc $results[$idx]};
      } elsif (uc $results[$idx] eq "UNCLASSIFIED" or uc $results[$idx] eq "NA") {
        $results[$idx] = "NA";
      } else {
        print "UNKNOWN CLADE! @results, $idx\n";
        $results[$idx] = "NA";
      }      
    }
    print OUTPUT "@results\n";    
  }
  close(OUTPUT);
  if (-e $temp_name) {
    `mv $temp_name $output_classification`;
  }
}

#Reroot a tree using reference package
sub reroot_tree {
  my $refpkg = rel2abs($_[0]);
  my $input_tree = rel2abs($_[1]);
  my $output_tree = rel2abs($_[2]);
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  my ($refname, $basename, $empty) = fileparse($refpkg);
  my ($inputname, $basename, $empty) = fileparse($input_tree);
  
  `cp $refpkg $temp_name/ -r`;
  `cp $input_tree $temp_name/$refname/$inputname`;
  `perl -p -i -e 's/"tree_file"\\s*:\\s+"[^"]+"/"tree": "$inputname"/g' $temp_name/$refname/CONTENTS.json`;
  `perl -p -i -e 's/"tree"\\s*:\\s+"[^"]+"/"tree": "$inputname"/g' $temp_name/$refname/CONTENTS.json`;
  `rppr-64-1.10.alpha reroot -c $temp_name/$refname -o $temp_name/temp.tree`;
  if (-e "$temp_name/temp.tree") {
    `mv $temp_name/temp.tree $output_tree`;
  }
  `rm $temp_name -rf`;
}

#Returns whether the json file contains placements or not
sub is_empty_jsonfile {
  my $input_json = rel2abs($_[0]);
  
  my $results = `grep placements $input_json`;
  Phylo::trim($results);
  if ($results =~ m/"placements": \[\],/) {
    return 1;
  }
  return 0;
}

#Generates classification from an input json file, but from full reference package
sub generate_classifications_2 {
  my $input_json = rel2abs($_[0]);
  my $refpkg = rel2abs($_[1]);
  my $taxonomy = rel2abs($_[2]);
  my $threshold = $_[3];
  my $output = rel2abs($_[4]);
  my $optional = $_[5];
  
  my %levels_map = ("species", 0, "genus", 1, "family", 2, "order", 3, "class", 4, "phylum", 5);
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$taxonomy", "lower");  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  `guppy classify -c $refpkg -o test.classification --out-dir $temp_name $input_json`;
  if (-e "$temp_name/test.classification") {
    my %classification = ();
    my %classification_best = ();
    open(INPUT, "<$temp_name/test.classification");
    my $line = <INPUT>;
    #Name, level, clade, threshold
    my @current_fragment = ("NA", "NA", "NA", "-1");
    while ($line = <INPUT>) {
      $line = Phylo::trim($line);      
      my @results = split(/\s+/, $line);
      if ($current_fragment[0] eq "NA") {
        $current_fragment[0] = $results[0];
      } elsif ($current_fragment[0] ne $results[0]) {
        #Fill out classification
        my @lineage = ("$current_fragment[0]","NA","NA","NA","NA","NA","NA");
        foreach my $clade (keys %levels_map) {
          my $name = $taxon_ref->{$current_fragment[2]}->[$key_ref->{$clade}];
          if ($name ne "") {
            $lineage[$levels_map{$clade}+1]=$name;
          }
        }
        $classification{$current_fragment[0]} = \@lineage;
        @current_fragment = ("NA", "NA", "NA", "-1");
      }
      if ($results[1] eq $results[2] and defined $levels_map{$results[1]} and $results[4] >= $threshold and ($current_fragment[3] < $results[4] or ($levels_map{$results[1]} < $levels_map{$current_fragment[1]}))) {
        @current_fragment = ($current_fragment[0], $results[1], $results[3],  $results[4]);
      }
    }
    close(INPUT);
    
    my @lineage = ("$current_fragment[0]","NA","NA","NA","NA","NA","NA");
    foreach my $clade (keys %levels_map) {
      my $name = $taxon_ref->{$current_fragment[2]}->[$key_ref->{$clade}];
      if ($name ne "") {
        $lineage[$levels_map{$clade}+1]=$name;
      }
    }
    $classification{$current_fragment[0]} = \@lineage;
    
    Phylo::write_multiple_mapping(\%classification, "$output", "\t", "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum");
    if (defined $optional) {
      my %seq = %{Phylo::read_fasta_file($optional, 0)};    
      open(OUTPUT, ">>$output");
      foreach my $key (sort(keys %seq)) {
        if (not defined $classification{$key}) {
          print OUTPUT "$key\tNA\tNA\tNA\tNA\tNA\tNA\n";
        }
      }
      close(OUTPUT);
    }
  }
  `rm $temp_name -rf`;
}


#Generates classification from an input json file
sub generate_classifications {
  my $input_json = rel2abs($_[0]);
  my $refpkg = rel2abs($_[1]);
  my $threshold = $_[2];
  my $output = rel2abs($_[3]);
  my $rename = $_[4];      
  my $optional = $_[5];  
  my $version = $_[6];
  
  my $temp_name = Phylo::get_temp_file();
  Phylo::make_directory($temp_name);
  
  my $result = `grep version $input_json`;
  Phylo::trim($result);
  
  if ($result =~ m/"version": 3/ or (defined $version and $version ne "")) {
		Phylo::make_db($input_json, $refpkg, "$temp_name/output.db", "", "new");
  } else {
		Phylo::make_db($input_json, $refpkg, "$temp_name/output.db");
  } 
  
  
  Phylo::generate_classification("$temp_name/output.db", "$temp_name/output.classification", $threshold);
  
  if (not -e "$temp_name/output.classification") {
		print "No classifications generated\n";
		exit;
  }
  my $class = "$temp_name/output.classification";
  if (defined $rename and -e $rename) {
    Meta::rename_classifications("$temp_name/output.classification", "$temp_name/renamed.classification", $rename);
    $class = "$temp_name/renamed.classification";
  }
  if (defined $optional) {
    my %seq = %{Phylo::read_fasta_file($optional, 0)};
    my %assignments = %{Phylo::read_mapping($class, "\\s+")};
    open(OUTPUT, ">>$class");
    foreach my $key (sort(keys %seq)) {
      if (not defined $assignments{$key}) {
        print OUTPUT "$key\tNA\tNA\tNA\tNA\tNA\tNA\n";
      }
    }
    close(OUTPUT);
  }
  if (-e $class) {
    `mv $class $output`;
  }
  `rm $temp_name -rf`;    
}


#Converts Megan output to a classification file
sub convert_megan_to_tab {
  my $input_file = rel2abs($_[0]);
  my $output_file = rel2abs($_[1]);
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping($taxonomy_file, "lower");
  my %taxonomy = %{$taxon_ref};
  my %level_map = %{$level_ref};
  my %key_map = %{$key_ref};
  
  my %assignments = %{Phylo::read_mapping($input_file, "\\s+")};  
  my ($filename, $basename, $empty) = fileparse($output_file);  
  Phylo::make_directory($basename);  
  
  my $temp_file = Phylo::get_temp_file() . ".output";
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


#Makes a reference package
sub make_reference_packages {
  my $name = $_[0];
  my $alignment = $_[1];
  my $tree = $_[2];
  my $stats = $_[3];
  my $output_dir = $_[4];
        
  #Now build reference package
  ` python2.7 /projects/sate7/tools/bin/taxit create  --package-name $output_dir/$name.refpkg --locus $name --author "Nam Nguyen <namphuon\@cs.utexas.edu>" --package-version 1.0 --tree-stats $stats --tree-file $tree --aln-fasta $alignment --taxonomy $taxonomy_file --seq-info $mapping_file`;  
}

sub reverse_dna {
  my $sequence = $_[0];
  my %map = ('A', 'T', 'a', 't', 'C', 'G', 'c', 'g',
             'T', 'A', 't', 'a', 'G', 'C', 'g', 'c', '-', '-');
  $sequence = reverse($sequence);
  my @chars = split(//, $sequence);
  foreach my $idx (0..(scalar @chars-1)) {
    if (not defined $map{$chars[$idx]}) {
      $chars[$idx] = 'N';
    }
    else {
      $chars[$idx] = $map{$chars[$idx]};
    }
  }
  local $" = "";
  return "@chars";
}

#Returns name of a fragment
sub get_name_fragment {
  my $frag = $_[0];
  my $name = undef;
  
  my @results = split(/\s+/, $frag);  
  
  return $results[0];
}

#Returns the direction of a fragment
sub get_direction {
  my $frag = $_[0];  
  my $direction = undef;
  $frag =~ m/KEY=([^}]+)}/;
  if (not defined $1) {
    return undef;
  }
  my @results = split(",", $1);
  return $results[1];
}

#Returns the source name of a fragment
sub get_source {
  my $frag = $_[0];  
  my $name = undef;
  if ($frag =~ m/SOURCE_1="([^"]+)"/) {
    $name = $1;
  }  
  if (not defined $name) {
    $frag =~ m/(.*)_(\d+)/;
    $name = $1;
  }
  return $name;
}

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
            #$output_mapping{$fragment}->[$miss] = "UNCLASSIFIED";
            $output_mapping{$fragment}->[$miss] = "NA";
          }
          last;
        }
        $lca--;
      }
    }
  }
  
  Phylo::write_multiple_mapping(\%output_mapping, $output_file, "\t", "fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum");    
  return;
}


#This function takes as input the fragment file and a clade and returns all sequences belonging to 
#that clade
sub filter_fragments {
  my $frag_file = $_[0];
  my $level = $_[1];
  my $clade = $_[2];
  my $mapping_file = $_[3];
  
  my $validation_directory = "$base_dir/data/validation_sets/";
  if (defined $_[4]) {
    $validation_directory = $_[4];
  }
  
  
  my $mapping_ref = undef;
  if (defined $mapping_file and -e $mapping_file) {
    $mapping_ref = Phylo::read_mapping($mapping_file, "\t");
  }
    
  my %frags = %{Phylo::read_fasta_file($frag_file, 0)};
  open(INPUT, "<$validation_directory/$level/$clade/sequences");
  my $line = <INPUT>;
  close(INPUT);
  my %sequences = map {$_ => $_ } split(/\s+/, $line);
  my %subfrags = ();
  foreach my $frag (keys %frags) {
    my $frag_name = $frag;
    if (defined $mapping_ref) {
      $frag_name = $mapping_ref->{$frag};
    }
    my $source = Meta::get_source($frag_name);
    if (defined $sequences{$source}) {
      $subfrags{$frag} = $frags{$frag};
    }
  }
  return \%subfrags;
}

#This function takes as input the alignment file and a clade and returns an induced alignment (alignment that does not contain members from the clade)
sub filter_sequences {
  my $alignment_file = $_[0];
  my $level = $_[1];
  my $clade = $_[2];
  my $validation_directory = "$base_dir/data/validation_sets/";
  if (defined $_[3]) {
    $validation_directory = $_[3];
  }
    
  my %alignment = %{Phylo::read_fasta_file($alignment_file, 0)};  
  open(INPUT, "<$validation_directory/$level/$clade/sequences");
  my $line = <INPUT>;
  close(INPUT);
  my %sequences = map {$_ => $_ } split(/\s+/, $line);  
  foreach my $sequence (keys %sequences) {
    delete $alignment{$sequence};
  }
  return \%alignment;
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
    
    my ($h1, $c1, $t1) = Meta::read_stats_file("$file");  
    my @t = ($h1, \%hit);  
    %hit = %{Phylo::sum_hash(\@t)};
    @t = ($c1, \%classified);  
    %classified = %{Phylo::sum_hash(\@t)};
    @t = ($t1, \%total);  
    %total = %{Phylo::sum_hash(\@t)};      
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
    print OUTPUT sprintf("$level,hits,%f\n", $hit{$level});
    print OUTPUT sprintf("$level,classified,%f\n", $classified{$level});
    print OUTPUT sprintf("$level,total,%f\n", $total{$level});
  }
  close(OUTPUT);  
}

#This function reads a stats file and prints it to a file
sub print_stats_file {
  my $stats_file = $_[0];
  my $output_file = $_[1];
  my @levels = ("species", "genus", "family", "order", "class", "phylum");
  
  my %hit = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %total = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  my %classify = ("species", 0, "genus", 0, "family", 0, "order", 0, "class", 0, "phylum", 0);
  foreach my $level (@levels) {
    my $p = `grep "precision" $stats_file | grep "$level" | awk {'print \$3'}`;    
    my $c = `grep "classified_count" $stats_file | grep "$level\t" | awk {'print \$3'}`;    
    my $t = `grep "total" $stats_file | grep "$level" | awk {'print \$3'}`;
    
    my $pc = $p/100.0*$c;
    $hit{$level}+=$pc;
    $total{$level}+=$t;
    $classify{$level}+=$c;    
  }
  
  open(OUTPUT,">$output_file");
  print OUTPUT "Level\t%Correct\t%Incorrect\t%Unclassified\n";
  foreach my $level (@levels) {
    my $correct = $hit{$level}/$total{$level};
    my $incorrect = ($classify{$level}-$hit{$level})/$total{$level};
    my $unclassified = ($total{$level}-$classify{$level})/$total{$level};    
    print OUTPUT sprintf("$level\t%0.2f\t%0.2f\t%0.2f\n", $correct, $incorrect, $unclassified);
  }
  close(OUTPUT);
}

sub convert_motu_to_classification {  
  #Read names
  open(INPUT,"\$WORK/metagen/data/taxonomy/names.dmp");  
  my %name_map = ();
  while (my $line = <INPUT>) {
    if ($line !~ m/scientific name/) {
      next;
    }
    my @results = map {Phylo::trim($_)} split(/\|/, $line);
    $name_map{$results[0]} = $name_map{$results[1]};    
  }
  close(INPUT);

  #Build taxonomy from nodes
  open(INPUT,"\$WORK/metagen/data/taxonomy/nodes.dmp");  
  my %nodes_map = ();  
  while (my $line = <INPUT>) {
    my @results = map {Phylo::trim($_)} split(/\|/, $line);
    my %node = ('rank', $results[2], 'id', $results[0], 'parent_id', $results[1], 'children', {});
    $nodes_map{$results[0]} = \%node;    
  }  
  close(INPUT);

  #Build children map  
  foreach my $node (keys %nodes_map) {
    my $map = $nodes_map{$nodes_map{$node}->{'parent_id'}}->{'children'};
    $map->{$node}=$nodes_map{$nodes_map{$node}};
  }  


  my $in = $_[0];
  my $out = $_[1];
  my $dataset = $_[2];
  my %replacement_map = ("269483","482957");
  my ($taxon_ref, $level_ref, $key_ref) = Phylo::read_taxonomy_mapping("$taxonomy_file", "lower");
  
  open(INPUT, "<$in");  
  my $start = 0;
  my %abundances = ();
  foreach my $level (@levels) {
    $abundances{$level}={};
  }
  while (my $line = <INPUT>) {
    if ($line =~ m/taxaid\s+motus\.processing\./) {
      $start = 1;
      $line = <INPUT>;
      next;
    } elsif ($start == 0) {
      next;
    } 
    my @results = split(/\s+/,$line);
    if ($results[1] == 0) {
      next;
    }
    if (defined $replacement_map{$results[0]}) {
      $results[0] = $replacement_map{$results[0]};
    }
    if (not defined $taxon_ref->{$results[0]}) {
      if (not defined $nodes_map{$results[0]}) {      
        print "Can't find $results[0] in node\n";
        $results[0] = 1;        
      } else {
        while (not defined $taxon_ref->{$results[0]}) {
          if (not defined $nodes_map{$results[0]}->{'parent_id'}) {
            $results[0] = 1;
          } else {
            $results[0]=$nodes_map{$results[0]}->{'parent_id'};          
          }
        }
      }       
      #exit();
    }
    my @lineage = @{$taxon_ref->{$results[0]}};
    foreach my $level(@levels) {
      my $clade = $lineage[$key_ref->{$level}];
      if ($clade eq "") {
        $abundances{$level}->{"unclassified"}+=$results[1];
      } else {
        $abundances{$level}->{$taxon_ref->{$clade}->[$key_ref->{'tax_name'}]}+=$results[1];
      } 
    }    
  }
  
  open(OUTPUT, ">$out");
  foreach my $level (@levels) {
    foreach my $clade (keys %{$abundances{$level}}) {      
      print OUTPUT "motu\t$clade\t$level\tWGS\t$dataset\tall\t$abundances{$level}->{$clade}\t$abundances{$level}->{$clade}\n";
    }    
  }
  close(OUTPUT);
}


1;


