package Kaviar;
use strict;
use warnings;
use Tabix;
use Carp;
use Data::Dumper;
##-------------------------------------------------------------------------##
##  File:
##      @(#) Kaviar.pm
##  Authors:
##      Gustavo Glusman <Gustavo@systemsbiology.org>
##      Denise Mauldin <dmauldin@systemsbiology.org>
##  Description:
##   Kaviar (.Known VARiants.) is a compilation of human SNVs collected 
##     from many and diverse sources, stressing accessibility and ease
##     of use. Kaviar answers a very specific question: What variants 
##     have been reported already for a given specific genomic location? 
##     For each SNV in a query, Kaviar reports the known variants and 
##     their source (e.g. in which population or individual genomes it
##     was observed), or the fact that no variants are known.
##     Where available, dbSNP identifiers are displayed. 
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2011,2012 Developed by
#* Gustavo Glusman
#* Denise Mauldin
#*
#* This work is licensed under the MIT License.  To view a copy of this
#* license, visit http://www.opensource.org/licenses/mit-license.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
#  ChangeLog:
#
#    $Log$
#
#
###############################################################################
#
=head1 NAME

  Kaviar.pm - Interface to the Known VARiants database.

=head1 DESCRIPTION

  The Kaviar module permits access to the Kaviar database, which contains 
  single nucleotide variation information from a variety of sources and
  populations.  This database is maintained and updated by the Institute
  for Systems Biology and can be downloaded from :
  http://db.systemsbiology.net/kaviar

  This module is used for local access to the Kaviar database.  The 
  KaviarWebService.pl in the bin directory of the Kaviar source installation
  can be used to access the ISB Kaviar database.

  The reference genome is always coded as heterozygous because there is only
  one allele in the reference genome.


=head1 SYNOPSIS

  Create a Kaviar object, then use the appropriate query function to query the
  databases that you are interested in.  Kaviar includes two databases, a SNP
  database and a Ranged database.  Functions beginning with knownVariant query
  the SNP database and functions starting with rangeVariant query the ranged database.
  Use the appropriate function depending on the output you'd like to receive and
  calculate statistics on the result using digest (hash) or decoct (array) 
  depending on which output type you requested.

=head1 BASIC USAGE

  There are a couple of very basic queries for the Kaviar module:

  1) knownVariant(<chromosome>, <position>)

     This takes a chromosome and a position and returns a hash suitable to 
     use in digest for information.  Queries the SNV database only.

  2) variantsInRange(<chromosome, <start>, <stop>)

     This takes in a chromsosome, a start, and a stop and returns an array of arrays

=head1 STANDARD OPTIONS FOR KAVIAR OBJECT

  start_match, start_overlap, end_match, and end_overlap are described in the query window section.

  -return_reference (DEFAULT: OFF)
      When a variant includes reference base pairs at the beginning (common in insertions), control
      whether Kaviar returns that reference base pair information.

  -exclude_insertions (DEFAULT: OFF)
      Skip all insertions.

  -expand_insertions (DEFAULT: OFF)
      Insertions are displayed as <base pair>+ in the return variation.  Ex: For a reference position,
      AT, a variant GCG is encountered.  If expand_insertions is OFF (default), this would return as
      position 1 => A vs G, position 2 => T vs C+.  If expand_insertions is ON, this would return as
      position 1 => A vs G, position 2 => T vs C, position 3 => G_INSERTION.  The G_INSERTION will not
      have reference information available (unless it happens to overlap with another variation that 
      includes reference information) and therefore statistics calculated on the G_INSERTION will be
      unreliable.

  -maxlimit (DEFAULT: 50000)
     The maximum number of queries to handle at once.

  -maxlength (DEFAULT: 100)
     The maximum length of a variant to return in allVariantPositions call.

  -warn (DEFAULT: OFF)
     Return informative warning messages about filtering on variant or person.

=head1 QUERY WINDOW DESIGNATIONS

  There are 6 different ways a variation can interact with a query window.
    - start_match - the end of the variation matches the start of the query window
    - end_match - the start of the variation matches the end of the query window
    - start_overlap - the variation overlaps the start of the query window
    - end_overlap - the variation overlaps the end of the query window
    - within - the variation is entirely within the query window
    - inside - the query window is entirely inside the variation

     start_match                                        end match
    |--------------|                               |-------------------|
  QUERY WINDOW     |-------------------------------|
                     |---------| range within  ^ SNP within
  start_overlap |---------|                   |----------|  end_overlap
           |--------------------------------------------| inside

  The Kaviar object defaults to :
	start_match OFF, end_match ON
	start_overlap OFF, end_overlap ON

  There is no way to turn inside and within OFF.

  The defaults mean that you will receive results for basepairs
  within the query window non-inclusive start and inclusive end.

  Suppose the database has three entries:
    1: 1295-1300
    2: 1300-1305
    3: 1304-1308

  If you DEFAULT query 1300 to 1303, then Kaviar will return variation #2 ONLY
    in the 1303 to 1303 window.
  If you turn start_match on and query 1300 to 1303, then Kaviar will return 
    variations #1 and #2, but only those basepairs in the 1300 to 1303 window.
  If you turn start_overlap on and query 1300 to 1303, then Kaviar will return
    variations #1 and #2, in the window 1295 to 1303.
  If you DEFAULT query 1301 to 1305, then Kaviar will return variations #2 and
    #3, but only those basepairs in the 1301 to 1305 window.
  If you turn end_match OFF and query 1301 to 1305, then Kaviar will return 
    variations #2 and #3, but only those basepairs in the 1301 to 1304 window.
  If you turn end_overlap on and query 1301 to 1305, then Kaviar will return
    variations #2 and #3, in the window 1301 to 1308.

=head1 KAVIAR VARIANT QUERIES

  There are three different ways to query the Kaviar database because Kaviar
    includes two different databases - a SNP database and a Ranged database.

    To query ONLY the SNP database:
      knownVariant  		- returns a hash of one hit
      knownVariants 		- returns an array of arrays of multiple hits
      knownVariantPositions 	- returns a hash of multiple hits separated by position

    To query ONLY the Ranged database:
      rangeVariants 		- returns an array of arrays of multiple hits
      rangeVariantPositions	- returns a hash of multiple hits separated by position

    To query BOTH databases:
      allVariants		- returns an array of arrays of multiple hits
      allVariantPositions	- returns a hash of multiple hits separated by position

  All of these methods require that you provide chromosome, start, and end parameters.

=head1 STANDARD OPTIONS FOR KAVIAR VARIANT QUERIES

  -people
    Provide an array reference of people Names or IDs (from $kav->listInfo)

  -variant
  

=head1 EXAMPLES

  Basic query of SNP database:
    use Kaviar;
    my $freeze = 'hg19';
    my $kav = new Kaviar(freeze => $freeze);
    die "Kaviar object not created correctly" unless (ref($kav) eq 'Kaviar');
    my $varRef = $kav->knownVariant(<chromosome>, <position>);

  Query of SNP database with return_reference and start_match:
    use Kaviar;
    my $freeze = 'hg19';
    my $kav = new Kaviar(freeze => $freeze, basedir => <basedir>, "warn" => 1, return_reference => 1);
    die "Kaviar object not created correctly" unless (ref($kav) eq 'Kaviar');
    my $varRef = $kav->knownVariant(chromosome => <chromosome>, start => <start>, end => <end>);

  Query of SNP database with warn and end_overlap:
    use Kaviar;
    my $freeze = 'hg19';
    my $kav = new Kaviar(freeze => $freeze, basedir => <basedir>, "warn" => 1, end_overlap => 1);
    die "Kaviar object not created correctly" unless (ref($kav) eq 'Kaviar');
    my $varRef = $kav->knownVariant(chromosome => <chromosome>, start => <start>, end => <end>);

  Query of Ranged database with expand_insertions:
    use Kaviar;
    my $freeze = 'hg19';
    my $kav = new Kaviar(freeze => $freeze, expand_insertions => 1);
    die "Kaviar object not created correctly" unless (ref($kav) eq 'Kaviar');
    my $varRef = $kav->rangeVariants(chromosome => <chromosome>, start => <start>, end => <end>);

  Query of all databases with variant specified:
    use Kaviar;
    my $freeze = 'hg19';
    my $kav = new Kaviar(freeze => $freeze, "warn" => 1);
    die "Kaviar object not created correctly" unless (ref($kav) eq 'Kaviar');
    my $varRef = $kav->allVariants(chromosome => <chromosome>, start => <start>, end => <end>, 
    				      variant => "A");



=head1 COPYRIGHT

  Copyright 2011,2012 Institute for Systems Biology

=head1 AUTHOR

  Gustavo Glusman <Gustavo@systemsbiology.org>
  Denise Mauldin <dmauldin@systemsbiology.org>

=cut

=head1 FUNCTIONS
=cut 

# Set $kaviar_root to the parent directory of this file

my $kaviar_root;
BEGIN {
  ($kaviar_root) = __FILE__ =~ m|^(\S+)/\S+?/\S+?$|;
}

# This belongs in $self, configurable in Kaviar.config
our @stat_display_order = qw (
  alleles
  refcall
  counts
  frequency
  total
  major
  minor
  maf
);

sub new {
	my $package = shift;
	
	my $obj = {};
	bless $obj, $package;
	my $status = $obj->initialize(@_);
	unless ((defined $status) && ($status eq "fail")) {
		return $obj;
	}
}

# Changes to initialize() October 2014
# - Two filesystem locations maintained: $basedir and $datadir. $datadir contains
#    subdirs for all values of $frz, each containing its Kaviar db files.
#    $basedir contains all other elements of Kaviar, including the Names,
#    Populations, and Projects files, and the software in $basedir/bin
# - No arguments/params required anymore.
# - New third argument and/or param, datadir. If not specified, basedir will be used.
# - Kaviar.config is expected to be in $basedir, not $datadir/$frz
# - Kaviar.config may contain any element allowable in params to $kav->initialize(),
#    including datadir. Values in params take priority.
# - Kaviar version is specified only in Kaviar.config, is not hardcoded in Kaviar.pm
sub initialize {
	my $self = shift;

	# set defaults
	$self->{'basedir'} = $kaviar_root;
	$self->{'warn'} = 0;
	$self->{'maxlimit'} = 50000;
	$self->{'maxlength'} = 100;
	$self->{'start_overlap'} = 0;
	$self->{'end_overlap'} = 1;
	$self->{'start_match'} = 0;
	$self->{'end_match'} = 1;
	$self->{'return_reference'} = 0; # return positions that only have reference basepairs
	$self->{'check_reference'} = 0; # print out debug parsing that looks for variants that have reference sequence
	$self->{'exclude_insertions'} = 0;
	$self->{'expand_insertions'} = 0;
	$self->{'tabix'} = "/tools/bin/tabix";
	$self->{'ref_tabix_dir'} = "$kaviar_root/tabixedRef";
	$self->{'freezes'} = "hg18,hg19,hg38";

	# Create a hash of allowed params that can be specified either
	# when calling initialize() or in Kaviar.config
	# Note: 'freeze' is used in params and Kaviar.config,
	# but 'frz' in in $frz and  $self->{'frz'}
	my %allowed_keys = map { $_ => 1} qw(warn freezes frz freeze basedir datadir start_match end_match start_overlap end_overlap return_reference check_reference exclude_insertions expand_insertions binsize sectionsize tabix ref_tabix_dir maxlength);

	my %params;
	# now store in %params any params passed to initialize()
	# Before Oct. 2014, params ($freeze, $basedir) were required .
	# Starting Oct. 2014, these may be omitted if specified in Kaviar.config.
	if (@_) {
	  # Are params being passed as a hash, or as ($freeze, $basedir)?
	  # If first item is not in %allowed_keys, we assume the latter.
	  if (defined $allowed_keys{$_[0]}) {
	    (%params) = @_;
	  } else {
	    my ($freeze, $basedir, $datadir) = @_;
	    $params{'freeze'} = $freeze if $freeze;
	    $params{'basedir'} = $basedir if $basedir;
	    $params{'datadir'} = $datadir if $datadir;
	  }
#	  die "Must provide a parameter hash to Kaviar.  my \$kav = new Kaviar(freeze => \"hg19\", basedir => \".\")" unless ref(\%params) eq 'HASH';
#	  die "Must provide a parameter hash to Kaviar.  my \$kav = new Kaviar(freeze => \"hg19\", basedir => \".\")" unless scalar (keys %params) >= 1;
#	  die "Must provide viar.  my \$kav = new Kaviar(freeze => \"hg19\", basedir => \".\")" unless scalar (keys %params) >= 1;
	}

	#print "parameters in initialize ".Dumper \%params if $self->{'warn'};
	
	$self->{'warn'} = 1 if ($params{'warn'} && $params{'warn'} == 1); # set warn early so that it'll print bad keys

	my $basedir = $params{'basedir'} || $self->{'basedir'};
	# Read values in optional config files.  Any file in cwd
	# supercedes any in Kaviar basedir.
	# In general, override defaults
	for my $config_loc ("$basedir", ".") {
	  if (open CONFIG, "$config_loc/Kaviar.config") {
	    while (<CONFIG>) {
	      chomp;
	      next if /^#/;
	      next if /^$/;
	      my($field, @values) = split /\t/;
	      next unless scalar @values;
	      # allow setting of anything that can be specified in params
	      # if not already specified in params
	      if ($field eq 'freeze') {
		$self->{'frz'} = $values[0];
	      } elsif ($field eq 'basedir') {
		print "basedir is not allowed to be set in config file.\n" if $self->{'warn'};
	      } elsif ($allowed_keys{$field}) {
		$self->{$field} = $values[0];
	      } elsif ($field eq 'version') {
		$self->{'version'} = $values[0];
	      } elsif ($field eq 'decoded') {
		$self->{'decoded'} = \@values;
	      } else {
		print "$field is not valid key in config file.\n" if $self->{'warn'};
	      }
	    }
	    close CONFIG;
	  }
	}

	# For any params passed to initialize(), store in %self,
	# overriding defaults and values in Kaviar.config
	foreach my $key (keys %params) {
	  if ($key eq 'freeze') {
	    $self->{'frz'} = $params{'freeze'};
	  } elsif ($allowed_keys{$key}){
	    $self->{$key} = $params{$key};
	  } else {
	    print "$key is not allowed to be set.\n" if $self->{'warn'};
	  }
	}
	$self->{'datadir'} = $basedir unless $self->{'datadir'};

	# Set these 3 variables because they will be used frequently
        $basedir = $self->{'basedir'};
        my $datadir = $self->{'datadir'};
	my $frz = $self->{'frz'};

	unless ($frz) {
	  print "Initializing Kaviar: Freeze must be provided in params or Kaviar.config\n";
	  return "fail";
	}

	unless ( -e "$datadir/$frz") {
		print "ERROR: can't find database for $frz in $datadir.\n";
		return "fail";
		# If any part of @INC is a Kaviar database, use it instead of 
		# the Kaviar provided by user. Use with caution.
		foreach my $line (@INC) {
			if ($line =~ /Kaviar/ && -e "$line/$frz") {
				$datadir = $basedir = $self->{'datadir'} = $self->{'basedir'} = $line;
				last;
			}
		}
	}
    
	my $freeze_dir = $self->{'dir'} = "$datadir/$frz";
	$self->{'frz'} = $frz;
	$self->{'dbs'}{'file'} ||= "$freeze_dir/Kaviar.db.s";
	$self->{'dbr'}{'file'} ||= "$freeze_dir/Kaviar.db.r";
	$self->{'dbg'}{'file'} ||= "$freeze_dir/Kaviar.db.g";
	#$self->{'dbn'}{'file'} ||= "$freeze_dir/Kaviar.db.n";  ### FOR BUILD 160113, disable position-specific AN
	$self->{'rsid_index'}{'file'} ||= "$freeze_dir/dbSNP_index";
	$self->{'binsize'} ||= 250;
	$self->{'sectionsize'} ||= 1e6;
	$self->{'decoded'} ||= ["null", 0..7, "tab", "newline", "verttab", 8, "carriagereturn", 9..42, "zero", 43..249];
	$self->{'alphabet'} ||= [qw/A C G T/];
	
	

	my @allIdentifiers;
	my(%name, %code, %sites, %complete, %genome_coverage);
	open CODES, "$freeze_dir/Kaviar.identifiers" || return "fail";
	while (<CODES>) {
		chomp;
		my($code, $name, $genome_coverage, $sites_s, $sites_r, $sites_g) = split /\t/;
		$sites_s ||= 0; $sites_r ||= 0; $sites_g ||= 0;
		$name{$code} = $name;
		$code{$name} = $code;
		$genome_coverage{$code} = $genome_coverage;
		$sites{'s'}{$code} = $sites_s; $sites{'r'}{$code} = $sites_r; $sites{'g'}{$code} = $sites_g;
		$complete{$code} = 1 if $sites{'s'}>=2e6;
		push @allIdentifiers, $code;
		$self->{'referenceCode'} = $code if $name =~ /^ref($|erence)/i;
		$self->{'dbsnpCode'} = $code if $name =~ /dbsnp/i;
		$self->{'ISBCode'} = $code if $name =~ /isb/i;
		$self->{'InovaCode'} = $code if $name =~ /inova/i;
		$self->{'exacCode'} = $code if $name =~ /63000exomes/i;
		$self->{'1000GenomesCode'} = $code if $name =~ /phase3/i;
	}
	close CODES;
	$self->{'name2code'} = \%code;
	$self->{'code2name'} = \%name;
	$self->{'genomeCoverage'} = \%genome_coverage;
	$self->{'allIdentifiers'} = \@allIdentifiers;
	my @allIdentifiersButRef = @allIdentifiers;
	shift @allIdentifiersButRef; #remove reference
	$self->{'allIdentifiersButRef'} = \@allIdentifiersButRef;
	$self->{'siteCounts'} = $sites{'s'};
	my $completeIndependent;

	$self->readProjects("$basedir/Projects") || return "fail";

	my %info;
	foreach my $file (qw/Names Populations/) {
		open NAMES, "$basedir/$file" || return "fail";
		$_ = <NAMES>;
		chomp;
		my @fields = split /\t/;
		while (<NAMES>) {
			next if /^#/;
			chomp;
			my(@values) = split /\t/;
			next unless defined $code{$values[0]};
			my $code = $code{$values[0]};
			
			foreach my $i (0..$#fields) {
				$info{$code}{$fields[$i]} = $values[$i];
			}
			$info{$code}{'type'} = lc($file);
			if ($file eq 'Populations') {
				delete $complete{$code};
			} elsif ($complete{$code} && !$info{$code}{'parents'}) {
				$completeIndependent++;
			} else {
				delete $complete{$code};
			}
			
			if ($info{$code}{'id'} =~ /^phase/) {
				$self->{'special'}{'1000Genomes-AN'} = $info{$code}{'AN'};
				# 1000Genomes used to be low power, but no more
				#$self->{'special'}{'1000Genomes-power'} = .7;
			} elsif ($info{$code}{'id'} =~ /^ESP/) {
				$self->{'special'}{'ESP-AN'} = $info{$code}{'AN'};
			}
			my $project = $info{$code}{'project'};
			$info{$code}{'Genome coverage'} =
			   $self->{'projectinfo'}{$project}{'Genome coverage'};
			$info{$code}{'Genotypes'} =
			   $self->{'projectinfo'}{$project}{'Genotypes'};
		}
		close NAMES;
	}

	my $wesAN=0;
	my $wgsAN=0;
	for my $project (keys %{$self->{'projectinfo'}}) {
	  my $AN = $self->{'projectinfo'}{$project}{'AN'};
	  my $coverage = $self->{'projectinfo'}{$project}{'Genome coverage'};
	  if ($coverage =~ /^WES$/i) {
	    $wesAN+= $AN;
	  } elsif ($coverage =~ /^WGS$/i) {
	    $wgsAN += $AN;
	  }
	}
	$self->{'wesAN'} = $wesAN;  # whole exomes, including populations
	$self->{'wgsAN'} = $wgsAN;  # whole genomes, including populations
	
	# Individual whole genome sequences
	$self->{'completeGenomes'} = \%complete;
	# Those independent from one another ("founders")
	$self->{'completeIndependent'} = $completeIndependent;
	
	$self->{'info'} = \%info;
	$self->{'loaded'} = {};
	$self->{'offset'} = {};
	$self->{'lastoffset'} = {};
	
	#foreach my $db (qw/dbs dbr dbg dbn rsid_index/) {  ### FOR BUILD 160113, disable position-specific AN
	foreach my $db (qw/dbs dbr dbg rsid_index/) {
		my $dbf = $self->{$db}{'file'};
		my $tbx = new Tabix ('-data', "$dbf.gz", '-index', "$dbf.gz.tbi");
		die "Couldn't create Tabix object\n" unless ref $tbx eq 'Tabix';
		$self->{$db}{'tbx'} = $tbx;
		#read chromosomes in index to avoid tabix crash
		my %chroms;
		unless (($db eq 'dbn') || ($db eq 'rsid_index')) {
		  open CHRMS, "$dbf.chroms";
		  while (<CHRMS>) {
		    chomp;
		    $chroms{$_} = 1 if $_;
		  }
		  close CHRMS;
		  $self->{$db}{'tbxchroms'} = \%chroms;
		}
	}
	
	my $exHet_file = "$basedir/exZ3_1000Genomes.v1.bed.gz";
	if (-e $exHet_file) {
	  $self->{'special'}{'exHet-ranges-file'} = $exHet_file;
	} else {
	  warn "Can't find exHet boundary file $exHet_file. No variants will be annotated as in an exHet region.";
	}
	# Takes about one second; speeds digest by 0.01s (bobama)
	($self->{'special'}{'exHet-binsize'}, $self->{'special'}{'exHet-regions'}) =
	  $self->loadAnnotationSegments($exHet_file);
	
	my $exome_file;
	if ($frz eq 'hg19') {
	  $exome_file = "$basedir/exome_calling_regions.v1.interval_list.gz";
	} else {
	  $exome_file = "$basedir/exome_calling_regions.v1.interval_list_$frz.gz";
	}
	if (-e $exome_file) {
	  $self->{'special'}{'exon-ranges-file'} = $exome_file;
	} else {
	  warn "Can't find exome boundary file $exome_file. All variants in whole genome data will be considered as though outside the exome.";
	}
	
}

#########################################################################################################
#########################################################################################################
# ACCESSOR METHODS
#########################################################################################################
#########################################################################################################


sub version {
	my($self) = @_;
	return ("$self->{'version'}");
}

sub AN {
  my ($self, $code) = @_;
  my $an = $self->{'info'}->{$code}{'AN'};
  $an = 2 unless defined $an;
  return $an;
}

sub posAN {
  my ($self, $chr, $pos) = @_;
  my $tbx = $self->{'dbn'}{'tbx'};
  my $res = $tbx->query($chr, $pos, $pos+1);
  return '' unless $res;
  my $line = $tbx->read($res);
  return '' unless $line;
  my @fields = split "\t", $line;
  my $ANcode = $fields[2];
  return $ANcode;
}

sub warn {
	my($self, $value) = @_;
	if ($value) {
		die "warn only accepts YES|1 or NO|0" if ($value ne 'YES' && $value ne 'NO' && $value ne '1' && $value ne '0');
		$self->{'warn'} = 0 if ($value eq 'NO' || $value eq 0);
		$self->{'warn'} = 1 if ($value eq 'YES' || $value eq 1);
	}
	return $self->{'warn'};
}

sub maxlimit{
	my($self, $value) = @_;
	if ($value) {
		$self->{'maxlimit'} = $value;
	}
	return $self->{'maxlimit'};
}

sub maxlength{
	my($self, $value) = @_;
	if ($value) {
		$self->{'maxlength'} = $value;
	}
	return $self->{'maxlength'};
}

=over 12

=item infoAbout

	$kav->infoAbout($code)

	For a code that describes a genome, return all of the information about that genome

=back

=cut
sub infoAbout {
	my($self, $code) = @_;
	
	return $self->{'info'}->{$code};
}

=over 12

=item datasetType

	$kav->datasetType($code)

	For a code that describes a genome, return the dataset type.  
	  Currently this is either 'names' or 'populations'.

=back

=cut
sub datasetType {
	my($self, $code) = @_;

	return $self->{'info'}->{$code}->{'type'};
}

=over 12

=item listInfo
  
	$kav->listInfo();

	Returns an ARRAY of the ID and Name of every person and population included in Kaviar

=cut
sub listInfo {
	my ($self) = @_;
	my @ids = $self->allIdentifiers();
	my @info;
	foreach my $id (@ids) {
		my $string = "ID: ".$self->infoAbout($id)->{'id'}.";";
		my @people = $self->whois($id);
		$string .= " Name: ".$people[0].";" if $self->whois($id);
		push(@info, $string);
	}
	return @info;
}

sub loadSection {
	die "Called deprecated Kaviar->loadSection() - please update your scripts.\n";
}

=back

=over 12

=item separateCodes

	$kav->separateCodes($codes, <detailed>)

	Takes in a string of codes -- 1-5 byte codes for sources, each followed immediately by
	  a 2-byte AC code if and only if that source is a population.
	Returns an array of codes if detailed is not set: $kav->separateCodes($codes)
	Returns a hash (codes hashed to AC) if detailed is set: $kav->separateCodes($codes, 1)
	This should only be used to provide entry to infoAbout() or datasetType().



=cut
sub separateCodes_old {
  my($self, $codes, $detailed) = @_;
  my(%splitcodes, @codelist);
  my $info = $self->{'info'};
  my $extralength = 3;
  while ($codes) {
    my $code;
    # Codes for data sources are 1-5 bytes. Check next 1 byte, 2 bytes, ... 5 bytes
    #  until we find a code, and peel it off.
    foreach my $l (1..5) {
      my $testcode = substr($codes, 0, $l);
      if (defined $info->{$testcode}) {
	$codes = substr($codes, $l);
	$code = $testcode;
	last;
      }
    }
    if (defined $code) {
      # If this is a population (AN>1), a 3-byte AC code should be next in the string.
      # Peel it off and store in hash.
      # (has to be > 1 to skip the reference, which causes substr errors)
      if ($info->{$code} && $info->{$code}{'AN'} && $info->{$code}{'AN'}>1) {
	if ($detailed) {
	  $splitcodes{$code} = $self->decodeValue3(substr($codes, 0, $extralength));
	} else {
	  push @codelist, $code;
	}
	next unless $codes; #removes the substr outside of range warning if codes is undef
	$codes = substr($codes, $extralength);
	# if not a population, there is no AC code to peel of.
      } else {
	if ($detailed) {
	  $splitcodes{$code} = "";
	} else {
	  push @codelist, $code;
	}
      }
    } else {
      print "code $code is not defined. Cannot interpret \"$codes\" (".length($codes).")\n" if $self->{'warn'};
      carp "Cannot interpret \"$codes\" (".length($codes).")\n";
      $codes='';
    }
  }
  return $detailed ? \%splitcodes : \@codelist;
}

=back

=item separateCodes2

        separateCodes, modified for genotype-centric kaviar to include encoded AN.

	$kav->separateCodes2($codes, <detailed>)

	Takes in a string of codes -- 1-5 byte codes for sources, each followed immediately by
	  a 2-byte AC code and a 2-byte AN code.
	Returns an array of codes if detailed is not set: $kav->separateCodes2($codes)
	Returns a hash if detailed is set: $kav->separateCodes($codes, 1)
	This is used by whois_detailed() and digest().
	Hash is deeper than that returned by original separateCodes.



=cut
sub separateCodes {
  my($self, $codes, $detailed) = @_;
  my(%splitcodes, @codelist);
  my $info = $self->{'info'};
  my $codelength = 3;
  $detailed = 0 unless defined $detailed;
  #print "Hello! In separateCodes Detailed is |$detailed|. \$codes is |$codes| (length " . length($codes) . ")\n";
  while ($codes) {
    #print " Processing next coded triplet in |$codes| (length " . length($codes) .")\n";
    my $code;
    my $data_source;
    # Codes for data sources are 1-5 bytes. Check next 1 byte, 2 bytes, ... 5 bytes
    #  until we find a code, and peel it off.
    foreach my $l (1..5) {
      my $testcode = substr($codes, 0, $l);
      if (defined $info->{$testcode}) {
	$codes = substr($codes, $l);
	$code = $testcode;
	$data_source = $self->{'info'}->{$code}->{'id'};
	#print "  data source code |$code| is for $data_source\n";
	last;
      }
    }
    if (defined $code) {  # code for data source is defined
      # If this is non-reference, a 3-byte AC code should be next in the string.
      # Peel it off and store in hash.
      if ($data_source !~ /reference/i) {
	if ($detailed) {
	  #print "     Now decoding codes |$codes| (length " . length($codes) .")\n";
	  $splitcodes{$code}->{'AC'} = $self->decodeValue3(substr($codes, 0, $codelength));
	  #print "      AC for $data_source is $splitcodes{$code}->{'AC'}\n";
	  #$splitcodes{$code}->{'AN'} = $self->decodeValue3(substr($codes, $codelength, $codelength));
	} else {
	  push @codelist, $code;
	}
	next unless $codes; #removes the substr outside of range warning if codes is undef
	#$codes = substr($codes, $codelength*2);  # used when we were storing AN here
	$codes = substr($codes, $codelength);
      } else {
	if ($detailed) {
	  $splitcodes{$code} = "";
	} else {
	  push @codelist, $code;
	}
      }
    } else {
      print "\$code is not defined. Cannot interpret codestring \"$codes\" (length ".length($codes).")\n" if $self->{'warn'};
      carp "Cannot interpret codestring \"$codes\" (length ".length($codes).")\n";
      $codes='';
    }
  }
  return $detailed ? \%splitcodes : \@codelist;
}

=back
=cut
####################################################################################################
####################################################################################################
# METHODS TO DECODE KAVIAR VALUES
####################################################################################################
####################################################################################################


=over 12

=item whois

	 $kav->whois($codes)

	For a string of codes that describe genomes and allele counts (typically for
          a specific variant at a specific position), return an array of descriptive 
	  names for those genomes.  The hash returned from a Kaviar call contains
	  two lists of codes in an array in the 'sources' key.  The first set of 
	  codes is the heterozygous codes.  The second set of codes is the homozygous
	  codes.  { 'sources' => ['hetcodes','homcodes'] }
	Parameters: A string of codes. 
	Returns: an array of names or ids. Does not return allele counts.

=back

=cut
sub whois {
	my($self, $codes) = @_;
	my @res;
	my @codes = @{$self->separateCodes($codes)};
	foreach my $code (@codes) {
		my $info = $self->{'info'}->{$code};
		my $name = $info->{'name'};
		my $id = $info->{'id'};
		#$name .= " ($id)" if defined $name && $name =~ /^anonymous/i;
		$name = '' if defined $name && $name =~ /^anonymous/i;
		push @res, $name || $id;
	}
	return @res;
}

# whois_detailed:
# Same as whois(), but appends allele count (AC) to each name.
sub whois_detailed_old {
	my($self, $codes) = @_;
	my @res;
	my %codes = %{$self->separateCodes($codes, 1)};
	foreach my $code (sort keys %codes) {
		my $info = $self->{'info'}->{$code};
		my $name = $info->{'name'};
		my $id = $info->{'id'};
		#$name .= " ($id)" if defined $name && $name =~ /^anonymous/i;
		$name = '' if defined $name && $name =~ /^anonymous/i;
		$name ||= $id;
		$name .= "=$codes{$code}" if $codes{$code};
		push @res, $name;
	}
	return @res;
}


# whois_detailed2:
# Same as whois_detailed(), but appends both AC and AN.
sub whois_detailed {
  my($self, $codes) = @_;
  my @res;
  my %codes = %{$self->separateCodes($codes, 1)};
  foreach my $code (sort keys %codes) {
    my $info = $self->{'info'}->{$code};
    my $name = $info->{'name'};
    my $id = $info->{'id'};
    my $AC = $codes{$code}->{'AC'};
    my $AN = $codes{$code}->{'AN'};
    $name = '' if defined $name && $name =~ /^anonymous/i;
    $name ||= $id;
    #$name .= "=$AC/$AN" if $AC && $AN;
    $name .= "=$AC" if $AC;
    push @res, $name;
  }
  return @res;
}


#########################################################################################################
#########################################################################################################
# SNP ONLY PARSING
#########################################################################################################
#########################################################################################################


=over 12

=item knownVariant

	knownVariant is an overloaded method that can take either a simple list of chromosome 
	and position or a hash reference to specify a more detailed query

	$kav->knownVariant($chrom, $position);
	$kav->knownVariant(chromosome => $chrom, position => $pos, start => $start, 
				 [var => $var, people => <array of people names or IDs>]);

	Takes a chromosome, a start position, an optional variant base pair (eg: ACGT), and
	  and optional array of people's names or IDs.
	Searches the SINGLE NUCLEOTIDE VARIANT DATABASE ONLY.  Returns ONLY ONE POSITION.
	If a variant is specified but not observed, the result is empty.
	Otherwise, it returns all variants at the specified position.

	Examples:

	1) chromosome and position (base 0)

		$kav->knownVariant("chr1", 10326)

		returns:
		{
			'sources' => {
				'C' => [<codes for C>],
				'T' => [<codes for T>],
			},
			'rsids' => ['rs112750067'],
			'chromosome' => "chr1",
			'start' => 10326,
			'end' => 10326,
		}

	2) chromosome, position (base 0), and variant

		$kav->knownVariant("chr1",10326,"C") returns: same as above
		$kav->knownVariant("chr1",10326,"G") returns: nothing

=back

=cut
sub knownVariant {
	my $self = shift;
	my %params;
	if (scalar @_ <= 2) {
		my ($chrom, $pos) = @_;
		$params{'chromosome'} = $chrom;
		$params{'position'} = $pos;
		die "Missing position in knownVariant call" if (!$pos);
	} else {
	  (%params) = @_;
	}
	die "Missing chromosome in knownVariant call" unless $params{'chromosome'};
	die "Missing start in knownVariant call" unless $params{'start'} || $params{'position'};
	print STDERR "\$kav->knownVariant does not support an end parameter.  Returns only the first variant at a position.\n" if $params{'end'} && $self->{'warn'};
	# end_overlap and start_overlap don't mean anything for a SNV query

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'} || $params{'position'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		$start_basezero = $start_basezero-1;
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
	}
	my $start_baseone = $start_basezero +1;
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);
	my $reqvar = $params{'variant'};
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my $chri = $self->cleanupChrom($chrom);
	my $db = 'dbs';
	my $decoded = $self->{'decoded'};
	return unless $self->{$db}{'tbxchroms'}{$chri};
	my $tbx = $self->{$db}{'tbx'};
	# even though tabix on the command line lets you call basezero-basezero, 
	# the tabix command line includes both sides - a position at bazezero will be returned and aposition at baseone
	# the tabix perl module includes only the right side - a position at basezero will NOT be returned
	# want to include and 'end' bp though because otherwise it will search the entire file
	print "KnownVariant ERROR: Start and end base pairs are in error.  start: $start_basezero end: $start_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $start_baseone) || ($start_basezero > $start_baseone) ) && $self->{'warn'};
	my $res = $tbx->query($chri, $start_basezero, $start_baseone);
	my %res;
	while (my $line = $tbx->read($res)) {
		my($c, $s, $rsids, @vars) = split /\t/, $line;
		if ($s==$start_baseone) {
			my %rsids;
			while ($rsids) {
				my $rsid = 
				$decoded->[ord(substr($rsids,0,1))]*15625000 +
				$decoded->[ord(substr($rsids,1,1))]*62500 +
				$decoded->[ord(substr($rsids,2,1))]*250 +
				$decoded->[ord(substr($rsids,3,1))];
				$rsids{"rs$rsid"} = 1;
				$rsids = substr($rsids, 4);
			}
			$res{'chromosome'} = $chrom;
			$res{'start'} = $start_basezero;
			$res{'end'} = $start_basezero;
			my @rsids = keys %rsids;
			$res{'rsids'} = \@rsids if scalar keys %rsids > 0;
			while (@vars) {
				my $var = shift @vars;
				if ($var eq '') { $var = '-'; }
				my $sources = shift @vars;
				# if people is defined, cull sources for those codes
				if (scalar @people > 0) {
					$sources = $self->cullCodes($sources, $people_codes_ref);
					next unless ($sources);
				}

				$res{'sources'}{$var}= $sources;
			}
			if ($reqvar && !$res{'sources'}{$reqvar}) {
				print "Did not find requested variant ($reqvar) in found variants (". join(", ", keys %{$res{'sources'}})."). Position ".($start_basezero)." filtered.\n" if $self->{'warn'};
				return;
			} else {
                return if (scalar @people > 0 && !$res{'sources'});
				return %res;
			}
		}
	}
	return;
}
=over 12

=item knownVariants


	$kav->knownVariants(chromosome => $chrom, position => $pos, start => $start, 
				  [end => $end, var => $var, people => \@people])

	Takes a chromosome, a start position, an optional end position, an optional variant 
	  base pair (eg: ACGT), and an optional array of people names or IDs.
	Searches the SINGLE NUCLEOTIDE VARIANT DATABASE ONLY.  
	If a variant is specified but not observed, the result is empty.
	Otherwise, it returns an array of all variants in the specified query window

	Examples:

	1) chromosome and position (base 0)

		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983)

		returns:
		[
		          [
		            13979,
		            13979,
		            {
		              'sources' => {
		                             'T' => [ <het codes>,<hom codes> ],
		                             'C' => [ <het codes>,<hom codes> ]
                		           },
		              'rsids' => [ 'rs201837439' ],
		              'chromosome' => 'chr5',
		              'end' => 13979,
		              'start' => 13979
		            }
		          ],
		          [
		            13983,
		            13983,
		            {
		              'sources' => {
		                             'T' => [ <het codes>,<hom codes> ],
		                             'G' => [ <het codes>,<hom codes> ]
		                           },
		              'rsids' => [ 'rs200281958' ],
		              'chromosome' => 'chr5',
		              'end' => 13983,
		              'start' => 13983
		            }
		          ]
	        ];


	2) chromosome, position (base 0), and variant

		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    variant => "T") 
			returns: same as above
		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    variant => "A") 
			returns: nothing
		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    variant => "G") 
			returns: 
			[
	        	  [
        		    13983,
		            13983,
		            {
		              'sources' => {
		                             'T' => [ <het codes>,<hom codes> ],
		                             'G' => [ <het codes>,<hom codes> ]
                		           },
		              'rsids' => [ 'rs200281958' ],
		              'chromosome' => 'chr5',
		              'end' => 13983,
		              'start' => 13983
		            }
		          ]
		        ];

	3) chromosome, position (base 0), and people

		Position 13979 contains NA10851 as a person and position 13983 does not.

		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    people => ['NA10851']) 
			returns:
			[
 		         [
		            13979,
		            13979,
		            {
		              'sources' => {
		                             'T' => [ <het codes>,<hom codes> ],
		                             'C' => [ <het codes>,<hom codes> ]
		                           },
		              'rsids' => [ 'rs201837439' ],
		              'chromosome' => 'chr5',
		              'end' => 13979,
		              'start' => 13979
		            }
		          ]
		        ];

	4) chromosome, position, variant, and people 

		Position 13979 contains NA10851 as a person and position 13983 does not.
		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    variant => "T", people => ['NA10851']) 
			returns: same as #3
		$kav->knownVariants(chromosome => "chr5", start => 13979, end => 13983, 
				    variant => "G", people => ['NA10851']) 
			returns: nothing

=back

=cut
sub knownVariants {
	my $self = shift;
	my (%params) = @_;
	die "Kaviar knownVariants requires hash parameterization.  \$kav->knownVariants(chromosome => \"1\", start => \"2938\").  Optional array of people may be provided using the 'people' parameter. \$kav->knownVariants(chromosome => \"1\", start => \"2938\", people => \\\@people).  Optional variant filter may be provided using the 'variant' parameter.  \$kav->knownVariants(chromosome => \"1\", start => \"2938\", variant => \"C\")" unless scalar (keys %params) >= 2;
	die "Missing chromosome in knownVariants call" unless $params{'chromosome'};
	die "Missing start in knownVariants call" unless $params{'start'} || $params{'position'};
	# end_overlap and start_overlap don't mean anything for a SNV query

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'} || $params{'position'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		$start_basezero = $start_basezero-1;
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
	}
	my $start_baseone = $start_basezero +1;
	my $end_basezero = $params{'end'} || $start_basezero;
	if ($self->{'end_match'} == 0) {
		$end_basezero = $end_basezero-1;
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
	}
	my $end_baseone = $end_basezero+1;
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);
	my $reqvar = $params{'variant'} || undef;
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my $chri = $self->cleanupChrom($chrom);
	my $db = 'dbs';
	my $decoded = $self->{'decoded'};
	return unless $self->{$db}{'tbxchroms'}{$chri};
	my $tbx = $self->{$db}{'tbx'};
	# even though tabix on the command line lets you call basezero-basezero, 
	# the tabix command line includes both sides - a position at bazezero will be returned and aposition at baseone
	# the tabix perl module includes only the right side - a position at basezero will NOT be returned
	print "knownVariants querying $chri $start_basezero $end_baseone\n" if $self->{'warn'};
	print "knownVariants ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
	my $res = $tbx->query($chri, $start_basezero, $end_baseone);
	my @results;
	while (my $line = $tbx->read($res)) {
		my %res;
		my($c, $s, $rsids, @vars) = split /\t/, $line;
		my @rsids;
		while ($rsids) {
			my $rsid = 
			$decoded->[ord(substr($rsids,0,1))]*15625000 +
			$decoded->[ord(substr($rsids,1,1))]*62500 +
			$decoded->[ord(substr($rsids,2,1))]*250 +
			$decoded->[ord(substr($rsids,3,1))];
			push(@rsids, "rs$rsid");
			$rsids = substr($rsids, 4);
		}
		while (@vars) {
			my $var = shift @vars;
			if ($var eq '') { $var = '-'; }
			my $sources = shift @vars;
			# if people is defined, cull sources for those codes
			if (scalar @people > 0) {
				$sources = $self->cullCodes($sources, $people_codes_ref);
				next unless ($sources);
			}

			$res{'sources'}{$var} = $sources;
		}
		next unless $res{'sources'}; # skip if there were no sources added to the result
		if ($reqvar && !$res{'sources'}{$reqvar}) {
			print "Did not find requested variant ($reqvar) in found variants (". join(", ", keys %{$res{'sources'}})."). Position ".($s-1)." filtered.\n" if $self->{'warn'};
			next;
		}
		$res{'chromosome'} = $chrom;
		$res{'start'} = $s -1;
		$res{'end'} = $s -1;
		$res{'rsids'} = \@rsids if scalar @rsids > 0;
		push(@results, [$s-1, $s-1, \%res]);
	}
	return \@results;
}

=over 12

=item variantsInRange

	$kav->variantsInRange($chrom, $start, $end)

	A simplified interface to knownVariants that doesn't require the hash input 
	for backwards compatibility
	Queries the SNP database only.

=back

=cut
sub variantsInRange {
	my ($self, $chrom, $start, $end) = @_;

	my $res = $self->knownVariants(chromosome => $chrom, start => $start, end => $end);
	return @$res;
}


=over 12

=item variantPositionsInRange

	$kav->variantPositionsInRange($chrom, $start, $end)

	Returns location and RSID information for all variants in a range, but does
	not return full information on those positions.  Queries SNP database only.
	The data structure is a hash of start -> end -> array of rsIDs.
	Limitation: only returns SNVs that have an rsID.
	Probably superceded by allVariants.

=back

=cut
sub variantPositionsInRange {
	my $self = shift;
	my %params;
    if (scalar @_ == 3) {
        my ($chrom, $start, $end) = @_;
        $params{'chromosome'} = $chrom;
        $params{'start'} = $start;
        $params{'end'} = $end;
	} else {
	  (%params) = @_;
	}
	die "Missing chromosome in variantPositionsInRange call" unless $params{'chromosome'};
	die "Missing start in variantPositionsInRange call" unless $params{'start'} || $params{'position'};
	die "variantPositionsInRange doesn't support variant or people parameters" if $params{'variant'} || $params{'people'};
	# end_overlap and start_overlap don't mean anything for a SNV query

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'} || $params{'position'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		$start_basezero = $start_basezero-1;
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
	}
	my $start_baseone = $start_basezero +1;
	my $end_basezero = $params{'end'} || $start_basezero;
	if ($self->{'end_match'} == 0) {
		$end_basezero = $end_basezero-1;
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
	}
	my $end_baseone = $end_basezero+1;
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);
	my $reqvar = $params{'variant'} || undef;
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my $chri = $self->cleanupChrom($chrom);
	my $db = 'dbs';
	my $decoded = $self->{'decoded'};
	return unless $self->{$db}{'tbxchroms'}{$chri};
	my $tbx = $self->{$db}{'tbx'};
	# even though tabix on the command line lets you call basezero-basezero, 
	# the tabix command line includes both sides - a position at bazezero will be returned and aposition at baseone
	# the tabix perl module includes only the right side - a position at basezero will NOT be returned
	print "variantPositionsInRange querying $chri $start_basezero $end_baseone\n" if $self->{'warn'};
	print "variantPositionsInRange ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
	my $res = $tbx->query($chri, $start_basezero, $end_baseone);
	my %res;
	while (my $line = $tbx->read($res)) {
		my($c, $s, $rsids, @vars) = split /\t/, $line;
		my @rsids;
		while ($rsids) {
			my $rsid = 
			$decoded->[ord(substr($rsids,0,1))]*15625000 +
			$decoded->[ord(substr($rsids,1,1))]*62500 +
			$decoded->[ord(substr($rsids,2,1))]*250 +
			$decoded->[ord(substr($rsids,3,1))];
			push(@rsids, "rs$rsid");
			$rsids = substr($rsids, 4);
		}
		$res{$s-1}{$s} = \@rsids if scalar @rsids > 0;
	}
	return %res;

}


=over 12

=item knownVariantPositions


	$kav->knownVariantPositions(chromosome => $chrom, position => $pos, 
				    start => $start, [end => $end, var => $var, 
				    people => \@people])

	Takes a chromosome, a start position, an optional end position, an optional 
	  variant base pair (eg: ACGT), and an optional array of people names or IDs.
	Searches the SINGLE NUCLEOTIDE VARIANT DATABASE ONLY.  
	If a variant is specified but not observed, the result is empty.
	Otherwise, it returns HASH of all variants in the specified query window

	For examples, see knownVariants.

=back

=cut

# this works with the method of filtering the reqVar in filterResults because there are no
# multiple basepair variants since this only reads the SNP file
# mostly here because it returns the hash version of results
sub knownVariantPositions {
	my $self = shift;
	my (%params) = @_;
	die "Kaviar knownVariantPositions requires hash parameterization.  \$kav->knownVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\").  Optional array of people may be provided using the 'people' parameter. \$kav->knownVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\", people => \\\@people)" unless scalar (keys %params) >= 3;
	die "Missing chromosome in knownVariantPositions call" unless $params{'chromosome'};
	die "Missing start in knownVariantPositions call" unless $params{'start'};
	die "Missing end in knownVariantPositions call" unless $params{'end'};
	die "end_overlap cannot be set to OFF for knownVariantPositions" if $self->{'end_overlap'} eq 0;

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		$start_basezero = $start_basezero-1;
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
	}
	my $end_basezero = $params{'end'};
	if ($self->{'end_match'} == 0) {
		$end_basezero = $end_basezero-1;
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
	}

	my $reqvar = $params{'variant'} || undef;
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);

	my $start_baseone = $start_basezero+1;
	my $end_baseone = $end_basezero+1;
	if ($self->{'end_match'} == 0) {
		# have to reduce baseone as well because tabix defaultly returns results
		# that overlap the end basepair, while it doesn't do that for start position
		$end_baseone = $end_baseone - 1;
	}

	my $chri = $self->cleanupChrom($chrom);
	my $decoded = $self->{'decoded'};
	my %res;
	foreach my $db (qw/dbs/) {
		next unless $self->{$db}{'tbxchroms'}{$chri};
		my $tbx = $self->{$db}{'tbx'};
	print "knownVariantPositions ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
		my $res = $tbx->query($chri, $start_basezero, $end_baseone);
		while (my $line = $tbx->read($res)) {
			my $ref = $self->processLineToPositions(\%res, $db, $line, $people_codes_ref, $reqvar);
			%res = %$ref;
		}
	}

	return \%res unless scalar keys %res > 0;
	my $filtered_results = $self->filterResults('HASH', $start_basezero, $end_baseone, \%res, $people_codes_ref);

	return $filtered_results;
}



#########################################################################################################
#########################################################################################################
# RANGE PARSING ONLY
#########################################################################################################
#########################################################################################################

=over 12

=item rangeVariants 

	$kav->rangeVariants(chromosome => <chromosome>, start => <start>, 
			    end => <end>, variant => <variant>, 
			    people => <array of people names or IDS>);

	Return all RANGED ONLY variants in the range (0-based)$start to (1-based)$end of $chrom.
	The output is a list of elements, each describing one genomic range as a 
	  three-element array.
	The first two elements in the array are the (1-based) beginning and end of 
	  the variant.
	The third element is a reference to a hash describing the observed variation, 
	  in the same format as that returned by knownVariant(). Use digest() to process 
	  these data structures.

	In other words:
	1) $kav->rangeVariants(chromosome => <chromosome>, start => <start>, end => <end>);
		will list all known variants in that range.
	2) $kav->rangeVariants(chromosome => <chromosome>, start => <start>, end => <end>, 
			   person => <array of people>);
		will list all known variants in those people in that range
	3) $kav->rangeVariants(chromosome => <chromosome>, start => <start>, end => <end>, 
	                  variant => <variant>)
		will list variants that match the supplied variant in that range
	4) $kav->rangeVariants(chromosome => <chromosome>, start => <start>, end => <end>, 
			  variant => <variant>, person => <array of people>)
		will list specified variants in those people in that range

	Results are by variant and do not combine multiple variants together

=back

=cut
sub rangeVariants {
	my ($self, %params) = @_;
	die "Kaviar AllVariants requires hash parameterization.  \$kav->rangeVariants(chromosome => \"1\", start => \"2938\", end => \"2940\").  Optional array of people may be provided using the 'people' parameter. \$kav->rangeVariants(chromosome => \"1\", start => \"2938\", end => \"2940\", people => \\\@people)" unless scalar (keys %params) >= 3;
	die "Missing chromosome in rangeVariants call" unless $params{'chromosome'};
	die "Missing start in rangeVariants call" unless $params{'start'};
	die "Missing end in rangeVariants call" unless $params{'end'};
	die "end_overlap cannot be set to OFF for rangeVariants" if $self->{'end_overlap'} eq 0;

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		$start_basezero = $start_basezero-1;
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
	}
	my $end_basezero = $params{'end'};
	if ($self->{'end_match'} == 0) {
		$end_basezero = $end_basezero-1;
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
	}

	my $reqvar = $params{'variant'} || undef;
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);

	my $start_baseone = $start_basezero+1;
	my $end_baseone = $end_basezero+1;
	if ($self->{'end_match'} == 0) {
		# have to reduce baseone as well because tabix defaultly returns results
		# that overlap the end basepair, while it doesn't do that for start position
		$end_baseone = $end_baseone - 1;
	}
	my $diff = $end_baseone - $start_baseone;

	my $chri = $self->cleanupChrom($chrom);
	my @res;
	my $decoded = $self->{'decoded'};
	next unless $self->{'dbr'}{'tbxchroms'}{$chri};
	my $tbx = $self->{'dbr'}{'tbx'};
	print "rangeVariants querying $chri $start_basezero $end_baseone\n" if $self->{'warn'};
	# use start_basezero because tabix doesn't return matches on the start position
	print "rangeVariants ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
	my $res = $tbx->query($chri, $start_basezero, $end_baseone);
	while (my $line = $tbx->read($res)) {
		# kaviar db results are baseone
		my($c, $s_baseone, $e_baseone, $rsids, @vars);
		($c, $s_baseone, $e_baseone, $rsids, @vars) = split /\t/, $line;
		my %single_res;
		my %rsids;
		while ($rsids) {
			my $rsid = 
			$decoded->[ord(substr($rsids,0,1))]*15625000 +
			$decoded->[ord(substr($rsids,1,1))]*62500 +
			$decoded->[ord(substr($rsids,2,1))]*250 +
			$decoded->[ord(substr($rsids,3,1))];
			$rsids{"rs$rsid"} =1;
			$rsids = substr($rsids, 4);
		}

		while (@vars) {
			my $var = shift @vars;
			if ($var eq '') { $var = '-'; }
			my $sources = shift @vars;
			# if people is defined, cull sources for those codes
			if (scalar @people > 0) {
				$sources = $self->cullCodes($sources, $people_codes_ref);
				next unless ($sources);
			}

			my $l_var = length($var); # if var was empty string and is now '-', $l_var will be 1. Problem?
			my $l_ref = ($e_baseone - $s_baseone) +1;

			# if it is a substitution followed by an insertion then the end bp will equal the start bp
			if ($l_var > $l_ref && $l_ref <= 1) {
				my $refCall = $self->referenceCall({ sources=> { $var => $sources}});
				if ($l_ref <= 1 && $refCall eq $var) {
					print "Warning: Skipping this variant $var at position $s_baseone to $e_baseone because it is the complement to a single base pair insertion.\n" if $self->{'warn'};
					next;
				}
				if ($self->{'exclude_insertions'}) {
					print "Warning: Skipping insertion variant $var at position $e_baseone because exclude_insertions is on\n" if $self->{'warn'};
					next;
				}
				$var.= '_INSERTION' if $l_ref == 0;
				$self->{'end_overlap'} = 1;
			}

			$single_res{'sources'}{$var} = $sources;
		}
		next unless $single_res{'sources'}; # skip if there were no sources added to the result
		my $s = $s_baseone - 1;  # convert back to basezero
		my $e = $e_baseone - 1;  # convert back to basezero

		if ($single_res{'sources'}{'-'}) {
			$single_res{'sources'}{'DELETION'} = $single_res{'sources'}{'-'};
			delete $single_res{'sources'}{'-'};
		}
		if ($reqvar && !$single_res{'sources'}{$reqvar}) {
			print "Did not find requested variant ($reqvar) in found variants (". join(", ", keys %{$single_res{'sources'}})."). Position ".$s." filtered.\n" if $self->{'warn'};
			next;
		}

		$single_res{'chromosome'} = $chrom;
		$single_res{'start'} = $s;
		$single_res{'end'} = $e;
		my @rsids = keys %rsids;
		$single_res{'rsids'} = \@rsids if scalar @rsids > 0;

		push @res, [$s, $e, \%single_res];
	}
	return \@res;
}

=over 12

=item rangeVariantPositions

	$kav->rangeVariantPositions(chromosome => <chromosome>, start => <start>, 
				    end => <end>, variant => <variant>, 
				    people => <array of people names or IDS>);

	Return all RANGED ONLY variants in the range (0-based)$start to (1-based)$end of $chrom.
	The output is a list of elements, each describing one genomic range as a 
	  three-element array.
	The first two elements in the array are the (1-based) beginning and end of 
	  the variant.
	The third element is a reference to a hash describing the observed variation, 
	  in the same format as that returned by knownVariant(). Use digest() to process 
	  these data structures.

	In other words:
	1) $kav->rangeVariantPositions(chromosome => <chromosome>, start => <start>, 
	 			       end => <end>);
		will list all known variants in that range.
	2) $kav->rangeVariantPositions(chromosome => <chromosome>, start => <start>, 
				       end => <end>, person => <array of people>);
		will list all known variants in those people in that range
	3) $kav->rangeVariantPositions(chromosome => <chromosome>, start => <start>, 
				       end => <end>, variant => <variant>)
		will list variants that match the supplied variant in that range
	4) $kav->rangeVariantPositions(chromosome => <chromosome>, start => <start>, 
				       end => <end>, variant => <variant>, 
				       person => <array of people>)
		will list specified variants in those people in that range

	Results are by variant and do not combine multiple variants together


=back

=cut
sub rangeVariantPositions {
	my $self = shift;
	my (%params) = @_;
	die "Kaviar AllVariants requires hash parameterization.  \$kav->rangeVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\").  Optional array of people may be provided using the 'people' parameter. \$kav->rangeVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\", people => \\\@people)" unless scalar (keys %params) >= 3;
	die "Missing chromosome in rangeVariantPositions call" unless $params{'chromosome'};
	die "Missing start in rangeVariantPositions call" unless $params{'start'};
	die "Missing end in rangeVariantPositions call" unless $params{'end'};
	die "end_overlap cannot be set to OFF for rangeVariantPositions" if $self->{'end_overlap'} eq 0;

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
		$start_basezero = $start_basezero-1;
	}
	my $end_basezero = $params{'end'};
	if ($self->{'end_match'} == 0) {
		$end_basezero = $end_basezero-1;
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
	}

	my $reqvar = $params{'variant'} || undef;
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);

	my $start_baseone = $start_basezero+1;
	my $end_baseone = $end_basezero+1;
	if ($self->{'end_match'} == 0) {
		# have to reduce baseone as well because tabix defaultly returns results
		# that overlap the end basepair, while it doesn't do that for start position
		$end_baseone = $end_baseone - 1;
	}

	my $chri = $self->cleanupChrom($chrom);
	my $decoded = $self->{'decoded'};
	my %res;
	foreach my $db (qw/dbr/) {
		next unless $self->{$db}{'tbxchroms'}{$chri};
		my $tbx = $self->{$db}{'tbx'};
	print "rangeVariantPositions ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
		my $res = $tbx->query($chri, $start_basezero, $end_baseone);
		while (my $line = $tbx->read($res)) {
			my $ref = $self->processLineToPositions(\%res, $db, $line, $people_codes_ref, $reqvar);
			%res = %$ref;
		}
	}
	return \%res unless scalar keys %res > 0;
	my $filtered_results = $self->filterResults('HASH', $start_basezero, $end_baseone, \%res, $people_codes_ref);

	return $filtered_results;
}


#########################################################################################################
#########################################################################################################
# ALL VARIANTS - SNP AND RANGE PARSING
#########################################################################################################
#########################################################################################################


# allVariants returns results by variant.
#
# Parameters:  A hash of parameters.  Requires 'chromosome', 'start', and 'end'.
# Optional Parameters:  An arrayref of people's names or IDs in the Kaviar database.
#			A 'type' parameter that can specify the HASH or ARRAY return.
# Returns: An arrayref of a triplet [start, stop, information] where information is
#           a hashref that can be provided to the decant subroutine to calculate
#           heterozygosity, allele frequency, and other statistics.
# Optional return: A hash of hashes of information about the position keyed on position 
#                    start and position end.  Identical formatting to allVariantPositions


=over 12

=item allVariants


	$kav->allVariants(chromosome => <chromosome>, start => <start>, end => <end>,
			  [variant => <variant>, people => <array of people names or IDs>);

	allVariants takes a hash of parameters describing a location in the genome
	and searches the SNP, ranged, and genotypes files of Kaviar for matches, then returns
	an arrayref of arrays, one for each VARIANT in or overlapping the query window.
        Kaviar settings of start_overlap, end_overlap, start_match, and end_match are
        ignored; all overlaps/matches are returned (TODO: correct this). Each
	array has start, end, and a hash of information about that VARIANT, which
	can be used to calculate statistics using the decoct method.

	Required Parameters:
		chromosome => <chromosome>
		start => <start position basezero>
		end => <end position basezero>

	Optional Parameters:
		variant => <string of variant ATGC>
		people => <arrayref of people names or IDs>
		type => 'HASH'

	Example calls:
		$kav->allVariants(chromosome => "5", 
				  start => "13740", 
				  end => "13744")
		$kav->allVariants(chromosome => "5", 
				  start => "13740", 
				  end => "13744", 
				  people => ["NA12877","Desmond Tutu"])
	
	Returns:
		An array triplet with basezero start, basezero end, and a hashref
                of information about the position
		[ 
		   [14255,14259,
		     { 
		   	sources => {
		 			'CCCAA' =>  codes,
					'CTCAA' =>  codes
				   },
			sourcesGT => {},
			chrom => 'chr5',
			start => '14255',
			end => '14256'
		     }
		   ],
		   [14256,14256,
		     { 
		     	sources => { 
					'C' => codes, 
					'T' => codes 
				   }, 
		     	sourcesGT => { 
					'C/C' => codes, 
					'C/T' => codes,
					'T/T' => codes 
				   }, 
		        rsids => [ rs112910437 ], 
			chrom => 'chr5', 
		        start => '14256',
		        end => '14256' 
                     }
		   ]
		]

=back

=cut

sub allVariants {
  my $self = shift;
  my (%params) = @_;
  die "Kaviar AllVariants requires hash parameterization.  \$kav->allVariants(chromosome => \"1\", start => \"2938\", end => \"2940\").  Optional array of people may be provided using the 'people' parameter. \$kav->allVariants(chromosome => \"1\", start => \"2938\", end => \"2940\", people => \\\@people)" unless scalar (keys %params) >= 3;
  die "Missing chromosome in allVariants call" unless $params{'chromosome'};
  die "Missing start in allVariants call" unless defined $params{'start'};
  die "Missing end in allVariants call" unless $params{'end'};
  my $chrom = $params{'chromosome'};
  my $start_basezero = $params{'start'};
  if ($self->{'start_overlap'} || $self->{'start_match'}) {
    print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
    $start_basezero = $start_basezero-1;
  }
  my $end_basezero = $params{'end'};
  my $reqvar = $params{'variant'};
  print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
  my @people = @{$params{'people'}} if $params{'people'};
  print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});

  my $start_baseone = $start_basezero+1;
  my $end_baseone = $end_basezero+1;
  if ($self->{'end_match'} == 0) {
    print "end_match set to off.\n" if $self->{'warn'};
    # move the window one to the left in order to avoid any positions that are end match
    $end_baseone = $end_baseone -1;
    $start_basezero = $start_basezero -1;
    if ($start_baseone == $end_baseone || $start_baseone > $end_baseone) {
      $start_basezero = $start_basezero-1;
    }
  }
  my $diff = $end_baseone - $start_baseone;
  die "Maximum limit on query size set to ".$self->maxlimit.".  Please bin your queries." if $diff > $self->maxlimit;
  my $chri = $self->cleanupChrom($chrom);
  my $decoded = $self->{'decoded'};
  my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);
  my %res;
  my @db;
  push @db, "dbs" unless $params{'no_dbs'};
  push @db, "dbr" unless $params{'no_dbr'};
  push @db, "dbg" unless $params{'no_dbg'};
  foreach my $db (@db) {
    next unless $self->{$db}{'tbxchroms'}{$chri};
    print "PROCESSING DB $db $self->{'frz'} $chri $start_basezero $end_baseone\n" if $self->{'warn'};
    my $index = ($db eq 'dbg') ? 'sourcesGT': 'sources';
    my $tbx = $self->{$db}{'tbx'};
    print "allVariants ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};

    # Get all Kaviar database lines  within the range desired.
    # Tabix query retrieves all lines that overlap the region.
    # Effectively, start_match, start_overlap, end_match, and end_overlap
    # are turned on.
    #print "Calling tabix on $chri, $start_basezero, $end_baseone\n";
    my $res = $tbx->query($chri, $start_basezero, $end_baseone);

    # Process each line returned
    while (my $line = $tbx->read($res)) {
      my ($chrom, $s_baseone, $e_baseone, $rsids, @vars);
      if (($db eq 'dbr') || ($db eq 'dbg')) {
	($chrom, $s_baseone, $e_baseone, $rsids, @vars) = split /\t/, $line;
      } elsif ($db eq 'dbs') {
	($chrom, $s_baseone, $rsids, @vars) = split(/\t/, $line);
	$e_baseone = $s_baseone;
      }
      #print " allVariants read tabix line $chrom $s_baseone $e_baseone @vars\n";

      my %single_res;

      # Decode the string describing the dbSNP rsid
      my %rsids;
      while ($rsids) {
	my $rsid = 
	$decoded->[ord(substr($rsids,0,1))]*15625000 +
	$decoded->[ord(substr($rsids,1,1))]*62500 +
	$decoded->[ord(substr($rsids,2,1))]*250 +
	$decoded->[ord(substr($rsids,3,1))];
	$rsids{"rs$rsid"} = 1;
	$rsids = substr($rsids, 4);
      }

      # Process each variant on this Kaviar database line
      while (@vars) {
	my $var = shift @vars;
	if ($var eq '') { $var = '-'; }

	# Gather the sources this var is seen in.
	my $sources = shift @vars;

	# if people is defined, cull sources for those codes
	if (scalar @people > 0) {
	  $sources = $self->cullCodes($sources, $people_codes_ref);
	  next unless ($sources);
	}

	# Check length of variant
	my $l_var = length($var);  # if var was empty string and is now '-', $l_var will be 1. Problem?
	if ($l_var > $self->{'maxlength'}) {
	  print "WARNING: Skipping this variant $var at position $s_baseone to $e_baseone becasue it is larger than the maxlength setting of $self->{'maxlength'}\n" if $self->{'warn'};
	  next;
	}

	# Handle some complex cases

	# Get the length of this variant as specified by its start/end.
	# If it's different from $l_var, we have an indel
	# 10/29/15: not sure this applies to dbg; didn't take time to check
	my $l_ref = ($e_baseone - $s_baseone) +1;

	# if it is a substitution followed by an insertion then the end bp will equal the start bp
	if ($l_var > $l_ref && $l_ref <= 1) {
	  my $refCall = $self->referenceCall({ sources=> { $var => $sources}});
	  if ($l_ref <= 1 && $refCall eq $var) {
	    print "Warning: Skipping this variant $var at position $s_baseone to $e_baseone because it is the complement to a single base pair insertion.\n" if $self->{'warn'};
	    next;
	  }
	  if ($self->{'exclude_insertions'}) {
	    print "Warning: Skipping insertion variant $var at position $e_baseone because exclude_insertions is on\n" if $self->{'warn'};
	    next;
	  }
	  $var.= '_INSERTION' if $l_ref == 0;
	  $self->{'end_overlap'} = 1; #what? why? TMF 11/15/14
	}
	$single_res{$index}{$var} = $sources;
      }

      next unless $single_res{$index}; # skip if there were no sources added to the result
      my $s = $s_baseone - 1;  # convert back to basezero
      my $e = $e_baseone - 1;  # convert back to basezero

      if ($single_res{$index}{'-'}) {
	$single_res{$index}{'DELETION'} = $single_res{$index}{'-'};
	delete $single_res{$index}{'-'};
      }

      # If the user is trying to confirm a specific variant,
      # and we didn't retrieve it, warn & don't return anything.
      # Possible TODO: if we did retrieve it, we are returning
      # all variants retrieved. Don't we want to return just the specific one?
      unless ($db eq 'dbg') {
	if ($reqvar && !$single_res{'sources'}{$reqvar}) {
	  print "Did not find requested variant ($reqvar) in found variants (". join(", ", keys %{$single_res{'sources'}})."). Position ".$s." filtered.\n" if $self->{'warn'};
	  next;
	}
      }

      # Fill in remaining fields for this variant
      $single_res{'chromosome'} = 'chr'.$chrom;
      $single_res{'start'} = $s;
      $single_res{'end'} = $e;
      $single_res{'rsids'} = \%rsids if scalar keys %rsids > 0;

      # If we already have some variants/genotypes for this start/end,
      # add this one on.
      if ($res{$s}{$e}) {
	my $old_result = $res{$s}{$e};
	foreach my $variant (keys %{$single_res{$index}}) {
	  if ($old_result->{$index}{$variant}) {
	    my $joined_res = $self->checkAddCodes($single_res{$index}{$variant}, $old_result->{$index}{$variant});
	    $old_result->{$s}{$e}{$variant} = $joined_res;
	  } else {
	    $old_result->{$index}{$variant} = $single_res{$index}{$variant};
	  }
	  if ($old_result->{'rsids'} && $single_res{'rsids'}) {
	    map {$old_result->{'rsids'}{$_}} keys %{$single_res{'rsids'}};
	  }
	}
	%single_res = %$old_result;
      }
      # Store the resulting, possibly compound, variant in %res
      $res{$s}{$e} = \%single_res;
    }
  }

  my @res;
  my %filtered_res;
  foreach my $start (sort {$a<=>$b} keys %res) {
    my %end_hash = %{$res{$start}};
    foreach my $end (sort {$a<=>$b} keys %end_hash) {
      if ($start < $params{'start'} && $end < $params{'start'}) {
	#print "skipping $start $end ".Dumper $res{$start}{$end};
	next;
      }
      if ($start > $params{'end'} && $end > $params{'end'}) {
	#print "skipping $start $end ".Dumper $res{$start}{$end};
	next;
      }
      if ($res{$start}{$end}{'rsids'}) {
	my @rsids = keys %{$res{$start}{$end}{'rsids'}};
	$res{$start}{$end}{'rsids'} = \@rsids;
      }
      $filtered_res{$start}{$end} = $res{$start}{$end};
      push(@res, [$start, $end, $res{$start}{$end}]);
    }
  }

  if ($params{'type'} && uc($params{'type'}) eq 'HASH') {
    return \%filtered_res;
  } else {
    return \@res;
  }
}

# allVariantPositions returns results by position in a hash.  This may result in a 
# deletion or substition that is misaligned being returned in the wrong position, but 
# users can use other queries to determine what the positioning is of the deletion 
# or substitution
#
# Parameters:  A hash of parameters.  Requires 'chromosome', 'start', and 'end'.
# Optional Parameters:  An arrayref of people's names or IDs in the Kaviar database.
# Returns: A hashref of information that can be provided to the decant subroutine to
#           calculate heterozygosity, allele frequency, and other statistics.

=over 12

=item allVariantPositions

	$kav->allVariantPositions(chromosome => <chromosome>, start => <start>, 
				  end => <end>, [variant => <variant>, 
				  people => <array of people names or IDs>);

	allVariantPositions takes a hash of parameters describing a location in 
	the genome and searches the SNP and ranged files of Kaviar for matches, 
	then returns a hashref of positions in the query window.  It attempts to
	integrate ranged variant information with snp variant information over
	the window of positions.  This can lead to odd results when dealing with 
	ranged variants that are prone to misalignment, so always check results.
	This function runs very slowly when used with very long variants.

	NOTE: Does not yet search the genotype files.

	Required Parameters:
		chromosome => <chromosome>
		start => <start position basezero>
		end => <end position basezero>

	Optional Parameters:
		variant => <variant>
		people => <arrayref of people names or IDs>

	Example calls:
		$kav->allVariantPositions(chromosome => "5", 
				  start => "13740", 
				  end => "13744")
		$kav->allVariantPositions(chromosome => "5", 
				  start => "13740", 
				  end => "13744", 
				  people => ["NA12877","Desmond Tutu"])
	
	Returns:
		A hash of hashes of information about the position keyed on position start and position end
		{ '13740' => {
			'13740' => { 
				  	sources => { 
						'T' => ['het codes','hom codes'], 
						'DELETION' => ['het codes', 'hom codes'] 
					     },	 
		        		  rsids => [ rs72498228 ], 
	 			          chrom => 'chr5', 
				          start => '13740', 
			  	  	  end => '13740' 
 		             	   }
			}
		}

=back

=cut

sub allVariantPositions {
	my $self = shift;
	my (%params) = @_;
	die "Kaviar AllVariantPositions requires hash parameterization.  \$kav->allVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\").  Optional array of people may be provided using the 'people' parameter. \$kav->allVariantPositions(chromosome => \"1\", start => \"2938\", end => \"2940\", people => \\\@people)" unless scalar (keys %params) >= 3;
	die "Missing chromosome in allVariantPositions call" unless $params{'chromosome'};
	die "Missing start in allVariantPositions call" unless $params{'start'};
	die "Missing end in allVariantPositions call" unless $params{'end'};

	my $chrom = $params{'chromosome'};
	my $start_basezero = $params{'start'};
	if ($self->{'start_overlap'} || $self->{'start_match'}) {
		print "Start match or overlap set to off, so won't match variants that match start basepair.\n" if $self->{'warn'};
		$start_basezero = $start_basezero-1;
	}
	my $end_basezero = $params{'end'};
	# if this is being called with a single position, we need to increment end to have it return results
	# because of how the Tabix.pm works.
	if ($start_basezero == $end_basezero) {
		$end_basezero++;
	}
	# do not return the end basepair or any information about it
	if ($self->{'end_match'} == 0 || $self->{'end_overlap'} == 0) {
		print "End match or overlap set to off, so won't match variants that match end basepair.\n" if $self->{'warn'};
		$end_basezero = $end_basezero-1;
	}
	my $reqvar = $params{'variant'};
	print STDERR "Warning: Filtering results by required variant: $reqvar\n" if ($reqvar && $self->{'warn'});
	my @people = @{$params{'people'}} if $params{'people'};
	print STDERR "Warning: Filtering results by people: ".join(", ", @people)."\n" if (scalar @people > 0 && $self->{'warn'});
	my $people_codes_ref = $self->getCodes(\@people) if (scalar @people > 0);

	my $start_baseone = $start_basezero+1;
	my $end_baseone = $end_basezero+1;
	my $diff = $end_baseone - $start_baseone;
	die "Maximum limit on query size set to ".$self->maxlimit.".  Please bin your queries." if $diff > $self->maxlimit;
	my $chri = $self->cleanupChrom($chrom);
	my $decoded = $self->{'decoded'};
	my %byPos;
	foreach my $db (qw/dbs dbr/) {
		next unless $self->{$db}{'tbxchroms'}{$chri};
		my $tbx = $self->{$db}{'tbx'};
		print "allVariantPositions querying $db $chri $start_basezero $end_baseone\n" if $self->{'warn'};
		print "allVariantPositions ERROR: Start and end base pairs are in error.  start: $start_basezero end: $end_baseone.  This query will return no results.  Change your query window parameters or base pairs.\n" if ( ($start_basezero == $end_baseone) || ($start_basezero > $end_baseone) ) && $self->{'warn'};
		my $res = $tbx->query($chri, $start_basezero, $end_baseone);
		while (my $line = $tbx->read($res)) {
			my $ref = $self->processLineToPositions(\%byPos, $db, $line, $people_codes_ref, $reqvar);
			%byPos = %$ref;
		}
	}
	return \%byPos unless scalar keys %byPos > 0;
	my $filtered_results = $self->filterResults('HASH',$start_basezero, $end_baseone, \%byPos, $people_codes_ref);
	return $filtered_results;

}


#########################################################################################################
#########################################################################################################
#
#  STATISTICS GENERATION METHODS
#
#########################################################################################################
#########################################################################################################

=over 12

=item digest

	$kav->digest(hashref or arrayref)

		Parse results that are in the format 
				        {                 sources => {
                                                      'T' => ['het codes','hom codes'],
                                                      'DELETION' => ['het codes', 'hom codes']
                                                   },
                                          rsids => [ rs72498228 ],
                                          chromosome => 'chr5',
                                          start => '13740',
                                          end => '13740'
                                        }

                         OR
                           [ start_position,
                             end_position,
                              {   sources => {
                                               'T' => ['het codes','hom codes'],
                                               'DELETION' => ['het codes', 'hom codes']
                                             },
                                  rsids => [ rs72498228 ],
                                  chromosome => 'chr5',
                                  start => '13740',
                                  end => '13740'
                               }
                  ]


	Parameters: Takes in a reference to a hash of position information
	Returns: A hashref of statistics
        Works for both alleles and genotypes.

=back

=cut
sub digest {
  my($self, @rest) = @_;
  my %res;
  if (scalar @rest == 1) {
    my ($start, $end, $hashRef) = @{$rest[0]};
    %res = %$hashRef;
  } elsif (scalar @rest == 3) {
    my ($start, $end, $hashRef) = @rest;
    %res = %$hashRef;
  } else {
    %res = @rest;
  }
  die "No information supplied to digest" unless scalar keys %res > 0;
  die "Incorrect information supplied to digest.  If the following says HASH, you need to dereference it before you pass it into digest. If the following is a hash that contains another hash, you need to provide the interior hash to digest. ".Dumper \%res unless $res{'sources'}; 
  if (ref(\%res) ne 'HASH') {
    die "Invalid value supplied to digest.  Digest processes a hash of information about a position.  If you want to provide an array of [start, end, info] then use decoct.";
  }
  my(%seen, @hist, %count, %id_count, %stat);
  my $adjust_for_phase3 = 1;
  my $refcode = $self->{'referenceCode'};
  my $complete = $self->{'completeGenomes'}; # indiv. complete genome seqs
  my $completeIndependent = $self->{'completeIndependent'};
  my $chrom = $res{'chromosome'};
  $chrom = substr($chrom, 3) if substr($chrom, 0, 3) eq 'chr';
  my $in_exome = $self->testInExome($chrom, $res{'start'}, $res{'end'});
  my $info = $self->{'info'};
  my $completeN = $completeIndependent;  # unrelated indiv. complete genome seqs
  my $refCall = '';
  my %alleles;
  my %AN;
  my %all_obs;
  my %platform_counts = ();
  my $total_platform_counts = 0;
  #print "Entering digest()\n";

  # For each allele/genotype observed at this position, tally various counts
  for my $calltype ('alleles', 'genotypes') {
    my $index = ($calltype eq 'alleles') ? 'sources' : 'sourcesGT';
    my %phase3_counts;
    my %exac_counts;
    foreach my $call (keys %{$res{$index}}) {
      # Get a hash of data sources the allele/genotype was seen in
      my %obs = %{$self->separateCodes($res{$index}{$call}, 1)};
      $all_obs{$calltype}{$call} = \%obs;
      $phase3_counts{$call} = 0;
      # For each data source
      for my $code (keys %obs) {

	# If data source is ref genome, store ref call and continue
	if ($code eq $refcode) {
	  $refCall = $call;
	  $stat{'refcall'} = $refCall;
	  next;
	}
	my ($person) = $self->whois($code);
	my $id =  $self->{'info'}->{$code}->{'id'};
	my $AC = $obs{$code}->{'AC'};
	#my $AN = $obs{$code}->{'AN'};  # we don't store this here anymore
	#$AN{$code} = $AN;
	next if $info->{$code}{'molecule type'} eq 'RNA';

	# If count is 0 (dbSNP), set it to ploidy (1)
	#$AC ||= $ploidy;
	my $project = $info->{$code}{'project'};
	$AC = 1 if $project eq 'dbSNP';
	$phase3_counts{$call} += $AC if $project eq '1000Genomes';
	$exac_counts{$call} += $AC if $project eq '63000exomes';
	# skip if AC=0 (e.g. some vars in 2000DanishExomes)
	# This actually messes up when the total AC ends up being zero:
	# the test below, "if ($sortedAlleles[1])", becomes false
	# and maf is not set.
	# We need a more sophisticated approach to handle AC=0
	# in Kaviar db.  08/27/15
	# Nov 2015: for now we are removing AC=0 from 2000DanishExomes
#      unless ($AC>0) {
#	delete $obs{$code};
#	next;
#      }
	# Increment the tally of observations in this data source
	$seen{$calltype}{$code} += $AC;

	# Set flag saying that we saw this allele in this data source
	$alleles{$calltype}{$code}{$call} = 1;

	# Increment the tally of observations of this allele
	$count{$calltype}{$call} += $AC;

	# development
	#$id_count{$person}{$call} += $AC;
	#print "==> person: $person call:$call count:$AC\n";
      }
    }

    # Go through all the calls again and adjust/count some stuff.
    foreach my $call (keys %{$res{$index}}) {
      unless ($call eq $refCall) {
	if ($calltype eq 'alleles') {
	  my %obs = %{$all_obs{'alleles'}{$call}};
      # Tally allele observations per platform
	for my $code (keys %obs) {
	    my $project = $info->{$code}{'project'};
	    my $platform = $self->{'projectinfo'}{$project}{'Depth'} || '';
	    if ($platform) {
	      $platform = 'CGI' if $platform =~ /cgi/i;
	      $platform = 'Illumina' if $platform =~ /illumina/i;
	      $platform_counts{$platform} += $obs{$code}->{'AC'};
	      $total_platform_counts += $obs{$code}->{'AC'};
	    }
	  }
	}
	# ExAC includes 1000Genomes; avoid double counting
	if ($adjust_for_phase3 && $in_exome) {
	  my $adjustment = 0;
	  $adjustment = $phase3_counts{$call} if defined $phase3_counts{$call};
	  $exac_counts{$call} = 0 unless defined $exac_counts{$call};
	  $adjustment = $exac_counts{$call} if $exac_counts{$call} < $phase3_counts{$call};
	  $count{$calltype}{$call} -= $adjustment;
	}
      }
    }
  }

  print STDERR "Warning: Reference information was not included in position, using unspecified 'Reference' instead of exact basepairs.\n" if !$refCall && $self->{'warn'};
  $refCall = 'Reference' if !$refCall;
  my $seenComplete = 0;
  my( $inWES, $inSNP);

  # For each data source in which any non-ref allele was seen
  foreach my $code (keys %{$seen{'alleles'}}) {

    next if $info->{$code}{'molecule type'} eq 'RNA';

    my $id = $info->{$code}{'id'};
    my $name = $self->{'code2name'}->{$code};
    my $genome_coverage = $info->{$code}{'Genome coverage'};
    $inWES = 1 if $genome_coverage =~ /WES/i;
    $inSNP = $code if $info->{$code}{'project'} eq 'dbSNP';

    #my $an = $info->{$code}{'AN'} || 2;
    # If the number of variant alleles is less than the total possible
    # and if this data source is not reference, increment the reference
    # tally by the difference.
#		if ($seen{$code}<$an && !$alleles{$code}{$refCall}) {
#		  $count{$refCall} += ($an-$seen{$code});
#		}
    # Increment the tally of indiv. genomes in which a non-ref
    # allele was seen at this position.
    $seenComplete++ if $complete->{$code};
  }

  ## If any alternate allele is seen only in dbSNP, we give it a count of one.
  ## If it is also seen elsewhere, we decrement the count because
  ## probably the dbSNP observation duplicates one of our other observations.
  if ($inSNP) {
    foreach my $call (keys %{$alleles{$inSNP}}) {
      $count{'alleles'}{$call}-- if $count{'alleles'}{$call} > 1;
    }
  }

  # On bobama this takes about 0.0003s (1s per 3500 calls) if pre-loading used.
  $stat{'inExHet'} = $self->testInAnnotationRegion(
    $self->{'special'}{'exHet-ranges-file'},
    $self->{'special'}{'exHet-binsize'},
    $self->{'special'}{'exHet-regions'},
    $chrom, $res{'start'}, $res{'end'});

  if ( $inWES ) {
    $stat{'inExome'} = 1;
  } else {
    if ($in_exome) {
      $stat{'inExome'} = 2;
    }
  }

  my %total;
  $total{'alleles'} = $self->{'wgsAN'};
  if ($in_exome) {
    $total{'alleles'} += $self->{'wesAN'};
    if ($adjust_for_phase3) {
      #correction b/c 1000G included in ExAC
      $total{'alleles'} -= $self->{'special'}{'1000Genomes-AN'};
    }
  }
  $total{'genotypes'} = int($total{'alleles'}/2); #this isn't quite right.

  # Use position-specific AN
  ### FOR BUILD 160113, disable position-specific AN
  if (0) {
  my $ANcode = $self->posAN($chrom, $res{'start'});
  if ($ANcode) {
    $total{'alleles'} = $self->decodeValue3(substr($ANcode, 0, 3));  #first value is total.
    $total{'genotypes'} = 0;
    my $i = 0;  # first code in ANcode is total AN; skip it
    for my $code ($self->allIdentifiersButRef()) {
      $i+=3;
      my $posAN = $self->decodeValue3(substr($ANcode, $i, 3));
      # don't count 1000Genomes; it is included in ExAC
      if ($code eq $self->{'1000GenomesCode'}) {
	$total{'alleles'} -= $posAN;
      } elsif ($info->{$code}{'Genotypes'}) {
	$total{'genotypes'} += $posAN;
      }
    }
    $total{'genotypes'} = int($total{'genotypes'}/2)
    #print "=>Final total AN is $total{'alleles'}\n";
  } else {
    print  STDERR "No entry in .n file for $chrom $res{'start'}; using max AN of $total{'alleles'}\n" if $self->{'warn'};
  }
  }

  $stat{'platform_specific'} = '';
  if ($total_platform_counts > 100) {
    my $threshold = $total_platform_counts * 0.95;
    for my $platform (keys %platform_counts) {
      if ($platform_counts{$platform} > $threshold) {
	$stat{'platform_specific'} = $platform;
	last;
      }
    }
  }

  my %sortedCalls;
  my %callFreq;
  foreach my $calltype ('alleles', 'genotypes') {
    my $counts_nonref=0;
    my $ref = ($calltype eq 'alleles') ? $refCall : "$refCall/$refCall";
    foreach my $call (keys %{$count{$calltype}}) {
      $counts_nonref += $count{$calltype}{$call} if $call ne $ref;
    }
    # not sure this is correct for genotypes.
    $count{$calltype}{$ref} = $total{$calltype} - $counts_nonref;
    @{$sortedCalls{$calltype}} =
      sort {$count{$calltype}{$b}<=>$count{$calltype}{$a}} keys %{$count{$calltype}};
    foreach my $call (@{$sortedCalls{$calltype}}) {
      next if $call =~ /\./;  #skip genotypes with '.'
      # why is total ever 0?
      $callFreq{$calltype}{$call} = $total{$calltype} ?
         $count{$calltype}{$call}/$total{$calltype} : 0;
      #print "frequency for $call is $count{$calltype}{$call} / $total, or $callFreq{$calltype}{$call}\n";
    }
  }

  ###WARNING: when MAF is 0.5, this essentially reports one of the two equal alleles as minor, entirely at random.
  ###Likewise if there are two equivalent minor alleles, one is reported at random.
  $stat{'alleles'} = $sortedCalls{'alleles'};
  $stat{'alleles_all'} = $stat{'calls'} = \%sortedCalls;
  $stat{'counts'} = $count{'alleles'};
  $stat{'counts_all'} = \%count;
  # development
  #$stat{'id_counts'} = \%id_count;
  $stat{'frequency'} = $callFreq{'alleles'};
  $stat{'frequency_all'} = \%callFreq;
  $stat{'total'} = $total{'alleles'};
  $stat{'total_all'} = \%total;
  $stat{'major'} = $sortedCalls{'alleles'}[0];
  $stat{'all_obs'} = \%all_obs;
  if ($sortedCalls{'alleles'}[1]) {
    my $minor_allele = $sortedCalls{'alleles'}[1];
    $stat{'minor'} = $minor_allele;
    if ($total{'alleles'}) {
      $stat{'maf'} = $count{'alleles'}{$minor_allele}/$total{'alleles'};
    } else {
      $stat{'maf'} = 0;
    }
  }
  return %stat;
}




# digest_lite: does only what's necessary for printVCF 11/20/14
# 06/08/15: modified to adapt to new data structures (see below)
#  following model of what was done with digest() several weeks ago
sub digest_lite {
  my($self, @rest) = @_;
  my %res;
  my ($start, $end, $hashRef) = @{$rest[0]};
  %res = %$hashRef;
  my(%seen, @hist, %count, %id_count, %stat);
  my $refcode = $self->{'referenceCode'};
  my $info = $self->{'info'};
  my $refCall;
  my %alleles;
  my %all_obs;

  # For each allele observed at this position, tally various counts
  foreach my $call (keys %{$res{'sources'}}) {
    # Get a hash of data sources the allele was seen in
    my %obs = %{$self->separateCodes($res{'sources'}{$call}, 1)};
    # For each data source
    # 06/08/15: %obs values are now hrefs, not scalars
    #  also, removed looping over HOM/HET
    for my $code (keys %obs) {
      # If data source is ref genome, store ref call and continue
      if ($code eq $refcode) {
	$stat{'refcall'} = $refCall = $call;
	next;
      }
      my $id =  $self->{'info'}->{$code}->{'id'};
      my $AC = $obs{$code}->{'AC'};
      $AC = 1 if $id =~ /dbSNP/i;
      # skip if AC=0 (e.g. some vars in 2000DanishExomes)
      # This actually messes up when the total AC ends up being zero:
      # the test below, "if ($sortedAlleles[1])", becomes false
      # and maf is not set.
      # We need a more sophisticated approach to handle AC=0
      # in Kaviar db.  08/27/15
#      unless ($AC>0) {
#	delete $obs{$code};
#	next;
#      }
      # skip if AC=0 (e.g. some vars in 2000DanishExomes)
      $seen{$code} += $AC;
      $alleles{$code}{$call} = 1;
      $count{$call} += $AC;
    }
    $all_obs{$call} = \%obs;
  }
  $refCall = 'Reference' if !$refCall;
  my( $inWES, $inSNP);
  foreach my $code (keys %seen) {
    my $id = $info->{$code}{'id'};
    my $genome_coverage = $info->{$code}{'Genome coverage'};
    $inWES = 1 if $genome_coverage =~ /WES/i;
    $inSNP = $code if $id =~ /dbSNP/i;
  }

  if ($inSNP) {
    foreach my $call (keys %{$alleles{$inSNP}}) {
      $count{$call}-- if $count{$call} > 1;
    }
  }

  my $chrom = $res{'chromosome'};
  $chrom = substr($chrom, 3) if substr($chrom, 0, 3) eq 'chr';

#  On bobama this takes about 0.0003s (1s per 3500 calls) if pre-loading used.
#  $stat{'inExHet'} = $self->testInAnnotationRegion(
#    $self->{'special'}{'exHet-ranges-file'},
#    $self->{'special'}{'exHet-binsize'},
#    $self->{'special'}{'exHet-regions'},
#    $chrom, $res{'start'}, $res{'end'});

  if ( $inWES ) {
    $stat{'inExome'} = 1;
  } else {
    if($self->testInExome($chrom, $res{'start'}, $res{'end'})) {
      $stat{'inExome'} = 2;
    }
  }

  ### TEMPORARY: get global AN for Kaviar
  my $total = $self->{'wgsAN'};
  if ($stat{'inExome'}) {
    $total += $self->{'wesAN'};
  }

  # Tally position-specific AN
  if (0) {  ### FOR BUILD 160113, disable position-specific AN
  my $ANcode = $self->posAN($chrom, $res{'start'});
  ### TODO: subtract 1000Genomes
  $total = $self->decodeValue3(substr($ANcode, 0, 3)) if $ANcode;  #first value is total.
  }

  my @sortedAlleles = sort {$count{$b}<=>$count{$a}} keys %count;
  my %alleleFreq;
  foreach my $allele (@sortedAlleles) {
    $alleleFreq{$allele} = $count{$allele}/$total;
  }
  $stat{'alleles'} = \@sortedAlleles;
  $stat{'counts'} = \%count;    # for digest_lite, omits reference
  $stat{'frequency'} = \%alleleFreq;
  $stat{'total'} = $total;
  $stat{'all_obs'} = \%all_obs;
  return %stat;
}

# Not recommended for testing a few sites, as it takes about half a second to load.
# Therefore, it is left to the caller application to call loadExomeSegments() if desired.
sub loadExomeSegments {
  my($self) = @_;
  my $file = $self->{'special'}{'exon-ranges-file'} || ''; 
  my ($binsize, $href) =
    $self->loadAnnotationSegments($file);
  $self->{'special'}{'ExAC-binsize'} = $binsize;;  
  return $self->{'special'}{'ExAC-regions'} = $href;

}

sub loadAnnotationSegments {
  # Pre-loads annotation segments for efficiency when testing many sites.

  my($self, $file) = @_;
  my %href;
  my $binsize = 10000; #arbitrary internal binsize for indexing ranges
  return unless -s $file > 0;
  #print "loading annotation segments in $file\n";
  open INFH, "gunzip -c $file |" or die "can't open file $file $!";
  while (<INFH>) {
    next if /[\@#]/;
    chomp;
    my($chrom, $start, $end, $rest) = split /\t/;
    foreach my $bin (int($start/$binsize)..int($end/$binsize)) {
      $href{$chrom}[$bin]{$start}{$end} = 1;
    }
  }
  close INFH;
  return $binsize, \%href;
}

sub testInExome {
  my($self, $chrom, $start, $end) = @_;
  my $href = $self->{'special'}{'ExAC-regions'};
  my $binsize = $self->{'special'}{'ExAC-binsize'};  
  my $file = $self->{'special'}{'exon-ranges-file'} || ''; 
  return $self->testInAnnotationRegion(
    $file, $binsize, $href,  $chrom, $start, $end);
}


sub testInAnnotationRegion {
  # Tests whether a given position is located within an annotation region
  # as listed in a file.
  # Uses region hash if pre-loaded. Otherwise, tabixes
  # into the file.

  my ($self, $file, $binsize, $href, $chrom, $start, $end) = @_;

  if (defined $href) {
    my $startbin = int($start/$binsize);
    return unless defined $href->{$chrom};
    return unless defined $href->{$chrom}[$startbin];
    foreach my $rstart (keys %{$href->{$chrom}[$startbin]}) {
      next if $rstart>$end;
      foreach my $rend (keys %{$href->{$chrom}[$startbin]{$rstart}}) {
	return $href->{$chrom}[$startbin]{$rstart}{$rend} if $rend>=$start;
      }
    }
  } else {
    my $in_region=0;
    return unless -s $file && -s "$file.tbi";
    open TBX, "$self->{'tabix'} $file $chrom:$start-$end |" or die "Can't open pipe on $file $chrom:$start-$end $!";
    while (<TBX>) {
      chomp;
      my(undef, $blockStart, $blockEnd, $rest) = split /\t/;
      $in_region = 1 if $blockStart<=$end && $blockEnd>=$start;
    }
    close TBX;
    return $in_region;
  }
  return 0;
}

# saving cool method name idea for later ;)
#sub decoct {}
#sub decant { }

#########################################################################################################
#########################################################################################################
#
#  UTILITY METHODS FOR USE BY KAVIAR CODE
# 
#########################################################################################################
#########################################################################################################
# returns hash of baseone positions
sub processLineToPositions {
	my ($self, $pos_ref, $db, $line, $people_codes_ref, $req_var) = @_;
	my ($chrom, $s_baseone, $e_baseone, $rsids, @vars);
	if ($db eq 'dbr') {
		($chrom, $s_baseone, $e_baseone, $rsids, @vars) = split /\t/, $line;
	} elsif ($db eq 'dbs') {
		($chrom, $s_baseone, $rsids, @vars) = split(/\t/, $line);
		$e_baseone = $s_baseone;
	} else {
		confess "Invalid database type supplied $db";
	}
	my $decoded = $self->{'decoded'};
	my %allRsids;
	while ($rsids) {
		my $rsid = 
			$decoded->[ord(substr($rsids,0,1))]*15625000 +
			$decoded->[ord(substr($rsids,1,1))]*62500 +
			$decoded->[ord(substr($rsids,2,1))]*250 +
			$decoded->[ord(substr($rsids,3,1))];
			$allRsids{"rs$rsid"}++;
			$rsids = substr($rsids, 4);
	}
	my @people;
	@people = @$people_codes_ref if $people_codes_ref && scalar @$people_codes_ref > 0;
    my $new_ref = $pos_ref;   
	my %byPos = %$new_ref;
	# in order to filter variants, we need to save the variant before it's been split up
	my %found_vars;
	while (@vars) {
		my $var = shift @vars;
		if ($var eq '') { $var = '-'; }
		my $sources = shift @vars;
		my $found_var = $var;
		# if people is defined, cull sources for those codes
		if (scalar @people > 0) {
			$sources = $self->cullCodes($sources, $people_codes_ref);
			next unless ($sources);
		}

		my @bps = split(/|/, $var);
		my $l_var = length($var);  # if var was empty string and is now '-', $l_var will be 1. Problem?
		if ($l_var > $self->{'maxlength'}) {
			print "WARNING: Skipping this variant $var at position $s_baseone to $e_baseone becasue it is larger than the maxlength setting of $self->{'maxlength'}\n" if $self->{'warn'};
			next;
		}
		my $l_ref = ($e_baseone - $s_baseone)+1;
		my $insertion_bp = $e_baseone;
		my $insertion_flag = 0;
		my $end = $e_baseone;
		if ($l_var > $l_ref && $l_ref <= 1) {
			my $refCall = $self->referenceCall({ sources=> { $var => $sources}});
			if ($l_ref <= 1 && $refCall eq $var) {
				print "Warning: Skipping this variant $var at position $s_baseone to $e_baseone because it is the complement to a single base pair insertion.\n" if $self->{'warn'};
				next;
			}
			if ($self->{'exclude_insertions'}) {
				print "Warning: Skipping insertion variant $var at position $e_baseone because exclude_insertions is on\n" if $self->{'warn'};
				next;
			}
			$insertion_flag = 1;
			$self->{'end_overlap'} = 1;
			$end = $e_baseone + $l_var - 1;
		}
		for (my $i = $s_baseone; $i <= $end; $i++) {
			# add the variant
			my $bp = $bps[$i-$s_baseone];
			$bp = '-' if (!$bp); # it is a deletion if there isn't a bp available - this makes substitutions left aligned
			my $curr_pos = $i;
			# this should be > and not >= because insertion_bp is set considering the reference
			if ($insertion_flag) {
				if ($curr_pos > $insertion_bp) {
					print "ERROR!!!!!! $bp shouldn't be an insertion!!!!\n" if ($bp eq '-');
					$bp .= '_INSERTION';
				} else {
					next;
				}
			}
			my %single_result;
			%single_result = %{$byPos{$curr_pos}{$curr_pos}} if $byPos{$curr_pos} && $byPos{$curr_pos}{$curr_pos};
			$single_result{'chromosome'} = 'chr'.$chrom;
			$single_result{'start'} = $curr_pos;
			$single_result{'end'} = $curr_pos;
			# need to add rsids to old ones
			if (scalar keys %allRsids > 0 ) {
				if (scalar keys %{$single_result{'rsids'}} > 0 ) {
					my %old_rsids = %{$single_result{'rsids'}};
					@old_rsids{keys %allRsids } = values %allRsids;
				$single_result{'rsids'} = \%old_rsids;
				} else {
					$single_result{'rsids'} = \%allRsids;
				}
			}
			if ($single_result{'sources'}{$bp}) {
				my $res = $self->checkAddCodes([$sources], $single_result{'sources'}{$bp});
				$single_result{'sources'}{$bp} = $res;
			} else {
				$single_result{'sources'}{$bp} = $sources;
			}

			$byPos{$curr_pos}{$curr_pos} = \%single_result;
		}
	} #end while vars

	if ($req_var) {
		if ($found_vars{$req_var}) {
			return \%byPos;
		} else {
			print "Did not find requested variant ($req_var) in found variants (". join(", ", keys %found_vars)."). Position ".($s_baseone-1)." filtered.\n" if $self->{'warn'};
			return $pos_ref;
		}
	} else {
		return \%byPos;
	}
} 

# take a hash of positions and filter by the kaviar object settings for start_overlap, start_match,
# end_overlap, end_match, and return_reference.  Return a hash or an array depending on the 
# result_type requested
sub filterResults {
	my ($self, $result_type, $start_basezero, $end_basezero, $resRef, $people_codes_ref) = @_;
	my %res = %$resRef;
	my %filtered_results;
	my @result;
	my @keys = keys %res;
	my @sorted_keys = sort {$a <=> $b} @keys;
	# correct pos from baseone to basezero to return results
	my $start = $start_basezero+1;
	$start = $sorted_keys[0] if $self->{'start_overlap'} || $self->{'start_match'};
	my $end = $end_basezero+1;
	# this searches through a lot of extra basepairs depending on how big the window is
	# could be improved?
	for (my $i = $start; $i <= $end; $i++) {
		my $pos = $i;
		if (!$res{$pos}) {
			next;
		}

		my $s = $res{$pos}{$pos}{'start'};
		my $s_basezero = $s -1;
		my $e = $res{$pos}{$pos}{'end'};	
		my $e_basezero = $e-1;
		my %localPos = %{$res{$pos}{$pos}};
		#print "localPos before ".Dumper \%localPos;
		$localPos{'start'} = $s_basezero;
		$localPos{'end'} = $e_basezero;
		#print "localPos after ".Dumper \%localPos;
		if ($localPos{'sources'}{'-'}) {
			$localPos{'sources'}{'DELETION'} = $localPos{'sources'}{'-'};
			delete $localPos{'sources'}{'-'};
		}

		# change the rsids from a hash to an array to be backwards compatible
		my @rsids = keys %{$localPos{'rsids'}};
		$localPos{'rsids'} = \@rsids;

		# check to see if return_reference is set, otherwise skip positions that have only one
		# allele and that allele includes reference in the heterozygous codes, unless people 
		# filtering is set, in which case there may be only one person requested and therefore
		# all results will have one allele
		if (scalar keys %{$localPos{'sources'}} <= 1) {
			# if there are people filtering provided then any position could be filtered
			# down to having one allele, so turn this check off in that case
			next if ( $people_codes_ref ); 
			my $test = $self->referenceCall(\%localPos);
			# test will return a ? if there is no reference in the codes
			# test will return a basepair if there is a reference in the codes
			if ($test ne '?' && !$self->{'return_reference'}) {
				print "Filtered out position $s_basezero because it had one allele ($test) that was the reference allele.\n" if $self->{'warn'};
				next;
			}
		}

		# is it valid to exclude part of a variant?
		# for example, if I query chr2 38217 38219 I get a variant:
		# 2  38216  38219  AAAA   AAG
		# Is it valid to exclude the basepair at 38216? (well one can argue
		# here that the AA shouldn't be present in the variant because it's
		# reference anyway, but you get the point).
		if ($self->{'start_overlap'} == 1 && $self->{'end_overlap'} == 1) {
			if ($result_type eq 'HASH') {
				$filtered_results{$s_basezero}{$e_basezero} = \%localPos;
			} else {
				push @result, [$s_basezero, $e_basezero, \%localPos];
			}
		} elsif ($self->{'start_overlap'} == 1) {
			if ($e_basezero <= $end_basezero) {
				if ($result_type eq 'HASH') {
					$filtered_results{$s_basezero}{$e_basezero} = \%localPos;
				} else {
					push @result, [$s_basezero, $e_basezero, \%localPos];
				}
			}
		} elsif ($self->{'end_overlap'} == 1) {
			if ($s_basezero >= $start_basezero) {

				if ($result_type eq 'HASH') {
					$filtered_results{$s_basezero}{$e_basezero} = \%localPos;
				} else {
					push @result, [$s_basezero, $e_basezero, \%localPos];
				}
			}
		} else {  # start_overlap = 0, end_overlap = 0
			if ($s_basezero >= $start_basezero && $e_basezero <= $end_basezero) {
				if ($result_type eq 'HASH') {
					$filtered_results{$s_basezero}{$e_basezero} = \%localPos;
				} else {
					push @result, [$s_basezero, $e_basezero, \%localPos];
				}
			}
		}
	}

	my $result;
	if ($result_type eq 'HASH') {
		$result = \%filtered_results;
	} else {
		$result = \@result;
	}
	return $result;
}

sub allIdentifiers {
	my($self) = @_;
	return @{$self->{'allIdentifiers'}};
}

sub allIdentifiersButRef {
	my($self) = @_;
	return @{$self->{'allIdentifiersButRef'}};
}

sub checkAddCodes {
	my ($self, $cur_ref, $new_ref) = @_;

	$cur_ref = $self->combineCodes($cur_ref, $new_ref);
	my @com = @$cur_ref;
	my $codes = join('', @com) || undef;

	return($codes);
}

sub combineCodes {
	my ($self, $cur_single, $new_single) = @_;

	my @cur_codes = @{$self->separateCodes($cur_single)};
	my %combined = map { $_ => 1 } @cur_codes if scalar @cur_codes > 0;
	my @new_codes = @{$self->separateCodes($new_single)};
	map { $combined{$_} = 1 } @new_codes if scalar @new_codes > 0;
	my @com;
	foreach my $key (sort keys %combined) {
		my $an = $self->{'info'}->{$key}{'AN'};
		# reference has a value of 1, need to skip that because it's not a 'true' population
		# according to how Kaviar is coded even though it has a value for AN
		if ($an && $an > 1) {
			my $encodedAn = $self->encodeValue3($an);
			push(@com, $key.$encodedAn);
		} else {
			push(@com, $key);
		}
	}

	return \@com;
}

# Decode a value encoded by encodeValue3()
sub decodeValue3 {
	my($self, $value) = @_;
	
	my $decoded = $self->{'decoded'};
	return $decoded->[ord(substr($value,0,1))]*62500+
		$decoded->[ord(substr($value,1,1))]*250+
		$decoded->[ord(substr($value,2,1))];
}

# Given an integer <= (250**3), return
# a three-byte code for the integer. For encoding AC (allele count).
sub encodeValue3 {
        my($self,$value) = @_;

        my $tmp = $value;
	my @encoded = (1..8, 12, 14..47, 49..255);
	my @encodedChr = map chr($encoded[$_]), (0..249);
        my $low = $encodedChr[$value % 250];
        $value = int($value/250);
	my $middle = $encodedChr[$value % 250];
        $value = int($value/250);
        if ($value>250) {
                print STDERR "FATAL encoding value: $tmp too large for 3 bytes\n";
                die "$tmp too large for three bytes\n";
        }
        return join("", $encodedChr[$value], $middle, $low);
}

# Decode a value encoded by encodeValue2(). Deprecated.
sub decodeValue2 {
	my($self, $value) = @_;
	
	my $decoded = $self->{'decoded'};
	return $decoded->[ord(substr($value,0,1))]*250+
		$decoded->[ord(substr($value,1,1))];
}

# Given an integer <= (250**2), return
# a two-byte code for the integer. For encoding AC (allele count).
# Deprecated.
sub encodeValue2 {
        my($self,$value) = @_;

        my $tmp = $value;
	my @encoded = (1..8, 12, 14..47, 49..255);
	my @encodedChr = map chr($encoded[$_]), (0..249);
        my $low = $encodedChr[$value % 250];
        $value = int($value/250);
        if ($value>250) {
                print STDERR "FATAL encoding value: $tmp too large for 2 bytes\n";
                die "$tmp too large for two bytes\n";
        }
        return join("", $encodedChr[$value], $low);
}

sub getCodes {
	my ($self, $people_ref) = @_;
	die "Incorrect data type supplied to getCodes.  Should be an array reference." unless ref($people_ref) eq 'ARRAY';
	my %names;
	if ($self->{'names'}) {
		%names =  %{$self->{'names'}};
	} else {
		my %info = %{$self->{'info'}};
		foreach my $code (keys %info) {
			$names{$info{$code}{'name'}} = $code if ($info{$code} && $info{$code}{'name'});
			$names{$info{$code}{'id'}} = $code;
		}
		$self->{'names'} = \%names;
	}

	my @codes;
	foreach my $person (@$people_ref) {
		if ($names{$person}) {
			push(@codes, $names{$person});
		} else {
			print "Could not find a code for $person as name or id.  Check that you are spelling the name or identifier correctly.\n";
		}
	}

	return \@codes;
}

# take in a string of codes and an array of filter codes and return only the codes that match the filter
sub cullCodes {
	my ($self, $found_ref, $filter_ref) = @_;
	my @found_codes = @{$self->separateCodes($found_ref)};
	my %filter_codes;
	map { $filter_codes{$_} = 1; } @$filter_ref;
    # add ref code to filter code so it doesn't get filtered out
    my $ref_code = $self->{'names'}{'reference'};
    $filter_codes{$ref_code} = 1;
	my @keep_codes;
	foreach my $code (@found_codes) {
		push(@keep_codes, $code) if $filter_codes{$code};
	}

    if (scalar @keep_codes > 0) {
        my $keep_string = join("", sort @keep_codes);
        return $keep_string;
    } else {
        return "";
    }
#	my @keep_encoded;
#	foreach my $key (sort @keep_codes) {
#		my $an = $self->{'info'}->{$key}{'AN'};
#        print "key $key an $an\n";
#		# reference has a value of 1, need to skip that because it's not a 'true' population
#		# according to how Kaviar is coded even though it has a value for AN
#		if ($an && $an > 1) {
#			my $encodedAn = $self->encodeValue2($an);
#			push(@keep_encoded, $key.$encodedAn);
#		} else {
#			push(@keep_encoded, $key);
#		}
#	}
#    print "keep encoded ".Dumper \@keep_encoded;
#	if (scalar @keep_encoded > 0){
#		my $reference_code = $self->{'names'}{'reference'};
#		push(@keep_encoded, $reference_code);
#		my $keep_string = join("", sort @keep_encoded);
#        print "keep string $keep_string\n";
#		return $keep_string;
#	} else {
#		return "";
#	}

}

sub referenceCall {
	
	my($self, $res) = @_;
	my $refcode = $self->{'referenceCode'};

	foreach my $call (keys %{$res->{'sources'}}) {
		# check the codes for the reference
		my $check_codes = $res->{'sources'}{$call}; 
		foreach my $code (@{$self->separateCodes($check_codes)}) {
			return $call if $code eq $refcode;
		}
	}
	return "?";
}

sub heterozygosity {
	
	### need to revisit
	
	my($self, %res) = @_;
	my(%seen, @hist, $total);
	my @alph = @{$self->{'alphabet'}};
	#my $refcode = $self->{'referenceCode'};
	foreach my $call (@alph) {
		foreach my $code (@{$self->separateCodes($res{'sources'}{$call})}) {
			$seen{$code}{$call} = 1;
		}
	}
	foreach my $code (keys %seen) {
		$hist[scalar keys %{$seen{$code}}]++;
		$total++;
	}
	$hist[0] = $total ? sprintf("%.4f", $hist[2]/$total) : -1;
	return \@hist;
}

sub cleanupChrom {
	my($self, $chrom) = @_;
	
	confess "No chromosome provided to cleanupChrom" unless $chrom;
	if ($chrom =~ /^chr/) {
		$chrom =~ s/^chromosome(?=.)//i;
		$chrom =~ s/^chrom(?=.)//i;
		$chrom =~ s/^chr(?=.)//i;
	}
	$chrom = uc $chrom;
	$chrom = "M" if $chrom eq "MT";
	return $chrom;
}

sub readProjects {
	my($self, $file) = @_;
	my(%pinfo);
	open PROJ, $file || return;
	$_ = <PROJ>;
	chomp;
	my @fields = split /\t/;
	while (<PROJ>) {
		next if /^#/;
		chomp;
		my(@values) = split /\t/;
		my $code = $values[0];
		foreach my $i (1..$#fields) {
			$pinfo{$code}{$fields[$i]} = $values[$i];
		}
	}
	close PROJ;
	$self->{'projectinfo'} = \%pinfo;
	return 1;
}

sub DESTROY {
	#nothing to do here anymore
}

sub variantsPerSource {
  my ($self, $verbose, $exclude_ISB, $exclude_Inova, $test) = @_;
  my $totalVariants=0;
  my %variantsPerSource;
  my %variantsPerProject;
  my $refcode = $self->{'referenceCode'};
  my $dbsnpcode = $self->{'dbsnpCode'};
  my $ISBcode = $self->{'ISBcode'} || '';
  my $Inovacode = $self->{'InovaCode'} || '';

  my @all_projects = keys %{$self->{'projectinfo'}};

  my %has;
  my %hasnt;

  #foreach my $db (qw/dbs dbr/) {
  foreach my $db (qw/dbs/) {          # only process SNVs for now
    my $file = $self->{$db}{'file'} . ".gz";
    open(INFH, "zcat $file|") || die "Can't pipe $file\n";
    print "Processing $db " if $verbose;
    my $i = 0;
    while (<INFH>) {
      next if /^#/;
      $i++;
      if ($verbose) {
	print "." unless $i % 100000;
	print "$i" unless $i % 1000000;
      }
      last if ($test  && ($i == 10000));
      chomp;
      my ($chrom, $start, $end, $rsids, @v);
      if ($db eq 'dbs') {
	($chrom, $start, $rsids, @v) = split("\t", $_);
      } elsif ($db eq 'dbr') {
	($chrom, $start, $end, $rsids, @v) = split("\t", $_);
      }
      #my $n_v = int((scalar @v - 2) / 3) + 1; @ what the hell?
      my $n_v = (scalar @v) / 2;
      my $uniqueSources="";
      my %complete_codehash=();
      my %projecthash = ();

      # For each allele described on this line, including reference
      for my $j (0 .. $n_v-1) {
	my $var = shift @v;
	my $codes = shift @v;
	my $codelist_href = $self->separateCodes($codes, 1);
	my %AC;

	my $this_var_is_ref=0;

	# Make a hash of codes and AC values
	# for the data sources this allele is seen in
	for my $code ( keys %{$codelist_href}) {
	  $this_var_is_ref =1 if $code eq $refcode;
	  last if $code eq $refcode;
	  next if ($code eq $ISBcode) && ($exclude_ISB);
	  next if ($code eq $Inovacode) && ($exclude_Inova);
	  $AC{$code} = $codelist_href->{$code}->{'AC'};
	}

	# If this allele is found in the reference,
	# skip it.
	next if $this_var_is_ref;

	my @codelist = keys %AC;

	# DEBUG
	##my $sources = $self->whois_list_to_string(@codelist);
	##print "$var $sources\n";

	# Add the sources and projects this
	# variant appears in to a hash of variant sources at this position
	for my $code (@codelist) {
	  $complete_codehash{$code}=$AC{$code};
	  my $project = $self->{'info'}->{$code}{'project'};
	  if (defined $projecthash{$project}) {
	    $projecthash{$project} += $AC{$code};
	  } else {
	    $projecthash{$project} = $AC{$code};
	  }
	}

      }
      my @complete_codelist = keys %complete_codehash;
      my @projectlist = keys %projecthash;

      # Increment the tallies of variants per source/project
      for my $code (@complete_codelist) {
	$variantsPerSource{$code}{'all'}++;
      }
      for my $project (@projectlist) {
	$variantsPerProject{$project}{'all'}++;
      }

      # If, by removing some data sources, we have no datasources
      # left for this position, go to the next one.
      next unless scalar @complete_codelist;

      $totalVariants++;


      # If list contains only dbSNP and one other source,
      # count this variant as unique to the other source so as not
      # to penalize the source for (presumably) contributing to dbSNP
      # Code below could be more elegant/efficient
#      if (((scalar @complete_codelist) == 2) &&
#	  (($complete_codelist[0] eq $dbsnpcode) ||
#	   ($complete_codelist[1] eq $dbsnpcode))) {
#	 my $other_code = $complete_codelist[0] eq $dbsnpcode ?
#	    $complete_codelist[1] : $complete_codelist[0];
#	 @complete_codelist = ($other_code);
#      }
#      # do the same for the project list ...
#      if (((scalar @projectlist) == 2) &&
#	  (($projectlist[0] =~ /dbsnp/i) ||
#	   ($projectlist[1] =~ /dbsnp/i))) {
#	 my $other_code = $projectlist[0] =~ /dbsnp/i ?
#	    $projectlist[1] : $projectlist[0];
#	 @projectlist = ($other_code);
#      }

      # Compile tallies for each project of how often it shares
      # a variant position with N total other projects, and, conversely,
      # how often it lacks a variant position that is present in N
      # total other projects.
      my $has = scalar @projectlist - 1;  # num projects that have variation here
                                          # subtract self
      for my $project (@all_projects) {
	if ($projecthash{$project}) {
	  $has{$project}[$has]++;
	} else {
	  $hasnt{$project}[$has]++;
	}
      }

      # If there is only one datasource that has a nonreference allele
      # at this position, count it as unique. If, also, AC==1, count as
      # singleton.
      if (scalar @complete_codelist == 1) {
	$variantsPerSource{$complete_codelist[0]}{'unique'}++;
	$variantsPerSource{$complete_codelist[0]}{'singleton'}++
	  if ($complete_codehash{$complete_codelist[0]} == 1);
      }
      if (scalar @projectlist == 1) {
	$variantsPerProject{$projectlist[0]}{'unique'}++;
	$variantsPerProject{$projectlist[0]}{'singleton'}++
	  if ($projecthash{$projectlist[0]} == 1);
      }



      #DEBUG
#      my $sources = $self->whois_list_to_string(@complete_codelist);
#      my $projects = join(",",@projectlist);

#      print "$chrom\t$start\t$sources\t$projects\n";
    }
    print "\n" if $verbose;
  }
  return \%variantsPerSource, \%variantsPerProject,
	 \%has, \%hasnt, $totalVariants;
}

sub whois_list_to_string{
	my($self, @codes) = @_;
	my @res;
	foreach my $code (@codes) {
		my $info = $self->{'info'}->{$code};
		my $name = $info->{'name'};
		my $id = $info->{'id'};
		#$name .= " ($id)" if defined $name && $name =~ /^anonymous/i;
		$name = '' if defined $name && $name =~ /^anonymous/i;
		push @res, $name || $id;
	}
	return join(",",@res);
}


# Get a specific base of a reference sequence using tabixed files
# created for gestalt. One-based. Use tabix to get an entire chrom
# at once. If used to get only one bin at a time, becomes slow for
# later bins. Return nothing if unsuccessful.
# Argument $chrom should be of form 1 2 ... 23 X Y M
sub refbase {
  my($self, $chrom, $pos, $len) = @_;
  $len = 1 unless defined $len;
  my $binsize = 10000;

  unless ($self->{'refbase_init'}) {
    $self->{'refbin'} = 0;
    $self->{'refchrom'} = '';
    $self->{'refbase_init'} = 1;
  }

  my $tbxfh;
  my $tabix_file = $self->{'ref_tabix_dir'} . "/" . $self->{'frz'} . ".gz";

  # Fetch new chromosome if needed
  my $bin = int(($pos-1)/$binsize+1e-6)+1;
  unless (($chrom eq $self->{'refchrom'}) && ($bin >= $self->{'refbin'})) {
    my $gotit = open $tbxfh, "tabix $tabix_file chr$chrom |";
    unless ($gotit) {
      print "can't open tabix on$tabix_file chr$chrom\n" if $self->{'warn'};
      return;
    }
    $self->{'refchrom'} = $chrom;
    $self->{'refbin'} = 0;
  } else {
    $tbxfh = $self->{'refbase_tbxfh'};
  }

  # Fetch new 10kb chunk of reference if needed
  while ($self->{'refbin'} < $bin) {
    my $line = <$tbxfh>;
    unless ($line) {
      print "couldn't read anything from tabix filehandle\n";
      return;
    }
    chomp $line;
    (undef, $self->{'refbin'}, $self->{'refchunk'}) = split /\t/, $line;
  }
  if ($self->{'refbin'} != $bin) {
    print "Couldn't get bin $bin for chr$chrom\n" if $self->{'warn'};
    return;
  }

  $self->{'refbase_tbxfh'} = $tbxfh;

  # Extract and return desired base
  return uc(substr($self->{'refchunk'}, ($pos-1)%$binsize, 1));
}

# Input: aref returned from allVariants()
# Output: VCF records (no headers)
sub printVCF {
  my ($self, $vars_aref, $bin_start, $bin_end) = @_;
  $bin_start = 0 unless $bin_start;
  my %out_lines = ();
  for my $var_aref (@$vars_aref) {
    my ($start, $end, $info_href) = @$var_aref;  # results are base 0
    my $chrom = $self->cleanupChrom($info_href->{'chromosome'});
    my $need_padding = 0;
    # Skip vars that begin before desired cutoff
    next if $start < $bin_start;
    my %stat = $self->digest_lite($var_aref);
    my $ref = $stat{'refcall'};
    $ref = ''  unless defined $ref;
    $need_padding = 1 unless $ref;
    my $freq_href = $stat{'frequency'};
    my $ac_href = $stat{'counts'};
    my (@af, @ac, @alt, @sources, @rsids);
    my $an=$stat{'total'};;
    for my $var (sort @{$stat{'alleles'}}) {
      next if $var eq $ref;
      next if $var eq 'Reference';
      my $printable_var = $var;
      $printable_var = '' if $printable_var eq 'DELETION';
      $printable_var =~ s/_INSERTION//;
      $need_padding = 1 unless $printable_var;
      push @alt, $printable_var;
      push @af, sprintf "%0.7f",$freq_href->{$var};
      push @ac, $ac_href->{$var};
      # Get sources for this var.
      my $codes =  $info_href->{'sources'}->{$var};
      my $rsids_aref = $info_href->{'rsids'};
      push @rsids, @{$rsids_aref} if defined $rsids_aref;;
      my @var_sources = ($self->whois($codes));
      # remove dbSNP from sources
      my %tmp_hash = map { $_ => 1 } @var_sources;
      delete $tmp_hash{'dbSNP'};
      @var_sources = sort keys %tmp_hash;
      push @sources, join "|", @var_sources;
    }

    # Convert start, end to base 1
    my $start_base1 = $start+1;
    my $end_base1 = $end+1;

    # Add a single reference base to ref & var
    # for positions that include indels. VCF spec requires this.
    # Adjust start/end accordingly.
    my $pad;
    if ($need_padding) {
      if ($start_base1 == 1) { #only occurs for chrM; pad at end
	$pad = $self->refbase($chrom,$end_base1+1);
	my @tmp = ();
	for my $alt (@alt) {
	  push @tmp, $alt . $pad;
	}
	@alt = @tmp;
	$ref = $ref . $pad;
	$end_base1++;
	$end++;
      } else { # usual case; pad at start
	$pad = $self->refbase($chrom,$start_base1-1);
	my @tmp = ();
	for my $alt (@alt) {
	  push @tmp, $pad . $alt;
	}
	@alt = @tmp;
	$ref = $pad . $ref;
	$start_base1--;
	$start--;
      }
    }

    # Prepare items for printing
    my $alt = join ",", @alt;
    my $sources = join ",", @sources;
    my %tmp_hash = map { $_ => 1} @rsids; # remove duplicates via hash
    @rsids = sort keys %tmp_hash;
    my $rsids = (join ";", @rsids) || ".";
    my $af = join ",", @af;
    my $ac = join ",", @ac;

    my $info = "AF=$af;AC=$ac;AN=$an";
    $info .= ";END=$end_base1" if $end != $start;
    $info .= ";DS=$sources" unless $sources =~ /^\,*$/;   # $sources is empty if var seen only in dbSNP
    my $qual = '.';
    my $filter = '.';

    my $line = join ("\t", $chrom, $start_base1, $rsids, $ref, $alt, $qual, $filter, $info) . "\n";
    push @{$out_lines{$start_base1}}, $line;
  }

  # The sequential ordering of variants is disturbed by the
  # padding of indels with reference bases, which causes those
  # variant positions to be decremented. So sort output lines by position.
  # If $bin_end provided, skip vars outside of bin.
  my $out = '';
  if ($bin_end) {
    for (my $pos = $bin_start; $pos <= $bin_end; $pos++) {
      for my $line (@{$out_lines{$pos}}) {
	$out .= $line;
      }
    }
  } else {
    for my $pos (sort {$a <=> $b} keys %out_lines) {
      for my $line (@{$out_lines{$pos}}) {
	$out .= $line;
      }
    }
  }

#  if ($debug) {
#    my $last_epoc = $epoc;
#    $epoc = time();
#    my $interval = $epoc - $last_epoc;
#    my $nvars = scalar @$vars_aref;
#    my $speed = int($nvars/($interval + 0.5));
#    print STDERR "chr$chrom bin $bin_start $nvars vars ${interval}s ($speed/s)\n";
#  }

  print $out;
}

# Input: aref returned from allGenotypes()
# Output: for each variant, HOM/HET frequency and overall variant frequency
sub printGTfreq {
  # Sources for which each genotype is explicitly called for every sample
  # 06/10/15: removed Wellderly; it seems to lack homozygous reference calls
  #  (or, perhaps, any homozygous calls at all!)
#    phase3-AHD #missing
#    phase3-AJM #missing
#    phase3-CHD #missing
#    phase3-GHN #missing
#    phase3-MAB #missing
#    phase3-MKK #missing
#    phase3-MRM #missing
#    phase3-RDH #missing

  my $min_samples = 1;

  my %complete_sources = map {$_ => 1} qw (
    phase3-ITU
    phase3-ESN
    phase3-STU
    phase3-BEB
    phase3-MSL
    phase3-ACB
    phase3-AHD
    phase3-AJM
    phase3-ASW
    phase3-CDX
    phase3-CEU
    phase3-CHB
    phase3-CHD
    phase3-CHS
    phase3-CLM
    phase3-FIN
    phase3-GBR
    phase3-GHN
    phase3-GIH
    phase3-GWD
    phase3-IBS
    phase3-JPT
    phase3-KHV
    phase3-LWK
    phase3-MAB
    phase3-MKK
    phase3-MRM
    phase3-MXL
    phase3-PEL
    phase3-PJL
    phase3-PUR
    phase3-RDH
    phase3-TSI
    phase3-YRI
    SSIP
    Malay
    ADNI
  );

  my %phase3_sources = map {$_ => 1} qw (
    phase3-ITU
    phase3-ESN
    phase3-STU
    phase3-BEB
    phase3-MSL
    phase3-ACB
    phase3-AHD
    phase3-AJM
    phase3-ASW
    phase3-CDX
    phase3-CEU
    phase3-CHB
    phase3-CHD
    phase3-CHS
    phase3-CLM
    phase3-FIN
    phase3-GBR
    phase3-GHN
    phase3-GIH
    phase3-GWD
    phase3-IBS
    phase3-JPT
    phase3-KHV
    phase3-LWK
    phase3-MAB
    phase3-MKK
    phase3-MRM
    phase3-MXL
    phase3-PEL
    phase3-PJL
    phase3-PUR
    phase3-RDH
    phase3-TSI
    phase3-YRI
  );

  my ($self, $gts_aref, $bin_start) = @_;
  $bin_start = 0 unless $bin_start;

  print join "\t", qw (chrom start %het frac_het sum %p3_het frac_p3_het p3_sum gts);
  print "\n";

  for my $gt_aref (@$gts_aref) {
    my ($start, $end, $info_href) = @$gt_aref;  # results are base 0
    my $chrom = $self->cleanupChrom($info_href->{'chromosome'});
    # Filter gts that begin before desired cutoff
    next if $start < $bin_start;

    my $rsids = defined $info_href->{'rsids'} ? join ",", @{$info_href->{'rsids'}} : '';
    my @gts;
    for my $gt (keys %{$info_href->{'sources'}}) {
      push @gts, $gt;
    }

    my %totals;
    my %phase3_totals;
    my %ids;
    for my $gt (@gts) {
      if ($gt eq '') { $gt = '-'; }
      my $sources = $info_href->{'sources'}->{$gt};
      my %obs = %{$self->separateCodes($sources, 1)};
      for my $code (keys %obs) {
	my $id =  $self->{'info'}->{$code}->{'id'};
	if (defined $complete_sources{$id}) {
	  my $hash_id = ($id =~ /phase3/) ? '1000G' : $id;
	  $ids{$hash_id} = 1;
	  if (ref $obs{$code} eq 'HASH') {
	    my $AC = $obs{$code}->{'AC'};
	    $totals{$gt} += $AC;
	    #my $AN = $obs{$code}->{'AN'};
	    $phase3_totals{$gt} += $AC  if defined $phase3_sources{$id};
	  }
	}
      }
    }

    my $sum = 0;
    my $phase3_sum = 0;
    my $nhet = 0;
    my $phase3_nhet = 0;
    for my $gt (@gts) {
      next unless (length $gt == 3 && defined $totals{$gt});
      my ($a, $b) = $gt =~ /(.).(.)/;
      $nhet += $totals{$gt} if ($a ne $b);
      $sum += $totals{$gt};
      next unless defined $phase3_totals{$gt};
      $phase3_nhet += $phase3_totals{$gt} if ($a ne $b);
      $phase3_sum += $phase3_totals{$gt};
    }

    if ($sum >= $min_samples) {
      my $fhet = sprintf "%0.5f", $nhet/$sum;
      my $percent_het = int($fhet*100);

      my $phase3_fhet = $phase3_sum ? sprintf "%0.5f", $phase3_nhet/$phase3_sum : 0;
      my $phase3_percent_het = int($phase3_fhet*100);

      # Convert start to base 1. (end==start as of June 2015)
      my $start_base1 = $start+1;

      my $ids = join ",", sort keys %ids;
      my $line = join ("\t", $chrom, $start_base1, $percent_het, $fhet, $sum, $phase3_percent_het, $phase3_fhet, $phase3_sum, $ids, %totals) . "\n";
      print $line;
    }
  }
}

sub get_rsid_pos {
  my ($self, $rsid) = @_;
  my ($id) = ($rsid =~ /^rs(\d+)$/);
  return '' unless $id;
  my $idx = substr($id, 0, 2);
  my $tbx = $self->{'rsid_index'}{'tbx'};
  my $res = $tbx->query($idx, $id-1, $id);
  return '' unless $res;
  my $line = $tbx->read($res);
  return '' unless $line;
  my (undef, undef, $chrom, $pos) = split "\t", $line;
  return ($chrom, $pos);
}

1;
