use strict;
use lib "$ENV{HOME}/perl/";
#   use SqlUtils;

# db connect
#my $DBH = SqlUtils::getDBH();

sub decodeCigar {
    
    my $cig = shift;
    my @list = split (/\d+/, $cig);
    shift(@list);
    my @list2 = split (/[MIDNSHP]/, $cig);
    my %hash = ();
    for ( my $i=0; $i<scalar(@list); $i++ ) {
        $hash{$list[$i]} += $list2[$i]
    }
    return $hash{'M'}; # the number of ungapped aligned matched & mismatched bases
}

sub getMismatch {
    my $info = shift;
   
    if ($info =~ /XM\:i\:(\d+)/) {
            return $1; # the number of mismatched bases
    }
    elsif ($info =~ /NM\:i\:(\d+)/) {
            return $1;
    }
    
    return -1;
}

#sub queryContig {
#    # this functions gets the contig IDs of the given BACT_ID
#    my ($contig_id) = @_; # the BACT_ID
#    my $sql = qq {
#        SELECT bact_id
#        FROM refContig
#        WHERE contig_id="$contig_id"
#    };
#    my $results = SqlUtils::query($DBH, $sql);
#    
#    return $results->[0][0]; # returns the bact ID
#}

my $usage = qq {
    $0 <sam> <minimum seq identity> <minimum M> <minimum A> <output option> <preserve header?>
    <sam>: a sam file (BWA v5.9 standard) or a bam file (with a suffix .sam or .bam)
    <minimum identity>: the minimum sequence identity (matched/read length)
    <minimum M>: the minimum percentage of M/L (M is total bases in ungapped alignment, L is the read length)
    <minimum A>: the percent of (Matched bases/M), this cutoff may be more appropriate for longer reads than <minimum identity>
    
    <output option>:
    1, print filtered SAM.
    2, print seq identity and alignment info. pileup2coverage.pl depends on this format.
    3, print both: sam to STDOUT and seq info to STDERR
    4, customized tab-delimited: BACT_ID, position, read length, seq ident, seq ident in aligned region, aligned length, mismatches. 
    
    <preserve header?>: 1, yes. 0, no.
    
    This script filters the sam file based on aligned read's sequence identity to the reference genome and generates a new sam.
    
    This script can also handle a BAM file.  It recognizes the file format based on the file's suffix: .sam or .bam
    
};

die $usage unless @ARGV==6;

my $sam = $ARGV[0];
my $ident = $ARGV[1];
my $M = $ARGV[2];
my $A = $ARGV[3];
my $option = $ARGV[4];
my $printHeader = $ARGV[5];

my $count = 0;
my %contig2bact = (); # a hash table to map contigs IDs to BACT IDs
if ($sam =~ /.sam$/){
    open SAM, "$sam" or die "Can't open $sam:$!\n";
}
else { # assuming that the input file is in BAM format
    if ($printHeader==1){
        open SAM, "samtools view -h $sam |" or die "Can't read $sam. Perhaps forgot to use Samtools?";
    }
    else {
        open SAM, "samtools view -h $sam |" or die "Can't read $sam, perhaps forgetting to use Samtools? :$!\n";
    }
}
while (<SAM>) {
    $count++;

    if ($_ =~ /^\@/){ # these headers are only available in <SAM> when the original file is in the sam format
        if ($printHeader == 1){
            print;
        }
        next;
    }
 
    
    # for each read line:
    my @temp = split /\t/;
    
    # check SAM format
    unless (scalar(@temp)>=9){ 
        print STDERR $_, " less than 9 columns\n";
        next;
    } # the read line in SAM format should have >= 9 columns
    
    my ($qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, @tabs) = split /\t/;
    
    # unmapped reads
    if ($rname eq "*" || $cigar eq "*"){
        #print STDERR $qname, " unmapped reads.  filtered\n";
        next;
    }
    
    # mapped reads
    my $mismatch = -1;
    for my $t (@tabs){ # search the tag for mismatch bases
        if ($t =~ /^[NX]M\:i/){
            $mismatch = getMismatch($t);
            last;
        }
    }
    
    my $contig = $rname;
    $contig =~ s/lcl\|//;
    #$contig2bact{$contig} = queryContig($contig) unless $contig2bact{$contig}; # fetch the BACT_ID from the database mbiome, (HMP genomes only)
    
    # error checking
    if ($mismatch == -1) {
        if ($cigar =~ /1S$/){ #BWA can't handle the cigar code properly when it ends with 1S?!
            # also these are 100% matched reads
            print join "\t", $contig, $pos, length($seq), 100, 100, length($seq), "0\n" if $option==4;
        }
        
        warn ("Error: Mismatch info missing: XM:i:\\d+.  This read has unknown seq identity. ");
        $mismatch = 0;
        warn ("Read ignored: $_");
        next;
    }
    
    # get BACT ID of the HMP genomes
    
    # calculate alignment related metrics
    my $alignLength = decodeCigar($cigar); # this returns the sum of the ungapped alignments
    my $seqIdent = 100*($alignLength-$mismatch)/length($seq); # the sequence identity over the total length
    my $alignIdent = 100*($alignLength-$mismatch)/$alignLength; # the sequence identity over the total aligned bases
    my $mPercent = 100*(decodeCigar($cigar)/length($seq)); # the percent of ungapped alignment of the total bases
    
    # apply the user-defined filters and output the final results to STDOUT
    if ($seqIdent>=$ident && $mPercent>=$M && $alignIdent>=$A) { # if the read passes the cutoff
        
        if ($option == 1 || $option == 3) {
            print;
            if ($option == 3) {
                print STDERR join "\t", $contig2bact{$contig}, $rname, $pos, $seqIdent, $alignIdent, $mPercent, $cigar, $mismatch, length($seq)."\n";
            }
        }
        elsif ($option == 4){
            print join "\t", $contig2bact{$contig}, $contig, $pos, length($seq), $seqIdent, $alignIdent, $alignLength, $mismatch."\n"; 
        }
    
        else {
            print join "\t", $contig2bact{$contig}, $rname, $pos, $seqIdent, $mPercent, $cigar, $mismatch, length($seq)."\n";
        }
    }
    else { # filtered reads are printed to STDERR
        # the fields are: BACT ID, READ ID, Seq identity, Ungapped Alignments, Aligned Identity, Seq length, "filtered" 
        #print STDERR join "\t", $contig2bact{$contig}, $rname, $seqIdent, $mPercent, $alignIdent, length($seq), "filtered\n";
    }
}
close SAM;
print STDERR $count." reads processed\n"; 
exit;
