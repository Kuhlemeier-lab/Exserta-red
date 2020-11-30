#!/usr/bin/perl
use strict;
use warnings;

my %contigs;



# Read in a list of contig names
open(IN, '<', "Pax_303_myb_list");

foreach my $line(<IN>) {
    chomp $line;
    $contigs{$line} = 1;
}
close IN;

# Read in the transcriptome file

my (%fastahash) = &arrays2hash(&get_fasta_names_and_seqs("./PeaxHIC303.pep.fasta"));
open (OUT, ">Pax_303_mybs_aa_seqs.fasta"); 
foreach my $prekey(keys %fastahash) {
    my @keycol = split(' ', $prekey);
    my $contig = $keycol[0];

# For each of the listed contig names above, pull out the sequence and print it to a file 
# in FASTA format
    if (exists $contigs{"$contig"}) {
	print OUT ">$contig\n$fastahash{$prekey}\n\n";
    }
}
close OUT;
print "Break here\n";





#####################################################
#arrays2hash(\@array1, \@array2)
#a subroutine to take two arrays of equal size and 
#convert them to a hash keys taken from first array,
#values from the second.
#####################################################

sub arrays2hash {
    use strict;
    use warnings;

    (my $keyarray, my $valuearray) = @_;
    if (scalar(@$keyarray) != scalar(@$valuearray)) {
        die "Arrays differ in size: Mismatched number of keys and values"; 
    }
    
    my %newhash = ( );
    
    @newhash{ @$keyarray } = @$valuearray;


    return (%newhash);

}

#####################################################
#END arrays2hash
#####################################################




#####################################################
#get_fasta_names_and_seqs("filename")
#returns matching arrays of fasta heders and seqs 
#given a fastname filename
#####################################################

sub get_fasta_names_and_seqs {
    use strict;
    use warnings;

    my ($inputfilename) = @_;
    my @fasta_names = ();
    my @fasta_seqs= ();

           
    unless ( open(FILEDATA, $inputfilename) ) {
        print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
        exit; #exit the program
    }    

    my @filedata = <FILEDATA>; #Read the lines of the file into an array
    close FILEDATA;
    
    my $seq_count = 0; #this will be used to keep track of the number of sequences
    foreach my $line (@filedata){
        chomp $line;
        if ($line =~ /^\s*$/) {next;} #ignore line if it is blank
        
        elsif ($line =~ /^>/) { #if the line is a header line (begins with ">")...
            if ($line =~ /^>.*[\w]+/){
                my $partialLine = substr ($&, 1);
                push (@fasta_names, $partialLine); #add that line to an array of fasta names
                push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
                ++$seq_count; #also increment our counter which keeps track of sequence number
            }
        }    
        
        else { #if the line's not blank or a header, add it to the current sequence 
            $fasta_seqs[$seq_count-1] .= $line;
        }
        
        $fasta_seqs[$seq_count-1] =~s/\s//g; #remove all whitespace from the current  sequence
    }
    
    return (\@fasta_names, \@fasta_seqs);

}

#####################################################
#end get_fasta_names_and_seqs()
#####################################################
