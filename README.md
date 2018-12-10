# PuppyTails
PuppyTails System Requirements
The PuppyTails program can be run on a desktop computer with the hardware listed below.  It can also been run on a Linux-based computational server with the appropriate version of Python.  The requirements below are minimum requirements for the desktop computer.  Program run time is much decreased when using a computational server.
Hardware Requirements
RAM: 5+ GB
CPU: Intel core i7-3770 CPU@3.40GHz

Software Requirements
Linux 3.13.0-46-generic (x86_64) OS: Ubuntu 14.04.2 LTS
Python 3.3.3

PuppyTails Documentation
Programs written by Douglas Porter

extract_inserts.py 
Identifies tails from every .fastq file in a directory. 

Usage:
python extract_inserts.py -d directory
The -d parameter sets the folder with the fastq files of interest. The program opens every .fastq file in the current directory if the -d parameter is not set.

If the -s option is not set, the program expects MS2-tagged tRNA ends.
The targeted RNA end can be given in a file passed with the -s filename argument.
If the targeted RNA is the yeast three-hybrid RNA, the –y3h mode option should be set.

The format of a tails file is as follows. For a tRNA, we expect, with X =some sequence: 
file format:
XCCA
XCC
XCA
XC
X
etc. 

Tails are searched for in the order given in the file. The X is the fixed sequence and must be present in all forms of the tail end. Since these are present in every tail output, they are not included in the output. The variations after the fixed sequence are the variable sequence.

The program will extract tRNA 3ʹ end sequences to determine the tRNA reporter 3ʹ end sequence, tail sequence, and random heptamer for each tRNA end.

It will then read tail sequences and remove PCR duplicates, using the random heptamer. It outputs a final .inserts file that contains the tRNA end, tail and number of tRNA ends (with different heptamers) containing that end type.

The .inserts file is sorted by the number of times each tail sequence is observed.

Expected output:
A,TTT,3
Meaning:
tRNA end,tail,number of unique heptamers.

fileInserts.py 

This is a class used by other PuppyTails programs.


analyzeTails.py 

Usage: 
python analyzeTails.py [options] 

Outputs a dataForGraph.txt file. If the control is subtracted, the dataForGraph.txt file is put in a controlSubtracted folder.
Automatically calls make figs.py after running.

Options:
-i is a folder of inserts files. By default, the ./inserts/ files will be searched.
-o will only use tails with a CCA end.
-n will normalize to dataset size (per 1 million).
-s will subtract a control.
-c is a control .inserts file to subtract.
-p will output values as percents.
-a will convert negative enrichment vs the control to zero. That is, there is no negative enrichment allowed.
-m will average the negative control files. If this option is used, the negative control files should be given as a list following the other inputs, rather than with the -c option.

Example usage: python analyzeTails.py -a -n -p -c -i folder of inserts files


make figs.py 
Creates figures from dataForGraph files. 

USAGE: python make figs.py -f data file 
Takes one argument, which is the data file to process. Outputs to a figs/ folder.





 
Total Percent Nucleotides Analyses
Perl script written by Christopher P. Lapointe

Software Requirements

OS: Windows 7, Windows 10

Perl command line


Total_Percent_nt_Analyses-vs1d.pl

Uses dataForGraph.txt files output from PuppyTails program as input

Uses SampleKey.out file to decode filenames/barcode numbers into sample names


Generates 2 files for each sample: 

Samplename.out 

File that has removed repeated lines from dataForGraph.txt file

Samplename-averages.out

File that summarizes total nucleotide addition metrics:

•	Sum Total unique tails	
•	Unweighted Average nucleotide incorporation (raw average of calculated percentages at each tail length)
•	Total of each nucleotide added	
•	Total number of nucleotides added	
•	Weighted Average nucleotide incorporation	sum of each nucleotide added divided total of all nucleotides added

Generates 2 files per set of data input, with metrics listed: 

UnweightedTailCompSummaryForAllSamples.out

•	Sum total unique tails	
•	Unweighted Average incorporation for each nucleotide
•	Average Tail Length	
•	Percent of reporter tRNAs tailed

WeightedTailCompSummaryForAllSamples.out

•	Sum Total unique tails
•	Total number of nucleotides added	
•	Weighted Average incorporation for each nucleotide
•	Average Tail Length	
•	Percent of reporter tRNAs tailed
