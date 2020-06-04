## Read me

In the source folder you can find python implementation of alignment of reads to reference genome sequence. The process of alignment is complex and demanding and divided into several steps each using a different algorithm. “Seed and Extend” method in combination with Burrows-Wheeler Transformation and FM index (Full-text index in Minute space) is used for indexed string search in a given text and string alignment is used Needleman-Wunsch algorithm.
Seed and Extend is a method that first performs “coarse” search by searching for patterns (seeds) extracted from a string that is searched for in the given text (BWT and FM index) and afterwards does fine alignments in the vicinity of the seeds and chooses the best scoring fine-grain alignment (Needleman-Wunsch).

In the files folder you can find the files which were used for testing. 

In the results folder you can find the log files, diagrams and also the sam file. 

## Modules

The project contains 7 modules:
seed_and_extend.py - the main module where all data is initialized and processed
util.py - a module where all parameters are read
FmIndex.py - a module where the creation of the FM index is implemented as well as all other operations with FM index 
FmCheckpoints.py - a module where creation and all operations with checkpoints are implemented 
bwt.py - module for creation of Burrows-Wheeler transformation and suffix array
NeedlemanWunsch.py - a module where global alignment is implemented
GIException - a class where a custom exception is implemented
AlignedReadDetails.py - a class that contains all read details
globalVariables.py - where all global variables are declared

## How to run

Run this program by running seed_and_extend.py script and passing all needed parameters:
	- path to fasta file containing reference genome
	- path to fastq file containing collection of reads
	- length of seed - a number greater than 0 that will be used in Seed and Extend algorithm and it represents the number of characters from the start of the read (pattern) that will be used in bwt algorithm for string searching
	- margin - a number between 0 and 3 (0 and 3 included) that represents how many more characters the reference genome will have compared to the rest of the read for string alignment
	- match - value for scoring matrix in case of a match
	- replacement - value for scoring matrix in case of a mismatch 
	- insertion - value for scoring matrix in case of an insertion or deletion

**Example:**

seed_and_extend.py ..\files\reference.fasta ..\files\reads.fastq 10 0 1 -2 -7

## Video presentation

Link to PowerPoint presentation: https://docs.google.com/presentation/d/18OK33_N61_ao7_ZU431dmy6bvtdFv7XwpUFq0SVmnQo/edit?usp=sharing

Link to video presentation: https://youtu.be/xY9z960uiX4
