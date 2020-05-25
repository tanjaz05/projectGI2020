## Read me

This repository contains python implementation of alignment of reads to reference genome sequence. The process of alignment is complex and demanding and divided into several steps each using a different algorithm. “Seed and Extend” method in combination with Burrows-Wheeler Transformation and FM index (Full-text index in Minute space) is used for indexed string search in a given text and for string alignment is used Needleman-Wunsch algorithm.
Seed and Extend is a method that first performs “coarse” search by searching for patterns (seeds) extracted from a string that is searched for in the given text (BWT and FM index) and afterwards does fine alignments in the vicinity of the seeds and chooses the best scoring fine-grain alignment (Needleman-Wunsch).

## Modules

The project contains 7 modules:
seed_and_extend.py - a main module where all all data is initialized and processed
util.py - a module where all parameters are read
FmIndex.py - a module where the creation of the FM index is implemented as well as all other operations with FM index 
FmCheckpoints.py - a module where creation and all operations with checkpoints are implemented 
bwt.py - module for creation of Burrows-Wheeler transformation and suffix array
NeedlemanWunsch.py - a module where global alignment is implemented
GIException - a class where a custom exception is implemented
AlignedReadDetails.py - a class that containes all read detailes
globalVariables.py - where all global variables are declared

## How to run

Run this program by running seed_and_extend.py script and passing all needed parameters:
-f path_to_fasta_file - file containing referent genome sequence
-q path_to_fastq_file - file containing a collection of reads which will be aligned to referent genome
-s seed - a number greater than 0 that will be used in Seed and Extend algorithm and it represents the number of characters from the start of the read (pattern) that will be used in bwt algorithm for string searching
-mg margin - a number between 0 and 3 (0 and 3 included) that represents how many more characters the reference genome will have compared to the rest of the read for string alignment
-m match - value for scoring matrix in case of a match
-r replacement - value for scoring matrix in case of a mismatch
-i insertion - value for scoring matrix in case of insertion or deletion

**example:**

seed_and_extend.py -f C:\Users\User\Downloads\genomska\sample.fa -q C:\Users\User\Downloads\genomska\sampletest.fq -s 10 -mg 0 -m 1 -i -7 -r -2

## Video presentation

Link to PowerPoint presentation: https://docs.google.com/presentation/d/18OK33_N61_ao7_ZU431dmy6bvtdFv7XwpUFq0SVmnQo/edit?usp=sharing

Link to video presentation: https://youtu.be/xY9z960uiX4




