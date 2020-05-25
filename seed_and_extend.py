import sys
from Bio import SeqIO
from FmIndex import FmIndex
from GIException import GIException
from AlignedReadDetails import AlignedReadDetails
from NeedlemanWunsch import NeedlemanWunsch
import util
import globalVariables


def align_read(x, y, occurences, is_reversed):
    list = []
    for index in occurences:
        start_position = index + globalVariables.seed_length
        read_length = len(y)
        dif_read_length = read_length + globalVariables.margin
        end_position = start_position + dif_read_length
        # print("fasta", x)
        # print("fastq", y)
        if (len(x) >= end_position):
            D, globalAlignmentScore = nw.globalAlignment(x[start_position:end_position], y, nw.scoringMatrix)
            # print(globalAlignmentScore)
            # print(D)
            a, transcript = nw.traceback(x[start_position: start_position + dif_read_length], y, D, nw.scoringMatrix)
            # print(a)
            # print(transcript)
            aligned_read = AlignedReadDetails(is_reversed, index, globalAlignmentScore, transcript)
        else:
            # print("Match was found too close to the end of referent sequnce, string can not be aligned")
            aligned_read = AlignedReadDetails(is_reversed, index, -sys.maxsize - 1, "")
        list.append(aligned_read)
        # print(transcript)
        # print(a)
    return list


def find_postion_and_align_read(fm, fasta_sequence, read, is_reversed):
    # fm = FmIndex(fasta_sequence)
    pattern = read[:globalVariables.seed_length]
    alignmetn_list = []
    if (fm.hasSubstring(pattern)):
        occurences = fm.occurrences(pattern)
        alignmetn_list = align_read(fasta_sequence, read_sequnce[globalVariables.seed_length:], occurences, is_reversed)
    return alignmetn_list


def process_whole_read(fm, fasta_sequence, read):
    reverse_read = read[len(read_sequnce)::-1]
    # print("reversed read", reverse_read)
    alignments_list = []
    # print("regular")
    alignments_list = alignments_list + find_postion_and_align_read(fm, fasta_sequence, read, False)
    # print("reverse")
    alignments_list = alignments_list + find_postion_and_align_read(fm, fasta_sequence, reverse_read, True)
    alignments_list.sort(key=lambda read: read.alignment_score, reverse=True)
    # print([(read.alignment_score, read.position, read.edit_transcript) for read in alignments_list])
    return alignments_list


util.read_command_line_arguments()

try:
    fasta_records = SeqIO.index(globalVariables.fastaFile, "fasta")  # expects one record in file(one genome seq)
    fastq_records = SeqIO.index(globalVariables.fastqFile, "fastq")
except Exception as ex:
    raise GIException("Failed to read sequence from file. Cause: " + str(ex))
# fasta_record = SeqIO.index(fastaFile, "fasta") can read multiple sequences, returns list of reads, supports none and multiple reads
fasta_record = fasta_records[list(fasta_records.keys())[0]]
result = dict()

if (len(fastq_records) == 0):
    print("There are no reads in stated fastq file: " + globalVariables.fastqFile)
else:
    fasta_sequence = fasta_record.seq
    fm = FmIndex(fasta_sequence)
    # global nw
    nw = NeedlemanWunsch(globalVariables.match, globalVariables.replacement, globalVariables.insertion)
    for read in list(fastq_records.keys()):
        read_sequnce = fastq_records[read].seq
        # print(read_sequnce)
        total_alignments = process_whole_read(fm, fasta_sequence, read_sequnce)
        # print("for loop", [(read.alignment_score, read.position, read.edit_transcript) for read in total_alignments])
        result[read] = total_alignments

    for key in result.keys():
        print(key)
        result_sort = sorted(result[key], key=lambda curr: curr.alignment_score, reverse=True)
        for read in result_sort:
            print(read.position, read.alignment_score, read.edit_transcript, read.is_revesed)
        # print(key, sizeof(result[key]))
