import sys
from Bio import SeqIO
from FmIndex import FmIndex
from GIException import GIException
from AlignedReadDetails import AlignedReadDetails
from NeedlemanWunsch import NeedlemanWunsch
import util
import globalVariables
import logging


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
# file_name = "/logs/test_logger_" + str(globalVariables.match) + "_" + str(globalVariables.replacement) + "_" + str(
#     globalVariables.insertion) + ".log"
file_name = "/logs/log_" + str(globalVariables.seed_length) + "_" + str(globalVariables.margin) + ".log"
logging.basicConfig(level=logging.DEBUG,
                    format='%(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=file_name,
                    filemode='w')
logger = logging.getLogger("")  # get the root logger
util.print_command_line_arguments(logger)
try:
    fasta_records = SeqIO.index(globalVariables.fastaFile, "fasta")  # expects one record in file(one genome seq)
    fastq_records = SeqIO.index(globalVariables.fastqFile, "fastq")
except Exception as ex:
    logger.error("Failed to read sequence from file. Cause: " + str(ex))
    raise GIException("Failed to read sequence from file. Cause: " + str(ex))
# fasta_record = SeqIO.index(fastaFile, "fasta") can read multiple sequences, returns list of reads, supports none and multiple reads
fasta_record = fasta_records[list(fasta_records.keys())[0]]
result = dict()

if (len(fastq_records) == 0):
    logger.info("There are no reads in stated fastq file: " + globalVariables.fastqFile)
    print("There are no reads in stated fastq file: " + globalVariables.fastqFile)
else:
    fasta_sequence = fasta_record.seq
    logger.info("fasta seq length " + str(len(fasta_sequence)))
    fm = FmIndex(fasta_sequence)
    # global nw
    nw = NeedlemanWunsch(globalVariables.match, globalVariables.replacement, globalVariables.insertion)
    logger.info("Total reads " + str(len(fastq_records)))
    for read in list(fastq_records.keys()):
        read_sequnce = fastq_records[read].seq
        # print(read_sequnce)
        total_alignments = process_whole_read(fm, fasta_sequence, read_sequnce)
        # print("for loop", [(read.alignment_score, read.position, read.edit_transcript) for read in total_alignments])
        result[read] = total_alignments
    number_of_mapped_reads = 0
    logger.info("key,position,alignmentScore,reversed")
    for key in result.keys():
        # logger.info(key)
        # print(key)
        result_sort = sorted(result[key], key=lambda curr: curr.alignment_score, reverse=True)
        if (len(result_sort) > 0):
            number_of_mapped_reads = number_of_mapped_reads+ 1
        for read in result_sort:
            logger.info(key +"," + str(read.position) +"," + str(read.alignment_score) + "," + read.edit_transcript + "," + str(read.is_revesed))
            # print(read.position, read.alignment_score, read.edit_transcript, read.is_revesed)
        # print(key, sizeof(result[key]))
    logger.info("Number of mapped reads " + str(number_of_mapped_reads))