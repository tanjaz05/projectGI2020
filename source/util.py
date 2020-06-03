from GIException import GIException
import globalVariables
import argparse


def set_fasta_file(value):
    globalVariables.fastaFile = value


def set_fastq_file(value):
    globalVariables.fastqFile = value


def set_seed(value):
    if value <= 0:
        raise GIException("Invalid parameter - seed length should be number greater than 0")
    globalVariables.seed_length = int(value)


def set_margin(value):
    globalVariables.margin = int(value)


def set_match(value):
    globalVariables.match = value


def set_replacement(value):
    globalVariables.replacement = value


def set_insertion_deletion(value):
    globalVariables.insertion = value


def print_command_line_arguments(logger):
    logger.info("Parameters successfully read:")
    logger.info("fasta file: " + globalVariables.fastaFile)
    logger.info("fastq file: " + globalVariables.fastqFile)
    logger.info("seed length: " + str(globalVariables.seed_length))
    logger.info("margin: " + str(globalVariables.margin))
    logger.info("Scoring matrix values:")
    logger.info("match: " + str(globalVariables.match))
    logger.info("mismatch: " + str(globalVariables.replacement))
    logger.info("insertion/deletion: " + str(globalVariables.insertion))


def init_arguments():
    parser = argparse.ArgumentParser("Performs alignment of reads to reference genome")
    parser.add_argument("f", help="path to fasta file containing reference genome", type=str)
    parser.add_argument("q", help="path to fastq file containing collection of reads", type=str)
    parser.add_argument("s", help="length of seed for seed and extend method", type=int)
    parser.add_argument("mg", help="how much longer the reference should be than the rest of the read", type=int, choices=[0,1,2,3])
    parser.add_argument("m", help="match value for scoring matrix", type=int)
    parser.add_argument("r", help="mismatch value for scoring matrix", type=int)
    parser.add_argument("i", help="value for insertion/deletion for scoring matrix", type=int)
    args = parser.parse_args()
    set_fasta_file(args.f)
    set_fastq_file(args.q)
    set_seed(args.s)
    set_margin(args.mg)
    set_match(args.m)
    set_replacement(args.r)
    set_insertion_deletion(args.i)
