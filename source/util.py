from GIException import GIException
import sys
import globalVariables


def set_fasta_file(value):
    globalVariables.fastaFile = value


def set_fastq_file(value):
    globalVariables.fastqFile = value


def set_seed(value):
    if (not value.isnumeric() or int(value) <= 0):
        raise GIException("Invalid parameter - seed length should be number greater than 0")
    globalVariables.seed_length = int(value)


def set_margin(value):
    if (not value.isnumeric() or int(value) < 0 or int(value) > 3):
        raise GIException("Invalid parameter - margin should be number between 0 and 3 (0 and 3 included)")
    globalVariables.margin = int(value)


def set_match(value):
    val = value
    sign = 1
    if (value[0] == "-"):
        val = value[1:]
        sign = -1
    if (not val.isnumeric()):
        raise GIException("Invalid parameter - match should be number")
    globalVariables.match = sign * int(val)


def set_replacement(value):
    val = value
    sign = 1
    if (value[0] == "-"):
        val = value[1:]
        sign = -1
    if (not val.isnumeric()):
        raise GIException("Invalid parameter - replacement should be number")
    globalVariables.replacement = sign * int(val)


def set_insertion_deletion(value):
    val = value
    sign = 1
    if (value[0] == "-"):
        val = value[1:]
        sign = -1
    if (not val.isnumeric()):
        raise GIException("Invalid parameter - insertion should be number")
    globalVariables.insertion = sign * int(val)


def invalid_parameter(value):
    raise GIException("Unknown parameter " + value)


def read_command_line_arguments():
    if (len(sys.argv) < 15):
        # logger.error("Missing parameters - required: fasta file path, fastq file path, seed length, margin, match. replace, insertion/deletions")
        raise GIException(
            "Missing parameters - required: fasta file path, fastq file path, seed length, margin, match. replace, insertion/deletions")
    switcher = {
        "-f": set_fasta_file,
        "-q": set_fastq_file,
        "-s": set_seed,
        "-mg": set_margin,
        "-m": set_match,
        "-r": set_replacement,
        "-i": set_insertion_deletion
    }
    for i in range(len(sys.argv)):
        if i % 2 == 1:
            set_parameter = switcher.get(sys.argv[i], invalid_parameter)
            arg = sys.argv[i] if set_parameter == invalid_parameter else sys.argv[i + 1];
            set_parameter(arg)
    # print_command_line_arguments(logger)


def print_command_line_arguments(logger):
    # print("Parameters successfully read:")
    logger.info("Parameters successfully read:")
    logger.info("fasta file: " + globalVariables.fastaFile)
    logger.info("fastq file: " + globalVariables.fastqFile)
    logger.info("seed length: " + str(globalVariables.seed_length))
    logger.info("margin: " + str(globalVariables.margin))
    logger.info("Scoring matrix values:")
    logger.info("match: " + str(globalVariables.match))
    logger.info("mismatch: " + str(globalVariables.replacement))
    logger.info("insertion/deletion: " + str(globalVariables.insertion))
