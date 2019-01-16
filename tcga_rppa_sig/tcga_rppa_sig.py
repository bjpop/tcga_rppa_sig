'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 17 Jan 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

The program reads one or more input FASTA files. For each file it computes a
variety of statistics, and then prints a summary of the statistics as output.
'''

from argparse import ArgumentParser
import csv
import sys
import logging
import pkg_resources


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_BAD_INPUT = 3
DEFAULT_VERBOSE = False
PROGRAM_NAME = "tcga_rppa_sig"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more FASTA files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--rppa',
        metavar='FILE',
        type=str,
        required=True,
        help='Filename of RPPA data')
    parser.add_argument(
        '--sigs',
        metavar='FILE',
        type=str,
        required=True,
        help='Filename of mutational signature data')
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    return parser.parse_args()


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))

# Sample IDs look like:
# TCGA-AO-A1KP-01A-21-A17I-20
# we only use the first 12 chars
def parse_rppa_sample_id(sample_id):
    if len(sample_id) >= 12:
        return sample_id[:12]
    else:
        exit_with_error("Cannot parse RPPA sample id: {}".format(sample_id), EXIT_BAD_INPUT)


def read_rppa(rppa_filename):
    # sample_id->cancer->sample_type->data
    results = {}
    rppa_sample_ids = set()
    rppa_cancer_types = set()
    rppa_sample_types = set()
    with open(rppa_filename) as rppa_file:
        reader = csv.DictReader(rppa_file)
        for row in reader:
            this_sample = parse_rppa_sample_id(row['Sample_ID'])
            rppa_sample_ids.add(this_sample)
            this_cancer = row['Cancer_Type']
            rppa_cancer_types.add(this_cancer)
            this_sample_type = row['Sample_Type']
            rppa_sample_types.add(this_sample_type)
            if this_sample not in results:
                results[this_sample] = {}
            if this_cancer not in results[this_sample]:
                results[this_sample][this_cancer] = {}
            results[this_sample][this_cancer][this_sample_type] = row
    return rppa_sample_ids, rppa_cancer_types, rppa_sample_types, results


# Sample IDs look like:
# TCGA.OR.A5J1 
# we join the fields to be consistent with the RPPA datah
# TCGA-OR-A5J1
# There are some non TCGA sample IDs in the data, we return them unchanged
def parse_sig_sample_id(sample_id):
    if sample_id.startswith('TCGA') and len(sample_id) == 12:
        return sample_id.replace('.', '-')
    else:
        #logging.info("Non TCGA sample ID in signature data: {}".format(sample_id))
        return sample_id 


def read_sigs(sigs_filename):
    # sample_id->data
    results = {}
    sigs_sample_ids = set()
    with open(sigs_filename) as sigs_file:
        reader = csv.DictReader(sigs_file, delimiter='\t')
        for row in reader:
            this_sample = parse_sig_sample_id(row['Tumor_Sample_Barcode'])
            sigs_sample_ids.add(this_sample)
    return sigs_sample_ids, results


def join_data(rppa, sigs):
    pass

def get_correlations(rppa_sigs):
    pass

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    rppa_sample_ids, rppa_cancer_types, rppa_sample_types, rppa = read_rppa(options.rppa)
    sigs_sample_ids, sigs = read_sigs(options.sigs)
    logging.info("Number of RRPA samples: {}".format(len(rppa_sample_ids)))
    logging.info("Number of SIGS samples: {}".format(len(sigs_sample_ids)))
    logging.info("Intersection of RPPA and SIGS samples: {}".format(len(sigs_sample_ids.intersection(rppa_sample_ids))))
    rppa_sigs = join_data(rppa, sigs)
    get_correlations(rppa_sigs)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
