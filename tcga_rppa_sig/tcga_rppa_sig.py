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
from collections import defaultdict
from itertools import product


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
    parser.add_argument(
        '--join',
        metavar='FILE',
        type=str,
        required=False,
        help='Options filename to output join of RPPA and SIG data')
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

def cancer_types_per_sample(sample_rppa):
    counts = [] 
    for sample in sample_rppa:
        num_cancers = len(sample_rppa[sample])
        counts.append(num_cancers)
    return min(counts), max(counts)


def read_rppa(rppa_filename):
    # sample_id->cancer->sample_type->data
    results = {}
    rppa_sample_ids = set()
    rppa_cancer_types = set()
    rppa_sample_types = defaultdict(int) 
    num_missing_protein_values = 0
    num_valid_protein_values = 0
    with open(rppa_filename) as rppa_file:
        reader = csv.DictReader(rppa_file)
        fieldnames = reader.fieldnames
        # fieldnames start with "Sample_ID,Cancer_Type,Sample_Type" then have protein names
        non_protein_fields = "Sample_ID,Cancer_Type,Sample_Type".split(',')
        protein_names = fieldnames[len(non_protein_fields):]
        for row in reader:
            this_sample = parse_rppa_sample_id(row['Sample_ID'])
            rppa_sample_ids.add(this_sample)
            this_cancer = row['Cancer_Type']
            rppa_cancer_types.add(this_cancer)
            this_sample_type = row['Sample_Type']
            rppa_sample_types[this_sample_type] += 1
            if this_sample not in results:
                results[this_sample] = {}
            if this_cancer not in results[this_sample]:
                results[this_sample][this_cancer] = {}
            if this_sample_type not in results[this_sample][this_cancer]:
                results[this_sample][this_cancer][this_sample_type] = {}
            for protein in protein_names:
                try:
                    this_protein_value = float(row[protein])
                    num_valid_protein_values += 1
                except:
                    this_protein_value = None 
                    num_missing_protein_values += 1
                results[this_sample][this_cancer][this_sample_type][protein] = this_protein_value
        logging.info("Number of missing protein values: {}".format(num_missing_protein_values))
        logging.info("Number of valid protein values: {}".format(num_valid_protein_values))
        return protein_names, rppa_sample_ids, rppa_cancer_types, rppa_sample_types, results


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
    # sample_id->field->value
    num_missing_sig_values = 0
    num_valid_sig_values = 0
    results = {}
    signature_associations = {}
    sigs_sample_ids = set()
    with open(sigs_filename) as sigs_file:
        reader = csv.DictReader(sigs_file, delimiter='\t')
        for row in reader:
            this_sample = parse_sig_sample_id(row['Tumor_Sample_Barcode'])
            sigs_sample_ids.add(this_sample)
            if this_sample not in results:
                results[this_sample] = {}
            this_signature = row['Signature']
            this_association = row['Association']
            if this_signature not in signature_associations:
                signature_associations[this_signature] = this_association
            try:
                this_signature_value = float(row['Contribution'])
                num_valid_sig_values += 1
            except:
                this_signature_value = None 
                num_missing_sig_values += 1
            if this_signature not in results[this_sample]:
                results[this_sample][this_signature] = this_signature_value
            else:
                exit_with_error("Duplicate signature {} for sample {}".format(this_signature, this_sample), EXIT_BAD_INPUT)
    logging.info("Number of missing signature values: {}".format(num_missing_sig_values))
    logging.info("Number of valid signature values: {}".format(num_valid_sig_values))
    return sigs_sample_ids, signature_associations, results


def join_on_sample(rppa, sigs):
    #results[this_sample][this_cancer][this_sample_type][protein] = this_protein_value
    # sample -> fields -> value
    results = {}
    for sample, cancer_types in rppa.items():
        if sample in sigs:
            for cancer_type, sample_types in cancer_types.items():
                for sample_type, proteins in sample_types.items():
                    # Only consider primaries at the moment
                    if sample_type == 'Primary':
                        if sample not in results:
                            results[sample] = {}
                        else:
                            exit_with_error("Duplicate sample in joined RPPA SIGs data: {}".format(sample), EXIT_BAD_INPUT)
                        if cancer_type not in results[sample]:
                            results[sample][cancer_type] = {'sigs': {}, 'proteins': {}}
                        sig_values = sigs[sample]
                        for s in sig_values:
                            results[sample][cancer_type]['sigs'][s] = sig_values[s]
                        for p in proteins:
                            results[sample][cancer_type]['proteins'][p] = proteins[p]
    return results

def output_join(output_filename, join_data, cancer_types, sig_names, protein_names):
    header = ["sample", "cancer"] + sig_names + protein_names
    with open(output_filename, "w") as out_file:
        print(",".join(header), file=out_file)
        for sample, cancer_types in join_data.items():
            for cancer_type in cancer_types:
                this_sigs = cancer_types[cancer_type]['sigs']
                this_proteins = cancer_types[cancer_type]['proteins']
                output_row = [sample, cancer_type]
                for sig in sig_names:
                    output_row.append(str(this_sigs[sig]))
                for protein in protein_names:
                    output_row.append(str(this_proteins[protein]))
                print(",".join(output_row), file=out_file)
        

def get_correlations(rppa_cancer_types, signature_names, protein_names, join_data):
    results = {}
    for cancer_type in rppa_cancer_types:
        results[cancer_type] = {}
        for s, p in product(signature_names, protein_names):
            results[cancer_type][(s, p)] = {'sigs': [], 'proteins': []}
    for sample, cancer_types in join_data.items():
        for cancer_type in cancer_types:
            this_sigs = cancer_types[cancer_type]['sigs']
            this_proteins = cancer_types[cancer_type]['proteins']
            if cancer_type in results:
                for s, p in product(signature_names, protein_names):
                    this_sig_val = this_sigs[s]
                    this_protein_val = this_proteins[p]
                    if this_sig_val is not None and this_protein_val is not None:
                        results[cancer_type][(s,p)]['sigs'].append(this_sig_val)
                        results[cancer_type][(s,p)]['proteins'].append(this_protein_val)
            else:
                exit_with_error("Unexpected cancer type in joined data: {}".format(cancer_type), EXIT_BAD_INPUT)
    return results


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    rppa_protein_names, rppa_sample_ids, rppa_cancer_types, rppa_sample_types, rppa = read_rppa(options.rppa)
    min_cancer_types_per_sample, max_cancer_types_per_sample = cancer_types_per_sample(rppa)
    logging.info("Minumum cancer types per sample: {}".format(min_cancer_types_per_sample))
    logging.info("Maximum cancer types per sample: {}".format(max_cancer_types_per_sample))
    logging.info("Number of RRPA samples: {}".format(len(rppa_sample_ids)))
    logging.info("RPPA sample types: {}".format(dict(rppa_sample_types)))
    logging.info("RRPA protein names: {}".format(rppa_protein_names))
    logging.info("Number of RRPA proteins: {}".format(len(rppa_protein_names)))
    sigs_sample_ids, sigs_associations, sigs = read_sigs(options.sigs)
    logging.info("Number of SIGS samples: {}".format(len(sigs_sample_ids)))
    logging.info("Intersection of RPPA and SIGS samples: {}".format(len(sigs_sample_ids.intersection(rppa_sample_ids))))
    join_data = join_on_sample(rppa, sigs)
    sig_names = sorted(sigs_associations.keys())
    if options.join:
        output_join(options.join, join_data, rppa_cancer_types, sig_names, rppa_protein_names)
    get_correlations(rppa_cancer_types, sig_names, rppa_protein_names, join_data)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
