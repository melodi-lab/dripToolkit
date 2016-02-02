#!/usr/bin/env python
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

from __future__ import with_statement

__authors__ = ['John Halloran <halloj3@uw.edu>' ]

import os
import re
import sys
import argparse
import cPickle as pickle
# import dripDigest as df
import multiprocessing
import shlex
import pyFiles.psm as psm
import struct
import glob
import itertools

from pyFiles.recalibrateCharge import recalibrate_charge_psms, write_drip_recal
from shutil import rmtree
from pyFiles.spectrum import MS2Spectrum
from pyFiles.peptide import Peptide, amino_acids_to_indices
from pyFiles.normalize import pipeline
from pyFiles.pfile.wrapper import PFile
# from pyFiles.args import make_dir_callback
from dripDigest import (check_arg_trueFalse,
                        parse_var_mods)
from pyFiles.dripEncoding import (create_drip_structure,
                                  create_drip_master,
                                  triangulate_drip,
                                  write_covar_file,
                                  create_pfile, 
                                  drip_peptide_sentence,
                                  drip_spectrum_sentence,
                                  make_master_parameters,
                                  make_master_parameters_lowres,
                                  interleave_b_y_ions,
                                  interleave_b_y_ions_var_mods,
                                  interleave_b_y_ions_lowres,
                                  interleave_b_y_ions_var_mods_lowres,
                                  filter_theoretical_peaks,
                                  filter_theoretical_peaks_lowres,
                                  load_drip_means)

from pyFiles.constants import allPeps
from pyFiles.shard_spectra import (calc_minMaxMz,
                                   pickle_candidate_spectra, 
                                   pickle_candidate_binarydb_spectra,
                                   candidate_spectra_generator,
                                   candidate_spectra_memeffic_generator,
                                   candidate_binarydb_spectra_generator)
from pyFiles.process_vitvals import make_ident_persid
from subprocess import call, check_output
from pyFiles.ioPsmFunctions import load_drip_psms

# set stdout and stderr for subprocess
# stdo = open(os.devnull, "w")

# stdo = sys.stdout
stdo = open(os.devnull, "w")
stde = sys.stdout
# stde = stdo

# stde = open('gmtk_err', "w")

# # set stdout for all default streams
# sys.stdout = open("dripSearch_output", "w")

def copyArgs(argsA, argsB):
    """ Copy previously selected arguments to current parsed arguments
    """
    argsA.max_length = argsB.max_length
    argsA.min_length = argsB.min_length
    argsA.max_mass = argsB.max_mass
    argsA.mods_spec = argsB.mods_spec
    argsA.cterm_peptide_mods_spec = argsB.cterm_peptide_mods_spec
    argsA.nterm_peptide_mods_spec = argsB.nterm_peptide_mods_spec
    argsA.max_mods = argsB.max_mods
    argsA.min_mods = argsB.min_mods
    argsA.decoy_format = argsB.decoy_format
    argsA.keep_terminal_aminos = argsB.keep_terminal_aminos
    argsA.seed = argsB.seed
    argsA.enzyme = argsB.enzyme
    argsA.custom_enzyme = argsB.custom_enzyme
    argsA.missed_cleavages = argsB.missed_cleavages
    argsA.digestion = argsB.digestion
    argsA.peptide_buffer = argsB.peptide_buffer

    argsA.monoisotopic_precursor = argsB.monoisotopic_precursor
    argsA.precursor_window = argsB.precursor_window
    argsA.precursor_window_type = argsB.precursor_window_type
    argsA.scan_id_list = argsB.scan_id_list
    argsA.charges = argsB.charges
    argsA.high_res_ms2 = argsB.high_res_ms2
    argsA.decoys = argsB.decoys
    argsA.num_threads = argsB.num_threads
    argsA.top_match = argsB.top_match
    argsA.beam = argsB.beam
    argsA.num_jobs = argsB.num_jobs
    argsA.cluster_mode = argsB.cluster_mode
    argsA.write_cluster_scripts = argsB.write_cluster_scripts
    argsA.random_wait = argsB.random_wait
    argsA.min_spectrum_length = argsB.min_spectrum_length
    argsA.max_spectrum_length = argsB.max_spectrum_length
    argsA.ident = argsB.ident
    # argsA.output_dir = argsB.output_dir
    # argsA.obs_dir = argsB.obs_dir
    argsA.shards = argsB.shards
    argsA.num_spectra = argsB.num_spectra
    argsA.ppm = argsB.ppm
    argsA.filter_ident = argsB.filter_ident
    argsA.max_obs_mass = argsB.max_obs_mass
    argsA.normalize = argsB.normalize
    argsA.filt_theo_peaks = argsB.filt_theo_peaks
    argsA.per_spectrum_mz_bound = argsB.per_spectrum_mz_bound
    argsA.mz_lb = argsB.mz_lb
    argsA.mz_ub = argsB.mz_ub
    argsA.structure_file = argsB.structure_file
    argsA.master_file = argsB.master_file
    argsA.collection_dir = argsB.collection_dir
    argsA.learned_means = argsB.learned_means
    argsA.learned_covars = argsB.learned_covars
    argsA.mean_file = argsB.mean_file
    # argsA.covar_file = argsB.covar_file
    argsA.gauss_file = argsB.gauss_file
    argsA.mixture_file = argsB.mixture_file
    argsA.collection_file = argsB.collection_file
    argsA.covar_name = argsB.covar_name
    argsA.mixture_name = argsB.mixture_name
    argsA.gaussian_component_name = argsB.gaussian_component_name
    argsA.collection_name = argsB.collection_name
    argsA.dpmf_name = argsB.dpmf_name
    # argsA.peptide_database_file = argsB.peptide_database_file
    # argsA.load_peptide_database_file = argsB.load_peptide_database_file
    # argsA.target_output_file = argsB.target_output_file
    # argsA.decoy_output_file = argsB.decoy_output_file
    argsA.logDir = argsB.logDir
    argsA.pepDBperSidList = argsB.pepDBperSidList
    argsA.recalibrate = argsB.recalibrate
    argsA.cluster_dir = argsB.cluster_dir

    # check input arguments, create necessary directories, 
    # set necessary global variables based on input arguments
    # constant string
    if args.max_obs_mass < 0:
        print "Supplied max-mass %d, must be greater than 0, exitting" % args.max_obs_mass
        exit(-1)
    # check arguments
    assert(args.shards > 0)
    if not args.output:
        print "No output file specified for PSMs, exitting"
        exit(8)

    # create obervation directories
    # vit vals output directory
    try:
        base = os.path.abspath(args.logDir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.logDir)
        exit(-1)

    # top level encode directory
    try:
        base = os.path.abspath(args.output_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.output_dir)
        exit(-1)
    # observation directory
    try:
        base = os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        exit(-1)
    # Gaussian collection directory
    try:
        base = os.path.abspath(args.collection_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.collection_dir)
        exit(-1)


    args.covar_file = os.path.join(os.path.abspath(args.collection_dir), args.covar_file)
    args.mean_file = os.path.join(os.path.abspath(args.collection_dir), args.mean_file)
    args.gauss_file = os.path.join(os.path.abspath(args.collection_dir), args.gauss_file)
    args.mixture_file = os.path.join(os.path.abspath(args.collection_dir), args.mixture_file)
    args.collection_file = os.path.join(os.path.abspath(args.collection_dir), args.collection_file)

    if not args.filter_ident:
        args.ident = ''

def parseInputOptions():
    parser = argparse.ArgumentParser(conflict_handler='resolve')
    ############## input and output options
    iFileGroup = parser.add_argument_group('iFileGroup', 'Necessary input files.')
    help_spectra = '<string> - The name of the file from which to parse fragmentation spectra, in ms2 format.'
    iFileGroup.add_argument('--spectra', type = str, action = 'store',
                            help = help_spectra)
    help_pepdb = '<string> - Output directory for dripDigest.'
    iFileGroup.add_argument('--digest-dir', type = str, action = 'store',
                            help = help_pepdb, default = 'dripDigest-output')
    ############## search parameters
    searchParamsGroup = parser.add_argument_group('searchParamsGroup', 'Search parameter options.')
    help_precursor_window = """<float> - Tolerance used for matching peptides to spectra.  Peptides must be within +/-'precursor-window' of the spectrum value. The precursor window units depend upon precursor-window-type. Default=3."""
    searchParamsGroup.add_argument('--precursor-window', type = float, action = 'store', default = 3.0, help = help_precursor_window)
    help_precursor_window_type = """<Da|ppm> - Specify the units for the window that is used to select peptides around the precursor mass location, either in Daltons (Da) or parts-per-million (ppm). Default=Da."""
    searchParamsGroup.add_argument('--precursor-window-type', type = str, action = 'store', default = 'Da', help = help_precursor_window_type)
    help_scan_id_list = """<string> - A file containing a list of scan IDs to search.  Default = <empty>."""
    searchParamsGroup.add_argument('--scan-id-list', type = str, action = 'store', default = '', help = help_scan_id_list)
    help_charges = """<comma-separated-integers|all> - precursor charges to search. To specify individual charges, list as comma-separated, e.g., 1,2,3 to search all charge 1, 2, or 3 spectra. Default=All."""
    searchParamsGroup.add_argument('--charges', type = str, action = 'store', default = 'All', help = help_charges)
    help_high_res_ms2 = """<T|F> - boolean, whether the search is over high-res ms2 (high-high) spectra. When this parameter is true, DRIP used the real valued masses of candidate peptides as its Gaussian means. For low-res ms2 (low-low or high-low), the observed m/z measures are much less accurate so these Gaussian means are learned using training data (see dripTrain). Default=False."""
    searchParamsGroup.add_argument('--high-res-ms2', type = str, action = 'store', default = 'false', help = help_high_res_ms2)
    help_decoys = '<T|F> - whether to create (shuffle target peptides) and search decoy peptides. Default = True'
    searchParamsGroup.add_argument('--decoys', type = str, action = 'store', 
                                   default = 'True', help = help_decoys)
    help_num_threads = '<integer> - the number of threads to run on a multithreaded CPU. If supplied value is greater than number of supported threads, defaults to the maximum number of supported threads minus one. Multithreading is not suppored for cluster use as this is typically handled by the cluster job manager. Default=1.'
    searchParamsGroup.add_argument('--num-threads', type = int, action = 'store', 
                                   default = 1, help = help_num_threads)
    help_top_match = '<integer> - The number of psms per spectrum written to the output files. Default=1.'
    searchParamsGroup.add_argument('--top-match', type = int, action = 'store', 
                                   default = 1, help = help_top_match)
    help_beam = """<int> - K-beam width to use to speed up inference. Default value of 0 means exact inference. Warning - identifications may be significantly poor if the beam width is too small, i.e., beam < 100. Default = 0."""
    searchParamsGroup.add_argument('--beam', type = int, action = 'store', default = 0, help = help_beam)
    ############## Cluster usage parameters
    clusterUsageGroup = parser.add_argument_group('clusterUsageGroup', 'Cluster data generation options.')
    help_num_cluster_jobs = '<integer> - the number of jobs to run in parallel. Default=1.'
    clusterUsageGroup.add_argument('--num-cluster-jobs', type = int, action = 'store', 
                                   dest = 'num_jobs',
                                   default = 1, help = help_num_cluster_jobs)
    help_cluster_mode = '<T|F> - evaluate dripSearch prepared data as jobs on a cluster.  Only set this to true once dripSearch has been run to prepare data for cluster use.  Default = False'
    clusterUsageGroup.add_argument('--cluster-mode', type = str, action = 'store', 
                                   default = 'False', help = help_cluster_mode)
    help_write_cluster_scripts = '<T|F> - write scripts to be submitted to cluster queue.  Only used when num-jobs > 1.  Job outputs will be written to log subdirectory in current directory. Default = True'
    clusterUsageGroup.add_argument('--write-cluster-scripts', type = str, action = 'store', 
                                   default = 'True', help = help_write_cluster_scripts)
    help_cluster_dir = '<string> - absolute path of directory to run cluster jobs. Default = /tmp'
    clusterUsageGroup.add_argument('--cluster-dir', type = str, action = 'store', 
                                   default = '/tmp', help = help_cluster_dir)
    help_random_wait = '<integer> - randomly wait up to specified number of seconds before accessing NFS. Default=10'
    clusterUsageGroup.add_argument('--random-wait', type = int, action = 'store', 
                                   default = 10, help = help_random_wait)
    help_merge_cluster_results = '<T|F> - merge dripSearch cluster results collected in local directory log.  Default = False'
    clusterUsageGroup.add_argument('--merge-cluster-results', type = str, action = 'store', 
                                   default = 'False', help = help_merge_cluster_results)
    ############ shard spectra options
    parser.add_argument('--min_spectrum_length', type = int, action = 'store',
                        dest = 'min_spectrum_length', default = 1)
    parser.add_argument('--max_spectrum_length', type = int, action = 'store',
                        dest = 'max_spectrum_length', default = 10000)
    parser.add_argument('--ident', type = str, action= 'store', 
                        default = '')
    parser.add_argument('--output-dir', type = str, default = 'encode')
    parser.add_argument('--obs-dir', type = str, action = 'store',
                        default = 'obs')
    parser.add_argument('--shards', type = int, action = 'store',
                        default = 1)
    parser.add_argument('--num_spectra', type = int, action = 'store')
    parser.add_argument('--ppm', action = 'store_true',
                        default = False)
    parser.add_argument('--filter_ident', action = 'store_true',
                        default = False,
                      help = "Filter sids by ident")
    ######### encoding options
    parser.add_argument("--max_obs_mass", type = int, default = 0)
    parser.add_argument('--normalize', dest = "normalize", type = str,
                        help = "Name of the spectrum preprocessing pipeline.", 
                        default = 'top300TightSequest')
    # help_filt_theo_peaks = "Filter theoretical peaks outside of observed spectra min and max"
    # parser.add_argument('--filt_theo_peaks', action = 'store_false', dest = 'filt_theo_peaks', 
    #                     default = True, help = help_filt_theo_peaks)
    parser.add_argument('--per_spectrum_mz_bound', action = 'store_true', 
                        default = False, help = "Calculate observed m/z bound per spectrum")
    parser.add_argument('--mz_lb', type = float, action = 'store',
                        dest = 'mz_lb', help = 'Lower m/z bound' ,default = 150.0)
    parser.add_argument('--mz_ub', type = float, action = 'store',
                        dest = 'mz_ub', help = 'Upper m/z bound' ,default = 2000.0)
    parser.add_argument('--structure-file', type = str, action = 'store',
                        default = 'drip.str',
                        help = "DRIP GMTK structure file")
    parser.add_argument('--master-file', type = str, action = 'store',
                        default = 'drip.mtr',
                        help = "DRIP GMTK master file")
    parser.add_argument('--collection-dir', action = 'store', dest = 'collection_dir',
                        type = str, help = 'Where to store Gaussian files', 
                        default = 'drip_collection')
    parser.add_argument('--learned-means', action = 'store', dest = 'learned_means',
                        type = str, help = 'Previously learned means', default = "riptideLearnedMeans.txt")
    parser.add_argument('--learned-covars', action = 'store', dest = 'learned_covars',
                        type = str, help = 'Previously learned covariances', default = "riptideLearnedCovars.txt")
    help_mean_file = """<string> - Where to write model means. Default = drip.means"""
    parser.add_argument('--mean-file', type = str, action = 'store', default = "drip_mz.txt",
                            help = help_mean_file)
    help_covar_file = """<string> - Where to write model covariances. Default = covar.txt"""
    parser.add_argument('--covar-file', type = str, action = 'store', default = "covar.txt",
                            help = help_covar_file)

    parser.add_argument('--gauss_file', action = 'store', dest = 'gauss_file',
                        type = str, help = 'Gaussians.', 
                        default = 'drip_Gaussian_components.txt')

    parser.add_argument('--mixture_file', action = 'store', dest = 'mixture_file',
                        type = str, help = 'Gaussian mixtures.',
                        default = 'drip_Gaussian_mixtures.txt')

    parser.add_argument('--collection_file', action = 'store', dest = 'collection_file',
                        type = str, help = 'Where to put collection file.',
                        default = 'drip_mz_Gaussians.txt')

    parser.add_argument('--covar_name', action = 'store', dest = 'covar_name',
                        type = str, help = 'Variance name.',
                        default = 'covar0')

    parser.add_argument('--mixture_name', action = 'store', dest = 'mixture_name',
                        type = str, help = 'Gaussian mixture name.',
                        default = 'mixture')

    parser.add_argument('--gaussian_component_name', action = 'store', dest = 'gaussian_component_name',
                        type = str, help = 'Gaussian mixture name.',
                        default = 'gc')

    parser.add_argument('--collection_name', action = 'store', dest = 'collection_name',
                        type = str, help = 'Gaussian mixture name.', 
                        default = 'drip_mz_Gaussians')

    parser.add_argument('--dpmf_name', action = 'store', dest = 'dpmf_name',
                        type = str, help = 'DPMF name.', 
                        default = 'unityDPMF')
    # output file options
    oFileGroup = parser.add_argument_group('oFileGroup', 'Output file options')
    oFileGroup.add_argument('--logDir', type = str, 
                      help = 'directory to collect output DRIP results.', default = 'log')
    oFileGroup.add_argument('--pepDBperSidList', type = str,
                      help = 'file listing number of peptides in database per spectrum')
    oFileGroup.add_argument('--output', type = str,
                      help = 'identification file name')

    return parser.parse_args()

def runVit(vitString, stdo, stde):
    return call(vitStr.split(), stdout = stdo, stderr = stde)

def process_args(args):
    """ Check whether relevant directories already exist (remove them if so), create
        necessary directories.  Adjust parameters.

        pre:
        - args has been created by calling parse_args()
        
        post:
        - Relevant directories will be first removed then created, boolean parameters
          will be adjusted to reflect their boolean command line strings
    """

    # set true or false strings to booleans
    args.high_res_ms2 = check_arg_trueFalse(args.high_res_ms2)
    args.cluster_mode = check_arg_trueFalse(args.cluster_mode)
    args.write_cluster_scripts = check_arg_trueFalse(args.write_cluster_scripts)
    args.merge_cluster_results = check_arg_trueFalse(args.merge_cluster_results)

    args.filt_theo_peaks = True

    if args.merge_cluster_results:
        chargeRecalibrate(args.logDir,
                          args.output + '.txt', 
                          args.top_match)
        exit(0)

    if args.cluster_mode:
        return 1

    # set absolute path of learned means and covariance
    args.learned_means = os.path.abspath(args.learned_means)
    args.learned_covars = os.path.abspath(args.learned_covars)

    # check input arguments, create necessary directories, 
    # set necessary global variables based on input arguments
    # constant string
    if args.max_obs_mass < 0:
        print "Supplied max-mass %d, must be greater than 0, exitting" % args.max_obs_mass
        exit(-1)
    # check arguments
    assert(args.shards > 0)
    if not args.output:
        print "No output file specified for PSMs, exitting"
        exit(8)

    base = os.path.abspath(args.digest_dir)
    if not os.path.exists(base) and not args.cluster_mode:
        print "Digest directory %s does not exist." % (args.digest_dir)
        print "Please run dripDigest first and specify the resulting directory with --digest-dir."
        exit(-1)
    else:
        args.digest_dir = os.path.abspath(args.digest_dir)

    # load dripDigest options used to digest the fasta file
    ddo = pickle.load(open(os.path.join(args.digest_dir, 'options.pickle')))
    # set parameters to those used by dripDigest
    args.max_length = ddo.max_length
    args.max_mass = ddo.max_mass
    args.min_length = ddo.min_length
    args.monoisotopic_precursor = ddo.monoisotopic_precursor
    args.mods_spec = ddo.mods_spec
    args.cterm_peptide_mods_spec = ddo.cterm_peptide_mods_spec
    args.nterm_peptide_mods_spec = ddo.nterm_peptide_mods_spec
    args.max_mods = ddo.max_mods
    args.min_mods = ddo.min_mods
    args.decoys = ddo.decoys
    args.decoy_format = ddo.decoy_format
    args.keep_terminal_aminos = ddo.keep_terminal_aminos
    args.seed = ddo.seed
    args.enzyme = ddo.enzyme
    args.custom_enzyme = ddo.custom_enzyme
    args.missed_cleavages = ddo.missed_cleavages
    args.digestion = ddo.digestion
    args.recalibrate = ddo.recalibrate
    args.peptide_buffer = ddo.peptide_buffer

    # create obervation directories
    # vit vals output directory
    try:
        base = os.path.abspath(args.logDir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.logDir)
        exit(-1)

    # top level encode directory
    try:
        base = os.path.abspath(args.output_dir)
        args.output_dir = base
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.output_dir)
        exit(-1)
    # observation directory
    try:
        base = os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        exit(-1)
    # Gaussian collection directory
    try:
        base = os.path.abspath(args.collection_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.collection_dir)
        exit(-1)

    args.covar_file = os.path.join(os.path.abspath(args.collection_dir), args.covar_file)

    if not args.filter_ident:
        args.ident = ''

    # check precursor mass type
    pmt = args.precursor_window_type.lower()
    if pmt == 'da':
        args.ppm = False
    else:
        args.ppm = True

    # make sure number of input threads does not exceed number of supported threads
    if args.num_threads > multiprocessing.cpu_count():
        args.num_threads = multiprocessing.cpu_count()
        # args.num_threads = max(multiprocessing.cpu_count()-1,1)

def make_drip_data_highres(args, data, 
                           mods, ntermMods, ctermMods,
                           varMods, ntermVarMods, ctermVarMods):
    """Generate test data .pfile. and create job scripts for cluster use (if num_jobs > 1).
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go.

       inputs:
       args - output of parsed input arguments (struct)
       data - instance of generator pyFiles.shard_spectra.candidate_spectra_generator, which is a dict with fields:
              "spectra" - spectra in shard
              "target" - target PSMs in shard
              "decoy" - decoy PSMs in shard (if args.decoys set to true)
              "minMz" - minimum mz in dataset
              "maxMz" - maximum mz in dataset
       mods - dict whose keys are amino acids and correspond to offsets
              for modifications denoted in args.mods_spec
       ntermMods - similar as above but for n-terminal mods and args.nterm-peptide-mods-spec
       ctermMods - similar to ntermMods but for the c-terminal


       outputs:
       sids - list of scan IDs for the generated data

       pre:
       - args has been created by parse_args(), directories have been created/checked for existence,
         relevant arguments have been processed (Booleans, mods, digesting enzyme, etc)
       - data has been created by candidate_spectra_generate() and contains the above mentioned fields

       post:
       - args.{mean_file, gauss_file, mixture_file, collection_file} will all be adjusted
       - args.max_mass will be updated to the size of the number of unique theoretical fragmentation locations (floating point if high-res ms2, integers if low-res ms2)
    """

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    # Load the pickled representation of this shard's data
    # prune multiples
    visited_spectra = set([])
    spectra = []
    # spectra = data['spectra']
    for s in data['spectra']:
        if s.spectrum_id not in visited_spectra:
            spectra.append(s)
            visited_spectra.add(s.spectrum_id)

    target = data['target']
    decoy = data['decoy']
    if args.recalibrate:
        recal_decoy = data['recal_decoy']

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # Write file names to
    print >> sys.stderr, 'Writing testing data to %s' % args.output_dir
    
    if(os.path.isfile(os.path.join(args.output_dir, 'inverted-pepDBperSid.txt'))):
        pep_db_per_spectrum = open(args.output_dir+'/inverted-pepDBperSid.txt', 'a')
    else:
        pep_db_per_spectrum = open(args.output_dir+'/inverted-pepDBperSid.txt', 'w')
        pep_db_per_spectrum.write('sid\tnumPeps\n')

    meanFile = args.mean_file
    gaussFile = args.gauss_file
    mixtureFile = args.mixture_file
    collectionFile = args.collection_file

    theo_spec_dict = {}
    ion_dict_per_sid = {}
    numBY_dict_per_sid = {}

    sids = []
    validcharges = args.charges

    for s in spectra:
        preprocess(s)
        ion_dict = {}
        for charge in validcharges:
            if (s.spectrum_id, charge) not in target:
                continue
            tl = target[(s.spectrum_id,charge)]
            dl = decoy[(s.spectrum_id,charge)]

            if args.recalibrate:
                recal_dl = recal_decoy[(s.spectrum_id,charge)]

            max_target_theo_peaks = 0
            max_decoy_theo_peaks = 0
            # check if we're filtering theoretical peaks outside observed m/z values
            if args.filt_theo_peaks:
                if args.per_spectrum_mz_bound:
                    minMz = s.mz[0]
                    maxMz = s.mz[-1]
                else:
                    minMz = args.mz_lb
                    maxMz = args.mz_ub

            # calculate maximum decoy and target theoretical spectra cardinalities
            # serialized data:
            #    peptide[0] = peptide mass (float)
            #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
            #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
            #    peptide[3] = nterm_flanking (character)
            #    peptide[4] = cterm_flanking (character)
            #    peptide[5] = binary string deoting variable modifications
            for tp in tl:
                p = tp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = tp[5][:len(p)]
                    theoSpecKey = p + varModSequence
                    bNy = interleave_b_y_ions_var_mods(Peptide(p), charge, 
                                                       mods, ntermMods, ctermMods,
                                                       varMods, ntermVarMods, ctermVarMods,
                                                       varModSequence)
                else:
                    theoSpecKey = p
                    bNy = interleave_b_y_ions(Peptide(p), charge, 
                                              mods, ntermMods, ctermMods)
                numBY_dict_per_sid[s.spectrum_id, theoSpecKey] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                theo_spec_dict[s.spectrum_id, theoSpecKey] = bNy

                for i in bNy:
                    ion_dict[i] = 1
            for dp in dl:
                d = dp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = dp[5][:len(d)]
                    theoSpecKey = d + varModSequence
                    bNy = interleave_b_y_ions_var_mods(Peptide(d), charge, 
                                                       mods, ntermMods, ctermMods,
                                                       varMods, ntermVarMods, ctermVarMods,
                                                       varModSequence)
                else:
                    theoSpecKey = d
                    bNy = interleave_b_y_ions_var_mods(Peptide(d), charge, 
                                                       mods, ntermMods, ctermMods)

                numBY_dict_per_sid[s.spectrum_id, theoSpecKey] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                theo_spec_dict[s.spectrum_id, theoSpecKey] = bNy
                for i in bNy:
                    ion_dict[i] = 1

            if args.recalibrate:
                for recal_dp in recal_dl:
                    recal_d = recal_dp[1].split('\x00')[0]
                    if varMods or ntermVarMods or ctermVarMods:
                        varModSequence = recal_dp[5][:len(recal_d)]
                        theoSpecKey = recal_d + varModSequence
                        bNy = interleave_b_y_ions_var_mods(Peptide(recal_d), charge, 
                                                           mods, ntermMods, ctermMods,
                                                           varMods, ntermVarMods, ctermVarMods,
                                                           varModSequence)
                    else:
                        theoSpecKey = recal_d
                        bNy = interleave_b_y_ions_var_mods(Peptide(recal_d), charge, 
                                                           mods, ntermMods, ctermMods)
                        
                    numBY_dict_per_sid[s.spectrum_id, theoSpecKey] = len(bNy)
                    if args.filt_theo_peaks:
                        filter_theoretical_peaks(bNy, minMz, maxMz)
                    theo_spec_dict[s.spectrum_id, theoSpecKey] = bNy
                    for i in bNy:
                        ion_dict[i] = 1

        ions = list(ion_dict.iterkeys())
        ions.sort()
        for i, ion in enumerate(ions):
            ion_dict[ion] = i

        if ion_dict:
            ion_dict_per_sid[s.spectrum_id] = ion_dict
            # make collection per spectrum
            args.mean_file = meanFile[:-4] + '-sid' + str(s.spectrum_id) + '.txt'
            args.gauss_file = gaussFile[:-4] + '-sid' + str(s.spectrum_id) + '.txt'
            args.mixture_file = mixtureFile[:-4] + '-sid' + str(s.spectrum_id) + '.txt'
            args.collection_file = collectionFile[:-4] + '-sid' + str(s.spectrum_id) + '.txt'
            make_master_parameters(args, ion_dict, ions)

    # serialize candidate peptide data with fields:
    # (1) - peptide type: 1 = target, 2 = decoy, 3 = recalibrate decoy
    # (2) - peptide string
    # (3) - num b- and y-ions
    # (4) - charge
    # (5) - amino acid flanking n-terminus
    # (6) - amino acid flanking c-terminus
    # (7) - string denoting variable modded amino acids
    s_bin = struct.Struct('I %ds I I s s %ds' % (args.max_length, args.max_length))

    for s in spectra:
        sid = s.spectrum_id

        if sid not in ion_dict_per_sid:
            continue

        # write peptide database to parse and identify GMTK segments later
        pepdb_list = open(args.output_dir+'/sid%d-pepDB.bin' % (s.spectrum_id), "wb")

        # # write peptide database to parse and identify GMTK segments later
        # pepdb_list = open(args.output_dir+'/sid%d-pepDB.txt' % (s.spectrum_id), "w")
        # pepdb_list.write("Kind\tPeptide\tNumBY\tCharge\tN_term\tC_term\tVar_mod_sequence\n")

        preprocess(s)

        ion_dict = ion_dict_per_sid[sid]

        num_candidates = 0
        for charge in validcharges:
            if (sid, charge) in target:
                num_candidates += len(target[sid, charge])
            if (sid, charge) in decoy:
                num_candidates += len(decoy[sid, charge])
            if args.recalibrate:
                if (s.spectrum_id, charge) in recal_decoy:
                    num_candidates += len(recal_decoy[s.spectrum_id, charge])

        peptide_pfile = create_pfile(pfile_dir,
                                     'pep-lengths-sid%d.pfile' % (sid), 
                                     0, 1)
            
        spectrum_pfile = create_pfile(pfile_dir,
                                      'spectrum-sid%d.pfile' % (sid),
                                      2,0)

        # write observation file; assume we've filtered spectra of interest before this,
        # so that we don't have to check whether we should create this observation file.
        # Though, for individual peptide observations, we still don't know what relevant 
        # charges are being searched
        drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)

        pep_num = 0 # total peptide decision trees being written
        # pep_dt = open(args.output_dir+'/temp_iterable.dts', "w")
        pep_dt = open(os.path.join(args.output_dir, 'sid%d-iterable.dts' % (sid)), "w")
        pep_dt.write('%d\n\n' % (num_candidates))

        ran = 0

        for charge in validcharges:
            if (sid, charge) not in target:
                continue
            ran = 1

            # Generate the peptide-spectrum matches.
            tl = target[(sid,charge)]
            dl = decoy[(sid,charge)]
            if args.recalibrate:
                recal_dl = recal_decoy[(sid, charge)]

            for tp in tl:
                pepType = 1
                p = tp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = tp[5][:len(p)]
                    theoSpecKey = p + varModSequence
                else:
                    varModSequence = ''.join(['0'] * len(p))
                    theoSpecKey = p
                bNy = theo_spec_dict[sid, theoSpecKey]
                bNy = [ion_dict[bOrY] for bOrY in bNy]
                drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                      pep_num, sid, args.max_obs_mass,
                                      peptide_pfile, True, len(bNy))

                curr_db_pep = (pepType, tp[1], numBY_dict_per_sid[sid, theoSpecKey], charge,
                               tp[3], tp[4], varModSequence)
                pepdb_list.write(s_bin.pack(*curr_db_pep))
                # pepdb_list.write("t\t%s\t%d\t%d\t%c\t%c\t%s\n" % (p,
                #                                                   numBY_dict_per_sid[sid, theoSpecKey], 
                #                                                   charge,
                #                                                   tp[3], tp[4], varModSequence))
                pep_num += 1

            for dp in dl:
                pepType = 2
                d = dp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = dp[5][:len(d)]
                    theoSpecKey = d + varModSequence
                else:
                    varModSequence = ''.join(['0'] * len(d))
                    theoSpecKey = d
                bNy = theo_spec_dict[sid, theoSpecKey]
                bNy = [ion_dict[bOrY] for bOrY in bNy]
                drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                      pep_num, sid, args.max_obs_mass,
                                      peptide_pfile, False, len(bNy))
                curr_db_pep = (pepType, dp[1], numBY_dict_per_sid[sid, theoSpecKey], charge,
                               dp[3], dp[4], varModSequence)
                pepdb_list.write(s_bin.pack(*curr_db_pep))
                # pepdb_list.write("d\t%s\t%d\t%d\t%c\t%c\t%s\n" % (d, 
                #                                                   numBY_dict_per_sid[sid, theoSpecKey], 
                #                                                   charge,
                #                                                   dp[3], dp[4], varModSequence))
                # write observation file
                # drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                pep_num += 1

            if args.recalibrate:
                for recal_dp in recal_dl:
                    pepType = 3
                    recal_d = recal_dp[1].split('\x00')[0]
                    if varMods or ntermVarMods or ctermVarMods:
                        varModSequence = recal_dp[5][:len(recal_d)]
                        theoSpecKey = recal_d + varModSequence
                    else:
                        varModSequence = ''.join(['0'] * len(recal_d))
                        theoSpecKey = recal_d
                    bNy = theo_spec_dict[sid, theoSpecKey]
                    bNy = [ion_dict[bOrY] for bOrY in bNy]
                    drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                          pep_num, sid, args.max_obs_mass,
                                          peptide_pfile, False, len(bNy))
                    curr_db_pep = (pepType, recal_dp[1], numBY_dict_per_sid[sid, theoSpecKey], charge,
                                   recal_dp[3], recal_dp[4], varModSequence)
                    pepdb_list.write(s_bin.pack(*curr_db_pep))
                    pep_num += 1

        # close streams for this spectrum
        pep_dt.close()
        # close peptide database stream for this spectrum
        pepdb_list.close()
        # compile dt using gmtkDTIndex
        call(['gmtkDTindex', '-decisionTreeFiles', os.path.join(args.output_dir,'sid%d-iterable.dts' % (s.spectrum_id))], 
             stdout = stdo, stderr = stde)
                
        # write pepDB info for later parsing
        pep_db_per_spectrum.write('%d\t%d\n' % (s.spectrum_id, pep_num))

        if ran: # guard against adding multiply charged spectra more than once
            sids.append(s)

    # close streams openend per ms2 file
    pep_db_per_spectrum.close()
    return sids

def make_drip_data_lowres(args, data, 
                          mods, ntermMods, ctermMods,
                          varMods, ntermVarMods, ctermVarMods):
    """Generate test data .pfile. and create job scripts for cluster use.
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go
    """

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    # Load the pickled representation of this shard's data
    # prune multiples
    visited_spectra = set([])
    spectra = []
    # spectra = data['spectra']
    for s in data['spectra']:
        if s.spectrum_id not in visited_spectra:
            spectra.append(s)
            visited_spectra.add(s.spectrum_id)
        
    target = data['target']
    decoy = data['decoy']
    if args.recalibrate:
        recal_decoy = data['recal_decoy']

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # Write file names to
    print >> sys.stderr, 'Writing testing data to %s' % args.output_dir
    
    if(os.path.isfile(os.path.join(args.output_dir, 'inverted-pepDBperSid.txt'))):
        pep_db_per_spectrum = open(args.output_dir+'/inverted-pepDBperSid.txt', 'a')
    else:
        pep_db_per_spectrum = open(args.output_dir+'/inverted-pepDBperSid.txt', 'w')
        pep_db_per_spectrum.write('sid\tnumPeps\n')

    dripMeans = load_drip_means(args.learned_means)

    meanFile = args.mean_file
    gaussFile = args.gauss_file
    mixtureFile = args.mixture_file
    collectionFile = args.collection_file

    make_master_parameters_lowres(args, dripMeans)

    # serialize candidate peptide data with fields:
    # (1) - peptide type: 1 = target, 2 = decoy, 3 = recalibrate decoy
    # (2) - peptide string
    # (3) - num b- and y-ions
    # (4) - charge
    # (5) - amino acid flanking n-terminus
    # (6) - amino acid flanking c-terminus
    # (7) - string denoting variable modded amino acids
    s_bin = struct.Struct('I %ds I I s s %ds' % (args.max_length, args.max_length))

    sids = []
    validcharges = args.charges
    for s in spectra:
        preprocess(s)
        ran = 0
        if args.filt_theo_peaks:
            if args.per_spectrum_mz_bound:
                minMz = s.mz[0]
                maxMz = s.mz[-1]
            else:
                minMz = args.mz_lb
                maxMz = args.mz_ub

        # write peptide database to parse and identify GMTK segments later
        pepdb_list = open(args.output_dir+'/sid%d-pepDB.bin' % (s.spectrum_id), "wb")

        # pepdb_list = open(args.output_dir+'/sid%d-pepDB.txt' % (s.spectrum_id), "w")
        # # pepdb_list.write("Kind\tPeptide\tNumBY\tCharge\n")
        # pepdb_list.write("Kind\tPeptide\tNumBY\tCharge\tN_term\tC_term\tVar_mod_sequence\n")

        # pfile for PSM observations: number of theoretical peaks
        peptide_pfile = create_pfile(pfile_dir,
                                     'pep-lengths-sid%d.pfile' % (s.spectrum_id), 
                                     0, 1)

        spectrum_pfile = create_pfile(pfile_dir,
                                      'spectrum-sid%d.pfile' % (s.spectrum_id),
                                      2,0)

        # write observation file; assume we've filtered spectra of interest before this,
        # so that we don't have to check whether we should create this observation file.
        # Though, for individual peptide observations, we still don't know what relevant 
        # charges are being searched
        drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)

        num_candidates = 0
        for charge in validcharges:
            if (s.spectrum_id, charge) in target:
                num_candidates += len(target[s.spectrum_id, charge])
            if (s.spectrum_id, charge) in decoy:
                num_candidates += len(decoy[s.spectrum_id, charge])
            if args.recalibrate:
                if (s.spectrum_id, charge) in recal_decoy:
                    num_candidates += len(recal_decoy[s.spectrum_id, charge])                

        pep_num = 0 # total peptide decision trees being written
        # pep_dt = open(args.output_dir+'/temp_iterable.dts', "w")
        pep_dt = open(os.path.join(args.output_dir,'sid%d-iterable.dts' % (s.spectrum_id)), "w")
        pep_dt.write('%d\n\n' % (num_candidates))

        for charge in validcharges:
            if (s.spectrum_id, charge) not in target:
                continue
            ran = 1
        # Generate the peptide-spectrum matches.
            tl = target[(s.spectrum_id,charge)]
            dl = decoy[(s.spectrum_id,charge)]

            if args.recalibrate:
                recal_dl = recal_decoy[(s.spectrum_id,charge)]

            # serialized data:
            #    peptide[0] = peptide mass (float)
            #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
            #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
            #    peptide[3] = nterm_flanking (character)
            #    peptide[4] = cterm_flanking (character)
            #    peptide[5] = binary string deoting variable modifications
            for tp in tl:
                pepType = 1
                p = tp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = tp[5][:len(p)]
                    theoSpecKey = p + varModSequence
                    bNy = interleave_b_y_ions_var_mods_lowres(Peptide(p), charge, 
                                                              mods, ntermMods, ctermMods,
                                                              varMods, ntermVarMods, ctermVarMods,
                                                              varModSequence)
                else:
                    varModSequence = ''.join(['0'] * len(p))
                    theoSpecKey = p
                    bNy = interleave_b_y_ions_lowres(Peptide(p), charge, 
                                                     mods, ntermMods, ctermMods)
                # numBY for DRIP features assumed all b-/y-ions, not just those
                # unfiltered per spectrum
                curr_db_pep = (pepType, tp[1], len(bNy), charge,
                               tp[3], tp[4], varModSequence)
                pepdb_list.write(s_bin.pack(*curr_db_pep))
                # pepdb_list.write("t\t%s\t%d\t%d\t%c\t%c\t%s\n" % (p,
                #                                                   len(bNy), 
                #                                                   charge,
                #                                                   tp[3], tp[4], varModSequence))
                # pepdb_list.write("t\t%s\t%d\t%d\n" % (p, len(bNy), charge))
                if args.filt_theo_peaks:
                    filter_theoretical_peaks_lowres(bNy, dripMeans,
                                                    minMz, maxMz)
                drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                      pep_num, s.spectrum_id, args.max_obs_mass,
                                      peptide_pfile, True, len(bNy))
                pep_num += 1

            for dp in dl:
                pepType = 2
                d = dp[1].split('\x00')[0]
                if varMods or ntermVarMods or ctermVarMods:
                    varModSequence = dp[5][:len(d)]
                    theoSpecKey = d + varModSequence
                    bNy = interleave_b_y_ions_var_mods_lowres(Peptide(d), charge, 
                                                              mods, ntermMods, ctermMods,
                                                              varMods, ntermVarMods, ctermVarMods,
                                                              varModSequence)
                else:
                    varModSequence = ''.join(['0'] * len(d))
                    theoSpecKey = d
                    bNy = interleave_b_y_ions_lowres(Peptide(d), charge, 
                                                     mods, ntermMods, ctermMods)
                # numBY for DRIP features assumed all b-/y-ions, not just those
                # unfiltered per spectrum
                curr_db_pep = (pepType, dp[1], len(bNy), charge,
                               dp[3], dp[4], varModSequence)
                pepdb_list.write(s_bin.pack(*curr_db_pep))

                # pepdb_list.write("d\t%s\t%d\t%d\t%c\t%c\t%s\n" % (d,
                #                                                   len(bNy), 
                #                                                   charge,
                #                                                   dp[3], dp[4], varModSequence))
                # pepdb_list.write("d\t%s\t%d\t%d\n" % (d, len(bNy), charge))
                if args.filt_theo_peaks:
                    filter_theoretical_peaks_lowres(bNy, dripMeans,
                                                    minMz, maxMz)
                drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                      pep_num, s.spectrum_id, args.max_obs_mass,
                                      peptide_pfile, False, len(bNy))
                # write observation file
                # drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                pep_num += 1

            if args.recalibrate:
                for recal_dp in recal_dl:
                    pepType = 3
                    recal_d = recal_dp[1].split('\x00')[0]
                    if varMods or ntermVarMods or ctermVarMods:
                        varModSequence = recal_dp[5][:len(recal_d)]
                        theoSpecKey = recal_d + varModSequence
                        bNy = interleave_b_y_ions_var_mods_lowres(Peptide(recal_d), charge, 
                                                                  mods, ntermMods, ctermMods,
                                                                  varMods, ntermVarMods, ctermVarMods,
                                                                  varModSequence)
                    else:
                        varModSequence = ''.join(['0'] * len(recal_d))
                        theoSpecKey = recal_d
                        bNy = interleave_b_y_ions_lowres(Peptide(recal_d), charge, 
                                                         mods, ntermMods, ctermMods)
                    # numBY for DRIP features assumed all b-/y-ions, not just those
                    # unfiltered per spectrum
                    curr_db_pep = (pepType, recal_dp[1], len(bNy), charge,
                                   recal_dp[3], recal_dp[4], varModSequence)
                    pepdb_list.write(s_bin.pack(*curr_db_pep))

                    if args.filt_theo_peaks:
                        filter_theoretical_peaks_lowres(bNy, dripMeans,
                                                        minMz, maxMz)
                    drip_peptide_sentence(pep_dt, theoSpecKey, bNy, 
                                          pep_num, s.spectrum_id, args.max_obs_mass,
                                          peptide_pfile, False, len(bNy))
                # write observation file
                # drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                    pep_num += 1

        pepdb_list.close()
        pep_dt.close()
        # compile dt using gmtkDTIndex
        call(['gmtkDTindex', '-decisionTreeFiles', os.path.join(args.output_dir,'sid%d-iterable.dts' % (s.spectrum_id))], 
             stdout = stdo, stderr = stde)
                
        # write pepDB info for later parsing
        pep_db_per_spectrum.write('%d\t%d\n' % (s.spectrum_id, pep_num))

        if ran: # guard against adding multiply charged spectra more than once
            sids.append(s)            

    # close streams openend per ms2 file
    pep_db_per_spectrum.close()
    return sids

def runDrip(args):    
    """ Run drip once per spectrum, collapsing all charge-varying candidates into a single GMTK call
    """
    # create constant gmtkViterbi command line string
    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -ckbeam " + str(args.beam) \
        + " -sdiffact1 rl -sdiffact2 rl -fdiffact2 rl" \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"

    # # formulate digestion regular expression
    # digest_re = df.create_digest_re(args.enzyme, args.custom_enzyme)
    # print "digest regular expression: %s" % digest_re
    # r = re.compile(digest_re)

    # parse modifications
    mods, var_mods = parse_var_mods(args.mods_spec, True)
    nterm_mods, nterm_var_mods = parse_var_mods(args.nterm_peptide_mods_spec, False)
    cterm_mods, cterm_var_mods = parse_var_mods(args.cterm_peptide_mods_spec, False)

    meanFile = args.mean_file
    gaussFile = args.gauss_file
    mixtureFile = args.mixture_file
    collectionFile = args.collection_file

    spec_eval = []
    # create GMTK observation files
    if args.num_jobs > 1:
        # pickle_candidate_spectra(args)

        # pickle all data and arguments for cluster use
        print "Cluster usage selected with %d jobs, disabling multithreading" % args.num_jobs
        args.num_threads = 1
        pickle_candidate_binarydb_spectra(args,
                                          mods, nterm_mods, cterm_mods,
                                          var_mods, nterm_var_mods, cterm_var_mods)
        exit(0)
    else:
        # iterate through what would be pickles and create pfiles
        # for data in candidate_spectra_generator(args, r,
        #                                         mods, 
        #                                         nterm_mods, cterm_mods):
        for data in candidate_binarydb_spectra_generator(args,
                                                         mods, nterm_mods, cterm_mods,
                                                         var_mods, nterm_var_mods, cterm_var_mods):
            args.mz_lb = data['minMz']
            args.mz_ub = data['maxMz']
            if args.high_res_ms2:
                spec_eval += make_drip_data_highres(args, data, 
                                                    mods, nterm_mods, cterm_mods,
                                                    var_mods, nterm_var_mods, cterm_var_mods)
            else:
                spec_eval += make_drip_data_lowres(args, data, 
                                                   mods, nterm_mods, cterm_mods,
                                                   var_mods, nterm_var_mods, cterm_var_mods)

    pepDBperSidList = os.path.join(args.output_dir, 'inverted-pepDBperSid.txt')
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # create structure and master files then triangulate
    try:
        create_drip_structure(args.high_res_ms2, args.structure_file, 
                              args.max_obs_mass)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(args.high_res_ms2, args.master_file, 
                           args.max_obs_mass,
                           "DRIP_MZ",
                           "drip_collection/covar.txt",
                           "DRIP_GAUSSIAN_COMPONENTS",
                           "DRIP_GAUSSIAN_MIXTURES",
                           "DRIP_MZ_GAUSSIANS")
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        # triangulate_drip(args.structure_file, args.master_file, stdo, stde)
        triangulate_drip(args.structure_file)
    except:
        print "Could not create triangulate structure file %s, exitting" % args.structure_file
        exit(-1)


    try:
        write_covar_file(args.high_res_ms2, args.covar_file, 
                         args.learned_covars)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    # run GMTK
    spec_dict = {} # create dictionary of spectra with sid as key, used for calculating DRIP features
    for s in spec_eval:
        sid = s.spectrum_id
        spec_dict[sid] = s
        currSid = str(sid)
        # first compile iterable dt
        dtFile = os.path.join(args.output_dir, 'sid%d-iterable.dts' % (sid))
        # call(['cp', dtFile, 'iterable.dts'])
        # call(['gmtkDTindex', '-decisionTreeFiles', 'iterable.dts'], 
        #      stdout = stdo, stderr = stde)
        if args.high_res_ms2:
            mean_file = meanFile[:-4] + '-sid' + currSid  + '.txt'
            gauss_file = gaussFile[:-4] + '-sid' + currSid  + '.txt'
            mixture_file = mixtureFile[:-4] + '-sid' + currSid + '.txt'
            collection_file = collectionFile[:-4] + '-sid' + currSid + '.txt'
        else:
            mean_file = meanFile
            gauss_file = gaussFile
            mixture_file = mixtureFile
            collection_file = collectionFile

        cppCommand = '\'-DITERABLE_DT=' + dtFile \
            + ' -DDRIP_MZ=' + mean_file \
            + ' -DDRIP_GAUSSIAN_COMPONENTS=' + gauss_file \
            + ' -DDRIP_GAUSSIAN_MIXTURES=' + mixture_file \
            + ' -DDRIP_MZ_GAUSSIANS=' + collection_file \
            + '\''

        # call gmtkViterbi
        # gmtkViterbi command line
        vitStr = vitStr0 + ' -vitValsFile ' + args.logDir + '/vitVals-sid' + currSid + '.txt' \
            + ' -of1 ' + pfile_dir + '/spectrum-sid' + currSid + '.pfile' \
            + ' -of2 ' + pfile_dir + '/pep-lengths-sid' + currSid + '.pfile' \
            + ' -cppCommand ' + cppCommand
        call(shlex.split(vitStr), stdout = stdo, stderr = stde)

    if args.recalibrate:
        psm.write_dripSearch_ident(args.output+ '-uncal.txt', args.logDir, 
                                   args.top_match, args.high_res_ms2, 
                                   spec_dict, meanFile,
                                   pepDBperSidList, args.output_dir,
                                   args.max_length, args.recalibrate,
                                   mods, nterm_mods, cterm_mods,
                                   var_mods, nterm_var_mods, cterm_var_mods)

        chargeRecalibrate(args.output+ '-uncal.txt',
                          args.output+ '.txt',
                          args.top_match)
    else:
        psm.write_dripSearch_ident(args.output+ '.txt', args.logDir, 
                                   args.top_match, args.high_res_ms2, 
                                   spec_dict, meanFile,
                                   pepDBperSidList, args.output_dir,
                                   args.max_length, args.recalibrate,
                                   mods, nterm_mods, cterm_mods,
                                   var_mods, nterm_var_mods, cterm_var_mods)

def runDripCluster(args):    
    """ Run drip for prepared pickle, collapsing all charge-varying candidates into a single GMTK call
        Assumed cluster usage
    """
    # create constant gmtkViterbi command line string
    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -ckbeam " + str(args.beam) \
        + " -sdiffact1 rl -sdiffact2 rl -fdiffact2 rl" \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"

    spec_eval = []
    # create GMTK observation files

    # load pickled data
    data = pickle.load(open(args.spectra))

    # load options used to pickle data
    copyArgs(args, data['args'])

    meanFile = args.mean_file
    gaussFile = args.gauss_file
    mixtureFile = args.mixture_file
    collectionFile = args.collection_file

    args.mz_lb = data['minMz']
    args.mz_ub = data['maxMz']
    mods = data['mods']
    nterm_mods = data['nterm_mods']
    cterm_mods = data['cterm_mods']
    var_mods = data['var_mods']
    nterm_var_mods = data['nterm_var_mods']
    cterm_var_mods = data['cterm_var_mods']

    if args.high_res_ms2:
        spec_eval += make_drip_data_highres(args, data, 
                                            mods, nterm_mods, cterm_mods,
                                            var_mods, nterm_var_mods, cterm_var_mods)
    else:
        spec_eval += make_drip_data_lowres(args, data, 
                                           mods, nterm_mods, cterm_mods,
                                           var_mods, nterm_var_mods, cterm_var_mods)

    pepDBperSidList = os.path.join(args.output_dir, 'inverted-pepDBperSid.txt')
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # create structure and master files then triangulate
    try:
        create_drip_structure(args.high_res_ms2, args.structure_file, 
                              args.max_obs_mass)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(args.high_res_ms2, args.master_file, 
                           args.max_obs_mass,
                           "DRIP_MZ",
                           "drip_collection/covar.txt",
                           "DRIP_GAUSSIAN_COMPONENTS",
                           "DRIP_GAUSSIAN_MIXTURES",
                           "DRIP_MZ_GAUSSIANS")
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        # triangulate_drip(args.structure_file, args.master_file, stdo, stde)
        triangulate_drip(args.structure_file)
    except:
        print "Could not create triangulate structure file %s, exitting" % args.structure_file
        exit(-1)


    try:
        write_covar_file(args.high_res_ms2, args.covar_file, 
                         args.learned_covars)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    # run GMTK
    spec_dict = {} # create dictionary of spectra with sid as key, used for calculating DRIP features
    for s in spec_eval:
        sid = s.spectrum_id
        spec_dict[sid] = s
        currSid = str(sid)
        # first compile iterable dt
        dtFile = os.path.join(args.output_dir, 'sid%d-iterable.dts' % (sid))
        # call(['cp', dtFile, 'iterable.dts'])
        # call(['gmtkDTindex', '-decisionTreeFiles', 'iterable.dts'], 
        #      stdout = stdo, stderr = stde)
        if args.high_res_ms2:
            mean_file = meanFile[:-4] + '-sid' + currSid  + '.txt'
            gauss_file = gaussFile[:-4] + '-sid' + currSid  + '.txt'
            mixture_file = mixtureFile[:-4] + '-sid' + currSid + '.txt'
            collection_file = collectionFile[:-4] + '-sid' + currSid + '.txt'
        else:
            mean_file = meanFile
            gauss_file = gaussFile
            mixture_file = mixtureFile
            collection_file = collectionFile

        cppCommand = '\'-DITERABLE_DT=' + dtFile \
            + ' -DDRIP_MZ=' + mean_file \
            + ' -DDRIP_GAUSSIAN_COMPONENTS=' + gauss_file \
            + ' -DDRIP_GAUSSIAN_MIXTURES=' + mixture_file \
            + ' -DDRIP_MZ_GAUSSIANS=' + collection_file \
            + '\''

        # call gmtkViterbi
        # gmtkViterbi command line
        vitStr = vitStr0 + ' -vitValsFile ' + args.logDir + '/vitVals-sid' + currSid + '.txt' \
            + ' -of1 ' + pfile_dir + '/spectrum-sid' + currSid + '.pfile' \
            + ' -of2 ' + pfile_dir + '/pep-lengths-sid' + currSid + '.pfile' \
            + ' -cppCommand ' + cppCommand
        call(shlex.split(vitStr), stdout = stdo, stderr = stde)

    # cluster use assumes multiple jobs, run recalibration over all resulting files
    psm.write_dripSearch_ident(args.output+ '.txt', args.logDir, 
                               args.top_match, args.high_res_ms2, 
                               spec_dict, meanFile,
                               pepDBperSidList, args.output_dir,
                               args.max_length, args.recalibrate,
                               mods, nterm_mods, cterm_mods,
                               var_mods, nterm_var_mods, cterm_var_mods)

def chargeRecalibrate(dripIdent,
                      output_file,
                      topPsms = 1, 
                      percentile = 0.99999999):

    if os.path.isdir(dripIdent):
        # dripIdents = [dripIdent + '/' + f for f in os.listdir(dripIdent)
        #               if any(f.endswith('txt'))]
        dripIdents = list(glob.glob(os.path.join(dripIdent, 'split*-ident.txt')))
    else:
        dripIdents = [dripIdent]

    training_decoys_by_charge = {}
    psms_by_charge = {}

    for ident in dripIdents:
        load_drip_psms(ident, training_decoys_by_charge, 
                       psms_by_charge)

    # check whether recalibration decoys were searched
    do_recalibration = True
    for c in psms_by_charge:
        if c not in training_decoys_by_charge:
            do_recalibration = False
            break

    if do_recalibration:
        recalibrate_charge_psms(output_file, psms_by_charge,
                                training_decoys_by_charge, 
                                percentile, topPsms)
    else:
        # write PSMs to output file
        targets = []
        decoys = []
        targetf = lambda r: r.kind == "t"
        decoyf = lambda r: r.kind == "d"
        for c in psms_by_charge:
            targets += itertools.ifilter(targetf, psms_by_charge[c])
            decoys += itertools.ifilter(decoyf, psms_by_charge[c])

        targets.sort(key = lambda r: r.scan)
        decoys.sort(key = lambda r: r.scan)

        assert len(targetf) == len(decoyf)
        write_drip_recal(output_file, targets, decoys)

if __name__ == '__main__':
    # read in options and process input arguments
    args = parseInputOptions()
    process_args(args)
    
    if args.cluster_mode:
        runDripCluster(args)
        exit(0)

    if args.num_threads <= 1 or args.num_jobs > 1:
        # more efficient to call runDrip and redefine it below when performing multithreading,
        # rather than defining it below only (per the restriction of using the multiprocessing module
        # within main) and setting a pool of size 1 for single-threading.  The major drawback is needing
        # make sure both implementations are up to speed
        runDrip(args)
        sys.exit(0)

    # otherwise, redefine runDrip in main() so that multiprocessing runs properly
    # create constant gmtkViterbi command line string

    ####################################################################################
    ##################################### start runDrip for multithreading
    ####################################################################################

    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -ckbeam " + str(args.beam) \
        + " -sdiffact1 rl -sdiffact2 rl -fdiffact2 rl" \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"

    # # formulate digestion regular expression
    # digest_re = df.create_digest_re(args.enzyme, args.custom_enzyme)
    # print "digest regular expression: %s" % digest_re
    # r = re.compile(digest_re)

    # parse modifications
    mods, var_mods = parse_var_mods(args.mods_spec, True)
    # print "mods:"
    # print mods
    nterm_mods, nterm_var_mods = parse_var_mods(args.nterm_peptide_mods_spec, False)
    # print "n-term mods:"
    # print nterm_mods
    cterm_mods, cterm_var_mods = parse_var_mods(args.cterm_peptide_mods_spec, False)
    # print "c-term mods:"
    # print cterm_mods

    meanFile = args.mean_file
    gaussFile = args.gauss_file
    mixtureFile = args.mixture_file
    collectionFile = args.collection_file

    spec_eval = []
    # create GMTK observation files
    if args.num_jobs > 1:
        # pickle_candidate_spectra(args)

        # pickle all data and arguments for cluster use
        print "Cluster usage selected with %d jobs, disabling multithreading" % args.num_jobs
        args.num_threads = 1
        pickle_candidate_binarydb_spectra(args,
                                          mods, nterm_mods, cterm_mods,
                                          var_mods, nterm_var_mods, cterm_var_mods)
        exit(0)
    else:
        # iterate through what would be pickles and create pfiles
        # for data in candidate_spectra_generator(args, r,
        #                                         mods, 
        #                                         nterm_mods, cterm_mods):
        for data in candidate_binarydb_spectra_generator(args,
                                                         mods, nterm_mods, cterm_mods,
                                                         var_mods, nterm_var_mods, cterm_var_mods):
            args.mz_lb = data['minMz']
            args.mz_ub = data['maxMz']
            if args.high_res_ms2:
                spec_eval += make_drip_data_highres(args, data, 
                                                    mods, nterm_mods, cterm_mods, 
                                                    var_mods, nterm_var_mods, cterm_var_mods)
            else:
                spec_eval += make_drip_data_lowres(args, data, 
                                                   mods, nterm_mods, cterm_mods,
                                                   var_mods, nterm_var_mods, cterm_var_mods)

    pepDBperSidList = os.path.join(args.output_dir, 'inverted-pepDBperSid.txt')
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # create structure and master files then triangulate
    try:
        create_drip_structure(args.high_res_ms2, args.structure_file, 
                              args.max_obs_mass)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(args.high_res_ms2, args.master_file, 
                           args.max_obs_mass,
                           "DRIP_MZ",
                           "drip_collection/covar.txt",
                           "DRIP_GAUSSIAN_COMPONENTS",
                           "DRIP_GAUSSIAN_MIXTURES",
                           "DRIP_MZ_GAUSSIANS")
        
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        triangulate_drip(args.structure_file)
    except:
        print "Could not create triangulate structure file %s, exitting" % args.structure_file
        exit(-1)


    try:
        write_covar_file(args.high_res_ms2, args.covar_file,
                         args.learned_covars)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    # set up pool for multithreading
    pool = multiprocessing.Pool(processes = args.num_threads)
    
    # add GMTK jobs to pool
    spec_dict = {} # create dictionary of spectra with sid as key, used for calculating DRIP features
    for s in spec_eval:
        sid = s.spectrum_id
        spec_dict[sid] = s
        currSid = str(sid)
        # first compile iterable dt
        dtFile = os.path.join(args.output_dir, 'sid%d-iterable.dts' % (sid))
        if args.high_res_ms2:
            mean_file = meanFile[:-4] + '-sid' + currSid  + '.txt'
            gauss_file = gaussFile[:-4] + '-sid' + currSid  + '.txt'
            mixture_file = mixtureFile[:-4] + '-sid' + currSid + '.txt'
            collection_file = collectionFile[:-4] + '-sid' + currSid + '.txt'
        else:
            mean_file = meanFile
            gauss_file = gaussFile
            mixture_file = mixtureFile
            collection_file = collectionFile

        cppCommand = '\'-DITERABLE_DT=' + dtFile \
            + ' -DDRIP_MZ=' + mean_file \
            + ' -DDRIP_GAUSSIAN_COMPONENTS=' + gauss_file \
            + ' -DDRIP_GAUSSIAN_MIXTURES=' + mixture_file \
            + ' -DDRIP_MZ_GAUSSIANS=' + collection_file \
            + '\''

        # call gmtkViterbi
        # gmtkViterbi command line
        vitStr = vitStr0 + ' -vitValsFile ' + args.logDir + '/vitVals-sid' + currSid + '.txt' \
            + ' -of1 ' + pfile_dir + '/spectrum-sid' + currSid + '.pfile' \
            + ' -of2 ' + pfile_dir + '/pep-lengths-sid' + currSid + '.pfile' \
            + ' -cppCommand ' + cppCommand
        # add to pool
        r = pool.apply_async(check_output, args=(shlex.split(vitStr),))
    pool.close()
    pool.join()

    ####################################################################################
    ##################################### end runDrip for multithreading
    ####################################################################################

    if args.recalibrate:
        psm.write_dripSearch_ident(args.output+ '-uncal.txt', args.logDir, 
                                   args.top_match, args.high_res_ms2, 
                                   spec_dict, meanFile,
                                   pepDBperSidList, args.output_dir,
                                   args.max_length, args.recalibrate,
                                   mods, nterm_mods, cterm_mods,
                                   var_mods, nterm_var_mods, cterm_var_mods)

        chargeRecalibrate(args.output+ '-uncal.txt',
                          args.output+ '.txt',
                          args.top_match)
    else:
        psm.write_dripSearch_ident(args.output+ '.txt', args.logDir, 
                                   args.top_match, args.high_res_ms2, 
                                   spec_dict, meanFile,
                                   pepDBperSidList, args.output_dir,
                                   args.max_length, args.recalibrate,
                                   mods, nterm_mods, cterm_mods,
                                   var_mods, nterm_var_mods, cterm_var_mods)


    if stdo:
        stdo.close()
    if stde:
        stde.close()
