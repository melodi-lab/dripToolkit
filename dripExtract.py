#!/usr/bin/env python
#
# Copyright 2015 <fill in later>

from __future__ import with_statement

__authors__ = ['John Halloran <halloj3@uw.edu>' ]

# todo: add count of number missed cleavages to features

import os
import re
import sys
import argparse
import cPickle as pickle
import pyFiles.digest_fasta as df
import multiprocessing
import shlex
import math
import pyFiles.psm as psm

from random import shuffle
from shutil import rmtree
from pyFiles.spectrum import MS2Spectrum
from pyFiles.peptide import Peptide, amino_acids_to_indices
from pyFiles.normalize import pipeline
from pyFiles.pfile.wrapper import PFile
# from pyFiles.args import make_dir_callback
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
                                  interleave_b_y_ions_lowres,
                                  filter_theoretical_peaks,
                                  filter_theoretical_peaks_lowres,
                                  load_drip_means)

from pyFiles.constants import allPeps
from pyFiles.ioPsmFunctions import load_psms, load_pin_file, load_pin_return_dict
from pyFiles.shard_spectra import (calc_minMaxMz,
                                   pickle_candidate_spectra, 
                                   candidate_spectra_generator,
                                   candidate_spectra_memeffic_generator,
                                   load_spectra_ret_dict,
                                   load_target_decoy_db_no_reshuffle)

from pyFiles.process_vitvals import parse_dripExtract, write_psm_ins_dels
from subprocess import call, check_output

debug=1

# set stdout and stderr for subprocess
# stdo = open(os.devnull, "w")
# # set stdout for all default streams
# sys.stdout = open("dripExtract_output", "w")

def process_args(args):
    """ Check whether relevant directories already exist (remove them if so), create
        necessary directories.  Adjust parameters.

        pre:
        - args has been created by calling parse_args()
        
        post:
        - Relevant directories will be first removed then created, boolean parameters
          will be adjusted to reflect their boolean command line strings
    """
    # check input arguments, create necessary directories, 
    # set necessary global variables based on input arguments
    # constant string
    if args.max_obs_mass < 0:
        print "Supplied max-mass %d, must be greater than 0, exitting" % args.max_obs_mass
        exit(-1)
    # check arguments
    assert(args.shards > 0)
    if not args.output_file or not args.fasta:
        parser.print_help()
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

    if not args.filter_ident:
        args.ident = ''

    # set true or false strings to booleans
    args.load_peptide_database_file = df.check_arg_trueFalse(args.load_peptide_database_file)
    args.write_pin = df.check_arg_trueFalse(args.write_pin)
    args.append_to_pin = df.check_arg_trueFalse(args.append_to_pin)
    args.decoys = df.check_arg_trueFalse(args.decoys)
    args.monoisotopic_precursor = df.check_arg_trueFalse(args.monoisotopic_precursor)
    args.high_res_ms2 = df.check_arg_trueFalse(args.high_res_ms2)
    args.randomize_ms2_spectra = df.check_arg_trueFalse(args.randomize_ms2_spectra)

    # check precursor mass type
    pmt = args.precursor_window_type.lower()
    if pmt == 'da':
        args.ppm = False
    else:
        args.ppm = True

    # make sure number of input threads does not exceed number of supported threads
    if args.num_threads > multiprocessing.cpu_count():
        args.num_threads = max(multiprocessing.cpu_count()-1,1)

def make_drip_data_highres(args, spectra, stdo, stde):
    """Generate test data .pfile. and create job scripts for cluster use (if num_jobs > 1).
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go.

       inputs:
       args - output of parsed input arguments (struct)

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
    # parse modifications
    mods = df.parse_mods(args.mods_spec, True)
    print "mods:"
    print mods
    ntermMods = df.parse_mods(args.nterm_peptide_mods_spec, False)
    print "n-term mods:"
    print ntermMods
    ctermMods = df.parse_mods(args.cterm_peptide_mods_spec, False)
    print "c-term mods:"
    print ctermMods

    if not args.append_to_pin:
        target,decoy,num_psms = load_psms(args.psm_file)
    else:
        target,decoy,num_psms = load_pin_file(args.psm_file)
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    # assume that we should randomize PSMs for multithreading purposes; only reason
    # why we are currently assuming this is that there is already a parameter for dripSearch
    # which signifies whether we should shuffle the data
    shuffle(sid_charges)

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    validcharges = args.charges

    ion_dict = {} # global dictionary for used fragment ions
    theo_spec_dict = {}
    numBY_dict_per_sid = {}
    # construct ion_dict
    for sid in spectra:
        s = spectra[sid]
        preprocess(s)
        for charge in validcharges:
            if (s.spectrum_id, charge) not in target:
                continue
            # check if we're filtering theoretical peaks outside observed m/z values
            if args.filt_theo_peaks:
                if args.per_spectrum_mz_bound:
                    minMz = s.mz[0]
                    maxMz = s.mz[-1]
                else:
                    minMz = args.mz_lb
                    maxMz = args.mz_ub

            # calculate maximum decoy and target theoretical spectra cardinalities
            for p in target[s.spectrum_id, charge]:
                pep = p.peptide
                bNy = interleave_b_y_ions(Peptide(pep), charge, mods,
                                          ntermMods, ctermMods)
                numBY_dict_per_sid[sid, pep] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                theo_spec_dict[s.spectrum_id, pep] = bNy

                for i in bNy:
                    ion_dict[i] = 1
            for d in decoy[s.spectrum_id, charge]:
                pep = d.peptide
                bNy = interleave_b_y_ions(Peptide(pep), charge, mods, 
                                          ntermMods, ctermMods)
                numBY_dict_per_sid[sid, pep] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                theo_spec_dict[s.spectrum_id, pep] = bNy
                for i in bNy:
                    ion_dict[i] = 1

    ions = list(ion_dict.iterkeys())
    ions.sort()
    for i, ion in enumerate(ions):
        ion_dict[ion] = i

    # make collection per spectrum
    make_master_parameters(args, ion_dict, ions)
    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,0)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    # write peptide database to parse and identify GMTK segments later
    pepdb_list = open(os.path.join(args.output_dir, 'pepDB.txt'), "w")
    pepdb_list.write("Kind\tSid\tPeptide\tNumBY\tCharge\n")

    spec_dict = {}
    pep_num = 0
    for sid, charge in sid_charges:
        if sid not in spec_dict:
            s = spectra[sid]
            preprocess(s)
            spec_dict[sid] = s
        else:
            s = spec_dict[sid]

        for p in target[sid,charge]:
            pep = p.peptide
            bNy = theo_spec_dict[s.spectrum_id, pep]
            bNy = [ion_dict[bOrY] for bOrY in bNy]
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy))
            drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
            pepdb_list.write("t\t%d\t%s\t%d\t%d\n" % (sid, 
                                                      pep, 
                                                      numBY_dict_per_sid[sid, pep],
                                                      charge))
            pep_num += 1

        if (sid,charge) in decoy:
            for d in decoy[sid,charge]:
                pep = d.peptide
                bNy = theo_spec_dict[s.spectrum_id, pep]
                bNy = [ion_dict[bOrY] for bOrY in bNy]
                drip_peptide_sentence(pep_dt, pep, bNy, 
                                      pep_num, s.spectrum_id, args.max_obs_mass,
                                      peptide_pfile, False, len(bNy))
                drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                pepdb_list.write("d\t%d\t%s\t%d\t%d\n" % (sid, 
                                                          pep, 
                                                          numBY_dict_per_sid[sid, pep],
                                                          charge))
                pep_num += 1

    # close streams for this spectrum
    pep_dt.close()
    pepdb_list.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

    return spec_dict, pep_num

def make_drip_data_lowres(args, spectra, stdo, stde):
    """Generate test data .pfile. and create job scripts for cluster use.
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go
    """

    # parse modifications
    mods = df.parse_mods(args.mods_spec, True)
    print "mods:"
    print mods
    ntermMods = df.parse_mods(args.nterm_peptide_mods_spec, False)
    print "n-term mods:"
    print ntermMods
    ctermMods = df.parse_mods(args.cterm_peptide_mods_spec, False)
    print "c-term mods:"
    print ctermMods

    # load means
    dripMeans = load_drip_means(args.learned_means)
    # make master file
    make_master_parameters_lowres(args, dripMeans)

    if not args.append_to_pin:
        target,decoy,num_psms = load_psms(args.psm_file)
    else:
        target,decoy,num_psms = load_pin_file(args.psm_file)
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    # assume that we should randomize PSMs for multithreading purposes; only reason
    # why we are currently assuming this is that there is already a parameter for dripSearch
    # which signifies whether we should shuffle the data
    shuffle(sid_charges)

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    validcharges = args.charges

    # write peptide database to parse and identify GMTK segments later
    pepdb_list = open(os.path.join(args.output_dir, 'pepDB.txt'), "w")
    pepdb_list.write("Kind\tSid\tPeptide\tNumBY\tCharge\n")

    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,0)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    spec_dict = {}
    pep_num = 0
    for sid, charge in sid_charges:
        if sid not in spec_dict:
            s = spectra[sid]
            preprocess(s)
            spec_dict[sid] = s
        else:
            s = spec_dict[sid]

        if args.filt_theo_peaks:
            if args.per_spectrum_mz_bound:
                minMz = s.mz[0]
                maxMz = s.mz[-1]
            else:
                minMz = args.mz_lb
                maxMz = args.mz_ub

        for p in target[sid,charge]:
            pep = p.peptide
            bNy = interleave_b_y_ions_lowres(Peptide(pep), charge, mods,
                                             ntermMods, ctermMods)
            pepdb_list.write("t\t%d\t%s\t%d\t%d\n" % (sid, pep, len(bNy), charge))
            # numBY for DRIP features assumes all b-/y-ions, not just those
            # unfiltered per spectrum
            if args.filt_theo_peaks:
                filter_theoretical_peaks_lowres(bNy, dripMeans,
                                                minMz, maxMz)
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy))
            drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
            pep_num += 1

        if (sid,charge) in decoy:
            for d in decoy[sid,charge]:
                pep = d.peptide
                bNy = interleave_b_y_ions_lowres(Peptide(pep), charge, mods, 
                                          ntermMods, ctermMods)
                pepdb_list.write("d\t%d\t%s\t%d\t%d\n" % (sid, pep, len(bNy), charge))
                # numBY for DRIP features assumes all b-/y-ions, not just those
                # unfiltered per spectrum
                if args.filt_theo_peaks:
                    filter_theoretical_peaks_lowres(bNy, dripMeans,
                                             minMz, maxMz)
                drip_peptide_sentence(pep_dt, pep, bNy, 
                                      pep_num, s.spectrum_id, args.max_obs_mass,
                                      peptide_pfile, False, len(bNy))
                drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                pep_num += 1

    # close streams for this spectrum
    pep_dt.close()
    pepdb_list.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

    return spec_dict, pep_num

def runDripExtract(args, stdo, stde):
    """ Run drip once per spectrum, collapsing all charge-varying candidates into a single GMTK call
    """
    # create constant gmtkViterbi command line string
    # don't need frame/segment difference actions since each PSM corresponds to a specific spectrum, 
    # so that there isn't much redudandancy to exploit
    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -fdiffact2 rl" \
        + " -inputMasterFile model.mtr -inputTrainableParameters trained.params -failOnZeroClique F"

    # for now, don't worry about checking whether peptide is in valid (i.e., present in the digested
    # set of peptide candidates given the protein database)

    # currently ignore ident file input for spectra filtering
    spectra, minMz, maxMz, validcharges, _ = load_spectra_ret_dict(args.spectra, args.charges)
    # update encountered charges
    args.charges = validcharges
    args.mz_lb = minMz
    args.mz_ub = maxMz

    # create GMTK observation files
    # add in support for cluster usage later; assume standalone with multithreading
    if args.high_res_ms2:
        spec_dict, num_psms = make_drip_data_highres(args, spectra, stdo, stde)
    else:
        spec_dict, num_psms = make_drip_data_lowres(args, spectra, stdo, stde)

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
                           "inverted_collection/covar.txt",
                           "DRIP_GAUSSIAN_COMPONENTS",
                           "DRIP_GAUSSIAN_MIXTURES",
                           "DRIP_MZ_GAUSSIANS")
    except:
        print "Could not create DRIP master file %s, exitting" % args.master_file
        exit(-1)

    try:
        triangulate_drip(args.structure_file, args.master_file)
    except:
        print "Could not create triangulate structure file %s, exitting" % args.structure_file
        exit(-1)

    print args.covar_file
    try:
        write_covar_file(args.high_res_ms2, args.covar_file,
                         args.learned_covars)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    # run GMTK
    dtFile = os.path.join(args.output_dir, 'iterable.dts')
    cppCommand = '\'-DITERABLE_DT=' + dtFile \
        + ' -DDRIP_MZ=' + args.mean_file \
        + ' -DDRIP_GAUSSIAN_COMPONENTS=' + args.gauss_file \
        + ' -DDRIP_GAUSSIAN_MIXTURES=' + args.mixture_file \
        + ' -DDRIP_MZ_GAUSSIANS=' + args.collection_file \
        + '\''

    # call gmtkViterbi
    # gmtkViterbi command line
    vitValsFile = os.path.join(args.logDir, 'vitVals.txt')
    vitStr = vitStr0 + ' -vitValsFile ' +  vitValsFile \
        + ' -of1 ' + pfile_dir + '/spectrum.pfile' \
        + ' -of2 ' + pfile_dir + '/pep-lengths.pfile' \
        + ' -cppCommand ' + cppCommand
    call(shlex.split(vitStr), stdout = stdo, stderr = stde)

    t,d = psm.parse_dripExtract(vitValsFile, os.path.join(args.output_dir, 'pepDB.txt'))

    return t,d, spec_dict

def write_output(targets, decoys, filename, meanFile, spec_dict):
    """ Write PSMs and features to file
    """

    dripMeans = load_drip_means(meanFile)

    try:
        identFid = open(filename, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    identFid.write('Kind\tSid\tFrames\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\n')

    for sid, charge in targets:
        s = spec_dict[sid]
        for psm in targets[sid,charge]:
            write_psm_ins_dels(psm, s, dripMeans, identFid)
        if (sid,charge) in decoys:
            for psm in decoys[sid,charge]:
                write_psm_ins_dels(psm, s, dripMeans, identFid)
    identFid.close()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(conflict_handler='resolve')
    ############## input and output options
    iFileGroup = parser.add_argument_group('iFileGroup', 'Necessary input files.')
    help_psm_file = """<string> - File containing peptide-spectrum matches (PSMs) in either format: tab-delimited, PIN, pepXML, mzIdentML."""
    iFileGroup.add_argument('--psm-file', type = str, action = 'store',
                            help = help_psm_file)
    help_spectra = """<string> - The name of the file from which to parse fragmentation spectra, in ms2 format."""
    iFileGroup.add_argument('--spectra', type = str, action = 'store',
                            help = help_spectra)
    help_pepdb = """<string> - Protein FASTA file."""
    iFileGroup.add_argument('--fasta', type = str, action = 'store',
                            help = help_pepdb)
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
    ############## Cluster usage parameters
    clusterUsageGroup = parser.add_argument_group('clusterUsageGroup', 'Cluster data generation options.')
    help_randomize_ms2_spectra = '<T|F> - whether to randomize order of spectra when splitting ms2 file. Default = True'
    clusterUsageGroup.add_argument('--randomize_ms2_spectra', type = str, action = 'store', 
                                   default = 'True', help = help_randomize_ms2_spectra)
    help_random_wait = '<integer> - randomly wait up to specified number of seconds before writing results back to NFS. Default=30'
    clusterUsageGroup.add_argument('--random_wait', type = int, action = 'store', 
                                   default = 30, help = help_random_wait)
    help_num_jobs = '<integer> - the number of jobs to run in parallel. Default=1.'
    clusterUsageGroup.add_argument('--num-jobs', type = int, action = 'store', 
                                   default = 1, help = help_num_jobs)
    ############## peptide properties
    peptidePropertiesGroup = parser.add_argument_group('peptidePropertiesGroup', 'Options for digested peptids.')
    max_length_help = """<integer> - The maximum length of peptides to consider. Default = 50."""
    peptidePropertiesGroup.add_argument('--max-length', type = int, action = 'store',
                                        default = 50, help = max_length_help)
    max_mass_help = """<foat> - The maximum mass (in Da) of peptides to consider. Default = 7200."""
    peptidePropertiesGroup.add_argument('--max-mass', type = float, action = 'store',
                                        default = 7200.0, help = max_mass_help)
    min_length_help = """<integer> - The minimum length of peptides to consider. Default = 6."""
    peptidePropertiesGroup.add_argument('--min-length', type = int, action = 'store',
                                        default = 6, help = min_length_help)
    min_mass_help = """<foat> - The minimum mass (in Da) of peptides to consider. Default = 7200."""
    peptidePropertiesGroup.add_argument('--min-mass', type = float, action = 'store',
                                        default = 200.0, help = min_mass_help)
    monoisotopic_precursor_help = """<T|F> - When computing the mass of a peptide, use monoisotopic masses rather than average masses. Default = true."""
    peptidePropertiesGroup.add_argument('--monoisotopic-precursor', type = str, action = 'store', 
                                        default = 'True', help = monoisotopic_precursor_help)
    ############## amino acid modifications
    aaModsGroup = parser.add_argument_group('aaModsGroup', 'Options for amino acid modifications.')
    aaModsGroup.add_argument('--mods-spec', type = str, action = 'store',
                             default = 'C+57.02146')
    aaModsGroup.add_argument('--cterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    aaModsGroup.add_argument('--nterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    aaModsGroup.add_argument('--max-mods', type = int, action = 'store',
                             default = 255)
    aaModsGroup.add_argument('--min-mods', type = int, action = 'store',
                             default = 0)
    ############## decoy database generation
    decoyDBGroup = parser.add_argument_group('decoyDBGroup', 'Options decoy database generation.')
    help_decoy_format = '<shuffle|peptide-reverse> - Include a decoy version of every peptide by shuffling or reversing the target sequence or protein. In shuffle or peptide-reverse mode, each peptide is either reversed or shuffled, leaving the N-terminal and C-terminal amino acids in place. Note that peptides appear multiple times in the target database are only shuffled once. In peptide-reverse mode, palindromic peptides are shuffled. Also, if a shuffled peptide produces an overlap with the target or decoy database, then the peptide is re-shuffled up to 100n times, where n is the length of the string to be shuffled. Note that, despite this repeated shuffling, homopolymers will appear in both the target and decoy database. The protein-reverse mode reverses the entire protein sequence, irrespective of the composite peptides. Default = shuffle.'
    decoyDBGroup.add_argument('--decoy-format', type = str, action = 'store', 
                              default = 'shuffle', help = help_decoy_format)
    help_keep_terminal_aminos = '<N|C|NC|none> - When creating decoy peptides using decoy-format=shuffle or decoy-format=peptide-reverse, this option specifies whether the N-terminal and C-terminal amino acids are kept in place or allowed to be shuffled or reversed. For a target peptide "EAMPK" with decoy-format=peptide-reverse, setting keep-terminal-aminos to "NC" will yield "EPMAK"; setting it to "C" will yield "PMAEK"; setting it to "N" will yield "EKPMA"; and setting it to "none" will yield "KPMAE". Default = NC.'
    decoyDBGroup.add_argument('--keep-terminal-aminos', type = str, action = 'store', 
                              default = 'NC', help = help_keep_terminal_aminos)
    help_seed = '<seed> - When given a unsigned integer value seeds the random number generator with that value. When given the string "time" seeds the random number generator with the system time. Default = 1.'
    decoyDBGroup.add_argument('--seed', type = str, action = 'store', 
                              default = '1', help = help_seed)
    ############ enzymatic digestion
    enzymeDigestGroup = parser.add_argument_group('enzymeDigestGroup', 'Options for the digesting enzyme.')
    enzyme_help = """"<no-enzyme|trypsin|trypsin/p|chymotrypsin|elastase|clostripain|cyanogen-bromide|iodosobenzoate|proline-endopeptidase|staph-protease|asp-n|lys-c|lys-n|arg-c|glu-c|pepsin-a|elastase-trypsin-chymotrypsin|custom-enzyme> - Specify the enzyme used to digest the proteins in silico. Available enzymes (with the corresponding digestion rules indicated in parentheses) include no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p ([RK]|[]), chymotrypsin ([FWYL]|{P}), elastase ([ALIV]|{P}), clostripain ([R]|[]), cyanogen-bromide ([M]|[]), iodosobenzoate ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease ([E]|[]), asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin ([ALIVKRWFY]|{P}). Specifying --enzyme no-enzyme yields a non-enzymatic digest. Warning: the resulting peptide database may be quite large. Default = trypsin."""
    enzymeDigestGroup.add_argument('--enzyme', type = str, action = 'store',
                                   default = 'trypsin', help = enzyme_help)
    custom_enzyme_help = """<string> - Specify rules for in silico digestion of protein sequences. Overrides the enzyme option. Two lists of residues are given enclosed in square brackets or curly braces and separated by a |. The first list contains residues required/prohibited before the cleavage site and the second list is residues after the cleavage site. If the residues are required for digestion, they are in square brackets, '[' and ']'. If the residues prevent digestion, then they are enclosed in curly braces, '{' and '}'. Use X to indicate all residues. For example, trypsin cuts after R or K but not before P which is represented as [RK]|{P}. AspN cuts after any residue but only before D which is represented as [X]|[D]. Default = <empty>."""
    enzymeDigestGroup.add_argument('--custom-enzyme', type = str, action = 'store',
                                   default = '', help = custom_enzyme_help)
    missed_cleavages_help = """<integer> - Maximum number of missed cleavages per peptide to allow in enzymatic digestion. Default = 0."""
    enzymeDigestGroup.add_argument('--missed-cleavages', type = int, action = 'store',
                                   default = 0, help = custom_enzyme_help)
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
    parser.add_argument('--shards', type = int, action = 'store')
    parser.add_argument('--num_spectra', type = int, action = 'store')
    parser.add_argument('--ppm', action = 'store_true',
                        default = False)
    parser.add_argument('--filter_ident', action = 'store_true',
                        default = False,
                      help = "Filter sids by ident")
    parser.add_argument('-z', '--use_gzcat', action = 'store_true',
                        dest = 'gzcat', default = False,
                        help = "Use gzcat to decompress .ms2.gz files, if avail.")
    ######### encoding options
    parser.add_argument("--max_obs_mass", type = int, default = 0)
    parser.add_argument('--normalize', dest = "normalize", type = str,
                        help = "Name of the spectrum preprocessing pipeline.", 
                        default = 'top300TightSequest')
    help_filt_theo_peaks = "Filter theoretical peaks outside of observed spectra min and max"
    parser.add_argument('--filt_theo_peaks', action = 'store_true', dest = 'filt_theo_peaks', 
                        default = False, help = help_filt_theo_peaks)
    parser.add_argument('--per_spectrum_mz_bound', action = 'store_true', 
                        default = False, help = "Calculate observed m/z bound per spectrum")
    parser.add_argument('--mz_lb', type = float, action = 'store',
                        dest = 'mz_lb', help = 'Lower m/z bound' ,default = 150.0)
    parser.add_argument('--mz_ub', type = float, action = 'store',
                        dest = 'mz_ub', help = 'Upper m/z bound' ,default = 2000.0)
    parser.add_argument('--structure-file', type = str, action = 'store',
                        default = 'model.str',
                        help = "DRIP GMTK structure file")
    parser.add_argument('--master-file', type = str, action = 'store',
                        default = 'model.mtr',
                        help = "DRIP GMTK master file")
    parser.add_argument('--collection-dir', action = 'store', dest = 'collection_dir',
                        type = str, help = 'Where to store Gaussian files')
    parser.add_argument('--learned-means', action = 'store', dest = 'learned_means',
                        type = str, help = 'Learned means', default = "riptideLearnedMeans.txt")
    parser.add_argument('--learned-covars', action = 'store', dest = 'learned_covars',
                        type = str, help = 'Learned covariances', default = "riptideLearnedCovars.txt")
    parser.add_argument('--mean-file', action = 'store', dest = 'mean_file',
                        type = str, help = 'Where to put mean file.')
    parser.add_argument('--covar-file', action = 'store',
                        type = str, help = 'Covariance file name', 
                        default = "covar.txt")
    parser.add_argument('--gauss_file', action = 'store', dest = 'gauss_file',
                        type = str, help = 'Gaussians.')

    parser.add_argument('--mixture_file', action = 'store', dest = 'mixture_file',
                        type = str, help = 'Gaussian mixtures.')

    parser.add_argument('--collection_file', action = 'store', dest = 'collection_file',
                        type = str, help = 'Where to put collection file.')

    parser.add_argument('--covar_name', action = 'store', dest = 'covar_name',
                        type = str, help = 'Variance name.')

    parser.add_argument('--mixture_name', action = 'store', dest = 'mixture_name',
                        type = str, help = 'Gaussian mixture name.')

    parser.add_argument('--gaussian_component_name', action = 'store', dest = 'gaussian_component_name',
                        type = str, help = 'Gaussian mixture name.')

    parser.add_argument('--collection_name', action = 'store', dest = 'collection_name',
                        type = str, help = 'Gaussian mixture name.')

    parser.add_argument('--dpmf_name', action = 'store', dest = 'dpmf_name',
                        type = str, help = 'DPMF name.')

    # add output file options as a separate group
    oFileGroup = parser.add_argument_group('oFileGroup', 'Output file options')
    help_write_pin = '<T|F> - Write output in percolator PIN format. Default = False'
    oFileGroup.add_argument('--write-pin', type = str, action = 'store', 
                            default = 'False', help = help_write_pin)
    help_append_to_pin = '<T|F> - Append DRIP features to the features of the input PIN file. Default = False'
    oFileGroup.add_argument('--append-to-pin', type = str, action = 'store', 
                            default = 'False', help = help_append_to_pin)
    help_output_file = '<string> - Output file for resulting peptides.'
    oFileGroup.add_argument('--output_file', type = str, action = 'store', 
                            default = 'peptides.txt', help = help_output_file)
    help_peptide_database_file = """<string> - Output file for resulting peptides. Assumed to be binary format (i.e., python pickle file) for quick future access. See target-output-file/decoy-output-file options for ASCII file equivalents. Default = None."""
    oFileGroup.add_argument('--peptide-database-file', type = str, action = 'store', 
                            default = '', help = help_peptide_database_file)
    help_load_peptide_database_file = """<T|F> - Load previously digested peptide database file. If the supplied --peptide-database-file does not exist, the protein database will be digested. If decoys is true and no decoys are contained in the pickle, decoys will be created on the fly. Default = F."""
    oFileGroup.add_argument('--load-peptide-database-file', type = str, action = 'store', 
                            default = 'False', help = help_load_peptide_database_file)
    oFileGroup.add_argument('--logDir', type = str, 
                      help = 'directory with output log files')
    oFileGroup.add_argument('--output', type = str,
                      help = 'ident file name')
    args = parser.parse_args()

    # process input arguments
    process_args(args)

    stde = open('gmtk_err', "w")

    # set absolute path of learned means and covariance
    args.learned_means = os.path.abspath(args.learned_means)
    args.learned_covars = os.path.abspath(args.learned_covars)

    if debug:
        stdo = sys.stdout
    else:
        stdo = open(os.devnull, "w")
        sys.stdout = open("dripExtract_output", "w")

    if args.num_threads <= 1:
        # more efficient to call runDripExtract and redefine it below when performing multithreading,
        # rather than defining it below only (per the restriction of using the multiprocessing module
        # within main) and setting a pool of size 1 for single-threading.  The major drawback is needing
        # make sure both implementations are up to speed
        targets, decoys, spec_dict = runDripExtract(args, stdo, stde)
    else:
        # otherwise, redefine runDrip in main() so that multiprocessing runs properly
        # create constant gmtkViterbi command line string

    ####################################################################################
    ##################################### start runDripExtract for multithreading
    ####################################################################################

        # create constant gmtkViterbi command line string
        # don't need frame/segment difference actions since each PSM corresponds to a specific spectrum, 
        # so that there isn't much redudandancy to exploit
        vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
            + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
            + " -fdiffact2 rl" \
            + " -inputMasterFile model.mtr -inputTrainableParameters trained.params -failOnZeroClique F"

        # for now, don't worry about checking whether peptide is in valid (i.e., present in the digested
        # set of peptide candidates given the protein database)

        # currently ignore ident file input for spectra filtering
        spectra, minMz, maxMz, validcharges, _ = load_spectra_ret_dict(args.spectra, args.charges)
        # update encountered charges
        args.charges = validcharges
        args.mz_lb = minMz
        args.mz_ub = maxMz

        # create GMTK observation files
        if args.high_res_ms2:
            spec_dict, num_psms = make_drip_data_highres(args, spectra, stdo, stde)
        else:
            spec_dict, num_psms = make_drip_data_lowres(args, spectra, stdo, stde)

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
                               "inverted_collection/covar.txt",
                               "DRIP_GAUSSIAN_COMPONENTS",
                               "DRIP_GAUSSIAN_MIXTURES",
                               "DRIP_MZ_GAUSSIANS")
        except:
            print "Could not create DRIP structure file %s, exitting" % args.structure_file
            exit(-1)

        try:
            triangulate_drip(args.structure_file, args.master_file)
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
        dtFile = os.path.join(args.output_dir, 'iterable.dts')
        cppCommand = '\'-DITERABLE_DT=' + dtFile \
            + ' -DDRIP_MZ=' + args.mean_file \
            + ' -DDRIP_GAUSSIAN_COMPONENTS=' + args.gauss_file \
            + ' -DDRIP_GAUSSIAN_MIXTURES=' + args.mixture_file \
            + ' -DDRIP_MZ_GAUSSIANS=' + args.collection_file \
            + '\''

        # set up pool for multithreading
        pool = multiprocessing.Pool(processes = args.num_threads)

        inc = int(math.ceil(float(num_psms) / float(args.num_threads)))
    
        dcdrng_start = 0
        dcdrng_end = inc-1
        v = []
        for thread in range(args.num_threads):
            outputFile = 'vitVals' + str(thread) + '.txt'
            # gmtkViterbi command line
            vitValsFile = os.path.join(args.logDir, outputFile)
            v.append(vitValsFile)
            vitStr = vitStr0 + ' -vitValsFile ' +  vitValsFile \
                + ' -of1 ' + pfile_dir + '/spectrum.pfile' \
                + ' -of2 ' + pfile_dir + '/pep-lengths.pfile' \
                + ' -cppCommand ' + cppCommand \
                + ' -dcdrng ' + str(dcdrng_start) + ':' + str(dcdrng_end)

            dcdrng_start += inc
            dcdrng_end += inc
            if dcdrng_end >= num_psms:
                dcdrng_end = num_psms - 1
            # add to pool
            r = pool.apply_async(check_output, args=(shlex.split(vitStr),))
        # close pool
        pool.close()
        # wait for processes to finish
        pool.join()
        # initialize target and decoy dictionaries, passed by reference to the parser
        targets = {}
        decoys = {}
        for vitValsFile in v:
            psm.parse_dripExtract_parallel(vitValsFile, 
                                           os.path.join(args.output_dir, 'pepDB.txt'), 
                                           targets, decoys)
    ####################################################################################
    ##################################### end runDripExtract for multithreading
    ####################################################################################

    if args.write_pin:
        if not args.append_to_pin: # load database features to write drip output
            # database information required
            # formulate digestion regular expression
            digest_re = df.create_digest_re(args.enzyme, args.custom_enzyme)
            print "digest regular expression: %s" % digest_re
            r = re.compile(digest_re)

            # parse modifications
            mods = df.parse_mods(args.mods_spec, True)
            print "mods:"
            print mods
            nterm_mods = df.parse_mods(args.nterm_peptide_mods_spec, False)
            print "n-term mods:"
            print nterm_mods
            cterm_mods = df.parse_mods(args.cterm_peptide_mods_spec, False)
            print "c-term mods:"
            print cterm_mods

            cleavage_req = df.enzyme_enzSet(args.enzyme, args.custom_enzyme)
            target_db, decoy_db = load_target_decoy_db_no_reshuffle(args, r, mods, nterm_mods, cterm_mods)
            psm.write_percolator_pin(targets, decoys, args.output, args.mean_file, spec_dict, 
                                     cleavage_req, target_db, decoy_db)
        else: # don't worry about database features, just append DRIP features to PIN file's features
            targets0,decoys0 = load_pin_return_dict(args.psm_file)
            # merge the loaded PSM features with DRIP
            psm.append_to_percolator_pin(targets, decoys, 
                                         targets0, decoys0,
                                         args.output, 
                                         args.mean_file, spec_dict)
    else:
        psm.write_output(targets, decoys, args.output, args.mean_file, spec_dict)

    if stdo:
        stdo.close()
    if stde:
        stde.close()
