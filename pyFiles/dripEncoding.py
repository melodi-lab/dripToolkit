# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

import re
import os
import math
import sys
import shlex

from shutil import rmtree
from pfile.wrapper import PFile
from dripBootParameters import init_strictOrbitrapRiptide_learned_means, writeCPTs
from pyFiles.constants import allPeps, mass_h, mass_h2o
from subprocess import call, check_output

def gaussPdf(x, mu, sig):
    Z = (2*math.pi*sig)**.5
    return math.exp(-(x-mu)**2/(2*sig)) / Z

############## 
class dripExtractParams(object):
    """ used by dripToolKit to plot PSMs for a collection
    """
    def __init__(self, psm_file = '',
                 spectra = '', charges = 'all',
                 mods = '', ntermMods = '', ctermMods = '',
                 highResMs2 = False, 
                 dripLearnedMeans = 'dripLearned.means',
                 dripLearnedCovars = 'dripLearned.covars'):

        collection_dir = 'drip_collection'
        collection_name = 'drip_mz_Gaussians'
        mean_file = 'drip_mz.txt'
        covar_file = 'covar.txt'
        gauss_file = 'drip_Gaussian_components.txt'
        mixture_file = 'drip_Gaussian_mixtures.txt' 
        collection_file = 'drip_mz_Gaussians.txt'
        gaussian_component_name = 'gc'
        mixture_name = 'mixture'
        dpmf_name = 'unityDPMF' 
        covar_name = 'covar0'
        structure_file = 'drip.str'
        master_file = 'drip.mtr'
        max_obs_mass = 2001

        self.high_res_ms2 = highResMs2
        self.psm_file = psm_file
        self.spectra = spectra
        self.charges = charges
        self.mods_spec = mods
        self.nterm_peptide_mods_spec = ntermMods
        self.cterm_peptide_mods_spec = ctermMods

        # self.output_dir = 'dtk_psmCollection'
        self.output_dir = 'encode'
        self.obs_dir = 'obs'
        self.logDir = 'log'
        self.normalize = 'top300TightSequest'
        self.mz_lb = 0.0
        self.mz_ub = 2000.0
        self.append_to_pin = False
        self.per_spectrum_mz_bound = False
        self.filt_theo_peaks = True

        self.covar_file = os.path.join(collection_dir, covar_file)
        self.mean_file = os.path.join(collection_dir, mean_file)
        self.gauss_file = os.path.join(collection_dir, gauss_file)
        self.mixture_file = os.path.join(collection_dir, mixture_file)
        self.collection_file = os.path.join(collection_dir, collection_file)
        self.collection_name = collection_name
        self.gaussian_component_name = gaussian_component_name
        self.mixture_name = mixture_name
        self.dpmf_name = dpmf_name
        self.covar_name = covar_name
        self.mixture_name = mixture_name
        self.structure_file = structure_file
        self.master_file = master_file
        self.max_obs_mass = max_obs_mass
        self.learned_means = dripLearnedMeans
        self.learned_covars = dripLearnedCovars

        # create obervation directories
        # vit vals output directory
        base = os.path.abspath(collection_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)

        base = os.path.abspath(self.output_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)

        base = os.path.join(os.path.abspath(self.output_dir), self.obs_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)

        base = os.path.abspath(self.logDir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)


class dripGaussianCollectionNames(object):
    """ used by dripToolKit python interactive interpreter
    """
    def __init__(self, collection_dir = 'drip_collection', 
                 collection_name = 'drip_mz_Gaussians',
                 mean_file = 'drip_mz.txt', 
                 covar_file = 'covar.txt',
                 gauss_file = 'drip_Gaussian_components.txt',
                 mixture_file = 'drip_Gaussian_mixtures.txt', 
                 collection_file = 'drip_mz_Gaussians.txt',
                 gaussian_component_name = 'gc', 
                 mixture_name = 'mixture',
                 dpmf_name = 'unityDPMF', 
                 covar_name = 'covar0',
                 structure_file = 'drip.str',
                 master_file = 'drip.mtr',
                 max_obs_mass = 2001):
        
        # create obervation directories
        # vit vals output directory
        base = os.path.abspath(collection_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)

        self.covar_file = os.path.join(collection_dir, covar_file)
        self.mean_file = os.path.join(collection_dir, mean_file)
        self.gauss_file = os.path.join(collection_dir, gauss_file)
        self.mixture_file = os.path.join(collection_dir, mixture_file)
        self.collection_file = os.path.join(collection_dir, collection_file)
        self.collection_name = collection_name
        self.gaussian_component_name = gaussian_component_name
        self.mixture_name = mixture_name
        self.dpmf_name = dpmf_name
        self.covar_name = covar_name
        self.mixture_name = mixture_name
        self.structure_file = structure_file
        self.master_file = master_file
        self.max_obs_mass = max_obs_mass

############## dripMean class
class dripMean(object):
    """ GMTK mean class definition.  For DRIP, means have dimension 1
    """

    def __init__(self, name = '',
                 val = 0.0, index = 0, dimension = 1):

        assert dimension == 1, "Mean %s with dimension $d, must be 1 for DRIP model.  Exitting" % (name, dimension)

        self.name = name
        self.val = val
        self.index = index

    def __cmp__(self,other):
        ####### just check names
        if self.name == other.name:
            return 1
        return 0

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        k = "%s:%f" % (self.name,self.val)
        return k

def create_empty_master(outputFile = "drip.mtr"):
    fid = open(outputFile, "w")
    print >> fid, """
"""
    fid.close()

def create_params_no_train(outputFile = "params.notrain"):
    fid = open(outputFile, "w")
    print >> fid, """DENSECPT *
DPMF unityDPMF
COVAR covar0
MEAN intensity_mean
"""
    fid.close()

def create_drip_master(highRes = False, mtrFile = "model.mtr", 
                       maxTermMass = "200002",     
                       meanFile = "drip_collection/drip_mz.txt",
                       covarFile = "drip_collection/covar.txt",
                       gaussComponents = "drip_collection/drip_Gaussian_components.txt",
                       gaussMixtures = "drip_collection/drip_Gaussian_mixtures.txt",
                       gaussCollection = "drip_collection/drip_mz_Gaussians.txt"):
    ciPenalty = "-2.525473431771295"
    deltaStates = "202"
    peakLength = "200"
    mzObs = "0"
    intensityObs = "1"
    peakLengthObs = "2"
    hrPenalty = "-4.5369118985307910"
    lrPenalty = "-6.831884549495013"

    try:
        fid = open(mtrFile, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % mtrFile
        exit(-1)

    print >> fid, """% Non-trainable parts of the model parameterization.
"""

    print >> fid, "MEAN_IN_FILE %s ascii\n" % meanFile
    print >> fid, "COVAR_IN_FILE %s ascii\n" % covarFile

    print >> fid, """DPMF_IN_FILE inline 1
0
unityDPMF
1
1.0

MEAN_IN_FILE inline 1
0
intensity_mean 1 1.5

MC_IN_FILE inline 1
0
1
0
gc_intensity
intensity_mean covar_intensity

MX_IN_FILE inline 1
0
1
gm_intensity
1
unityDPMF
gc_intensity
"""

    print >> fid, "MC_IN_FILE %s ascii\n" % gaussComponents

    print >> fid, "MX_IN_FILE %s ascii\n" % gaussMixtures

    print >> fid, "NAME_COLLECTION_IN_FILE %s ascii\n" % gaussCollection

    print >> fid, """DT_IN_FILE inline 3
%STRUCTURE:
%dt_num
%dt_name
%num_parents
%query_parent splits split0 split1 ... spliti ... splitn
%-1 branch

0
increment_fragment_mass
ITERABLE_DT
% iterable.dts

1
switch_delta_distribution
% parents: p0=num_peaks_minus_one(0), p1=peak_position(-1), p2=insertion(-1)
3
2 2 0 default
  -1 { max(0, p0-p1) } %n-1 deletes left, switch to  appropriate distribution
  -1 { 0 } %do not increment

2
increment_position
2
   -1 { min(p0 + p1, 199) }


DETERMINISTIC_CPT_IN_FILE inline 5
%CPT_num
%CPT_name
%num_parents
%parent_cardinalities
%self_cardinality
%dt_mapping_name
"""
    print >> fid, "0\nend_reached_cpt\n0\n%s\ninternal:alwaysZero\n" % deltaStates

    print >> fid, "1\nincrement_position_cpt\n2\n%s %s\n%s\nincrement_position\n" % (deltaStates, 
                                                                                                   peakLength,
                                                                                                   peakLength)
    
    print >> fid, "2\nincrement_fragment_mass_cpt\n1\n%s\n%s\nincrement_fragment_mass\n" % (peakLength,
                                                                                            maxTermMass)

    print >> fid, """3
not_insertion_cpt
0
2
internal:alwaysZero
"""

    print >> fid, "4\nprologue_position_cpt\n1\n%s\n%s\ninternal:copyParent\n" % (deltaStates,
                                                                                  peakLength)

    fid.close()

def create_drip_constants():
    with open('constants.inc', 'w') as out:
        print >> out, """#define MAX_TERM_MASS 200002
#define PEAK_LENGTH 200
#define DELTA_STATES 202

%%%%%%%% observations
#define MZ_OBS 0
#define INTENSITY_OBS 1
#define PEAK_LENGTH_OBS 2

%%%%%%%% insertion penalties
#define HR_PEN -4.5369118985307910
#define LR_PEN -6.831884549495013
"""

def create_drip_structure(highRes = False, strFile = "model.str", 
                          maxTermMass = "200002", forcedAlignment = False, 
                          training = False,
                          highResMzWidth = 0.1):
    writeCPTs()
    create_drip_constants()
    if not training:
        ciPenalty = "-2.525473431771295"
    else:
        ciPenalty = "-4.4444460178209342"
    deltaStates = "202"
    peakLength = "200"
    if not forcedAlignment:
        mzObs = "0"
        intensityObs = "1"
        peakLengthObs = "2"
    else:
        mzObs = "0"
        intensityObs = "1"
        insertionObs = "2"
        peakLengthObs = "3"

    # default insertion penalties
    hrPenalty = -4.5369118985307910
    lrPenalty = -6.831884549495013

    # compute high-resolution insertion penalty
    if highRes:
        sig = (highResMzWidth / 8.0) ** 2 # 99.9 percent of Gaussian mass, i.e., 8 standard deviations, lies within highResMzWidth
        mu = highResMzWidth / 2.0
        hrPenalty = math.log(gaussPdf(0.0, mu, sig))

    try:
        fid = open(strFile, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % strFile
        exit(-1)

    # build General string for variables present in each frame
    # prologue specific variables
    prologueString="""   variable: DELTA { %Bernoulli random variable; if 1, we move down the theoretical spectrum, else stay at current theoretical peak.  If Insertion(-1)=1, Delta cannot be 1
"""
    
    prologueString +="     type: discrete hidden cardinality %s;\n" % deltaStates
    
    prologueString +="""     switchingparents: NUM_PEAKS_MINUS_ONE(0) using mapping ("internal:copyParent");
     conditionalparents: nil using DeterministicCPT("end_reached_cpt")
			| nil using DenseCPT("delta_cpt")
"""
    
    for i in range(1,199):
        prologueString += "			| nil using DenseCPT(\"delta_jump%d_cpt\")\n" % i

    prologueString += """			| nil using DenseCPT("delta_jump199_cpt");
   }

   variable: PEAK_POSITION { %Index of the theoretical peak currently being accessed
"""

    prologueString +="     type: discrete hidden cardinality %s;\n" %  peakLength


    prologueString +="""     conditionalparents: DELTA(0) using DeterministicCPT("prologue_position_cpt");
   }
"""

    # chunkEpilogue specific variables
    chunkEpilogueString="""   variable: DELTA { %Bernoulli random variable; if 1, we move down the theoretical spectrum, else stay at current theoretical peak.  If Insertion(-1)=1, Delta cannot be 1
"""
    
    chunkEpilogueString +="     type: discrete hidden cardinality %s;\n" % deltaStates
    
    chunkEpilogueString +="""     switchingparents: NUM_PEAKS_MINUS_ONE(0), PEAK_POSITION(-1), INSERTION(-1) using mapping ("switch_delta_distribution");
     conditionalparents: nil using DeterministicCPT("end_reached_cpt")
			| nil using DenseCPT("delta_cpt")
"""
   
    for i in range(1,200):
        chunkEpilogueString += "			| nil using DenseCPT(\"delta_jump%d_cpt\")\n" % i

    chunkEpilogueString += """			| nil using DenseCPT("delta_jump200_cpt");
   }

   variable: PEAK_POSITION { %Index of the theoretical peak currently being accessed
"""

    chunkEpilogueString +="     type: discrete hidden cardinality %s;\n" %  peakLength


    chunkEpilogueString +="""     conditionalparents: DELTA(0), PEAK_POSITION(-1)
     			 using DeterministicCPT("increment_position_cpt");
   }
"""

    epilogueString ="""   variable: SPECTRUM_LENGTH_OFFSET { %Arbitrary constant
      type: discrete observed value 1 cardinality 2;
     conditionalparents: nil using DenseCPT("score_offset");
    }
"""

    frameString="""    variable: NUM_PEAKS_MINUS_ONE {
"""

    frameString += "      type: discrete observed %s:%s cardinality %s;\n" % (peakLengthObs, peakLengthObs, peakLength)

    frameString +="""      conditionalparents: nil using DenseCPT("internal:UnityScore");
    }
"""

    frameString +="""    variable: FRAGMENT_MASS { %m/z location of a peak in the theoretical spectrum
"""

    frameString += "      type: discrete hidden cardinality %s;\n" % maxTermMass

    if not forcedAlignment:
        frameString +="""      conditionalparents: PEAK_POSITION(0) using DeterministicCPT("increment_fragment_mass_cpt");
    }

    variable: INSERTION { %Bernoulli random variable, determines whether we score observed peak with small var gaussian (Insertion=0) or constant penalty (Insertion=1)
      type: discrete hidden cardinality 2;
      conditionalparents: nil using DenseCPT("insertion_cpt");
   }

    variable: CMZ { %m/z spectrum observation; floating point value
"""
    else:
        frameString +="""      conditionalparents: PEAK_POSITION(0) using DeterministicCPT("increment_fragment_mass_cpt");
    }
"""

        frameString +="""    variable: INSERTION { %Bernoulli random variable, determines whether we score observed peak with small var gaussian (Insertion=0) or constant penalty (Insertion=1)
"""
        frameString+="      type: discrete observed %s:%s cardinality 2;\n" % (insertionObs,insertionObs)
        frameString +="""      conditionalparents: nil using DenseCPT("insertion_cpt");
   }

    variable: CMZ { %m/z spectrum observation; floating point value
"""

    frameString += "      type: continuous observed %s:%s;\n" % (mzObs, mzObs)

    frameString += """      weight:  scale 1 
"""

    if highRes:
        frameString += "	      | penalty %.16f;\n" % hrPenalty
    else:
        frameString += "	      | penalty %.16f;\n" % lrPenalty

    frameString += """      switchingparents: INSERTION(0) using mapping("internal:copyParent");
      conditionalparents: FRAGMENT_MASS(0) using mixture collection("drip_mz_Gaussians") mapping("internal:copyParent")
      			  | nil using mixture("internal:UnityScore");
     }

     variable: CI {
"""

    frameString += "       type: continuous observed %s:%s;\n" % (intensityObs, intensityObs)

    frameString += """       switchingparents: INSERTION(0) using mapping("internal:copyParent");
       weight: scale 1
"""

    frameString += """      	      | penalty %s;\n""" % ciPenalty

    frameString +="""       conditionalparents: nil using mixture("gm_intensity")
       		  	  | nil using mixture("internal:UnityScore");
     }
"""

    print >> fid, """GRAPHICAL_MODEL inverted_model
"""
    
    for frame in range(3):
        print >> fid, "frame: %d {" % frame
        if frame==0:
            print >> fid, prologueString
        else:
            print >> fid, chunkEpilogueString

        print >> fid, frameString

        if frame==2:
            print >> fid, epilogueString
        print >> fid, "}"

    print >> fid,  "chunk 1:1"
    fid.close()

def triangulate_drip(strFile = "model.str", mtrFile = "model.mtr"):
    triFile = strFile + '.trifile'

    triStr = 'gmtkTriangulate' \
        + ' -strFile ' + strFile \
        + ' -inputMasterFile ' + mtrFile
    try:
        os.remove(triFile)
    except OSError:
        pass
    check_output(shlex.split(triStr))
    # call(shlex.split(triStr), 
    #      stdout = stdo, stderr = stde)
    # call(['gmtkTriangulate', '-strFile', strFile, '-inputMasterFile', mtrFile] , 
    #      stdout = stdo, stderr = stde)

def write_covar_file(highRes = False, covarFile = 'covar.txt', 
                     learnedCovars = '', riptidePrior = True,
                     highResMzWidth = 0.1):
    hrMzCovar = 1.5625000000e-04
    lrMzCovar = 1.5625000000e-02

    if highRes:
        hrMzCovar = (highResMzWidth / 8.0) ** 2 # 99.9 percent of Gaussian mass, i.e., 8 standard deviations, lies within highResMzWidth

    if riptidePrior:
        intensityCovar = "8.7238661945e-02"
    else:
        intensityCovar = "0.12"

    if learnedCovars:
        try:
            covar_file = open(learnedCovars)
            log = covar_file.read()
            covar_file.close()
            covar_pattern = '[^%]?P<covarName>\S+) 1 (?P<covarVal>\S+)'
            found = False
            for match in re.finditer(covar_pattern, log):
                if match.group('covarName') == 'covar_intensity':
                    intensityCovar = match.group('covarVal')
                    found = True
            if not found:
                if riptidePrior:
                    intensityCovar = "8.7238661945e-02"
                    # intensityCovar = "0.0872357562"
                else:
                    intensityCovar = "0.12"
        except:
            if riptidePrior:
                intensityCovar = "8.7238661945e-02"
                # intensityCovar = "0.0872357562"
            else:
                intensityCovar = "0.12"
    try:
        fid = open(covarFile, 'w')
    except IOError:
        print "Could not open file %s for writing, exitting" % covarFile
        exit(-1)
        
    fid.write("2\n0\n")
    if highRes:
        fid.write("covar0 1 %e\n" % hrMzCovar)
    else:
        fid.write("covar0 1 %e\n" % lrMzCovar)

    fid.write("1\ncovar_intensity 1 %s" % intensityCovar)

    fid.close()

def create_pfile(directory, fn, num_flots, num_ints):
    return PFile(num_flots, num_ints, os.path.join(directory, fn))

def drip_peptide_sentence(pep_dt, peptide, bNy, pep_num, spectra_id, max_mass,
                          peptide_pfile, istarget, max_candidate_theo_peaks):
    """Add a sentence/segment to the peptide ascii file.

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Effects:
        Adds a sentence to 'pep_dt'.

    """
    if(istarget):
        pep_dt.write('%d\nsid%dtargetseq%s\n1\n0 %d ' % (pep_num, spectra_id, peptide, len(bNy)))
    else:
        pep_dt.write('%d\nsid%ddecoyseq%s\n1\n0 %d ' % (pep_num, spectra_id, peptide, len(bNy)))

    for i in range(0, len(bNy)-1): 
        pep_dt.write('%d ' % i)
    pep_dt.write('default\n')

    for ion in bNy:
        pep_dt.write('\t-1 { %d }\n' % ion)
    pep_dt.write('\n')

    #write number of peaks and number of spurious peaks to pfile
    arglist = []
    arglist.append(max_candidate_theo_peaks)

    peptide_pfile.add_frame(*arglist)
    peptide_pfile.end_segment()

def drip_spectrum_sentence(pf, mz, intensity):
    """Add a sentence/segment to the spectrum PFile.

    Arguments:
       pf: The PFile object, instance of pfile.wrapper.PFile.
       bins: Iterable of floats, representing the binned spectrum.

    Effects:
       Adds a sentence to 'pf'. Every call must be paired with a call
       to peptide_sentence.

    """
    
    for m, i in zip(mz, intensity):
        arglist = []
        arglist.append(m)
        arglist.append(i)
        pf.add_frame(*arglist)

    pf.end_segment()

def training_spectrum_sentence(pf, mz, intensity, dripPsm):
    """Add a sentence/segment to the spectrum PFile.

    Arguments:
       pf: The PFile object, instance of pfile.wrapper.PFile.
       bins: Iterable of floats, representing the binned spectrum.

    Effects:
       Adds a sentence to 'pf'. Every call must be paired with a call
       to peptide_sentence.

    """
    
    assert len(mz)==len(dripPsm.insertion_sequence), "Scan %d, PSM %s number of insertions does not equal number of frames, exitting" % (dripPsm.scan, dripPsm.peptide)

    for m, i, ins in zip(mz, intensity, dripPsm.insertion_sequence):
        arglist = []
        arglist.append(m)
        arglist.append(i)
        arglist.append(ins)
        pf.add_frame(*arglist)

    pf.end_segment()

def make_master_parameters(options, fragment_ion_dict, ions):
    """Mean file and collection file."""
    mean_file = open('%s' % options.mean_file, "w")
    gauss_file = open('%s' % options.gauss_file, "w")
    mixture_file = open('%s' % options.mixture_file, "w")
    collection_file = open('%s' % options.collection_file, "w")

    # collection file comment first
    collection_file.write('%number_of_collections collection_num name num_elements element_names\n')

    if len(fragment_ion_dict) > options.max_obs_mass:
        options.max_obs_mass = len(fragment_ion_dict)
    mean_file.write('%d\n' % (len(fragment_ion_dict)))
    gauss_file.write('%d\n\n' % (len(fragment_ion_dict)))
    mixture_file.write('%d\n\n' % (len(fragment_ion_dict)))
    collection_file.write('1\n\n0\n%s\n%d\n' % (options.collection_name, len(fragment_ion_dict)))
    # comments describing file parameters
    mean_file.write('%mean_number mean_name dimensionality mean_values\n\n')
    gauss_file.write('%component_number dimensionality type name mean_name covariance_name\n\n')
    mixture_file.write('%mixture_number dimensionality name num_components dpmf components\n\n')

    for ion in ions:
        i = fragment_ion_dict[ion]
        mean_file.write('%d mean%d 1 %.8f\n' % (i, i, ion))
        gauss_file.write('%d 1 0 %s%d mean%d %s\n' % (i, options.gaussian_component_name, i, i, options.covar_name))
        mixture_file.write('%d 1 %s%d 1 %s %s%d\n' % (i, options.mixture_name, i, options.dpmf_name, options.gaussian_component_name, i))
        collection_file.write('%s%d\n' % (options.mixture_name, i))
        
    mean_file.close()
    gauss_file.close()
    mixture_file.close()
    collection_file.close()

def make_master_parameters_lowres(options, dripMeans):
    """Mean file and collection file."""
    mean_file = open('%s' % options.mean_file, "w")
    gauss_file = open('%s' % options.gauss_file, "w")
    mixture_file = open('%s' % options.mixture_file, "w")
    collection_file = open('%s' % options.collection_file, "w")

    ions = sorted(list(dripMeans.iterkeys()))

    # collection file comment first
    collection_file.write('%number_of_collections collection_num name num_elements element_names\n')

    if len(dripMeans) > options.max_obs_mass:
        options.max_obs_mass = len(dripMeans)
    mean_file.write('%d\n' % (len(dripMeans)))
    gauss_file.write('%d\n\n' % (len(dripMeans)))
    mixture_file.write('%d\n\n' % (len(dripMeans)))
    collection_file.write('1\n\n0\n%s\n%d\n' % (options.collection_name, len(dripMeans)))
    # comments describing file parameters
    mean_file.write('%mean_number mean_name dimensionality mean_values\n\n')
    gauss_file.write('%component_number dimensionality type name mean_name covariance_name\n\n')
    mixture_file.write('%mixture_number dimensionality name num_components dpmf components\n\n')

    # for i, ind in enumerate(ions):
    #     ion = dripMeans[ind]
    #     mean_file.write('%d mean%d 1 %.8f\n' % (i, i, ion))
    #     gauss_file.write('%d 1 0 %s%d mean%d %s\n' % (i, options.gaussian_component_name, i, i, options.covar_name))
    #     mixture_file.write('%d 1 %s%d 1 %s %s%d\n' % (i, options.mixture_name, i, options.dpmf_name, options.gaussian_component_name, i))
    #     collection_file.write('%s%d\n' % (options.mixture_name, i))

    for i, ind in enumerate(ions):
        ion = dripMeans[ind]
        meanName = 'mean' + str(ind)
        gcName = options.gaussian_component_name + str(ind)
        mixtureName = options.mixture_name + str(ind)
        mean_file.write('%d %s 1 %.8f\n' % (i, meanName, ion))
        gauss_file.write('%d 1 0 %s %s %s\n' % (i, gcName, meanName, options.covar_name))
        mixture_file.write('%d 1 %s 1 %s %s\n' % (i, mixtureName, options.dpmf_name, gcName))
        collection_file.write('%s\n' % (mixtureName))
        
    mean_file.close()
    gauss_file.close()
    mixture_file.close()
    collection_file.close()

def return_b_y_ions(peptide, charge, mods = {},
                    ntermMods = {}, ctermMods = {}, 
                    ion_to_index_map = {}):
    """Create sets of b- and y-ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    ion_to_index_map = mapping from ion values (floats, drip mean values) to indices in the global collection of means.  A python dictionary with keys ions/drip mean values and values indices in the global collection of means

    Returns:
        Dictionaries of b- and y-ions where the keys are indices of the
        interleaved spectra and the elements are the floating point values

    """
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    bions = {}
    yion = {}
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = cf*mass_h
        offset = ntermOffset
        for b, aa in zip(ntm[1:-1], peptide.seq[:-1]):
            if aa in mods:
                offset += mods[aa]
            ion = (b+offset+bOffset)/cf
            bions[ion_to_index_map[ion]] = ion

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = mass_h2o + cf*mass_h
        offset = ctermOffset
        for y, aa in zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:])):
            if aa in mods:
                offset += mods[aa]
            ion = (y+offset+yOffset)/cf
            yions[ion_to_index_map[ion]] = ion

    return bions, yions

def return_b_y_ions_lowres(peptide, charge, mods = {},
                           ntermMods = {}, ctermMods = {},
                           ion_to_index_map = {}):
    """Create vector of interleaved b and y ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    bions = {}
    yions = {}
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = int(round(cf*mass_h))
        offset = ntermOffset
        numRes = 1
        for b, aa in zip(ntm[1:-1], peptide.seq[:-1]):
            if aa in mods:
                offset += mods[aa]
            ion = float(b+int(round(offset))+bOffset)/cf
            bions[int(round(ion))] = (numRes, c)
            numRes += 1

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = int(round(mass_h2o + cf*mass_h))
        offset = ctermOffset
        numRes = 1
        for y, aa in zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:])):
            if aa in mods:
                offset += mods[aa]
            ion = (y+int(round(offset))+yOffset)/cf
            yions[int(round(ion))] = (numRes, c)
            numRes += 1

    return bions, yions

def interleave_b_y_ions(peptide, charge, 
                        mods = {}, ntermMods = {}, ctermMods = {}):
    """Create vector of interleaved b and y ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    Returns:
        Vector of interleaved b and y ions

    """
    # (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    # ntm_fragments = []
    # ctm_fragments = []
    # ntm_charges = []
    # ctm_charges = []
    # bNy = []
    # for c in range(1,charge):
    #     ntm_fragments += [int(round(float(b+c)/float(c))) for b in ntm[1:-1]]
    #     ctm_fragments += [int(round(float(y+18+c)/float(c))) for y in ctm[1:-1]]
    #         # vector of charges
    #     bNy += [c for b in ntm[1:-1]]
    #     bNy += [c for y in ctm[1:-1]]

    # return sorted(list(set(bNy)))

    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    bNy = []
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = cf*mass_h
        offset = ntermOffset
        for b, aa in zip(ntm[1:-1], peptide.seq[:-1]):
            if aa in mods:
                offset += mods[aa]
            bNy.append((b+offset+bOffset)/cf)

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = mass_h2o + cf*mass_h
        offset = ctermOffset
        for y, aa in zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:])):
            if aa in mods:
                offset += mods[aa]
            bNy.append((y+offset+yOffset)/cf)

    return sorted(list(set(bNy)))

def interleave_b_y_ions_var_mods(peptide, charge, 
                                 mods = {}, ntermMods = {}, ctermMods = {},
                                 varMods = {}, ntermVarMods = {}, ctermVarMods = {},
                                 varModSequence = []):
    """Create vector of interleaved b and y ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    Returns:
        Vector of interleaved b and y ions

    """
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    bNy = []
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    elif p[0] in ntermVarMods:
        if varModSequence[0] == '2': # denotes an nterm variable modification
            ntermOffset = ntermVarMods[p[0]][1]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]
    elif p[-1] in ctermVarMods:
        if varModSequence[-1] == '3': # denotes a cterm variable modification
            ctermOffset = ctermVarMods[p[-1]][1]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = cf*mass_h
        offset = ntermOffset
        for ind, (b, aa) in enumerate(zip(ntm[1:-1], peptide.seq[:-1])):
            if aa in mods:
                offset += mods[aa]
            elif aa in varMods:
                if varModSequence[ind]=='1':
                    offset += varMods[aa][1]
            bNy.append((b+offset+bOffset)/cf)

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = mass_h2o + cf*mass_h
        offset = ctermOffset
        for ind, (y, aa) in enumerate(zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:]))):
            if aa in mods:
                offset += mods[aa]
            elif aa in varMods:
                if varModSequence[ind]=='1':
                    offset += varMods[aa][1]
            bNy.append((y+offset+yOffset)/cf)

    return sorted(list(set(bNy)))

def interleave_b_y_ions_lowres(peptide, charge, mods = {},
                               ntermMods = {}, ctermMods = {}):
    """Create vector of interleaved b and y ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    bNy = []
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = int(round(cf*mass_h))
        offset = ntermOffset
        for b, aa in zip(ntm[1:-1], peptide.seq[:-1]):
            if aa in mods:
                offset += mods[aa]
            bNy.append(int(round(float(b+int(round(offset))+bOffset)/cf)))

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = int(round(mass_h2o + cf*mass_h))
        offset = ctermOffset
        for y, aa in zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:])):
            if aa in mods:
                offset += mods[aa]
            bNy.append(int(round((y+int(round(offset))+yOffset)/cf)))

    return sorted(list(set(bNy)))

def interleave_b_y_ions_var_mods_lowres(peptide, charge, 
                                        mods = {}, ntermMods = {}, ctermMods = {},
                                        varMods = {}, ntermVarMods = {}, ctermVarMods = {},
                                        varModSequence = []):
    """Create vector of interleaved b and y ions

    Arguments:
    Peptide = current peptide, instance of object protein.Peptide
    charge = current charge, calculate quantities of mass / charge
    ntermStatMods = dictionary of static modifications specific to the n-terminus
    ctermStatMods = dictionary of static modifications specific to the c-terminus

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    bNy = []
    ntermOffset = 0
    ctermOffset = 0

    p = peptide.seq
    # check n-/c-term amino acids for modifications
    if p[0] in ntermMods:
        ntermOffset = ntermMods[p[0]]
    elif p[0] in ntermVarMods:
        if varModSequence[0] == '2': # denotes an nterm variable modification
            ntermOffset = ntermVarMods[p[0]][1]
    if p[-1] in ctermMods:
        ctermOffset = ctermMods[p[0]]
    elif p[-1] in ctermVarMods:
        if varModSequence[-1] == '3': # denotes a cterm variable modification
            ctermOffset = ctermVarMods[p[-1]][1]

    for c in range(1,max(charge,2)):
        cf = float(c)
        bOffset = int(round(cf*mass_h))
        offset = ntermOffset
        for ind, (b, aa) in enumerate(zip(ntm[1:-1], peptide.seq[:-1])):
            if aa in mods:
                offset += mods[aa]
            elif aa in varMods:
                if varModSequence[ind]=='1':
                    offset += varMods[aa][1]
            bNy.append(int(round(float(b+int(round(offset))+bOffset)/cf)))

    for c in range(1,max(charge,2)):
        cf = float(c)
        yOffset = int(round(mass_h2o + cf*mass_h))
        offset = ctermOffset
        for ind, (y, aa) in enumerate(zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:]))):
            if aa in mods:
                offset += mods[aa]
            elif aa in varMods:
                if varModSequence[ind]=='1':
                    offset += varMods[aa][1]
            bNy.append(int(round((y+int(round(offset))+yOffset)/cf)))

    return sorted(list(set(bNy)))

def filter_theoretical_peaks_lowres(theo_spectrum, 
                                    drip_means, minMz, maxMz, 
                                    binWidth = 1.0005079):
    """ Filter theoretical peaks outside of the observed minimum and maximum m/z range.
        theo_spectrum = vector of theoretical peaks(integer valued)
        drip_means = dictionary of DRIP learned means
        minMz = minimum observed m/z value
        maxMz = maximum observed m/z value
        pre-conditions:
        -theo_spectrum has been constructed(note that it is passed by reference)
        -dictionary of learned DRIP means, drip_means, has been loaded
        -minMz and maxMz have been calculated
        post-conditions:
        -theo_spectrum has set of theoretical peaks X removed such that, 
        for theoretical peak x \in X, x < minMz-binWidth/2 or x > maxMz+binWidth/2
    """
    # first remove theoretical peaks with learned mean value, x, such that
    # x < minMz - binWidth/2
    lb = minMz - binWidth/2.0
    filt_lb_ind = -1
    for i, m in enumerate(theo_spectrum):
        if m >= lb:
            filt_lb_ind = i
            break
    # make sure there are actually elements to delete
    if filt_lb_ind > -1:
        if filt_lb_ind > 0:
            del theo_spectrum[:filt_lb_ind]

    # next remove theoretical peaks with learned mean value, x, such that
    # x > maxMz + binWidth/2
    ub = maxMz + binWidth/2.0
    filt_ub_ind = -1
    for i, m in enumerate(theo_spectrum[::-1]):
        if m <= ub:
            filt_ub_ind = i
            break
    # make sure we are actually 
    if filt_ub_ind > -1:
        del theo_spectrum[(len(theo_spectrum)-filt_ub_ind):len(theo_spectrum)]
    ## print calculated bound indices
    # print "lb=%d,ub=%d" % (filt_lb_ind, filt_ub_ind)

def filter_theoretical_peaks(theo_spectrum, minMz, maxMz, binWidth = 0.1):
    """ Filter theoretical peaks outside of the observed minimum and maximum m/z range.
        theo_spectrum = vector of theoretical peaks(integer valued)
        drip_means = dictionary of DRIP learned means
        minMz = minimum observed m/z value
        maxMz = maximum observed m/z value
        pre-conditions:
        -theo_spectrum has been constructed(note that it is passed by reference)
        -dictionary of learned DRIP means, drip_means, has been loaded
        -minMz and maxMz have been calculated
        post-conditions:
        -theo_spectrum has set of theoretical peaks X removed such that, 
        for theoretical peak x \in X, x < minMz-binWidth/2 or x > maxMz+binWidth/2
    """
    # first remove theoretical peaks, x, such that
    # x < minMz - binWidth/2
    lb = minMz - binWidth/2.0
    filt_lb_ind = -1
    for i, m in enumerate(theo_spectrum):
        if m > lb:
            filt_lb_ind = i
            break
    # make sure there are actually elements to delete
    if filt_lb_ind > -1:
        if filt_lb_ind > 0:
            del theo_spectrum[:filt_lb_ind]
        # else:
        #     del theo_spectrum[filt_lb_ind]

    # next remove theoretical peaks, x, such that
    # x > maxMz + binWidth/2
    ub = maxMz + binWidth/2.0
    filt_ub_ind = -1
    for i, m in enumerate(theo_spectrum[::-1]):
        if m < ub:
            filt_ub_ind = i
            break
    # make sure we are actually 
    if filt_ub_ind > -1:
        del theo_spectrum[(len(theo_spectrum)-filt_ub_ind):len(theo_spectrum)]

def load_drip_means(master_file):
    """ Load DRIP learned means from master file.
        Return dictionary whose keys are the mean indices
    """
    try:
        mean_file = open(master_file)
        log = mean_file.read()
        mean_file.close()
        mean_pattern = '[^%]mean(?P<meanInd>\d+) 1 (?P<meanVal>\S+)'

        means = {}
        for match in re.finditer(mean_pattern, log):
            means[int(match.group('meanInd'))] = float(match.group('meanVal'))
    except IOError:
        print "Could not load means from %s, loading default means" % master_file
        means = init_strictOrbitrapRiptide_learned_means()

    return means
