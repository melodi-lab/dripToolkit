#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

__authors__ = [ 'John T. Halloran <halloj3@uw.edu>' ]

"""DRIP Peptide-Spectrum Matches (PSMs)
"""

import string
import csv
import re
import struct

# import matplotlib
# matplotlib.use('Agg')
# import pylab
import math

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
    plottingActive = True
except ImportError:
    plottingActive = False


from pyFiles.peptide import Peptide
from pyFiles.dripEncoding import load_drip_means, return_b_y_ions, return_b_y_ions_lowres, return_b_y_ions_var_mods, return_b_y_ions_lowres_var_mods

validPeps = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')-set('JOBZUX')
class PSM(object):
    """ Simple PSM class.  Far fewer fields then DRIP PSM class, meant for
    general PSMs which may be read from files and fields with low memory footprint,
     i.e., fields do not contain any lists
     
     other field contains all other fields stored in an output file
    """

    def __init__(self, sequence = '', score = float("-inf"), 
                 scan = -1, kind = 't', charge = 0, other = {},
                 protein = '', specId = ''):
        # we don't pass in the observed spectrum during initilization since we'll
        # likely be creating a PSM instance for all candidates per spectrum
        # Also, it makes sense to do this for the observed spectrum and not the 
        # other sequences (insertions, fragments) since the other sequences are
        # peptide dependent and, when finding the top N PSMs per spectrum,
        # the observed spectrum does not change

        # todo: decide whether we should just do away with scan being a member, but
        # do to the above this becomes tricky
        self.peptide = sequence
        self.scan = scan
        self.score = score
        self.kind = kind
        self.charge = charge
        self.other  = other 
        self.protein = protein # percolator file specific field
        self.specId = specId # percolator file specific field
        self.left_flanking_aa = ''
        self.right_flanking_aa = ''

        # best to make this default behavior so that, if simply a sequence of amino acids
        # is specified, we support adding that field in and, if flanking information is included,
        # we support addding flanking information.  The latter is most commonly encountered in
        # program outputs.
        # todo: Add support for reading modifications from an input file
        if len(sequence.split('.')) > 1: # flanking information included, split string
            s  = sequence.split('.')
            # should be 3 strings after the split
            # do some checking to make sure flanking amino acids are valid
            if s[0] in validPeps:
                self.left_flanking_aa = s[0]
            else:
                self.left_flanking_aa = '-'
            if s[-1] in validPeps:
                self.right_flanking_aa = s[-1]
            else:
                self.right_flanking_aa = '-'
            self.peptide = s[1]

    def __hash__(self):
        return hash((self.scan, self.peptide))

    def __str__(self):
        return "scan%d-%s" % (self.scan, self.peptide)

class dripPSM(object):
    """ DRIP Peptide-Spectrum Match class.
    """

    def __init__(self, sequence = '', score = float("-inf"), 
                 scan = -1, kind = 't', charge = 0, numBY = -1, 
                 fragments = [], insertions = [], usedFragments = {},
                 flanking_nterm = '', flanking_cterm = '', var_mod_sequence = '',
                 protein = -1):
        # we don't pass in the observed spectrum during initilization since we'll
        # likely be creating a PSM instance for all candidates per spectrum
        # Also, it makes sense to do this for the observed spectrum and not the 
        # other sequences (insertions, fragments) since the other sequences are
        # peptide dependent and, when finding the top N PSMs per spectrum,
        # the observed spectrum does not change

        # todo: decide whether we should just do away with scan being a member, but
        # do to the above this becomes tricky
        self.peptide = sequence
        self.spectrum = None
        self.scan = scan
        self.score = score
        self.kind = kind
        self.charge = charge
        self.num_ions = numBY
        self.ion_sequence = fragments
        self.insertion_sequence = insertions
        self.used_ions = usedFragments
        self.bions = {} # only used by dripToolKit
        self.yions = {} # only used by dripToolKit
        self.flanking_nterm = flanking_nterm # used by dripSearch, dripExtract
        self.flanking_cterm = flanking_cterm # used by dripSearch, dripExtract
        self.var_mod_sequence = var_mod_sequence # currently only used by dripSearch
        self.protein = protein

        if len(self.insertion_sequence) != len(self.ion_sequence):
            raise ValueError('Sequence of fragments not equal to number of insertions')

    def __cmp__(self,other):
        if self.scan == other.scan:
            if self.score < other.score:
                return -1
            else:
                return 1
        else:
            return 0

    def __hash__(self):
        return hash((self.scan, self.peptide))

    def __str__(self):
        return "scan%d-%s" % (self.scan, self.peptide)

    @property
    def num_frames(self):
        """Number of DRIP frames for this PSM"""
        return len(self.insertion_sequence)

    @property
    def num_used_ions(self):
        """Number of non-deleted theoretical peaks"""
        return len(self.used_ions)

    @property
    def length(self):
        """Number of amino acids in the peptide. Read-only computed property."""
        return len(self.peptide)

    @property
    def num_dels(self):
        """Number of amino acids in the peptide. Read-only computed property."""
        return self.num_ions - self.num_used_ions

    def add_obs_spectrum(self, spectrum):
        if len(spectrum.mz) != self.num_frames:
            raise ValueError('PSM sid %d: Sequence of %d insertions not equal to number of %d observed frames for scan %d' % (self.scan,
                                                                                                                              self.num_frames,
                                                                                                                              len(spectrum.mz),
                                                                                                                              spectrum.spectrum_id))
        self.spectrum = spectrum

    def calc_by_sets(self, c,
                     mods = {}, ntermMods = {}, ctermMods = {},
                     highResMs2 = False, ion_to_index_map = {}, 
                     varMods = {}, varNtermMods = {}, varCtermMods = {},
                     varModSequence = ''):
        """ Used by dripToolKit for plotting purposes; 
            the sequences of b- and y-ions must be recomputed to tell
            which set each fragment ion belongs to
        """
        if varMods or ntermVarMods or ctermVarMods:
            assert varModSequence, "Variable modifications enyme options specified, but string indicating which amino acids were var mods not supplied.  Exitting"
            if highResMs2:
                bions, yions = return_b_y_ions_var_mods(Peptide(self.peptide), c, 
                                                        mods, ntermMods, ctermMods,
                                                        ion_to_index_map,
                                                        varMods, varNtermMods, varCtermMods,
                                                        varModSequence)
            else:
                bions, yions = return_b_y_ions_lowres_var_mods(Peptide(self.peptide), c, 
                                                               mods, ntermMods, ctermMods,
                                                               ion_to_index_map,
                                                               varMods, varNtermMods, varCtermMods,
                                                               varModSequence)
        else:
            if highResMs2:
                bions, yions = return_b_y_ions(Peptide(self.peptide), c, mods,
                                               ntermMods, ctermMods,
                                               ion_to_index_map)
            else:
                bions, yions = return_b_y_ions_lowres(Peptide(self.peptide), c, mods,
                                                      ntermMods, ctermMods,
                                                      ion_to_index_map)

        self.bions = bions
        self.yions = yions

    def calculate_drip_features(self, drip_means):
        """Calculate all PSM features from the decoded DRIP Viterbi path for 
           post-processor (i.e., Percolator) use
        """
        sum_scored_intensities = 0.0
        sum_scored_mz_dists = 0.0
        num_ins = 0
        num_frames = self.num_frames

        num_dels = self.num_dels
        num_non_dels = self.num_used_ions

        for intensity,mz,ins,fragment in zip(self.spectrum.intensity,
                                             self.spectrum.mz,
                                             self.insertion_sequence,
                                             self.ion_sequence):
            if not ins:
                sum_scored_intensities += intensity
                sum_scored_mz_dists += abs(mz - drip_means[fragment])
            else:
                num_ins += 1

        # calculate number of non-insertions
        num_non_ins = self.num_frames - num_ins

        return num_ins, num_dels, num_non_ins, num_non_dels, \
            sum_scored_intensities, sum_scored_mz_dists

    if plottingActive:
        def plot_drip_viterbi(self, plotname,
                          xtick_increments = 100):
            """ Plot PSM's DRIP Viterbi path
            """
            peak_colors = 'kbr'
            s = self.spectrum
            # x axis tick labels
            lower_mz = int(math.floor(min(s.mz)/xtick_increments)*xtick_increments)
            upper_mz = int(math.ceil(max(s.mz)/xtick_increments + 1)*xtick_increments)
            ticks = range(lower_mz, upper_mz, xtick_increments)

            pylab.figure()

            pylab.xlabel('m/z', fontsize=20)
            pylab.spectral()

            # unit-normalize
            max_intensity = max(s.intensity)
            intensity = [intensity/max_intensity for intensity in s.intensity]
            s.intensity = intensity

            eps = 7
            bheight = 0.1
            yheight = 0.1

            f = pylab.gca()
            pylab.ylabel('intensity', fontsize = 20)
            pylab.spectral()

            inserted_mz = []
            inserted_intensity = []
            scored_bions_mz = []
            scored_bions_intensity = []
            scored_bions_numRes = []
            scored_yions_mz = []
            scored_yions_intensity = []
            scored_yions_numRes = []

            for i, (mz, intensity, ins, ion) in enumerate(zip(self.spectrum.mz, self.spectrum.intensity, 
                                                              self.insertion_sequence, self.ion_sequence)):
                if ins:
                    inserted_mz.append(mz)
                    inserted_intensity.append(intensity)
                else:
                    if ion in self.yions:
                        scored_yions_mz.append(mz)
                        scored_yions_intensity.append(intensity)
                        scored_yions_numRes.append(self.yions[ion])
                    else:
                        assert ion in self.bions, "Non-deleted ion %d not in theoretical spectrum of PSM, exitting" % ion
                        scored_bions_mz.append(mz)
                        scored_bions_intensity.append(intensity)
                        scored_bions_numRes.append(self.bions[ion])

            annotate_shift = 0.1

            # plot spectrum
            pylab.vlines(inserted_mz, 0, inserted_intensity, color = '0.3' , linestyles='solid', hold = True)

            if scored_bions_mz: # make sure scored b-ions is not empty
                pylab.vlines(scored_bions_mz, 0, scored_bions_intensity, color = 'b', linestyles = 'solid', hold = True)
                for i, mz, intensity in zip(scored_bions_numRes, scored_bions_mz, scored_bions_intensity):
                    ion_str = 'b%d' % i[0]
                    for c in range(i[1]):
                        ion_str += '+'
                    pylab.annotate('%s' % ion_str, xy=(mz, intensity), xytext=(mz-eps, intensity+0.01), color='b')
            
            if scored_yions_mz: # make sure scored y-ions is not empty
                pylab.vlines(scored_yions_mz, 0, scored_yions_intensity, color = 'r', linestyles = 'solid', hold = True)
                for i, mz, intensity in zip(scored_yions_numRes, scored_yions_mz, scored_yions_intensity):
                    ion_str = 'y%d' % i[0]
                    for c in range(i[1]):
                        ion_str += '+'
                    pylab.annotate('%s' % ion_str, xy=(mz, intensity), xytext=(mz-eps, intensity+0.01), color='r')

            if scored_bions_mz:
                if scored_yions_mz:
                    pylab.legend(('insertions', 'b-ions', 'y-ions'), loc = 'center left')
                else:
                    pylab.legend(('insertions', 'b-ions'), loc = 'center left')
            else:
                if scored_yions_mz:
                    pylab.legend(('insertions', 'y-ions'), loc = 'center left')
                else:
                    pylab.legend(('insertions'), loc = 'center left')

            pylab.ylim(0, 1.1)
            pylab.xlim(ticks[0], ticks[-1])
            pylab.xticks(ticks)
            # pylab.setp(f.get_xticklabels(), visible=False)

            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(12.5,6.5)    
            pylab.savefig(plotname, dpi=100)

    else:
        def plot_drip_viterbi(self, plotname,
                          xtick_increments = 100):
            """ Plot PSM's DRIP Viterbi path
            """
            print "Failed to import matplotlib, which is necessary to generate figures.  Exitting"
            exit(-1)

############### read candidate peptides written to binary
def read_binary_peptide_tuples(s, f):
    """ Read binary records from file
        inputs:
        s - compiled struct, instance of struct.Struct
        f - input file stream, assumed open for binary reading
    """
    chunks = iter(lambda: f.read(s.size), b'')
    return (s.unpack(chunk) for chunk in chunks)
        
############### function to calculate DRIP PSMs and write to output file
def write_dripPSM_to_pin(fid, psm, drip_means, fields,
                         dm, absDm, psm_mass,
                         enzN, enzC,
                         peptide_str, parent_protein,
                         one_hot_charge):
    """ header: (1)SpecId
                (2)Label
                (3)ScanNr
                (4)score
                (5)Mass
                (6)PepLen
                (7)dm
                (8)absdM
                (9)Insertions
                (10)Deletions
                (11)NumObsPeaksScored
                (12)NumTheoPeaksUsed
                (13)SumObsIntensities
                (14)SumScoredMzDist
                (15)enzN
                (16)enzC
                (17)Charge1
                (18)Charge2
                (19)Charge3
                ...
                (16+N)ChargeN
                (16+N+1)Peptide
                (16+N+2)Protein
    """
    num_ins, num_dels, num_non_ins, num_non_dels, sum_scored_intensities, sum_scored_mz_dists = psm.calculate_drip_features(drip_means)
    # create SpecId string
    if psm.kind == 't':
        s = 'target-' + 'scan' + str(psm.scan) + 'charge' + str(psm.charge)
    else:
        s = 'decoy-' + 'scan' + str(psm.scan) + 'charge' + str(psm.charge)
    for f in fields[:-1]:
        if f=='SpecId':
            fid.write("%s\t" % s)
        elif f=='Label':
            if psm.kind == 't':
                fid.write("1\t")
            else:
                fid.write("-1\t")
        elif f=='ScanNr':
                fid.write("%d\t" % psm.scan)
        elif f=='score':
                fid.write("%f\t" % psm.score)
        elif f=='Mass':
                fid.write("%f\t" % psm_mass)
        elif f=='PepLen':
                fid.write("%d\t" % len(psm.peptide))
        elif f=='dm':
                fid.write("%f\t" % dm)
        elif f=='absdM':
                fid.write("%f\t" % absDm)
        elif f=='Insertions':
                fid.write("%d\t" % num_ins)
        elif f=='Deletions':
                fid.write("%d\t" % num_dels)
        elif f=='NumObsPeaksScored':
                fid.write("%d\t" % num_non_ins)
        elif f=='NumTheoPeaksUsed':
                fid.write("%d\t" % num_non_dels)
        elif f=='SumObsIntensities':
                fid.write("%f\t" % sum_scored_intensities)
        elif f=='SumScoredMzDist':
                fid.write("%f\t" % sum_scored_mz_dists)
        elif f=='enzN':
                fid.write("%d\t" % enzN)
        elif f=='enzC':
                fid.write("%d\t" % enzC)
        elif f=='Peptide':
                fid.write("%s\t" % peptide_str)
        else: # should only be charge fields left
            try:
                fid.write("%d\t" % one_hot_charge[f])
            except KeyError:
                print "Field %s not valid pin field, exitting" % f
    # protein field is last
    fid.write("%s\n" % (parent_protein))
    return 1

def write_appended_dripPSM_to_pin(fid, psm, psm0,
                                  drip_means, 
                                  fields, otherFields, 
                                  peptide_str, parent_protein):
    """ header: (1)SpecId
                (2)Label
                (3)ScanNr
                (4)Insertions
                (5)Deletions
                (6)NumObsPeaksScored
                (7)NumTheoPeaksUsed
                (8)SumObsIntensities
                (9)SumScoredMzDist
                (10)User input PIN field 1
                (11)User input PIN field 2
                ...
                (9+N)User input PIN field N
                (9+N+1)Peptide
                (9+N+2)Protein
    """
    num_ins, num_dels, num_non_ins, num_non_dels, sum_scored_intensities, sum_scored_mz_dists = psm.calculate_drip_features(drip_means)
    # create SpecId string
    if psm.kind == 't':
        s = 'target-' + 'scan' + str(psm.scan) + 'charge' + str(psm.charge)
    else:
        s = 'decoy-' + 'scan' + str(psm.scan) + 'charge' + str(psm.charge)
    for f in fields:
        if f=='SpecId':
            fid.write("%s\t" % s)
        elif f=='Label':
            if psm.kind == 't':
                fid.write("1\t")
            else:
                fid.write("-1\t")
        elif f=='ScanNr':
                fid.write("%d\t" % psm.scan)
        elif f=='Insertions':
                fid.write("%d\t" % num_ins)
        elif f=='Deletions':
                fid.write("%d\t" % num_dels)
        elif f=='NumObsPeaksScored':
                fid.write("%d\t" % num_non_ins)
        elif f=='NumTheoPeaksUsed':
                fid.write("%d\t" % num_non_dels)
        elif f=='SumObsIntensities':
                fid.write("%f\t" % sum_scored_intensities)
        elif f=='SumScoredMzDist':
                fid.write("%f\t" % sum_scored_mz_dists)
        else:
            print "Non-standard DRIP field for appending to PIN files encountered, exitting"
            exit(-1)
    # write out other features from input PIN
    if otherFields:
        try:
            other = psm0.other
        except:
            print "Could not access PSM other member, exitting"
            exit(-1)
    for f in otherFields:
        fid.write("%s\t" % other[f])

    # protein field is last
    fid.write("%s\t%s\n" % (peptide_str, parent_protein))
    return 1

def write_dripPSM_to_ident(fid, psm, drip_means):
    """ header: (1)Kind
                (2)Scan
                (3)Frames
                (4)Score
                (5)Peptide
                (6)Obs_Inserts
                (7)Theo_Deletes
                (8)Obs_peaks_scored
                (9)Theo_peaks_used
                (10)Sum_obs_intensities
                (11)Sum_scored_mz_dist
                (12)Charge
    """
    num_ins, num_dels, num_non_ins, num_non_dels, sum_scored_intensities, sum_scored_mz_dists = psm.calculate_drip_features(drip_means)
    try:
        fid.write('%c\t%d\t%d\t%f\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\n' % (psm.kind, 
                                                                        psm.scan,
                                                                        psm.num_frames,
                                                                        psm.score,
                                                                        psm.peptide,
                                                                        num_ins,
                                                                        num_dels,
                                                                        num_non_ins, 
                                                                        num_non_dels,
                                                                        sum_scored_intensities,
                                                                        sum_scored_mz_dists,
                                                                        psm.charge))
    except IOError:
        print "Could not write to ident stream, exitting"
        exit(-1)

    return 1

def write_dripPSM_to_ident_var_mods(fid, psm, drip_means,
                                    mods, nterm_mods, cterm_mods,
                                    var_mods, nterm_var_mods, cterm_var_mods, 
                                    isVarMods = 0):
    """ old header: (1)Kind
                    (2)Scan
                    (3)Frames
                    (4)Score
                    (5)Peptide
                    (6)Obs_Inserts
                    (7)Theo_Deletes
                    (8)Obs_peaks_scored
                    (9)Theo_peaks_used
                    (10)Sum_obs_intensities
                    (11)Sum_scored_mz_dist
                    (12)Charge
                    (13)Flanking_nterm
                    (14)Flanking_ncterm

       header as of 1/4/2016: (1)Kind
                              (2)Scan
                              (3)Score
                              (4)Peptide
                              (5)Obs_Inserts
                              (6)Theo_Deletes
                              (7)Obs_peaks_scored
                              (8)Theo_peaks_used
                              (9)Sum_obs_intensities
                              (10)Sum_scored_mz_dist
                              (11)Charge
                              (12)Flanking_nterm
                              (13)Flanking_ncterm
                              (14)Protein_id

       header as of 3/15/2016: (1)Kind
                              (2)Scan
                              (3)Score
                              (4)Peptide
                              (5)Obs_Inserts
                              (6)Theo_Deletes
                              (7)Obs_peaks_scored
                              (8)Theo_peaks_used
                              (9)Sum_obs_intensities
                              (10)Sum_scored_mz_dist
                              (11)Charge
                              (12)Flanking_nterm
                              (13)Flanking_ncterm
                              (14)Protein_id
                              (15)Var_mod_seq (only apppears if variable mods selected)
                    
    """
    num_ins, num_dels, num_non_ins, num_non_dels, sum_scored_intensities, sum_scored_mz_dists = psm.calculate_drip_features(drip_means)
    try:
        var_mod_string = psm.var_mod_sequence.split('\x00')[0]

        c = psm.peptide[0]
        vm = psm.var_mod_sequence[0]
        pep_str = [c]

        # check n-term and c-term
        if c in nterm_mods:
            pep_str += str('[%1.0e]' % nterm_mods[c])
        elif vm == '2':
            pep_str += str('[%1.0e]' % nterm_var_mods[c][1])
        elif c in mods:
            pep_str += str('[%1.0e]' % mods[c])

        for c, vm in zip(psm.peptide[1:-1], psm.var_mod_sequence[1:-1]):
            pep_str += c
            if c in mods:
                pep_str += str('[%1.0e]' % mods[c])
            elif vm == '1': # variable mods
                pep_str += str('[%1.0e]' % var_mods[c][1])

        c = psm.peptide[-1]
        vm = psm.var_mod_sequence[-1]
        pep_str += [c]
        # check n-term and c-term
        if c in cterm_mods:
            pep_str += str('[%1.0e]' % cterm_mods[c])
        elif vm == '3':
            pep_str += str('[%1.0e]' % cterm_var_mods[c][1])
        elif c in mods:
            pep_str += str('[%1.0e]' % mods[c])

        if not isVarMods:
            fid.write('%c\t%d\t%f\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%c\t%c\t%s\n' % (psm.kind, 
                                                                                    psm.scan,
                                                                                    psm.score,
                                                                                    ''.join(pep_str),
                                                                                    num_ins,
                                                                                    num_dels,
                                                                                    num_non_ins, 
                                                                                    num_non_dels,
                                                                                    sum_scored_intensities,
                                                                                    sum_scored_mz_dists,
                                                                                    psm.charge,
                                                                                    psm.flanking_nterm,
                                                                                    psm.flanking_cterm,
                                                                                    psm.kind + str(psm.protein)))
        else:
            fid.write('%c\t%d\t%f\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%c\t%c\t%s\t%s\n' % (psm.kind, 
                                                                                        psm.scan,
                                                                                        psm.score,
                                                                                        ''.join(pep_str),
                                                                                        num_ins,
                                                                                        num_dels,
                                                                                        num_non_ins, 
                                                                                        num_non_dels,
                                                                                        sum_scored_intensities,
                                                                                        sum_scored_mz_dists,
                                                                                        psm.charge,
                                                                                        psm.flanking_nterm,
                                                                                        psm.flanking_cterm,
                                                                                        psm.kind + str(psm.protein),
                                                                                        var_mod_string))

    except IOError:
        print "Could not write to ident stream, exitting"
        exit(-1)

    return 1

############### dripSearch processing
def parse_drip_segments_binDb(s, filename, dripMeans, 
                              topMatch, sid, 
                              currSpec, numPeps, 
                              pepDBlist, recalibrate = False):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    try:
        f = open(filename, 'r')
    except IOError:
        return [(None, float("-inf")), (None, float("-inf"))]
    log = f.read()
    f.close()

    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    target_psms = {}
    decoy_psms = {}
    if recalibrate:
        recal_decoy_psms = {}

    pepLookUp = open(pepDBlist, "rb")
    # (0) - peptide type: 1 = target, 2 = decoy, 3 = recalibrate decoy
    # (1) - peptide string
    # (2) - num b- and y-ions
    # (3) - charge
    # (4) - amino acid flanking n-terminus
    # (5) - amino acid flanking c-terminus
    # (6) - string denoting variable modded amino acids
    # (7) - integer denoting protein number
    pepDB = list(read_binary_peptide_tuples(s, pepLookUp))
    pepLookUp.close()
    
    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        ck = pepDB[curr_segment][0]
        if ck == 1:
            curr_kind = 't'
        elif ck == 2:
            curr_kind = 'd'
        elif ck == 3:
            curr_kind = 'r'
        else:
            print "Encountered peptide type outside of [1,3], exitting"
            exit(-1)
        curr_pep_seq = pepDB[curr_segment][1].split('\x00')[0]
        curr_numby = pepDB[curr_segment][2]
        curr_charge = pepDB[curr_segment][3]
        curr_flanking_nterm = pepDB[curr_segment][4]
        curr_flanking_cterm = pepDB[curr_segment][5]
        curr_var_mod_sequence = pepDB[curr_segment][6]
        curr_protein_id = pepDB[curr_segment][7]
            
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        if curr_frames != len(currSpec.mz):
            print "Spectrum %d number of frames does not match up with number of observations, exitting" % sid
            exit(-1)

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # log all PSM info
        currPsm = dripPSM(curr_pep_seq, curr_score,
                          sid, curr_kind, curr_charge, curr_numby,
                          fragments, ins_sequence, used_peaks,
                          curr_flanking_nterm, curr_flanking_cterm, curr_var_mod_sequence,
                          curr_protein_id)
        if(curr_kind=='t'):
            if curr_charge not in target_psms:
                target_psms[curr_charge] = []
            target_psms[curr_charge].append(currPsm)
        elif(curr_kind=='d'):
            if curr_charge not in decoy_psms:
                decoy_psms[curr_charge] = []
            decoy_psms[curr_charge].append(currPsm)
        elif(curr_kind=='r'):
            if recalibrate:
                if curr_charge not in recal_decoy_psms:
                    recal_decoy_psms[curr_charge] = []
                recal_decoy_psms[curr_charge].append(currPsm)
            else:
                print "Encountered recalibration decoy but recalibrate set to false, exitting"
                exit(-1)
                           
        count += 1

    if count != numPeps: #didn't read all spectra scores, return empty scores
        # raise ValueError(" Number of decoded peptides not equal to number of encoded peptides")
        # don't want to exit since we may still completed runs not yet considered
        return [dripPSM('', float("-inf"),sid, 't'), dripPSM('', float("-inf"),sid, 'd')]

    sid_td_psms = []
    for charge in target_psms:
        # return topMatch number of targets and decoys
        target_psms[charge].sort(key = lambda x: x.score, reverse = True)
        decoy_psms[charge].sort(key = lambda x: x.score, reverse = True)
        if recalibrate:
            recal_decoy_psms[charge].sort(key = lambda x: x.score, reverse = True)
        for i in range(min(topMatch, len(target_psms[charge]))):
            currPsm = target_psms[charge][i]
            currPsm.add_obs_spectrum(currSpec)
            # currPsm.spectrum = currSpec
            sid_td_psms.append(currPsm)

        for i in range(min(topMatch, len(decoy_psms[charge]))):
            currPsm = decoy_psms[charge][i]
            # currPsm.spectrum = currSpec
            currPsm.add_obs_spectrum(currSpec)
            sid_td_psms.append(currPsm)

        if recalibrate:
            for i in range(min(topMatch, len(recal_decoy_psms[charge]))):
                currPsm = recal_decoy_psms[charge][i]
                currPsm.add_obs_spectrum(currSpec)
                sid_td_psms.append(currPsm)
                
    return sid_td_psms

def parse_drip_segments(filename, dripMeans, 
                        topMatch, sid, 
                        currSpec, numPeps, 
                        pepDBlist):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    try:
        f = open(filename, 'r')
    except IOError:
        return [(None, float("-inf")), (None, float("-inf"))]
    log = f.read()
    f.close()

    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    target_psms = {}
    decoy_psms = {}

    pepLookUp = open(pepDBlist, "r")
    # header: (1)Kind (2)Peptide (3) NumBY (4) Charge
    pepDB = [pepRow for pepRow in csv.DictReader(pepLookUp, delimiter = '\t')]
    pepLookUp.close()

    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        curr_kind = pepDB[curr_segment]['Kind']
        curr_pep_seq = pepDB[curr_segment]['Peptide']

        try:
            curr_numby = int(pepDB[curr_segment]['NumBY'])
        except ValueError:
            print "Trouble trying to convert numBy %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['NumBY'])

        try:
            curr_charge = int(pepDB[curr_segment]['Charge'])
        except ValueError:
            print "Trouble trying to convert charge %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['Peptide'])
            
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        if curr_frames != len(currSpec.mz):
            print "Spectrum %d number of frames does not match up with number of observations, exitting" % sid
            exit(-1)

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # log all PSM info
        currPsm = dripPSM(curr_pep_seq, curr_score,
                          sid, curr_kind, curr_charge, curr_numby,
                          fragments, ins_sequence, used_peaks)
        if(curr_kind=='t'):
            if curr_charge not in target_psms:
                target_psms[curr_charge] = []
            target_psms[curr_charge].append(currPsm)
        else:
            if curr_charge not in decoy_psms:
                decoy_psms[curr_charge] = []
            decoy_psms[curr_charge].append(currPsm)
                           
        count += 1

    if count != numPeps: #didn't read all spectra scores, return empty scores
        # raise ValueError(" Number of decoded peptides not equal to number of encoded peptides")
        # don't want to exit since we may still completed runs not yet considered
        return [dripPSM('', float("-inf"),sid, 't'), dripPSM('', float("-inf"),sid, 'd')]

    sid_td_psms = []
    for charge in target_psms:
        # return topMatch number of targets and decoys
        target_psms[charge].sort(key = lambda x: x.score, reverse = True)
        decoy_psms[charge].sort(key = lambda x: x.score, reverse = True)
        for i in range(min(topMatch, len(target_psms[charge]))):
            currPsm = target_psms[charge][i]
            currPsm.add_obs_spectrum(currSpec)
            # currPsm.spectrum = currSpec
            sid_td_psms.append(currPsm)

        for i in range(min(topMatch, len(decoy_psms[charge]))):
            currPsm = decoy_psms[charge][i]
            # currPsm.spectrum = currSpec
            currPsm.add_obs_spectrum(currSpec)
            sid_td_psms.append(currPsm)

    return sid_td_psms

def parse_dripsearch_persidBin_generator(logDir, topMatch, highResMs2, 
                                         spec_dict, meanFile,
                                         pepDBperSidList, output_dir,
                                         max_length, recalibrate = False):
    """ Parse all DRIP Viterbi paths for a dataset
    """

    # serialized candidate peptide data have fields:
    # (1) - peptide type: 1 = target, 2 = decoy, 3 = recalibrate decoy
    # (2) - peptide string
    # (3) - num b- and y-ions
    # (4) - charge
    # (5) - amino acid flanking n-terminus
    # (6) - amino acid flanking c-terminus
    # (7) - string denoting variable modded amino acids
    # (8) - integer denoting protein number
    s = struct.Struct('I %ds I I s s %ds I' % (max_length, max_length))

    # means are constant for low-res
    if not highResMs2:
        dripMeans = load_drip_means(meanFile)

    sid_numPeps_lookup = open(pepDBperSidList, "r")
    reader = csv.DictReader(sid_numPeps_lookup, delimiter = '\t')

    target_peptide = None
    decoy_peptide = None
    for row in reader:
        sid = int(row['sid'])
        currSpec = spec_dict[sid]
        numPeps = int(row['numPeps'])
        testOutputFile = '%s/vitVals-sid%d.txt' % (logDir, sid)
        pepDBlist = '%s/sid%d-pepDB.bin' % (output_dir, sid)

        if highResMs2:
            # load this spectrum's collection of means
            mean_file = meanFile[:-4] + '-sid' + str(sid) + '.txt'
            dripMeans = load_drip_means(mean_file)

        td = parse_drip_segments_binDb(s, testOutputFile, dripMeans, 
                                       topMatch, sid, 
                                       currSpec, numPeps, 
                                       pepDBlist, recalibrate)

        yield td, dripMeans

def parse_dripsearch_persid_generator(logDir, topMatch, highResMs2, 
                                      spec_dict, meanFile,
                                      pepDBperSidList, output_dir):
    """ Parse all DRIP Viterbi paths for a dataset
    """
    # means are constant for low-res
    if not highResMs2:
        dripMeans = load_drip_means(meanFile)

    sid_numPeps_lookup = open(pepDBperSidList, "r")
    reader = csv.DictReader(sid_numPeps_lookup, delimiter = '\t')

    target_peptide = None
    decoy_peptide = None
    for row in reader:
        sid = int(row['sid'])
        currSpec = spec_dict[sid]
        numPeps = int(row['numPeps'])
        testOutputFile = '%s/vitVals-sid%d.txt' % (logDir, sid)
        pepDBlist = '%s/sid%d-pepDB.txt' % (output_dir, sid)

        if highResMs2:
            # load this spectrum's collection of means
            mean_file = meanFile[:-4] + '-sid' + str(sid) + '.txt'
            dripMeans = load_drip_means(mean_file)

        td = parse_drip_segments(testOutputFile, dripMeans, 
                                 topMatch, sid, 
                                 currSpec, numPeps, 
                                 pepDBlist)

        yield td, dripMeans

def write_dripSearch_ident(output, logDir, topMatch, 
                           highResMs2, spec_dict, meanFile,
                           pepDBperSidList, output_dir,
                           max_length, recalibrate,
                           mods, nterm_mods, cterm_mods,
                           var_mods, nterm_var_mods, cterm_var_mods):
    try:
        identFid = open(output, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    isVarMods = len(var_mods) + len(nterm_var_mods) + len(cterm_var_mods)

    # identFid.write('Kind\tScan\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\tFlanking_nterm\tFlanking_cterm\tProtein_id\n')

    if isVarMods:
        identFid.write('Kind\tScan\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\tFlanking_nterm\tFlanking_cterm\tProtein_id\tVar_mod_seq\n')
    else:
        identFid.write('Kind\tScan\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\tFlanking_nterm\tFlanking_cterm\tProtein_id\n')    

    for td, drip_means in parse_dripsearch_persidBin_generator(logDir, topMatch, 
                                                               highResMs2, spec_dict, meanFile,
                                                               pepDBperSidList, output_dir,
                                                               max_length, recalibrate):
        for psm in td:
            write_dripPSM_to_ident_var_mods(identFid, psm, drip_means,
                                            mods, nterm_mods, cterm_mods,
                                            var_mods, nterm_var_mods, cterm_var_mods,
                                            isVarMods)

    identFid.close()
    return 1

################ dripExtract processing
#### should just rename the below san _parallel and make passing target_psms and decoy_psms
#### by reference the default 
#### Todo: do the above after testing the parallel implementation of dripExtract
def parse_dripExtract_parallel(filename, pepDBlist, target_psms, decoy_psms):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    try:
        f = open(filename, 'r')
    except IOError:
        print "Could not open vitValsFile %s for reading, exitting" % filename
        exit(-1)
    log = f.read()
    f.close()

    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    pepLookUp = open(pepDBlist, "r")
    # header: (1)Kind (2)Peptide (3) NumBY (4) Charge
    pepDB = [pepRow for pepRow in csv.DictReader(pepLookUp, delimiter = '\t')]
    pepLookUp.close()
    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        sid = int(pepDB[curr_segment]['Sid'])
        curr_kind = pepDB[curr_segment]['Kind']
        curr_pep_seq = pepDB[curr_segment]['Peptide']

        # print "segment %d, sid %d" % (curr_segment, sid)

        try:
            curr_numby = int(pepDB[curr_segment]['NumBY'])
        except ValueError:
            print "Could not convert numBy %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], 
                                                                                          sid, 
                                                                                          pepDB[curr_segment]['NumBY'])

        try:
            curr_charge = int(pepDB[curr_segment]['Charge'])
        except ValueError:
            print "Could not convert charge %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], 
                                                                                           sid, 
                                                                                           pepDB[curr_segment]['Peptide'])
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # print "%d: %d, len(ins_sequence)=%d" % (sid, curr_frames, len(ins_sequence))

        currPsm = dripPSM(curr_pep_seq, curr_score,
                          sid, curr_kind, curr_charge, curr_numby,
                          fragments, ins_sequence, used_peaks)
        if(curr_kind=='t'):
            if (sid,curr_charge) not in target_psms:
                target_psms[sid,curr_charge] = []
            target_psms[sid,curr_charge].append(currPsm)
        else:
            if (sid,curr_charge) not in decoy_psms:
                decoy_psms[sid,curr_charge] = []
            decoy_psms[sid,curr_charge].append(currPsm)
                           
        count += 1

def parse_dripExtract(filename, pepDBlist):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    try:
        f = open(filename, 'r')
    except IOError:
        print "Could not open vitValsFile %s for reading, exitting" % filename
        exit(-1)
    log = f.read()
    f.close()

    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    target_psms = {}
    decoy_psms = {}

    pepLookUp = open(pepDBlist, "r")
    # header: (1)Kind (2)Peptide (3) NumBY (4) Charge
    pepDB = [pepRow for pepRow in csv.DictReader(pepLookUp, delimiter = '\t')]
    pepLookUp.close()
    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        sid = int(pepDB[curr_segment]['Sid'])
        curr_kind = pepDB[curr_segment]['Kind']
        curr_pep_seq = pepDB[curr_segment]['Peptide']

        # print "segment %d, sid %d" % (curr_segment, sid)

        try:
            curr_numby = int(pepDB[curr_segment]['NumBY'])
        except ValueError:
            print "Could not convert numBy %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], 
                                                                                          sid, 
                                                                                          pepDB[curr_segment]['NumBY'])

        try:
            curr_charge = int(pepDB[curr_segment]['Charge'])
        except ValueError:
            print "Could not convert charge %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], 
                                                                                           sid, 
                                                                                           pepDB[curr_segment]['Peptide'])
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # print "%d: %d, len(ins_sequence)=%d" % (sid, curr_frames, len(ins_sequence))

        currPsm = dripPSM(curr_pep_seq, curr_score,
                          sid, curr_kind, curr_charge, curr_numby,
                          fragments, ins_sequence, used_peaks)
        if(curr_kind=='t'):
            if (sid,curr_charge) not in target_psms:
                target_psms[sid,curr_charge] = []
            target_psms[sid,curr_charge].append(currPsm)
        else:
            if (sid,curr_charge) not in decoy_psms:
                decoy_psms[sid,curr_charge] = []
            decoy_psms[sid,curr_charge].append(currPsm)
                           
        count += 1

    return target_psms, decoy_psms

def write_output(targets, decoys, filename, meanFile, spec_dict):
    """ Write PSMs and features to file

    Note: since this assumes a static set of means, this is not a general function for
    writing output files given a set of targets and decoys.  That is, if dripSearch was run in 
    high-res mode, then there are several sets of means specific to each searched spectrum and these
    are necessary to calculate the DRIP features.  Any other run of dripSearch and all runs of dripExtract
    may use this function to write their output since the set of means does not change per spectrum.

    Todo: we could get rid of this mean dependency, but we would need to recalculate each PSM's theoretical spectrum,
    including all the theoretical spectrum preprocessing to get it right.  We could do this efficiently by only 
    calculating these quantities (i.e., the DRIP features) for the top-match PSMs per spectrum.
    """

    dripMeans = load_drip_means(meanFile)

    try:
        identFid = open(filename, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    identFid.write('Kind\tScan\tFrames\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\n')

    sidCharges = targets.keys()
    sidCharges.sort(key = lambda r: r[0])

    # for sid, charge in targets:
    for sid, charge in sidCharges:
        s = spec_dict[sid]
        for psm in targets[sid,charge]:
            psm.add_obs_spectrum(s)
            write_dripPSM_to_ident(identFid, psm, dripMeans)
        if (sid,charge) in decoys:
            for psm in decoys[sid,charge]:
                psm.add_obs_spectrum(s)
                write_dripPSM_to_ident(identFid, psm, dripMeans)
    identFid.close()

def write_percolator_pin_pepdb(targets, decoys, filename, meanFile, spec_dict,
                               cleavage_req, target_db, decoy_db):
    """ Write PSMs and features to file

    Note: since this assumes a static set of means, this is not a general function for
    writing output files given a set of targets and decoys.  That is, if dripSearch was run in 
    high-res mode, then there are several sets of means specific to each searched spectrum and these
    are necessary to calculate the DRIP features.  Any other run of dripSearch and all runs of dripExtract
    may use this function to write their output since the set of means does not change per spectrum.

    We could get rid of this mean dependency, but we would need to recalculate each PSM's theoretical spectrum,
    including all the theoretical spectrum preprocessing to get it right.  We could do this efficiently by only 
    calculating these quantities (i.e., the DRIP features) for the top-match PSMs per spectrum.
    """

    dripMeans = load_drip_means(meanFile)

    # fields = ["SpecId","Label","ScanNr","lnrSp","deltLCn","deltCn","score","Sp","IonFrac","Mass","PepLen","Charge1","Charge2","Charge3","enzN","enzC","enzInt","lnNumSP","dm","absdM", "insertions", "deletions", "peaksScoredA", "theoPeaksUsedA", "Peptide","Proteins"]

    try:
        identFid = open(filename, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    # find maximum PSM charge
    max_charge = 0
    for sid, charge in targets:
        max_charge = max(max_charge,max([psm.charge for psm in targets[sid,charge]]))
        if (sid,charge) in decoys:
            max_charge = max(max_charge,max([psm.charge for psm in decoys[sid,charge]]))

    assert max_charge > 0, "Maximum charge for encountered PSMs is zero, exitting"

    # Create list of fields to write out
    fields = ["SpecId","Label","ScanNr",
              "score","Mass","PepLen",
              "dm","absdM", "Insertions","Deletions", 
              "NumObsPeaksScored", "NumTheoPeaksUsed",
              "SumObsIntensities", "SumScoredMzDist",
              "enzN", "enzC"]
    for c in range(1,max_charge+1):
        fields.append("Charge" + str(c))
    fields += ["Peptide","Proteins"]    

    for f in fields[:-1]:
        identFid.write('%s\t' % f)
    identFid.write('%s\n' % fields[-1])

    left_of_cleavage = cleavage_req[0]
    right_of_cleavage = cleavage_req[1]
    print left_of_cleavage
    print right_of_cleavage

    for sid, charge in targets:
        s = spec_dict[sid]
        precursor_mass = {} # dict whose keys are precursor charge, entries are precursor mass
        for c, m in s.charge_lines:
            precursor_mass[c] = m
        for psm in targets[sid,charge]:
            # add spectrum to field of psm
            psm.add_obs_spectrum(s)

            # create one hot charge vector
            one_hot_charge = {}
            for c in range(1,max_charge+1):
                c_key = 'Charge' + str(c)
                if c == charge:
                    one_hot_charge[c_key] = 1
                else:
                    one_hot_charge[c_key] = 0

            # digestedPep class
            try:
                p = target_db[psm.peptide]
            except KeyError:
                print "Target PSM %s not in supplied target database, exitting" % psm.peptide
                exit(-1)

            # calculate dm
            dm = p.peptideMass - precursor_mass[charge]
            absDm = abs(dm)
            pep_seq = psm.peptide
            enzN = 0
            enzC = 0
            if p.ntermFlanking:
                if (p.ntermFlanking in left_of_cleavage) and (pep_seq[0] in right_of_cleavage):
                    enzN = 1
            if p.ctermFlanking:
                if (pep_seq[-1] in left_of_cleavage) and (p.ctermFlanking in right_of_cleavage):
                    enzC = 1

            # form peptide for output
            if p.ntermFlanking:
                peptide_str = p.ntermFlanking
            else:
                peptide_str = '-'
            peptide_str += '.' + psm.peptide + '.'
            if p.ctermFlanking:
                peptide_str += p.ctermFlanking
            else:
                peptide_str += '-'

            write_dripPSM_to_pin(identFid, psm, 
                                 dripMeans, fields,
                                 dm, absDm, p.peptideMass,
                                 enzN, enzC,
                                 peptide_str, p.proteinName,
                                 one_hot_charge)
        if (sid,charge) in decoys:
            for psm in decoys[sid,charge]:
                # add spectrum to field of psm
                psm.add_obs_spectrum(s)

                # create one hot charge vector
                one_hot_charge = {}
                for c in range(1,max_charge+1):
                    c_key = 'Charge' + str(c)
                    if c == charge:
                        one_hot_charge[c_key] = 1
                    else:
                        one_hot_charge[c_key] = 0

                # digestedPep class
                try:
                    p = decoy_db[psm.peptide]
                except KeyError:
                    print "Decoy PSM %s not in supplied decoy database, exitting" % psm.peptide
                    exit(-1)

                # calculate dm
                dm = p.peptideMass - precursor_mass[charge]
                absDm = abs(dm)
                pep_seq = psm.peptide
                enzN = 0
                enzC = 0
                if p.ntermFlanking:
                    if (p.ntermFlanking in left_of_cleavage) and (pep_seq[0] in right_of_cleavage):
                        enzN = 1
                if p.ctermFlanking:
                    if (pep_seq[-1] in left_of_cleavage)  and (p.ctermFlanking in right_of_cleavage):
                        enzC = 1

                if p.ntermFlanking:
                    peptide_str = p.ntermFlanking
                else:
                    peptide_str = '-'
                peptide_str += '.' + psm.peptide + '.'
                if p.ctermFlanking:
                    peptide_str += p.ctermFlanking
                else:
                    peptide_str += '-'

                write_dripPSM_to_pin(identFid, psm, 
                                     dripMeans, fields,
                                     dm, absDm, p.peptideMass,
                                     enzN, enzC,
                                     peptide_str, p.proteinName,
                                     one_hot_charge)
    identFid.close()

def append_to_percolator_pin(targets, decoys, 
                             targets0, decoys0,
                             filename, meanFile, spec_dict):
    """ Write PSMs and features to file

    Note: since this assumes a static set of means, this is not a general function for
    writing output files given a set of targets and decoys.  That is, if dripSearch was run in 
    high-res mode, then there are several sets of means specific to each searched spectrum and these
    are necessary to calculate the DRIP features.  Any other run of dripSearch and all runs of dripExtract
    may use this function to write their output since the set of means does not change per spectrum.

    Todo: we could get rid of this mean dependency, but we would need to recalculate each PSM's theoretical spectrum,
    including all the theoretical spectrum preprocessing to get it right.  We could do this efficiently by only 
    calculating these quantities (i.e., the DRIP features) for the top-match PSMs per spectrum.
    """

    dripMeans = load_drip_means(meanFile)

    try:
        identFid = open(filename, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    # Create list of fields to write out
    fields = ["SpecId","Label","ScanNr",
              "Insertions","Deletions", 
              "NumObsPeaksScored", "NumTheoPeaksUsed",
              "SumObsIntensities", "SumScoredMzDist"]
    endFields = ['Peptide', 'Proteins']
    # don't write out these fields.  For instance, when the field is already a mandatory Percolator field
    # so that inclusion means double counting
    nullFields = ['Protein_id'] 
    otherFields = []
    # Get extra fields from PIN file
    for p in targets0:
        otherFields = list((set(targets0[p].other.iterkeys()) - set(fields)) - set(endFields) - set(nullFields))
        break

    for f in fields:
        identFid.write('%s\t' % f)
    for f in otherFields:
        identFid.write('%s\t' % f)
    # print end fields
    identFid.write('Peptide\tProteins\n')

    for sid, charge in targets:
        s = spec_dict[sid]
        for psm in targets[sid,charge]:
            # add spectrum to field of psm to calculate sum of intensities of non-inserted peaks
            psm.add_obs_spectrum(s)
            # hash key is (scan, peptide_string)
            psm0 = targets0[psm.scan, psm.peptide] # since read in from PIN file, contains flanking info

            peptide_str = psm0.left_flanking_aa + '.' \
                + psm0.peptide + '.' \
                + psm0.right_flanking_aa

            write_appended_dripPSM_to_pin(identFid, psm, psm0,
                                          dripMeans, fields,
                                          otherFields,
                                          peptide_str, psm0.protein)
        if (sid,charge) in decoys:
            for psm in decoys[sid,charge]:
                # add spectrum to field of psm
                psm.add_obs_spectrum(s)
                # PSM hash key is (scan, peptide_string)
                psm0 = decoys0[psm.scan, psm.peptide]

                peptide_str = psm0.left_flanking_aa + '.' \
                    + psm0.peptide + '.' \
                    + psm0.right_flanking_aa

                write_appended_dripPSM_to_pin(identFid, psm, psm0,
                                              dripMeans, fields,
                                              otherFields,
                                              peptide_str, psm0.protein)
    identFid.close()

def write_percolator_pin_pepdb(targets, decoys, filename, meanFile, spec_dict,
                               cleavage_req, target_db, decoy_db):
    """ Write PSMs and features to file

    Note: since this assumes a static set of means, this is not a general function for
    writing output files given a set of targets and decoys.  That is, if dripSearch was run in 
    high-res mode, then there are several sets of means specific to each searched spectrum and these
    are necessary to calculate the DRIP features.  Any other run of dripSearch and all runs of dripExtract
    may use this function to write their output since the set of means does not change per spectrum.

    We could get rid of this mean dependency, but we would need to recalculate each PSM's theoretical spectrum,
    including all the theoretical spectrum preprocessing to get it right.  We could do this efficiently by only 
    calculating these quantities (i.e., the DRIP features) for the top-match PSMs per spectrum.
    """

    dripMeans = load_drip_means(meanFile)

    # fields = ["SpecId","Label","ScanNr","lnrSp","deltLCn","deltCn","score","Sp","IonFrac","Mass","PepLen","Charge1","Charge2","Charge3","enzN","enzC","enzInt","lnNumSP","dm","absdM", "insertions", "deletions", "peaksScoredA", "theoPeaksUsedA", "Peptide","Proteins"]

    try:
        identFid = open(filename, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    # find maximum PSM charge
    max_charge = 0
    for sid, charge in targets:
        max_charge = max(max_charge,max([psm.charge for psm in targets[sid,charge]]))
        if (sid,charge) in decoys:
            max_charge = max(max_charge,max([psm.charge for psm in decoys[sid,charge]]))

    assert max_charge > 0, "Maximum charge for encountered PSMs is zero, exitting"

    # Create list of fields to write out
    fields = ["SpecId","Label","ScanNr",
              "score","Mass","PepLen",
              "dm","absdM", "Insertions","Deletions", 
              "NumObsPeaksScored", "NumTheoPeaksUsed",
              "SumObsIntensities", "SumScoredMzDist",
              "enzN", "enzC"]
    for c in range(1,max_charge+1):
        fields.append("Charge" + str(c))
    fields += ["Peptide","Proteins"]    

    for f in fields[:-1]:
        identFid.write('%s\t' % f)
    identFid.write('%s\n' % fields[-1])

    left_of_cleavage = cleavage_req[0]
    right_of_cleavage = cleavage_req[1]
    print left_of_cleavage
    print right_of_cleavage

    for sid, charge in targets:
        s = spec_dict[sid]
        precursor_mass = {} # dict whose keys are precursor charge, entries are precursor mass
        for c, m in s.charge_lines:
            precursor_mass[c] = m
        for psm in targets[sid,charge]:
            # add spectrum to field of psm
            psm.add_obs_spectrum(s)

            # create one hot charge vector
            one_hot_charge = {}
            for c in range(1,max_charge+1):
                c_key = 'Charge' + str(c)
                if c == charge:
                    one_hot_charge[c_key] = 1
                else:
                    one_hot_charge[c_key] = 0

            # digestedPep class
            try:
                p = target_db[psm.peptide]
            except KeyError:
                print "Target PSM %s not in supplied target database, exitting" % psm.peptide
                exit(-1)

            # calculate dm
            dm = p.peptideMass - precursor_mass[charge]
            absDm = abs(dm)
            pep_seq = psm.peptide
            enzN = 0
            enzC = 0
            if p.ntermFlanking:
                if (p.ntermFlanking in left_of_cleavage) and (pep_seq[0] in right_of_cleavage):
                    enzN = 1
            if p.ctermFlanking:
                if (pep_seq[-1] in left_of_cleavage) and (p.ctermFlanking in right_of_cleavage):
                    enzC = 1

            # form peptide for output
            if p.ntermFlanking:
                peptide_str = p.ntermFlanking
            else:
                peptide_str = '-'
            peptide_str += '.' + psm.peptide + '.'
            if p.ctermFlanking:
                peptide_str += p.ctermFlanking
            else:
                peptide_str += '-'

            write_dripPSM_to_pin(identFid, psm, 
                                 dripMeans, fields,
                                 dm, absDm, p.peptideMass,
                                 enzN, enzC,
                                 peptide_str, p.proteinName,
                                 one_hot_charge)
        if (sid,charge) in decoys:
            for psm in decoys[sid,charge]:
                # add spectrum to field of psm
                psm.add_obs_spectrum(s)

                # create one hot charge vector
                one_hot_charge = {}
                for c in range(1,max_charge+1):
                    c_key = 'Charge' + str(c)
                    if c == charge:
                        one_hot_charge[c_key] = 1
                    else:
                        one_hot_charge[c_key] = 0

                # digestedPep class
                try:
                    p = decoy_db[psm.peptide]
                except KeyError:
                    print "Decoy PSM %s not in supplied decoy database, exitting" % psm.peptide
                    exit(-1)

                # calculate dm
                dm = p.peptideMass - precursor_mass[charge]
                absDm = abs(dm)
                pep_seq = psm.peptide
                enzN = 0
                enzC = 0
                if p.ntermFlanking:
                    if (p.ntermFlanking in left_of_cleavage) and (pep_seq[0] in right_of_cleavage):
                        enzN = 1
                if p.ctermFlanking:
                    if (pep_seq[-1] in left_of_cleavage)  and (p.ctermFlanking in right_of_cleavage):
                        enzC = 1

                if p.ntermFlanking:
                    peptide_str = p.ntermFlanking
                else:
                    peptide_str = '-'
                peptide_str += '.' + psm.peptide + '.'
                if p.ctermFlanking:
                    peptide_str += p.ctermFlanking
                else:
                    peptide_str += '-'

                write_dripPSM_to_pin(identFid, psm, 
                                     dripMeans, fields,
                                     dm, absDm, p.peptideMass,
                                     enzN, enzC,
                                     peptide_str, p.proteinName,
                                     one_hot_charge)
    identFid.close()
