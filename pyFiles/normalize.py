#!/usr/bin/env python
# 
# MS/MS spectra pre-processing routines.
#
# Written by John Halloran <halloj3@uw.washington.edu>, Ajit Singh <ajit.ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

from __future__ import with_statement

__authors__ = [ 'John T. Halloran <halloj3@uw.edu>', 'Ajit Singh <ajit@ee.washington.edu>' ]

import sys
from math import sqrt

_bin_width_mono = 1.0005079
_bin_offset = 0.68
_bin_width_average = 1.0011413

def integerize(value, bin_size,bo):
    return int(value / bin_size + 1 - bo)

def normalize_intensity(spectrum, tbl = 'average'):
    """Preprocessing of spectra which mimics Crux's processing for XCorr.

    In Crux preprocessing of spectra consists of the following algorithm:

    For each peak (mz, intensity):
        mz = round(mz)
        intensity = sqrt(intensity)
        normalize_each_region
    FastSEQUEST transform.

    We replicate everything except the FastSEQUEST transform. The spectrum
    object is altered in-place.

    See:
        Crux, scorer.cpp::create_intensity_array_observed

    """
    pre_mz = spectrum.precursor_mz
    experimental_mass_cutoff = spectrum.precursor_mz + 50
    bin_width = _bin_width_average if tbl == 'average' else _bin_width_mono
    max_peak = max(p for p in spectrum.mz if p < experimental_mass_cutoff)
    region_selector = int(max_peak / 10)
    max_intensity_overall = -77
    max_intensity_per_region = [0] * 10

    # initialize the intensity array. The output array will be larger than the
    # input, but the m/z range should be about the same.
    max_mz = 512
    if experimental_mass_cutoff > 512:
        a = int(experimental_mass_cutoff/1024)
        b = experimental_mass_cutoff - (1024 * a)
        max_mz = a * 1024
        if b > 0:
            max_mz = max_mz + 1024
    observed = [0] * max_mz

    for (mz, h) in zip(spectrum.mz, spectrum.intensity):
        # skip peaks larger than experimental mass.
        if mz > experimental_mass_cutoff:
            continue

        # skip peaks within precursor ion mz +/- 15 units.
        if mz < pre_mz + 15 and mz > pre_mz - 15:
            continue

        # map peak location to bin
        x = integerize(mz, bin_width, _bin_offset)
        region = int(x / region_selector)
        if region > 9:
            continue

        # sqrt-transform intensities
        y = sqrt(h)
        max_intensity_overall = max(y, max_intensity_overall)

        if observed[x] < y:
            observed[x] = y
            max_intensity_per_region[region] = max(y, max_intensity_per_region[region])

    # normalize each of the 10 regions to max intensity of 50
    region_idx = 0
    max_intensity = max_intensity_per_region[region_idx]
    for (idx, h) in enumerate(observed):
        if idx >= region_selector * (region_idx + 1) and region_idx < 9:
            region_idx = region_idx + 1
            max_intensity = max_intensity_per_region[region_idx]

        # Only normalize if there are peaks in this region. Moreover,
        # drop peaks with intensity less than 1/20 of the overall max
        # intensity. This is for compatability with SEQUEST.
        if max_intensity != 0 and observed[idx] > 0.05 * max_intensity_overall:
            observed[idx] = (observed[idx] / max_intensity) * 50
        else:
            observed[idx] = 0

        if idx > 10 * region_selector:
            break

    # new spectrum is relatively sparse: i.e., many entries of the intensity
    # vector are zero. Taking the average intensity over any m/z-window is not
    # a sensible operation.
    spectrum.mz = range(len(observed))
    spectrum.intensity = observed

    # Check postconditions
    assert(len(spectrum.mz) == len(spectrum.intensity))
    assert(all(h >= 0 and h <= 50 for h in spectrum.intensity))

def hardIntensity(spectrum):
    spectrum.rank_normalize()
    spectrum.threshold_intensity()

def log0(spectrum):
    spectrum.log_normalize()

def normLog0(spectrum):
    spectrum.normlog_normalize()

def rank_normLog0(spectrum):
    spectrum.rank_normalize()
    spectrum.normLog()

def rank_sqrt(spectrum):
    spectrum.rank_normalize()
    spectrum.sqrt()

def rank_linComp(spectrum):
    spectrum.rank_normalize()
    spectrum.linearComp_normalize()

def rank_linCompStrict(spectrum):
    spectrum.rank_normalize()
    spectrum.linearComp_normalize(a=0.5,b=0.2)

def linComp(spectrum):
    spectrum.max_normalize()
    spectrum.linearComp_normalize()

def filter0(spectrum):
    spectrum.max_normalize()

def linCompStrict(spectrum):
    spectrum.linearComp_normalize(a=0.5,b=0.2)

def norm_exp(spectrum):
    spectrum.max_normalize()
    spectrum.expNorm()

def rank_exp(spectrum):
    spectrum.rank_normalize()
    spectrum.expNorm()

def norm_expStrict(spectrum):
    spectrum.max_normalize()
    spectrum.expNorm(2)

def rank_expStrict(spectrum):
    spectrum.rank_normalize()
    spectrum.expNorm(2)

def lin1(spectrum):
    spectrum.linear_normalize(1, 1)

def lin10(spectrum):
    spectrum.linear_normalize(10, 1)

def lin20(spectrum):
    spectrum.linear_normalize(20, 1)

def expSkew1(spectrum):
    spectrum.negexp_normalize()

def expSkew2(spectrum):
    spectrum.negexp_normalize(1, -2)

def expSkew3(spectrum):
    spectrum.negexp_normalize(1, -3)

def exp2Skew3(spectrum):
    spectrum.negexp_normalize(2, -3)

def exp10Skew3(spectrum):
    spectrum.negexp_normalize(10, -3)

def lin100(spectrum):
    spectrum.linear_normalize(100, 1)

def log1(spectrum):
    """Pipeline for log-normalization.

    Removes the bottom 1/10 of the points, in terms of intensity.
    Replaces the intensity by the log(intensity) for non-zero points.

    Arguments:
       spectrum: Instance of MS2Spectrum. Mutated by the routine.

    """
    spectrum.remove_low_peaks(0.10)
    spectrum.log_normalize()

def log2(spectrum):
    spectrum.remove_low_peaks(0.20)
    spectrum.log_normalize()

def squareroot(spectrum):
    spectrum.sqrt_normalize()

def rankTop400(spectrum):
# note: the order of these two matters
    spectrum.most_intense_n(400)
    spectrum.rank_normalize()

def rankTop350(spectrum):
# note: the order of these two matters
    spectrum.most_intense_n(350)
    spectrum.rank_normalize()

def rankTop300(spectrum):
# note: the order of these two matters
    spectrum.most_intense_n(300)
    spectrum.rank_normalize()

def rankTop200(spectrum):
# note: the order of these two matters
#    spectrum.most_intense_200()
    spectrum.most_intense_n(200)
    spectrum.rank_normalize()

def rankTop100(spectrum):
# note: the order of these two matters
    spectrum.most_intense_n(100)
    spectrum.rank_normalize()

def rank0(spectrum):
    spectrum.rank_normalize()

def rank05(spectrum):
    """Pipeline for rank normalization."""
    spectrum.remove_low_peaks(0.05)
    spectrum.rank_normalize()

def rank1(spectrum):
    """Pipeline for rank normalization."""
    spectrum.remove_low_peaks(0.10)
    spectrum.rank_normalize()

def filter10(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.10)

def filter20(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.20)

def filter30(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.30)

def filter40(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.40)

def filter50(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.40)

def filter60(spectrum):
    """pipeline for grass filtering."""
    spectrum.remove_low_peaks(0.40)

def region10Normalize(spectrum):
    spectrum.local_region_normalize(10)

def sequestNormalize(spectrum):
    spectrum.sqrt_normalize()
    spectrum.region_normalize(10, 0, 2000)

def chop5percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.05)
#    spectrum.remove_low_peaks(0.05)
    spectrum.rank_normalize()

def chop5percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.05)

def chop4percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.04)
#    spectrum.remove_low_peaks(0.04)
    spectrum.rank_normalize()

def chop4percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.04)

def chop3percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.03)
#    spectrum.remove_low_peaks(0.03)
    spectrum.rank_normalize()

def chop3percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.03)

def chop2percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.02)
#    spectrum.remove_low_peaks(0.02)
    spectrum.rank_normalize()

def chop2percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.02)

def chop1percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.01)
#    spectrum.remove_low_peaks(0.01)
    spectrum.rank_normalize()

def chop1percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.01)

def chop10percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.10)
#    spectrum.remove_low_peaks(0.05)
    spectrum.rank_normalize()

def chop10percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.10)

def chop20percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.20)
#    spectrum.remove_low_peaks(0.05)
    spectrum.rank_normalize()

def chop20percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.20)

def chop30percent_rank(spectrum):
    spectrum.remove_peaks_below_fraction(0.30)
#    spectrum.remove_low_peaks(0.05)
    spectrum.rank_normalize()

def chop30percent(spectrum):
    spectrum.remove_peaks_below_fraction(0.30)

def strictSequestNormalize(spectrum):
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 15.0)

def top325TightSequest(spectrum):
    spectrum.most_intense_n(325)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top325TightLooseSequest(spectrum):
    spectrum.most_intense_n(325)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, precursor_tol = 1.5)

def top350TightSequest(spectrum):
    spectrum.most_intense_n(350)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top350TightLooseSequest(spectrum):
    spectrum.most_intense_n(350)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, precursor_tol = 1.5)

def top400TightSequest(spectrum):
    spectrum.most_intense_n(400)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top400TightLooseSequest(spectrum):
    spectrum.most_intense_n(400)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, precursor_tol = 1.5)

def top300StrictSequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 15.0)

def top300TightSequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top300TightLooseSequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, precursor_tol = 1.5)

def top200TightSequest(spectrum):
    spectrum.most_intense_n(200)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top100TightSequest(spectrum):
    spectrum.most_intense_n(100)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def top200StrictSequest(spectrum):
    spectrum.most_intense_n(200)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 15.0)

def top200TightSequest(spectrum):
    spectrum.most_intense_n(200)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, precursor_tol = 1.5)

def filterMZ200_top300Sequest(spectrum):
    spectrum.filter_lower_mz(200.0)
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def filterMZ_top300Sequest(spectrum):
    spectrum.filter_lower_mz(57.5)
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top400Sequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top300Sequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top200Sequest(spectrum):
    spectrum.most_intense_n(200)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top150Sequest(spectrum):
    spectrum.most_intense_n(150)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top100Sequest(spectrum):
    spectrum.most_intense_n(100)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top300FastSequestTransform(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)
    spectrum.fast_sequest_intensity()

def top300TrueSequestCh2(spectrum):
    spectrum.sqrt()
    spectrum.most_intense_n_thresh(300,0.05,2.0)
    spectrum.region_normalize_maxMz(10, 0, 2000, 0.0, 15.0, 2.0)

def top300TrueSequestCh3(spectrum):
    spectrum.sqrt()
    spectrum.most_intense_n_thresh(300,0.05,3.0)
    spectrum.region_normalize_maxMz(10, 0, 2000, 0.0, 15.0, 3.0)

def top300TrueSequestNoPreFiltCh2(spectrum):
    spectrum.sqrt()
    spectrum.most_intense_n_thresh(300,0.05,2.0)
    spectrum.region_normalize_maxMz(10, 0, 2000, 0.0, 0.0, 2.0)

def top300TrueSequestNoPreFiltCh3(spectrum):
    spectrum.sqrt()
    spectrum.most_intense_n_thresh(300,0.05,3.0)
    spectrum.region_normalize_maxMz(10, 0, 2000, 0.0, 0.0, 3.0)

def top300PerSpectrumSequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.spectrum_region_normalize(10, 0.05, 0.0)

def top300LowerPerSpectrumSequest(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.lower_spectrum_region_normalize(10, 0.05, 0.0)

def top400Sequest(spectrum):
    spectrum.most_intense_n(400)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top500Sequest(spectrum):
    spectrum.most_intense_n(500)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def top600Sequest(spectrum):
    spectrum.most_intense_n(600)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)


def sequest(spectrum):
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def strictSequestNormalizePrecursorFiltOff(spectrum):
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.05, 0.0)

def sequestFiltOff(spectrum):
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, 0.0)

def top300SequestFiltOff(spectrum):
    spectrum.most_intense_n(300)
    spectrum.sqrt()
    spectrum.region_normalize(10, 0, 2000, 0.0, 0.0)

def strictSequestUnnormalize(spectrum):
    spectrum.region_unnormalize(10, 0, 2000, 0.05)

def rank2(spectrum):
    spectrum.remove_low_peaks(0.20)
    spectrum.rank_normalize()

def rank3(spectrum):
    spectrum.remove_low_peaks(0.30)
    spectrum.rank_normalize()

def rank4(spectrum):
    spectrum.remove_low_peaks(0.40)
    spectrum.rank_normalize()

def rank5(spectrum):
    spectrum.remove_low_peaks(0.50)
    spectrum.rank_normalize()

def rank6(spectrum):
    spectrum.remove_low_peaks(0.60)
    spectrum.rank_normalize()

def logmax0(spectrum):
    spectrum.log_normalize()
    spectrum.max_normalize()

def logmax1(spectrum):
    spectrum.remove_low_peaks(0.10)
    spectrum.log_normalize()
    spectrum.max_normalize()

def logmax2(spectrum):
    spectrum.remove_low_peaks(0.20)
    spectrum.log_normalize()
    spectrum.max_normalize()

def maxheight(spectrum):
    spectrum.clamp_intensities(1.0)

def filter10AndClamp(spectrum):
    spectrum.remove_low_peaks(0.10)    
    spectrum.clamp_intensities(1.0)

def filter20AndClamp(spectrum):
    spectrum.remove_low_peaks(0.20)    
    spectrum.clamp_intensities(1.0)

def filter30AndClamp(spectrum):
    spectrum.remove_low_peaks(0.30)    
    spectrum.clamp_intensities(1.0)

def filter40AndClamp(spectrum):
    spectrum.remove_low_peaks(0.40)    
    spectrum.clamp_intensities(1.0)

def filter50AndClamp(spectrum):
    spectrum.remove_low_peaks(0.50)    
    spectrum.clamp_intensities(1.0)

def top400(spectrum):
    spectrum.most_intense_n(400)

def top300(spectrum):
    spectrum.most_intense_n(300)

def top200(spectrum):
    spectrum.most_intense_200()

def top150(spectrum):
    spectrum.most_intense_n(150)

def top100(spectrum):
    spectrum.most_intense_n(100)

def top50(spectrum):
    spectrum.most_intense_n(50)

def top400Normalize(spectrum):
    spectrum.most_intense_n(400)
    spectrum.max_normalize()

def top300Normalize(spectrum):
    spectrum.most_intense_n(300)
    spectrum.max_normalize()

def top200Normalize(spectrum):
    spectrum.most_intense_200()
    spectrum.max_normalize()

def top150Normalize(spectrum):
    spectrum.most_intense_n(150)
    spectrum.max_normalize()

def top100Normalize(spectrum):
    spectrum.most_intense_n(100)
    spectrum.max_normalize()

def top50Normalize(spectrum):
    spectrum.most_intense_n(50)
    spectrum.max_normalize()

def top300NormalizeSqrt(spectrum):
    spectrum.most_intense_n(300)
    spectrum.max_normalize()
    spectrum.sqrt()

def top200NormalizeSqrt(spectrum):
    spectrum.most_intense_200()
    spectrum.max_normalize()
    spectrum.sqrt()

def top150NormalizeSqrt(spectrum):
    spectrum.most_intense_n(150)
    spectrum.max_normalize()
    spectrum.sqrt()

def top100NormalizeSqrt(spectrum):
    spectrum.most_intense_n(100)
    spectrum.max_normalize()
    spectrum.sqrt()

def top50NormalizeSqrt(spectrum):
    spectrum.most_intense_n(50)
    spectrum.max_normalize()
    spectrum.sqrt()

def null(spectrum):
    pass

def pipeline(name):
    """Return a function that processes a spectrum: e.g.,

    preprocess = pipeline('log1')
    preprocess(spectrum) % applies the log1 pipeline

    """
    # print "*** Preprocessing spectra using %s" % name
    # print >> sys.stderr, "*** Preprocessing spectra using %s" % name
    return eval(name)
