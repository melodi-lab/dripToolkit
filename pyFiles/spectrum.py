#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>, Ajit Singh <ajit.ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

import gzip
import math
import mmap
import numpy
import os
import random
import re
import subprocess
import itertools

class MS2Spectrum(object):

    def __init__(self):

        self.spectrum_id = -1    # Spectrum / scan id
        self.precursor_mz = -1   # Mass report on S line (MS1)
        self.charge_lines = [ ]  # List of pairs of (charge, M+H mass)
        self.mz = [ ]            # Vector of m/z values
        self.intensity = [ ]     # Vector of intensities

    def __str__(self):
        repn = 'Spectrum id: %d' % self.spectrum_id
        return repn

    def __hash__(self):
        return hash(self.spectrum_id)

    def __cmp__(self, other):
        if self.spectrum_id < other.spectrum_id:
            return -1
        elif self.spectrum_id > other.spectrum_id:
            return 1
        else:
            return 0

    @property
    def length(self):
        return len(self.mz)

    @property
    def charges(self):
        return tuple(r[0] for r in self.charge_lines)

    def mass_at_charge(self, charge):
        for c, m in self.charge_lines:
            if c == charge:
                return m
        raise ValueError('Charge %d not present %s' % (charge,
                                                       str(self.charge_lines)))

    def write(self, f):
        f.write('S\t%d\t%d\t%.4f\n' % (self.spectrum_id, self.spectrum_id,
                                       self.precursor_mz))
        for c, m in self.charge_lines:
            f.write('Z\t%d\t%.4f\n' % (c, m))
        for v, i in zip(self.mz, self.intensity):
            f.write('%.2f %.13f\n' % (v, i))

    def threshold_intensity(self):
        """ Assume rank normalization first, so all values are unique in [0,1]
        """
        values = self.intensity
        inc = float(1/3)
        self.intensity = [ 0.0 if v <= inc else 1.0 if v <= 2*inc else 2.0 for v in values ]
        # for i in range(len(values)):
        #     if(values[i] <= inc ):
        #         values[i] = 0.0
        #     elif(values[i] <= 2*inc):
        #         values[i] = 1.0
        #     else:
        #         values[i] = 2.0

        # self.intensity = values

    def normlog_normalize(self):
        """Normalize, than Log normalize the non-zero intensities of the spectrum.

        Log-normalization reduces the dynamic range of the intensities,
        which can be useful since the range in intensities across spectra
        are not calibrated to a common scale.
        """
        imax = max(self.intensity) # max values, normalize to map intensities -> [0, 1]
        values = [i/imax for i in self.intensity]
        # find nonzero minimum
        imin = 1 #make sure it is always defined
        for i in sorted(values): #sort, find first nonzero values
            if (i > 0):
                imin = math.log(i)
                break
        if(imin < 0): 
            imin = -1*imin
        self.intensity =  [ (math.log(v)+imin)/imin if v > 0 else 0 for v in values ]

    def normLog(self):
        """Normalize, than Log normalize the non-zero intensities of the spectrum.

        Log-normalization reduces the dynamic range of the intensities,
        which can be useful since the range in intensities across spectra
        are not calibrated to a common scale.
        """
        values = self.intensity
        # find nonzero minimum
        imin = 1 #make sure it is always defined
        for i in sorted(values): #sort, find first nonzero values
            if (i > 0):
                imin = math.log(i)
                break
        if(imin < 0): 
            imin = -1*imin
        self.intensity =  [ (math.log(v)+imin)/imin if v > 0 else 0 for v in values ]

    def expNorm(self, a = 1):
        values = [ math.exp(v) for v in self.intensity ]
        imax = max(values)
        imin = min(values)
        self.intensity = [ (x-imin)/(imax-imin) for x in values ]

    def negexp_normalize(self, b = 1, skew = -1):
        """Log normalize the non-zero intensities of the spectrum.

        Log-normalization reduces the dynamic range of the intensities,
        which can be useful since the range in intensities across spectra
        are not calibrated to a common scale.
        """
        values = self.intensity
        imax = max(values)
        self.intensity =  [ b*math.exp(skew*x/imax)+(1-b*math.exp(skew)) for x in values ]

    def linear_normalize(self, b = 10, offset = 0):
        """Log normalize the non-zero intensities of the spectrum.

        Fit the normalized data x to a line y = b(1-x)+offset.
        Offset is used to shift vertically, ie set a lower bound at 1 or
        other suitable values.
        """
        values = self.intensity
        imax = max(values)
        self.intensity =  [(b*(1-peak/imax)+offset) for peak in values ]

    def linearComp_normalize(self, a = 0.7, b = 0.2):
        """Log normalize the non-zero intensities of the spectrum.

        Composition of two linear functions
        Effectively, scale weight of lower peaks with the first line and
        higher-intensity peaks with the second (strictness of scoring
        a high peak).
        """
        values = self.intensity
        imax = max(values)
        self.intensity = [ a/b*x if x < a else ((1-a)/(1-b)*(x-b)+a) for x in values ]

    def log_normalize(self):
        """Log normalize the non-zero intensities of the spectrum.

        Log-normalization reduces the dynamic range of the intensities,
        which can be useful since the range in intensities across spectra
        are not calibrated to a common scale.
        """
        values = self.intensity
        self.intensity =  [ math.log(v) if v > 0 else 0 for v in values ]

    def max_normalize(self):
        """Max normalize, so that intensities range in [0,1].
        """
        values = self.intensity
        max_val = max(values)
        self.intensity = [ v/max_val for v in values ]

    def spectrum_region_normalize(self, regions, thresh, precursor_tol = 0):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """
        min_mz = min(self.mz)
        max_mz = max(self.mz)
        values = self.intensity
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
#                bin_max = sum(values[start:stop])
                # normalize bins to unity
                bin_max = max(values[start:stop])
                self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak
        precursor_lb = self.precursor_mz - precursor_tol
        precursor_ub = self.precursor_mz + precursor_tol
        for mz, intensity in zip(self.mz, self.intensity):
            if intensity > thresh:
                if precursor_tol != 0.0:
                    if mz > precursor_ub or mz < precursor_lb:
                        processed_mz.append(mz)
                        processed_intensity.append(intensity)
                else:
                    processed_mz.append(mz)
                    processed_intensity.append(intensity)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def lower_spectrum_region_normalize(self, regions, thresh, precursor_tol = 0):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """
        min_mz = 0.0
        max_mz = max(self.mz)
        values = self.intensity
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
#                bin_max = sum(values[start:stop])
                # normalize bins to unity
                bin_max = max(values[start:stop])
                self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak
        precursor_lb = self.precursor_mz - precursor_tol
        precursor_ub = self.precursor_mz + precursor_tol
        for mz, intensity in zip(self.mz, self.intensity):
            if intensity > thresh:
                if precursor_tol != 0.0:
                    if mz > precursor_ub or mz < precursor_lb:
                        processed_mz.append(mz)
                        processed_intensity.append(intensity)
                else:
                    processed_mz.append(mz)
                    processed_intensity.append(intensity)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def region_normalize(self, regions, min_mz, max_mz, thresh, precursor_tol = 0):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """
        values = self.intensity[:]
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
#                bin_max = sum(values[start:stop])
                # normalize bins to unity
                bin_max = max(values[start:stop])
                self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak
        precursor_lb = self.precursor_mz - precursor_tol
        precursor_ub = self.precursor_mz + precursor_tol
        for mz, intensity in zip(self.mz, self.intensity):
            if intensity > thresh:
                if precursor_tol != 0.0:
                    if mz > precursor_ub or mz < precursor_lb:
                        processed_mz.append(mz)
                        processed_intensity.append(intensity)
                else:
                    processed_mz.append(mz)
                    processed_intensity.append(intensity)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def region_normalize_unnorm(self, regions, min_mz, max_mz, thresh, precursor_tol = 0):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """
        values = [math.sqrt(i) for i in self.intensity]
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
                # normalize bins to unity
                bin_max = max(values[start:stop])
                self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        precursor_lb = self.precursor_mz - precursor_tol
        precursor_ub = self.precursor_mz + precursor_tol
        for mz, intensity, i in zip(self.mz, values, self.intensity):
            if intensity > thresh:
                if precursor_tol != 0.0:
                    if mz > precursor_ub or mz < precursor_lb:
                        processed_mz.append(mz)
                        processed_intensity.append(i)
                else:
                    processed_mz.append(mz)
                    processed_intensity.append(i)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def fast_sequest_intensity(self, max_mass = 2000, bo=0.68, bw=1.0005079):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """

        max_mass = int(max(self.mz)+0.5) # do this so we don't have to worry about going beyond bin ranges
        binned_intensity = [0.0]*(max_mass+1)
        for mz, i in zip(self.mz,self.intensity):
            curr_bin = int(math.floor((mz+bo)/bw))
            if curr_bin < max_mass:
                binned_intensity[curr_bin] = max(binned_intensity[curr_bin], i)

        fst_intensity = [0.0]*(max_mass+1)
        for i, ci in enumerate(binned_intensity):
            fst_intensity[i] = ci - sum([binned_intensity[tau] for tau in range(max(1, i-75), 1+min(max_mass-1, i+75))])/(151.0)

        for i, mz in enumerate(self.mz):
            curr_bin = int(math.floor((mz+bo)/bw))
            if curr_bin < max_mass:
                self.intensity[i] = fst_intensity[curr_bin]
            else:
                self.intensity[i] = min(fst_intensity)

    def region_normalize_maxMz(self, regions, min_mz, max_mz, thresh, precursor_tol = 0,charge=2.0):
        """Max-normalize to number of N evenly-spaced regions, pruning
        peaks which are >= 5% of maximum peak intensity
        """
        values = self.intensity[:]
        renormVals = [0.0]*len(values)
        max_mz = max(self.mz)
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    renormVals[last_el] = 1.0
                    # self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
#                bin_max = sum(values[start:stop])
                # normalize bins to unity
                bin_max = max(values[start:stop])
                # print "start=%d,stop=%d" % (start,stop)
                for index,peak in zip(range(start,stop),values[start:stop]):
                    # print "%d, %f" % (index,self.intensity[index])
                    renormVals[index] = peak/bin_max
                    # self.intensity[index] = peak/bin_max
                    # print "%d: %f" % (index,peak/bin_max)
                # self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak
        precursor_lb = self.precursor_mz - precursor_tol
        precursor_ub = self.precursor_mz + precursor_tol
        for mz, intensity in zip(self.mz, renormVals):
            if intensity > thresh:
                if precursor_tol != 0.0:
                    if mz > precursor_ub or mz < precursor_lb:
                        processed_mz.append(mz)
                        processed_intensity.append(intensity)
                else:
                    processed_mz.append(mz)
                    processed_intensity.append(intensity)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def region_unnormalize(self, regions, min_mz, max_mz, thresh):
        """Prune peaks which are >= 5% of maximum peak intensity over 
        N evenly-spaced regions
        """
        values = [val for val in self.intensity]
        processed_mz = []
        processed_intensity = []
        inc = float(math.ceil((max_mz-min_mz)/float(regions+1)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        # normalize N evenly spaced regions
        for right_edge in bins:
            if stop > last_el:
                if(boundary == 0): # last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #broke at last element, check where last element should belong to
                if(self.mz[stop]<=right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge):
                    boundary = 1
            if(stop != start):
#                bin_max = sum(values[start:stop])
                # normalize bins to unity
                bin_max = max(values[start:stop])
                self.intensity[start:stop] = [math.sqrt(peak/bin_max) for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak
        for mz, intensity, og_intensity in zip(self.mz, self.intensity, values):
            if intensity > thresh:
                processed_mz.append(mz)
                processed_intensity.append(og_intensity)
        
        self.mz = processed_mz
        self.intensity = processed_intensity

    def local_region_normalize(self, regions):
        """Max-normalize to number of N evenly-spaced regions.
        """
        values = self.intensity
        min_mz = min(self.mz)
        max_mz = max(self.mz)
        inc = float(math.ceil((max_mz-min_mz)/float(regions+1)))
        start = 0
        stop = 0
        bins = numpy.arange(min_mz+inc, max_mz+inc, inc)
        last_el = len(values)-1
        boundary = 0

        for right_edge in bins:
            binned = 0
            if stop > last_el: #boundary check
                if(boundary == 0): #last element is in a bin by itself
                    self.intensity[last_el] = 1.0
                break
            while(self.mz[stop] < right_edge): 
#                print ("%f: (%d, %f), (%d,%f)" % (right_edge, start, self.mz[start], stop, self.mz[stop])) #print elements in current bin
                stop = stop + 1
                if stop > last_el: 
                    break
            if(stop == last_el): #check where last element should belong to
                if(self.mz[stop]==right_edge): #last element belongs in current bin
                    boundary = 1
                    stop = stop+1
            elif(stop > last_el):
                if(self.mz[stop-1]<=right_edge): #last element need not be checked again
                    boundary = 1
            if(stop != start):
                bin_max = sum(values[start:stop])
                self.intensity[start:stop] = [peak/bin_max for peak in values[start:stop]]
                start = stop

        # print("%d %d" % (start, stop))
        # for peak in self.intensity: print peak

    def sqrt_normalize(self):
        """Square root all values.
        """
        imax = max(self.intensity)
        values = [ i/imax for i in self.intensity ]
        self.intensity = [ math.sqrt(float(v)) for v in values ]

    def sqrt(self):
        """Square root all values.
        """
        values = self.intensity
        self.intensity = [ math.sqrt(float(v)) for v in values ]

    def rank_normalize(self):
        """Rank normalize the intensities of the spectrum.

        Based on the relative intensity used in Wan et al. (2005)
        implementation of PepHMM, except we invert the ranks (like in
        Riptide). Note that this is a calibrated scale, the intensities
        of all spectra are in [0.0, 1.0], and the maximum intensity after
        normalization is *always* 1.0. The smallest intensity after
        normalization is always 1/length(self.intensity).
        """
        pairs = zip(self.mz, self.intensity)

        # High peaks have high rank. In PepHMM, low peaks have high rank.
        pairs.sort(lambda x, y: cmp(x[1], y[1]))
        assert(pairs[0][1] <= pairs[-1][1])

        pairs = [ (mz, float(idx + 1)/len(pairs))
                  for idx, (mz, _) in enumerate(pairs) ]
        pairs.sort(lambda x, y: cmp(x[0], y[0]))
        self.mz, self.intensity = zip(*pairs)

    def remove_peaks_below_fraction(self, fraction):
        """Remove points lower than fraction of the largest peak intensity.

        Arguments:
           fraction: Intensity percentage such that if a peak's intensity
           is below this value relative to the maximum peak intensity, it
           is deleted

        """
        assert(fraction < 1.0 and fraction >= 0.0)

        mz_keep = []
        intensity_keep = []
        max_intensity = max(self.intensity)
        for mz, intensity in zip(self.mz, self.intensity):
            if intensity/max_intensity > fraction:
                mz_keep.append(mz)
                intensity_keep.append(intensity)

        self.mz = mz_keep
        self.intensity = intensity_keep

    def remove_low_peaks(self, fraction):
        """Remove the lowest intensity points from the spectrum.

        Arguments:
           fraction: The fraction of peaks to remove. The actual fraction
              removed might be slightly higher, since we use ceil.

        """
        assert(fraction < 1.0 and fraction >= 0.0)
        pairs = zip(self.mz, self.intensity)
        pairs.sort(lambda x, y: cmp(x[1], y[1])) # sort by incr intensity.
        assert(pairs[0][1] <= pairs[-1][1])

        m = int(min(math.ceil(fraction * len(pairs)), len(pairs)))
        del pairs[0:m]
        pairs.sort(lambda x, y: cmp(x[0], y[0]))
        self.mz, self.intensity = zip(*pairs)

    def most_intense_200(self):
        """Remove all but the 200 most intense peaks.

        """
        if(len(self.mz) > 200):
            pairs = zip(self.mz, self.intensity)
            pairs.sort(lambda x, y: cmp(x[1], y[1])) # sort by incr intensity.
            del pairs[0:-200]
#            del pairs[0:(200-len(self.mz))] # delete bottom 200 peaks
            pairs.sort(lambda x, y: cmp(x[0], y[0]))
            self.mz, self.intensity = zip(*pairs)

    def most_intense_n(self, n):
        """Remove all but the n most intense peaks.

        """
        # remove all zero peaks
        nonzero_x = []
        nonzero_y = []
        for x,y in zip(self.mz,self.intensity):
            if y > 0.0:
                nonzero_x.append(x)
                nonzero_y.append(y)

        self.mz = nonzero_x
        self.intensity = nonzero_y

        if(len(self.mz) > n):
            pairs = zip(self.mz, self.intensity)
            pairs.sort(lambda x, y: cmp(x[1], y[1])) # sort by incr intensity.
            del pairs[0:-n]
#            del pairs[0:(n-len(self.mz))] # delete bottom n peaks
            pairs.sort(lambda x, y: cmp(x[0], y[0]))
            self.mz, self.intensity = zip(*pairs)

    def most_intense_n_thresh(self, n,thresh,charge):
        """Remove all but the n most intense peaks.

        """
        
        mz = []
        intensity = []
        experimental_mass_cutoff = max(self.mz)*charge+50.0 # filter peaks above this range
        intensity_thresh = max(self.intensity)*thresh # filter peaks below this intensity

        for cmz,ci in zip(self.mz,self.intensity):
            if cmz < experimental_mass_cutoff and ci > intensity_thresh:
                mz.append(cmz)
                intensity.append(ci)

        if(len(mz) > n):
            pairs = zip(mz, intensity)
            pairs.sort(lambda x, y: cmp(x[1], y[1])) # sort by incr intensity.
            del pairs[0:-n]
            pairs.sort(lambda x, y: cmp(x[0], y[0]))
            # set pointers to filtered spectrum
            self.mz, self.intensity = zip(*pairs)
        else:
            self.mz = mz
            self.intensity = intensity

    def filter_lower_mz(self, mz_cutoff = 200):
        """Remove all peaks such that mz < mz_cutoff.

        """
        if(min(self.mz) <= mz_cutoff):
            del_ind = 0
            for ind, mz in enumerate(self.mz):
                if mz > mz_cutoff: # peaks are sorted, we can exit now
                    break
                else:
                    del_ind += 1

            # delete peaks with mz < mz_cutoff
            del self.mz[:del_ind]
            del self.intensity[:del_ind]
            
    def clamp_intensities(self, value):
        self.intensity = [ value ] * len(self.intensity)

    def sort_points(self):
        """Sort the points in order of increasing m/z.

        Especially useful when plotting spectra.
        """
        pairs = zip(self.mz, self.intensity)
        pairs.sort(lambda x, y: cmp(x[0], y[0])) # sort by increasing mz.
        assert(pairs[0][1] <= pairs[-1][1])
        self.mz, self.intensity = zip(*pairs)


def MS2Iterator(filename, has_gzcat = False):

    if os.path.splitext(filename)[1] == '.gz':
        if has_gzcat:
            f = subprocess.Popen(['gzcat', filename], stdout = subprocess.PIPE)
            filebuf = f.stdout.read()
        else:
            f = gzip.open(filename)
            filebuf = f.read()
    else:
        fileobj = open(filename)
        filebuf = mmap.mmap(fileobj.fileno(), os.path.getsize(filename),
                            access = mmap.ACCESS_READ)

    pat = re.compile('(H.*(\r\n|\n))*^S\s+(?P<id>(\d+))\s+(?P=id)\s+(?P<premz>\S+)(\r\n|\n)(^(H|I|D).*(\r\n|\n))*(?P<charges>(Z\s+(\d+)\s+(\S+)(\r\n|\n))+)(^(I|D).*(\r\n|\n))*(?P<payload>((\S+)\s(\S+)(\r\n|\n))+)', re.MULTILINE)
    chg = re.compile('Z\s+(?P<charge>(\d))\s+(?P<mass>(\S+))')
    for m in pat.finditer(filebuf):
        s = MS2Spectrum()
        s.spectrum_id = int(m.group('id'))
        s.precursor_mz = float(m.group('premz'))
        s.charge_lines = list( (int(c.group('charge')), float(c.group('mass')))
                              for c in chg.finditer(m.group('charges')) )
        s.charge_lines.sort() # sort by increasing charge

        mz = [ ]
        intensity = [ ]
        for line in m.group('payload').split('\n'):
            tokens = line.split()
            if tokens:
                mz.append(float(tokens[0]))
                intensity.append(float(tokens[1]))
        s.mz = mz
        s.intensity = intensity
        yield s

def SampleSpectra(spectra, n):

    reservoir = list(itertools.islice(spectra, n))
    if len(reservoir) == n:
        for number, s in enumerate(spectra):
            index = random.randint(0, number)
            if index < n:
                reservoir[index] = s

    return iter(reservoir)

def MS2IteratorSample(filename, n, has_gzcat = False):

    spectra = MS2Iterator(filename, has_gzcat)
    return SampleSpectra(spectra, n)
