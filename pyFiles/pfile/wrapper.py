#!/usr/bin/env python
#
# Written by Ajit Singh <ajit@ee.washington.edu>

import os
import sys
import types
import libpfile as lib

class PFile(object):

    def __init__(self, nf, ni, f, doswap = None):
        """Constructor

        Creates a PFile with the name specified in fn. Each segment contains
        nf floats, followed by ni integers.

        Arguments:
            nf: Number of floats in a segment.
            ni: Number of integers (labels) in a segment.
            fn: Name of the pfile to create, no extension is forced.
            doswap: If set, it will force a particular byte order in the
                generated PFile. Useful if you're writing a PFile on one
                platform for use on another platform. If None, use whatever
                sys.byteorder returns.

        """

        if not doswap:
            if sys.byteorder == 'little':
                doswap = 1
            elif sys.byteorder == 'big':
                doswap = 0
            else:
                raise Exception("Could not infer byteorder.")

        index = 1
        if type(f) == types.FileType:
            self.f = f
        elif type(f) == types.StringType:
            self.f = open(f, 'w')
        else:
            raise Exception("Bad filename argument: %s" % str(f))
        self.nf = nf
        self.ni = ni
        self.doswap = doswap
        self.pf = lib.OutFtrLabStream_PFile(0, '', self.f, nf, ni, index,
                                            doswap)

        # Create buffers for translating Python lists of floats or ints to
        # float* and unsigneed int*
        self.buf_floats = lib.new_doubleArray(self.nf)
        self.buf_ints = lib.new_uintArray(self.ni)

    def __del__(self):
        """Destructor.

        TODO(ajit): Calling pfile.fclose of self.pf causes a segmentation fault.
        Determine where the file is really being deleted (it may only be on
        exit, or deletion of the class).

        """
        del self.pf
        lib.delete_doubleArray(self.buf_floats)
        lib.delete_uintArray(self.buf_ints)

    @property
    def name(self):
        return self.f.name

    def check_frame(self, *args):
        if len(args) != self.nf + self.ni:
            raise Exception("Wrong length %d vs. %d" % (len(args),
                                                        self.nf + self.ni))
        for i in xrange(0, self.nf, 1):
            if not type(args[i]) == types.FloatType:
                raise Exception("Wrong type arg[%d]: wanted float, got %s" % (
                    i, str(type(args[i]))))

        for i in xrange(self.nf, self.nf+self.ni, 1):
            if not type(args[i]) == types.IntType:
                raise Exception("Wrong type arg[%d]: wanted int, got %s" % (
                    i, str(type(args[i]))))

    def add_frame(self, *args):
        for i in xrange(0, self.nf, 1):
            lib.doubleArray_setitem(self.buf_floats, i, args[i])
        for i in xrange(self.nf, self.nf + self.ni, 1):
            lib.uintArray_setitem(self.buf_ints, i-self.nf, args[i])
        self.pf.write_ftrslabs(1, self.buf_floats, self.buf_ints)

    def add_segment(self, nframes, floats, ints):
        """Copy a whole sentence in one shot.

        Can be useful in reducing the Python -> C++ overhead required to
        generate one sentence: e.g., creating one list for all the floats
        in a sentence, instead of one list per-frame. You do not need to
        call end_segment after using this function.

        TODO(ajit): It's not clear whether the segment ID is actually used
        anywhere. The code in pfile.cc:doneseg does not appear to use the
        segment ID, and ignoring it doesn't seem to cause any problems.

        Arguments:
            nframes: Number of frames in the sentence.
            floats: Iterable with all of the floats in the sentence. First,
                all the floats in frame 0, then frame 1, etc.
            ints: Iterable with all the integers in the sentence.
        """
        pass

    def end_segment(self, i = None):
        if not i:
            i = lib.SEGID_UNKNOWN
        self.pf.doneseg(i)
