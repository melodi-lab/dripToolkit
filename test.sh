#!/bin/bash
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

# adjust the number of threads to be used for all tests in the variable below
NUMTHREADS=4

#################################################################
################## Low-resolution MS2 searches
#################################################################

############################################
######### train and search
############################################
function trainTest {
    echo "Training DRIP"
    python -OO dripTrain.py \
    	--psm-library data/riptideTrainingData/strict-orbitrap.psm \
    	--spectra data/riptideTrainingData/strict-orbitrap.ms2 \
    	--output-mean-file dripLearned.means \
    	--output-covar-file dripLearned.covars \
    	--mods-spec 'C+57.0214'

    echo "Digesting protein database"
    python -OO dripDigest.py \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    echo "Searching spectra"
    time python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
	--num-threads $NUMTHREADS \
    	--top-match 1 \
	--charges 2 \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output
}

############################################
######### test using beam pruning
############################################
function testBeam {
    BEAM=75
    echo "Digesting protein database"
    python -OO dripDigest.py \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    time python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
	--num-threads $NUMTHREADS \
    	--top-match 1 \
    	--beam $BEAM \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output-beam$BEAM
}

###################################################
######### Extract features for output of trainTest,
######### write output to PIN file format
###################################################
function dripExtractLowRes {
    python -OO dripExtract.py \
	--write-pin true \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
	--psm-file dripSearch-test-output.txt \
	--num-threads $NUMTHREADS \
	--mods-spec 'C+57.0214' \
	--spectra data/test.ms2 \
	--output dripExtract-test-output.txt
}

############################################
######### train and search ch3 PSMs
############################################
function trainTestCh3 {
    echo "Training DRIP"
    python -OO dripTrain.py \
    	--psm-library data/riptideTrainingData/strict-orbitrap-ch3.psm \
    	--spectra data/riptideTrainingData/strict-orbitrap-ch3.ms2 \
    	--output-mean-file dripLearned-ch3.means \
    	--output-covar-file dripLearned-ch3.covars \
    	--mods-spec 'C+57.0214'

    echo "Digesting protein database"
    python -OO dripDigest.py \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    echo "Searching spectra"
    time python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned-ch3.means \
    	--learned-covars dripLearned-ch3.covars \
    	--charges 3 \
	--num-threads $NUMTHREADS \
    	--top-match 1 \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output
}

#########################################
######### train, search, and find max PSM
######### over different charge states
#########################################
function trainTestRecalibrate {
    echo "Training DRIP"
    python -OO dripTrain.py \
    	--psm-library data/riptideTrainingData/strict-orbitrap.psm \
    	--spectra data/riptideTrainingData/strict-orbitrap.ms2 \
    	--output-mean-file dripLearned.means \
    	--output-covar-file dripLearned.covars \
    	--mods-spec 'C+57.0214'

    echo "Digesting protein database"
    python -OO dripDigest.py \
    	--recalibrate True \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    echo "Searching spectra"
    python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
	--num-threads $NUMTHREADS \
    	--top-match 1 \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output
}

############################################
######### search high-res MS2 spectra
############################################
function dripSearchHighres {
    # digest directory
    ./dripDigest.py  \
    	--fasta /s0/halloj3/samplePlasmodium/sample-malaria/plasmo_Pfalciparum3D7_NCBI.fasta \
    	--min-length 7 \
    	--custom-enzyme '[K]|[X]' \
    	--mods-spec 'C+57.0214,K+229.16293' \
    	--nterm-peptide-mods-spec 'X+229.16293' \
    	--monoisotopic-precursor true \
    	--recalibrate True \
    	--decoys True

    python -OO dripSearch.py \
	--digest-dir 'dripDigest-output' \
	--precursor-window 50 \
	--num-threads $NUMTHREADS \
	--high-res-ms2 true \
	--precursor-window-type 'ppm' \
	--precursor-filter 'True' \
	--spectra data/malariaTest.ms2 \
	--output dripSearch-malariaTest-output
}

############################################
######### search high-res MS2 spectra with 
######### variable mods
############################################
function dripSearchHighresVarMods {
    # digest directory
    ./dripDigest.py  \
    	--fasta /s0/halloj3/samplePlasmodium/sample-malaria/plasmo_Pfalciparum3D7_NCBI.fasta \
    	--min-length 7 \
    	--custom-enzyme '[K]|[X]' \
    	--mods-spec '3M+15.9949,C+57.0214,K+229.16293' \
    	--nterm-peptide-mods-spec 'X+229.16293' \
    	--monoisotopic-precursor true \
    	--recalibrate True \
    	--decoys True

    python -OO dripSearch.py \
	--digest-dir 'dripDigest-output' \
	--precursor-window 50 \
	--num-threads $NUMTHREADS \
	--high-res-ms2 true \
	--precursor-window-type 'ppm' \
	--precursor-filter 'True' \
	--spectra data/malariaTest.ms2 \
	--output dripSearch-malariaTestVarmods-output
}

############################################
######### extract features for high-res MS2,
######### output of dripSearchHighresVarMods
############################################
function dripExtractHighResVarMods {
    python -OO dripExtract.py \
	--append-to-pin false \
	--high-res-ms2 true \
	--precursor-filter 'True' \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
	--psm-file dripSearch-malariaTest-output.txt \
	--num-threads $NUMTHREADS \
    	--mods-spec '3M+15.9949,C+57.0214,K+229.16293' \
    	--nterm-peptide-mods-spec 'X+229.16293' \
	--spectra data/malariaTest.ms2 \
	--output dripExtract-malariaTestVarmods-output.txt
}

#############################################
######### Split data for cluster use, run 
######### individual jobs to simulate cluster
######### use, collect and merge results
#############################################
function clusterTest {
    if [ ! $DRIPTOOLKIT ]
    then
	echo "Please set DRIPTOOLKIT environment variable a directory containing the built toolkit."
    else

	echo "Digesting protein database"
	python -OO dripDigest.py \
    	    --min-length 6 \
    	    --fasta data/yeast.fasta \
    	    --enzyme 'trypsin/p' \
    	    --monoisotopic-precursor true \
    	    --missed-cleavages 0 \
    	    --digestion 'full-digest'

	echo "Creating cluster jobs"
	python -OO dripSearch.py \
    	    --digest-dir 'dripDigest-output' \
    	    --precursor-window 3.0 \
    	    --learned-means dripLearned.means \
    	    --learned-covars dripLearned.covars \
	    --num-threads $NUMTHREADS \
	    --num-cluster-jobs 4 \
    	    --top-match 1 \
    	    --spectra data/test.ms2 \
    	    --output dripSearch-clusterTest-output

	echo "Running created jobs (locally)"
	for j in encode/*.sh
	do
    	    echo $j
    	    ./$j
	done

	echo "Merging cluster results"
	python -OO dripSearch.py --merge-cluster-results True \
    	    --logDir log \
    	    --output dripSearch-clusterTest-output
    fi
}

# available examples (see function above for description): 
# trainTest
# testBeam
# dripExtractLowRes (run trainTest first)
# trainTestCh3
# trainTestRecalibrate
# dripSearchHighres
# dripSearchHighresVarMods
# dripExtractHighResVarMods (run dripSearchHighresVarMods first)
# clusterTest

# traing and test low-res MS2
trainTest

# run several tests
runTests=( dripSearchHighres dripSearchHighresVarMods \
    dripExtractHighResVarMods )

# loop through array of tests
for dripTest in ${runTests[@]}
do
    echo $dripTest
    $dripTest
done