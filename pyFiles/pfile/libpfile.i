// SWIG wrapper around _pfile.so
// 
// Written by Ajit Singh <ajit@ee.washington.edu>
//
// The %{ ... %} code are declarations, and are not parsed by SWIG. general.h
// is included because we need some of the typedefs. The stuff after the
// declarations are exposed to Python.

%module libpfile
%{
#include "pfile.h"
#include "general.h"
%}

// Allows you to pass Python file objects to C++ functions that expect FILE*
// automatically.
//
// WARNING: PyFile_AsFile has been removed from Python3, so this typemap will
// fail if you try compiling the wrapper against Python 3. If this is a problem,
// rewrite your Python code to use pfile.fopen instead of fopen.

%typemap(in) FILE* {
    if ( PyFile_Check($input) ){
        $1 = PyFile_AsFile($input);
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name must be a file type.");
        return NULL;
    }
}


%include "pfile.h"
typedef unsigned int UInt32;
typedef int Int32;

%include "carrays.i"
%array_functions(float, doubleArray);
%array_functions(UInt32, uintArray);

FILE* fopen(const char* filename, const char* mode);
int fclose(FILE* stream);


