/* This file is automatically generated. DO NOT EDIT! */

#ifndef _sf_file_h
#define _sf_file_h


#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>


#include <stdio.h>


#include "_bool.h"
#include "c99.h"


#define SF_MAX_DIM 9


typedef struct sf_File *sf_file;


typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX, SF_SHORT, SF_DOUBLE, SF_LONG} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;


void sf_file_error(bool err);
/*< set error on opening files >*/


sf_file sf_input (/*@null@*/ const char* tag);
/*< Create an input file structure >*/


sf_file sf_output (/*@null@*/ const char* tag);
/*< Create an output file structure.
  ---
  Should do output after the first call to sf_input. >*/


sf_datatype sf_gettype (sf_file file);
/*< return file type >*/


sf_dataform sf_getform (sf_file file);
/*< return file form >*/


size_t sf_esize(sf_file file);
/*< return element size >*/


void sf_settype (sf_file file, sf_datatype type);
/*< set file type >*/


void sf_setpars (sf_file file);
/*< change parameters to those from the command line >*/


void sf_expandpars (sf_file file);
/*< add parameters from the command line >*/


size_t sf_bufsiz(sf_file file);
/*< return buffer size for efficient I/O >*/


void sf_setform (sf_file file, sf_dataform form);
/*< set file form >*/


void sf_setformat (sf_file file, const char* format);
/*< Set file format.
  ---
  format has a form "form_type", i.e. native_float, ascii_int, etc.
  >*/


void sf_fileclose (sf_file file);
/*< close a file and free allocated space >*/


bool sf_histint (sf_file file, const char* key,/*@out@*/ int* par);
/*< read an int parameter from file >*/


bool sf_histints (sf_file file, const char* key,/*@out@*/ int* par,size_t n);
/*< read an int array of size n parameter from file >*/


bool sf_histlargeint (sf_file file, const char* key,/*@out@*/ off_t* par);
/*< read a sf_largeint parameter from file >*/


bool sf_histfloat (sf_file file, const char* key,/*@out@*/ float* par);
/*< read a float parameter from file >*/


bool sf_histdouble (sf_file file, const char* key,/*@out@*/ double* par);
/*< read a float parameter from file >*/


bool sf_histfloats (sf_file file, const char* key,
		    /*@out@*/ float* par,size_t n);
/*< read a float array of size n parameter from file >*/


bool sf_histbool (sf_file file, const char* key,/*@out@*/ bool* par);
/*< read a bool parameter from file >*/


bool sf_histbools (sf_file file, const char* key,
		   /*@out@*/ bool* par, size_t n);
/*< read a bool array of size n parameter from file >*/


char* sf_histstring (sf_file file, const char* key);
/*< read a string parameter from file (returns NULL on failure) >*/


void sf_fileflush (sf_file file, sf_file src);
/*< outputs parameter to a file (initially from source src)
  ---
  Prepares file for writing binary data >*/


void sf_putint (sf_file file, const char* key, int par);
/*< put an int parameter to a file >*/


void sf_putints (sf_file file, const char* key, const int* par, size_t n);
/*< put an int array of size n parameter to a file >*/


void sf_putlargeint (sf_file file, const char* key, off_t par);
/*< put a sf_largeint parameter to a file >*/


void sf_putfloat (sf_file file, const char* key,float par);
/*< put a float parameter to a file >*/


void sf_putfloats (sf_file file, const char* key, const float* par, size_t n);
/*< put a float array of size n parameter to a file >*/


void sf_putstring (sf_file file, const char* key,const char* par);
/*< put a string parameter to a file >*/


void sf_putline (sf_file file, const char* line);
/*< put a string line to a file >*/


void sf_setaformat (const char* format /* number format (.i.e "%5g") */, 
		    int line /* numbers in line */ );
/*< Set format for ascii output >*/


void sf_complexwrite (sf_complex* arr, size_t size, sf_file file);
/*< write a complex array arr[size] to file >*/


void sf_complexread (/*@out@*/ sf_complex* arr, size_t size, sf_file file);
/*< read a complex array arr[size] from file >*/


void sf_charwrite (char* arr, size_t size, sf_file file);
/*< write a char array arr[size] to file >*/


void sf_ucharwrite (unsigned char* arr, size_t size, sf_file file);
/*< write an unsigned char array arr[size] to file >*/


void sf_charread (/*@out@*/ char* arr, size_t size, sf_file file);
/*< read a char array arr[size] from file >*/


int sf_try_charread(const char* test, sf_file file);
/*< check if you can read test word >*/


int sf_try_charread2 (/*@out@*/ char* arr, size_t size, sf_file file);
/*< try to read size bytes.  return number bytes read >*/


void sf_ucharread (/*@out@*/ unsigned char* arr, size_t size, sf_file file);
/*< read a uchar array arr[size] from file >*/


void sf_intwrite (int* arr, size_t size, sf_file file);
/*< write an int array arr[size] to file >*/


void sf_intread (/*@out@*/ int* arr, size_t size, sf_file file);
/*< read an int array arr[size] from file >*/


void sf_shortread (/*@out@*/ short* arr, size_t size, sf_file file);
/*< read a short array arr[size] from file >*/


void sf_longread (/*@out@*/ off_t* arr, size_t size, sf_file file);
/*< read a long array arr[size] from file >*/


void sf_shortwrite (short* arr, size_t size, sf_file file);
/*< write a short array arr[size] to file >*/


void sf_floatwrite (float* arr, size_t size, sf_file file);
/*< write a float array arr[size] to file >*/


void sf_floatread (/*@out@*/ float* arr, size_t size, sf_file file);
/*< read a float array arr[size] from file >*/


off_t sf_bytes (sf_file file);
/*< Count the file data size (in bytes) >*/


off_t sf_tell (sf_file file);
/*< Find position in file >*/


FILE *sf_tempfile(char** dataname, const char* mode);
/*< Create a temporary file with a unique name >*/


void sf_seek (sf_file file, off_t offset, int whence);
/*< Seek to a position in file. Follows fseek convention. >*/


FILE* sf_filestream (sf_file file);
/*< Returns file descriptor to a stream >*/


void sf_unpipe (sf_file file, off_t size);
/*< Redirect a pipe input to a direct access file >*/


void sf_close(void);
/*< Remove temporary files >*/

#endif
