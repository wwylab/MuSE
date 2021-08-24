

#include <assert.h>
#include "PBLocalFile.h"
#include <limits.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

size_t PBLocalFile::pbread( void * ptr, size_t size, size_t count ){
    return fread(ptr, size, count, file);
}

size_t PBLocalFile::pbwrite(const void * ptr, size_t size, size_t count ){
    sum += size * count;
    size_t chunk = fwrite(ptr, size, count, file);
    if (chunk != count){
        fprintf( stderr, "Cannot write to file %s. \n", filename);
        exit(EXIT_FAILURE); 
    }
    return chunk;    
}

off_t PBLocalFile::pbseek(long int offset, int whence){
    int success = fseek(file, offset, whence);
    if (success != 0){
        if (ferror(file))
        {
            fprintf( stderr, "fseek failed on %s file. \n", filename);
            exit(EXIT_FAILURE);
        }
    }
    
    return ftell(file);
}

PBLocalFile::PBLocalFile(const char * _filename, const char * mode){
    filename = _filename;
    file = fopen(filename, mode);
    if (file == NULL){
        fprintf( stderr, "%s file cannot be open. \n", filename);
        failed = true;
    }
}

PBLocalFile::~PBLocalFile(){
    if (file != NULL){
        fclose(file);
    }
}


char * PBLocalFile::realPath(){
    return   realpath(filename, NULL);;
}
