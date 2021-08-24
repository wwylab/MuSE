

#ifndef __PBLOCALFILE_H__
#define __PBLOCALFILE_H__
#include <stdio.h>
class PBLocalFile
{
    private:
        FILE* file;
        const char * filename;
        size_t sum = 0;
        bool failed = false;
    public:
        bool isFailed() { return failed;}
        PBLocalFile(const char * _filename, const char * mode );
        ~PBLocalFile();
        size_t pbread( void * ptr, size_t size, size_t count );
        size_t pbwrite(const void * ptr, size_t size, size_t count );
        off_t pbseek(long int offset, int whence);
        char * realPath();

  
};


#endif
