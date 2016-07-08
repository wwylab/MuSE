/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#define BGZF_BLOCK_SIZE 0x10000 // 64k

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8

typedef struct {
	char open_mode;  // 'r' or 'w'
	int16_t compress_level;
	int16_t errcode;
	int cache_size;
	int block_length;
	int block_offset;
	int64_t block_address;
	void *uncompressed_block;
	void *compressed_block;
	void *cache; // a pointer to a hash table
//	void *fp; // "fp" is used in tabix-0.2.6 bgzf.h; actual file handler; FILE* on writing; FILE* or knetFile* on reading
	FILE* file; // "file" is used in samtools-0.1.18 bgzf.h; actual file handler; FILE* on writing; FILE* or knetFile* on reading
	
	int16_t owned_file;
	int file_descriptor;
	int uncompressed_block_size;
	int compressed_block_size;
	const char* error;
} BGZF;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

	/******************
	 * Basic routines *
	 ******************/

	/**
	 * Open an existing file descriptor for reading or writing.
	 *
	 * @param fd    file descriptor
	 * @param mode  mode matching /[rwu0-9]+/: 'r' for reading, 'w' for writing and a digit specifies
	 *              the zlib compression level; if both 'r' and 'w' are present, 'w' is ignored.
	 * @return      BGZF file handler; 0 on error
	 */
	BGZF* bgzf_dopen(int fd, const char *mode);

	/*
	 * Open an existing file descriptor for reading or writing.
	 * Mode must be either "r" or "w".
	 * A subsequent bgzf_close will not close the file descriptor.
	 * Returns null on error.
	 */
	BGZF* bgzf_fdopen(int fd, const char* __restrict mode);

	/*
	 * Open the specified file for reading or writing.
	 * Mode must be either "r" or "w".
	 * Returns null on error.
	 */
	BGZF* bgzf_open(const char* path, const char* __restrict mode);

	/*
	 * Close the BGZ file and free all associated resources.
	 * Does not close the underlying file descriptor if created with bgzf_fdopen.
	 * Returns zero on success, -1 on error.
	 */
	int bgzf_close(BGZF* fp);

	/*
	 * Read up to length bytes from the file storing into data.
	 * Returns the number of bytes actually read.
	 * Returns zero on end of file.
	 * Returns -1 on error.
	 */
	int bgzf_read(BGZF* fp, void* data, int length);

	/*
	 * Write length bytes from data to the file.
	 * Returns the number of bytes written.
	 * Returns -1 on error.
	 */
	int bgzf_write(BGZF* fp, const void* data, int length);

	/*
	 * Return a virtual file pointer to the current location in the file.
	 * No interpetation of the value should be made, other than a subsequent
	 * call to bgzf_seek can be used to position the file at the same point.
	 * Return value is non-negative on success.
	 * Returns -1 on error.
	 */
	#define bgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))

	/*
	 * Set the file to read from the location specified by pos, which must
	 * be a value previously returned by bgzf_tell for this file (but not
	 * necessarily one returned by this file handle).
	 * The where argument must be SEEK_SET.
	 * Seeking on a file opened for write is not supported.
	 * Returns zero on success, -1 on error.
	 */
	int64_t bgzf_seek(BGZF* fp, int64_t pos, int where);

	/**
	 * Check if the BGZF end-of-file (EOF) marker is present
	 *
	 * @param fp    BGZF file handler opened for reading
	 * @return      1 if EOF is present; 0 if not or on I/O error
	 */
	int bgzf_check_EOF(BGZF *fp);

	/**
	 * Write the data in the buffer to the file.
	 */
	int bgzf_flush(BGZF* fp);

	int bgzf_check_bgzf(const char *fn);

	/**
	 * Check if a file is in the BGZF format
	 *
	 * @param fn    file name
	 * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
	 */
	int bgzf_is_bgzf(const char *fn);

	/*********************
	 * Advanced routines *
	 *********************/

	/*
	 * Set the cache size. Zero to disable. By default, caching is
	 * disabled. The recommended cache size for frequent random access is
	 * about 8M bytes.
	 *
	 * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
	 *
	 * @param fp    BGZF file handler
	 * @param size  size of cache in bytes; 0 to disable caching (default)
	 */
	void bgzf_set_cache_size(BGZF *fp, int cache_size);

	/**
	 * Flush the file if the remaining buffer size is smaller than _size_
	 */
	int bgzf_flush_try(BGZF *fp, int size);

	/**
	 * Read one byte from a BGZF file. It is faster than bgzf_read()
	 * @param fp     BGZF file handler
	 * @return       byte read; -1 on end-of-file or error
	 */
	int bgzf_getc(BGZF *fp);

	/**
	 * Read one line from a BGZF file. It is faster than bgzf_getc()
	 *
	 * @param fp     BGZF file handler
	 * @param delim  delimitor
	 * @param str    string to write to; must be initialized
	 * @return       length of the string; 0 on end-of-file; negative on error
	 */
	int bgzf_getline(BGZF *fp, int delim, kstring_t *str);

	/**
	 * Read the next BGZF block.
	 */
	int bgzf_read_block(BGZF *fp);

#ifdef __cplusplus
}
#endif

#endif
