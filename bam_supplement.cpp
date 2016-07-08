//
//  bam_supplement.c
//  SamTest
//
//  Created by Daniel on 1/28/13.
//  Copyright (c) 2013 Daniel. All rights reserved.
//

#include <ctype.h>
#include <assert.h>
#include "bam.h"
#include "khash.h"
#include "ksort.h"
#include "kprobaln.h"

/*------------------------------------------------------------------------------------------------*/
// bam_aux.c

KHASH_MAP_INIT_STR(s, int)

void bam_destroy_header_hash(bam_header_t *header) {
	if(header->hash)
		kh_destroy(s, (khash_t(s)*)header->hash);
}

void bam_init_header_hash(bam_header_t *header) {
	if(header->hash == 0) {
		int ret, i;
		khiter_t iter;
		khash_t(s) *h;
		header->hash = h = kh_init(s);
		for(i = 0; i < header->n_targets; ++i) {
			iter = kh_put(s, h, header->target_name[i], &ret);
			kh_value(h, iter) = i;
		}
	}
}

int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *beg, int *end) {
	char *s;
	int i, l, k, name_end;
	khiter_t iter;
	khash_t(s) *h;
	
	bam_init_header_hash(header);
	h = (khash_t(s)*)header->hash;
	
	*ref_id = *beg = *end = -1;
	name_end = l = strlen(str);
	s = (char*)malloc(l+1);
	// remove space
	for(i = k = 0; i < l; ++i) {
		if(!isspace(str[i]))
			s[k++] = str[i];
	}
	s[k] = 0;
	l = k;
	// determine the sequence name
	for(i = l - 1; i >= 0; --i) {
		if(s[i] == ':')
			break; // look for colon from the end
	}
	if(i >= 0)
		name_end = i;
	if(name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for(i = name_end + 1; i < l; ++i) {
			if(s[i] == '-')
				++n_hyphen;
			else if(!isdigit(s[i]) && s[i] != ',')
				break;
		}
		if(i < l || n_hyphen > 1)
			name_end = l; // malformated region string; then take str as the name
		s[name_end] = 0;
		iter = kh_get(s, h, s);
		if(iter == kh_end(h)) { // cannot find the sequence name
			iter = kh_get(s, h, str); // try str as the name
			if(iter == kh_end(h)) {
				if(bam_verbose >= 2)
					fprintf(stderr, "[%s] fail to determine the sequence name.\n", __func__);
				free(s);
				return -1;
			}
			else
				s[name_end] = ':', name_end = l;
		}
	}
	else
		iter = kh_get(s, h, str);
	*ref_id = kh_val(h, iter);
	// parse the interval
	if(name_end < l) {
		for(i = k = name_end + 1; i < l; ++i) {
			if(s[i] != ',')
				s[k++] = s[i];
		}
		s[k] = 0;
		*beg = atoi(s + name_end + 1);
		for(i = name_end + 1; i != k; ++i) {
			if(s[i] == '-')
				break;
		}
		*end = i < k ? atoi(s + i + 1) : 1<<29;
		if(*beg > 0)
			--*beg;
	}
	else
		*beg = 0, *end = 1<<29;
	free(s);
	return *beg <= *end ? 0 : -1;
}

#define __skip_tag(s) do { \
int type = toupper(*(s)); \
++(s); \
if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
else (s) += bam_aux_type2size(type); \
} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]) {
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while(s < b->data + b->data_len) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if(x == y)
			return s;
		__skip_tag(s);
	}
	return 0;
}

int bam_aux_del(bam1_t *b, uint8_t *s) {
	// s MUST BE returned by bam_aux_get()
	uint8_t *p, *aux;
	aux = bam1_aux(b);
	p   = s - 2;
	__skip_tag(s);
	memmove(p, s, b->l_aux - (s - aux));
	b->data_len -= s - p;
	b->l_aux    -= s - p;
	return 0;
}

void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data) {
	int ori_len  = b->data_len;
	b->data_len += 3 + len;
	b->l_aux    += 3 + len;
	if(b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	b->data[ori_len]     = tag[0];
	b->data[ori_len + 1] = tag[1];
	b->data[ori_len + 2] = type;
	memcpy(b->data + ori_len + 3, data, len);
}

int32_t bam_aux2i(const uint8_t *s) {
	int type;
	if(s == 0) return 0;
	type = *s++;
	if(type == 'c') return (int32_t)*(int8_t*)s;
	else if(type == 'C') return (int32_t)*(uint8_t*)s;
	else if(type == 's') return (int32_t)*(int16_t*)s;
	else if(type == 'S') return (int32_t)*(uint16_t*)s;
	else if(type == 'i' || type == 'I') return *(int32_t*)s;
	else return 0;
}

float bam_aux2f(const uint8_t *s) {
	int type;
	type = *s++;
	if(s == 0) return 0.0;
	if(type == 'f') return *(float*)s;
	else return 0.0;
}

double bam_aux2d(const uint8_t *s) {
	int type;
	type = *s++;
	if(s == 0) return 0.0;
	if(type == 'd') return *(double*)s;
	else return 0.0;
}

char bam_aux2A(const uint8_t *s) {
	int type;
	type = *s++;
	if(s == 0) return 0;
	if(type == 'A') return *(char*)s;
	else return 0;
}

char *bam_aux2Z(const uint8_t *s) {
	int type;
	type = *s++;
	if(s == 0) return 0;
	if(type == 'Z' || type == 'H') return (char*)s;
	else return 0;
}

/*------------------------------------------------------------------------------------------------*/
// bam_index.c

#define BAM_LIDX_SHIFT 14    // 1<<14 is the size of minimum bin.
#define BAM_MAX_BIN    37450 // =(8^6-1)/7+1

//typedef struct {
//	uint64_t u, v;
//} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(off, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;

KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct __bam_index_t {
	int32_t n;
	uint64_t n_no_coor; // unmapped reads without coordinate
	khash_t(i) **index;
	bam_lidx_t *index2;
};

void bam_index_destroy(bam_index_t *idx) {
	khint_t k;
	int i;
	if(idx == 0)
		return;
	for(i = 0; i < idx->n; ++i) {
		khash_t(i) *index  = idx->index[i];
		bam_lidx_t *index2 = idx->index2 + i;
		for(k = kh_begin(index); k != kh_end(index); ++k) {
			if(kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}
	free(idx->index);
	free(idx->index2);
	free(idx);
}

static bam_index_t *bam_index_load_core(FILE *fp) {
	int i;
	char magic[4];
	bam_index_t *idx;
	if(fp == 0) {
		fprintf(stderr, "[bam_index_load_core] fail to load index.\n");
		return 0;
	}
	fread(magic, 1, 4, fp);
	if(strncmp(magic, "BAI\1", 4)) {
		fprintf(stderr, "[bam_index_load] wrong magic number.\n");
		fclose(fp);
		return 0;
	}
	idx = (bam_index_t*)calloc(1, sizeof(bam_index_t));
	fread(&idx->n, 4, 1, fp);
	if(bam_is_be)
		bam_swap_endian_4p(&idx->n);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (bam_lidx_t*)calloc(idx->n, sizeof(bam_lidx_t));
	for(i = 0; i < idx->n; ++i) {
		khash_t(i) *index;
		bam_lidx_t *index2 = idx->index2 + i;
		uint32_t key, size;
		khint_t k;
		int j, ret;
		bam_binlist_t *p;
		index = idx->index[i] = kh_init(i);
		// load binning index
		fread(&size, 4, 1, fp);
		if(bam_is_be)
			bam_swap_endian_4p(&size);
		for(j = 0; j < (int)size; ++j) {
			fread(&key, 4, 1, fp);
			if(bam_is_be)
				bam_swap_endian_4p(&key);
			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			fread(&p->n, 4, 1, fp);
			if(bam_is_be)
				bam_swap_endian_4p(&p->n);
			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m * 16);
			fread(p->list, 16, p->n, fp);
			if(bam_is_be) {
				int x;
				for(x = 0; x < p->n; ++x) {
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}
		// load linear index
		fread(&index2->n, 4, 1, fp);
		if(bam_is_be)
			bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		fread(index2->offset, index2->n, 8, fp);
		if(bam_is_be) {
			for(j = 0; j < index2->n; ++j)
				bam_swap_endian_8p(&index2->offset[j]);
		}
	}
	if(fread(&idx->n_no_coor, 8, 1, fp) == 0)
		idx->n_no_coor = 0;
	if(bam_is_be)
		bam_swap_endian_8p(&idx->n_no_coor);
	return idx;
}

bam_index_t *bam_index_load_local(const char *_fn) {
	FILE *fp;
	char *fnidx, *fn;
	
	fn = strdup(_fn);
	fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcpy(fnidx, fn);
	strcat(fnidx, ".bai");
	fp = fopen(fnidx, "rb");
	if(fp == 0) { // try "{base}.bai"
		char *s = strstr(fn, "bam");
		if(s == fn + strlen(fn) - 3) {
			strcpy(fnidx, fn);
			fnidx[strlen(fn)-1] = 'i';
			fp = fopen(fnidx, "rb");
		}
	}
	free(fnidx);
	free(fn);
	if(fp) {
		bam_index_t *idx = bam_index_load_core(fp);
		fclose(fp);
		return idx;
	}
	else
		return 0;
}

bam_index_t *bam_index_load(const char *fn) {
	bam_index_t *idx;
	idx = bam_index_load_local(fn);
	if(idx == 0)
		fprintf(stderr, "[bam_index_load] fail to load BAM index.\n");
	return idx;
}

//struct __bam_iter_t {
//	int from_first; // read from the first record; no random access
//	int tid, beg, end, n_off, i, finished;
//	uint64_t curr_off;
//	pair64_t *off;
//};

void bam_iter_destroy(bam_iter_t iter) {
	if(iter) {
		free(iter->off);
		free(iter);
	}
}

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[BAM_MAX_BIN]) {
	int i = 0, k;
	if(beg >= end)
		return 0;
	if(end >= 1u<<29)
		end = 1u<<29;
	--end;
	list[i++] = 0;
	for(k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for(k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for(k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for(k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for(k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}

bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end) {
	uint16_t *bins;
	int i, n_bins, n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;
	bam_iter_t iter = 0;
	
	if(beg < 0)
		beg = 0;
	if(end < beg)
		return 0;
	// initialize iter
	iter = (bam_iter_t)calloc(1, sizeof(struct __bam_iter_t));
	iter->tid = tid, iter->beg = beg, iter->end = end;
	iter->i = -1;
	//
	bins   = (uint16_t*)calloc(BAM_MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index  = idx->index[tid];
	if(idx->index2[tid].n > 0) {
		min_off = (beg>>BAM_LIDX_SHIFT >= idx->index2[tid].n) ? idx->index2[tid].offset[idx->index2[tid].n-1] : idx->index2[tid].offset[beg>>BAM_LIDX_SHIFT];
		if(min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
			int n = beg>>BAM_LIDX_SHIFT;
			if(n > idx->index2[tid].n)
				n = idx->index2[tid].n;
			for(i = n - 1; i >= 0; --i) {
				if(idx->index2[tid].offset[i] != 0)
					break;
			}
			if(i >= 0)
				min_off = idx->index2[tid].offset[i];
		}
	}
	else
		min_off = 0; // tabix 0.1.2 may produce such index files
	for(i = n_off = 0; i < n_bins; ++i) {
		if((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}
	if(n_off == 0) {
		free(bins);
		return iter;
	}
	off = (pair64_t*)calloc(n_off, 16);
	for(i = n_off = 0; i < n_bins; ++i) {
		if((k = kh_get(i, index, bins[i])) != kh_end(index)) {
			int j;
			bam_binlist_t *p = &kh_value(index, k);
			for(j = 0; j < p->n; ++j) {
				if(p->list[j].v > min_off)
					off[n_off++] = p->list[j];
			}
		}
	}
	free(bins);
	if(n_off == 0) {
		free(off);
		return iter;
	}
	{
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		int l;
		ks_introsort(off, n_off, off);
		// resolve completely contained adjacent blocks
		for(i = 1, l = 0; i < n_off; ++i) {
			if(off[l].v < off[i].v)
				off[++l] = off[i];
		}
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for(i = 1; i < n_off; ++i) {
			if(off[i-1].v >= off[i].u)
				off[i-1].v = off[i].u;
		}
		{
			// merge adjacent blocks
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
			for(i = 1, l = 0; i < n_off; ++i) {
#ifdef BAM_TRUE_OFFSET
				if(off[l].v + BAM_MIN_CHUNK_GAP > off[i].u)
					off[l].v = off[i].v;
#else
				if(off[l].v>>16 == off[i].u>>16)
					off[l].v = off[i].v;
#endif
				else
					off[++l] = off[i];
			}
			n_off = l + 1;
#endif
		}
		bam_destroy1(b);
	}
	iter->n_off = n_off;
	iter->off   = off;
	return iter;
}

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b) {
	uint32_t rbeg = b->core.pos;
	uint32_t rend = b->core.n_cigar ? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
	return (rend > beg && rbeg < end);
}

int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b) {
	int ret;
	if(iter && iter->finished)
		return -1;
	if(iter == 0 || iter->from_first) {
		ret = bam_read1(fp, b);
		if(ret < 0 && iter)
			iter->finished = 1;
		return ret;
	}
	if(iter->off == 0)
		return -1;
	for(;;) {
		if(iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if(iter->i == iter->n_off - 1) { // no more chunks
				ret = -1;
				break;
			}
			if(iter->i >= 0) // otherwise bug
				assert(iter->curr_off == iter->off[iter->i].v);
			if(iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
				bam_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bam_tell(fp);
			}
			++iter->i;
		}
		if((ret = bam_read1(fp, b)) >= 0) {
			iter->curr_off = bam_tell(fp);
			if(b->core.tid != iter->tid || b->core.pos >= iter->end) { // no need to proceed
				ret = bam_validate1(NULL, b) ? -1 : -5; // determine whether end of region or error
				break;
			}
			else if (is_overlap(iter->beg, iter->end, b))
				return ret;
		}
		else
			break; // end of file or error
	}
	iter->finished = 1;
	return ret;
}

/*------------------------------------------------------------------------------------------------*/
// bam_pileup.c

typedef struct {
	int k, x, y, end;
} cstate_t;

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg, end;
	cstate_t s;
	struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
	int cnt, n, max;
	lbnode_t **buf;
} mempool_t;

struct __bam_plp_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	int32_t tid, pos, max_tid, max_pos;
	int is_eof, flag_mask, max_plp, error, maxcnt;
	bam_pileup1_t *plp;
	// for the "auto" interface only
	bam1_t *b;
	bam_plp_auto_f func;
	void *data;
};

struct __bam_mplp_t {
	int n;
	uint64_t min, *pos;
	bam_plp_t *iter;
	int *n_plp;
	const bam_pileup1_t **plp;
};

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

static inline lbnode_t *mp_alloc(mempool_t *mp) {
	++mp->cnt;
	if(mp->n == 0)
		return (lbnode_t*)calloc(1, sizeof(lbnode_t));
	else
		return mp->buf[--mp->n];
}

static mempool_t *mp_init() {
	mempool_t *mp;
	mp = (mempool_t*)calloc(1, sizeof(mempool_t));
	return mp;
}

static void mp_destroy(mempool_t *mp) {
	int k;
	for(k = 0; k < mp->n; ++k) {
		free(mp->buf[k]->b.data);
		free(mp->buf[k]);
	}
	free(mp->buf);
	free(mp);
}

static inline void mp_free(mempool_t *mp, lbnode_t *p) {
	--mp->cnt;
	p->next = 0; // clear lbnode_t::next here
	if(mp->n == mp->max) {
		mp->max = mp->max ? mp->max<<1 : 256;
		mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
	}
	mp->buf[mp->n++] = p;
}

static inline int resolve_cigar2(bam_pileup1_t *p, uint32_t pos, cstate_t *s) {
	// s->k: the index of the CIGAR operator that has just been processed.
	// s->x: the reference coordinate of the start of s->k
	// s->y: the query coordiante of the start of s->k

	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam1_cigar(b);
	int k, is_head = 0;
	// determine the current CIGAR operation
	// fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam1_qname(b), pos, s->end, s->k, s->x, s->y);
	if(s->k == -1) { // never processed
		is_head = 1;
		if(c->n_cigar == 1) { // just one operation, save a loop
			if(_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF)
				s->k = 0, s->x = c->pos, s->y = 0;
		}
		else { // find the first match or deletion
			for(k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
				int op = _cop(cigar[k]);
				int l  = _cln(cigar[k]);
				if(op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF)
					break;
				else if(op == BAM_CREF_SKIP)
					s->x += l;
				else if(op == BAM_CINS || op == BAM_CSOFT_CLIP)
					s->y += l;
			}
			assert(k < c->n_cigar);
			s->k = k;
		}
	}
	else { // the read has been processed before
		int op, l = _cln(cigar[s->k]);
		if(pos - s->x >= l) { // jump to the next operation
			assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
			op = _cop(cigar[s->k+1]);
			if(op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
				if(_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF)
					s->y += l;
				s->x += l;
				++s->k;
			}
			else { // find the next M/D/N/=/X
				if(_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF)
					s->y += l;
				s->x += l;
				for(k = s->k + 1; k < c->n_cigar; ++k) {
					op = _cop(cigar[k]), l = _cln(cigar[k]);
					if(op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF)
						break;
					else if(op == BAM_CINS || op == BAM_CSOFT_CLIP)
						s->y += l;
				}
				s->k = k;
			}
			assert(s->k < c->n_cigar); // otherwise a bug
		} // else, do nothing
	}
	{ // collect pileup information
		int op, l;
		op = _cop(cigar[s->k]);
		l  = _cln(cigar[s->k]);
		p->is_del = p->indel = p->is_refskip = 0;
		if(s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
			int op2 = _cop(cigar[s->k+1]);
			int l2  = _cln(cigar[s->k+1]);
			if(op2 == BAM_CDEL)
				p->indel = -(int)l2;
			else if(op2 == BAM_CINS)
				p->indel = l2;
			else if(op2 == BAM_CPAD && s->k + 2 < c->n_cigar) { // no working for adjacent padding
				int l3 = 0;
				for(k = s->k + 2; k < c->n_cigar; ++k) {
					op2 = _cop(cigar[k]);
					l2  = _cln(cigar[k]);
					if(op2 == BAM_CINS)
						l3 += l2;
					else if(op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF)
						break;
				}
				if(l3 > 0)
					p->indel = l3;
			}
		}
		if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			p->qpos = s->y + (pos - s->x);
		}
		else if(op == BAM_CDEL || op == BAM_CREF_SKIP) {
			p->is_del = 1;
			p->qpos = s->y; // FIXME: distinguish D and N!!!!!
			p->is_refskip = (op == BAM_CREF_SKIP);
		} // cannot be other operations; otherwise a bug
		p->is_head = (pos == c->pos);
		p->is_tail = (pos == s->end);
	}
	return 1;
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp) {
	if(iter->error) {
		*_n_plp = -1;
		return 0;
	}
	*_n_plp = 0;
	if(iter->is_eof && iter->head->next == 0)
		return 0;
	while(iter->is_eof || iter->max_tid > iter->tid || (iter->max_tid == iter->tid && iter->max_pos > iter->pos)) {
		int n_plp = 0;
		lbnode_t *p, *q;
		// write iter->plp at iter->pos
		iter->dummy->next = iter->head;
		for(p = iter->head, q = iter->dummy; p->next; q = p, p = p->next) {
			if(p->b.core.tid < iter->tid || (p->b.core.tid == iter->tid && p->end <= iter->pos)) { // then remove
				q->next = p->next;
				mp_free(iter->mp, p);
				p = q;
			}
			else if(p->b.core.tid == iter->tid && p->beg <= iter->pos) { // here: p->end > pos; then add to pileup
				if(n_plp == iter->max_plp) { // then double the capacity
					iter->max_plp = iter->max_plp ? iter->max_plp<<1 : 256;
					iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
				}
				iter->plp[n_plp].b = &p->b;
				if(resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s))
					++n_plp; // actually always true...
			}
		}
		iter->head = iter->dummy->next; // dummy->next may be changed
		*_n_plp = n_plp;
		*_tid   = iter->tid;
		*_pos   = iter->pos;
		// update iter->tid and iter->pos
		if(iter->head->next) {
			if(iter->tid > iter->head->b.core.tid) {
				fprintf(stderr, "[%s] unsorted input. Pileup aborts.\n", __func__);
				iter->error = 1;
				*_n_plp = -1;
				return 0;
			}
		}
		if(iter->tid < iter->head->b.core.tid) { // come to a new reference sequence
			iter->tid = iter->head->b.core.tid;
			iter->pos = iter->head->beg; // jump to the next reference
		}
		else if(iter->pos < iter->head->beg) { // here: tid == head->b.core.tid
			iter->pos = iter->head->beg; // jump to the next position
		}
		else
			++iter->pos; // scan contiguously
		// return
		if(n_plp)
			return iter->plp;
		if(iter->is_eof && iter->head->next == 0)
			break;
	}
	return 0;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b) {
	if(iter->error)
		return -1;
	if(b) {
		if(b->core.tid < 0)
			return 0;
		if(b->core.flag & iter->flag_mask)
			return 0;
		if(iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt)
			return 0;
		bam_copy1(&iter->tail->b, b);
		iter->tail->beg   = b->core.pos;
		iter->tail->end   = bam_calend(&b->core, bam1_cigar(b));
		iter->tail->s     = g_cstate_null;
		iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
		if(b->core.tid < iter->max_tid) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (chromosomes out of order)\n");
			iter->error = 1;
			return -1;
		}
		if((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (reads out of order)\n");
			iter->error = 1;
			return -1;
		}
		iter->max_tid = b->core.tid;
		iter->max_pos = iter->tail->beg;
		if(iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
			iter->tail->next = mp_alloc(iter->mp);
			iter->tail = iter->tail->next;
		}
	}
	else
		iter->is_eof = 1;
	return 0;
}

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp) {
	const bam_pileup1_t *plp;
	if(iter->func == 0 || iter->error) {
		*_n_plp = -1;
		return 0;
	}
	if((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
		return plp;
	else { // no pileup line can be obtained; read alignments
		*_n_plp = 0;
		if(iter->is_eof)
			return 0;
		while(iter->func(iter->data, iter->b) >= 0) {
			if(bam_plp_push(iter, iter->b) < 0) {
				*_n_plp = -1;
				return 0;
			}
			if((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
				return plp;
			// otherwise no pileup line can be returned; read the next alignment.
		}
		bam_plp_push(iter, 0);
		if((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
			return plp;
		return 0;
	}
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp) {
	int i, ret = 0;
	uint64_t new_min = (uint64_t)-1;
	for(i = 0; i < iter->n; ++i) {
		if(iter->pos[i] == iter->min) {
			int tid, pos;
			iter->plp[i] = bam_plp_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
			iter->pos[i] = (uint64_t)tid<<32 | pos;
		}
		if(iter->plp[i] && iter->pos[i] < new_min)
			new_min = iter->pos[i];
	}
	iter->min = new_min;
	if(new_min == (uint64_t)-1)
		return 0;
	*_tid = new_min>>32;
	*_pos = (uint32_t)new_min;
	for(i = 0; i < iter->n; ++i) {
		if(iter->pos[i] == iter->min) { // FIXME: valgrind reports "uninitialised value(s) at this line"
			n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
			++ret;
		}
		else
			n_plp[i] = 0, plp[i] = 0;
	}
	return ret;
}

void bam_plp_destroy(bam_plp_t iter) {
	mp_free(iter->mp, iter->dummy);
	mp_free(iter->mp, iter->head);
	if(iter->mp->cnt != 0)
		fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
	mp_destroy(iter->mp);
	if(iter->b)
		bam_destroy1(iter->b);
	free(iter->plp);
	free(iter);
}

void bam_mplp_destroy(bam_mplp_t iter) {
	int i;
	for(i = 0; i < iter->n; ++i)
		bam_plp_destroy(iter->iter[i]);
	free(iter->iter);
	free(iter->pos);
	free(iter->n_plp);
	free(iter->plp);
	free(iter);
}

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data) {
	bam_plp_t iter;
	iter = (bam_plp_t)calloc(1, sizeof(struct __bam_plp_t));
	iter->mp        = mp_init();
	iter->head      = iter->tail = mp_alloc(iter->mp);
	iter->dummy     = mp_alloc(iter->mp);
	iter->max_tid   = iter->max_pos = -1;
	iter->flag_mask = BAM_DEF_MASK;
	iter->maxcnt    = 8000;
	if(func) {
		iter->func = func;
		iter->data = data;
		iter->b    = bam_init1();
	}
	return iter;
}

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data) {
	int i;
	bam_mplp_t iter;
	iter = (bam_mplp_t)calloc(1, sizeof(struct __bam_mplp_t));
	iter->pos   = (uint64_t*)calloc(n, 8);
	iter->n_plp = (int*)calloc(n, sizeof(int));
	iter->plp   = (const bam_pileup1_t**)calloc(n, sizeof(void*));
	iter->iter  = (bam_plp_t*)calloc(n, sizeof(void*));
	iter->n     = n;
	iter->min   = (uint64_t)-1;
	for(i = 0; i < n; ++i) {
		iter->iter[i] = bam_plp_init(func, data[i]);
		iter->pos[i]  = iter->min;
	}
	return iter;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt) {
	int i;
	for(i = 0; i < iter->n; ++i)
		iter->iter[i]->maxcnt = maxcnt;
}

/*------------------------------------------------------------------------------------------------*/
// bam_import.c

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};


/*------------------------------------------------------------------------------------------------*/
// bam_md.c

char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

/*
int bam_prob_realn_core(bam1_t *b, const char *ref, int flag) {
	int k, i, bw, x, y, yb, ye, xb, xe, apply_baq = flag&1, extend_baq = flag>>1&1;
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	kpa_par_t conf = kpa_par_def;
	uint8_t *bq = 0, *zq = 0, *qual = bam1_qual(b);
	if((c->flag & BAM_FUNMAP) || b->core.l_qseq == 0)
		return -1; // do nothing
	// test if BQ or ZQ is present
	if((bq = bam_aux_get(b, "BQ")) != 0)
		++bq;
	if((zq = bam_aux_get(b, "ZQ")) != 0 && *zq == 'Z')
		++zq;
	if(bq && zq) { // remove the ZQ tag
		bam_aux_del(b, zq-1);
		zq = 0;
	}
	if(bq || zq) {
		if((apply_baq && zq) || (!apply_baq && bq))
			return -3; // in both cases, do nothing
		if(bq && apply_baq) { // then convert BQ to ZQ
			for(i = 0; i < c->l_qseq; ++i)
				qual[i] = qual[i] + 64 < bq[i] ? 0 : qual[i] - ((int)bq[i] - 64);
			*(bq - 3) = 'Z';
		}
		else if(zq && !apply_baq) { // then convert ZQ to BQ
			for(i = 0; i < c->l_qseq; ++i)
				qual[i] += (int)zq[i] - 64;
			*(zq - 3) = 'B';
		}
		return 0;
	}
	// find the start and end of the alignment
	x = c->pos, y = 0, yb = ye = xb = xe = -1;
	for(k = 0; k < c->n_cigar; ++k) {
		int op, l;
		op = cigar[k]&0xf;
		l  = cigar[k]>>4;
		if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			if(yb < 0)
				yb = y;
			if(xb < 0)
				xb = x;
			ye = y + l;
			xe = x + l;
			x += l;
			y += l;
		}
		else if(op == BAM_CSOFT_CLIP || op == BAM_CINS)
			y += l;
		else if(op == BAM_CDEL)
			x += l;
		else if(op == BAM_CREF_SKIP)
			return -1; // do nothing if there is a reference skip
	}
	// set bandwidth and the start and the end
	bw = 7;
	if(abs((xe - xb) - (ye - yb)) > bw)
		bw = abs((xe - xb) - (ye - yb)) + 3;
	conf.bw = bw;
	xb -= yb + bw/2;
	if(xb < 0)
		xb = 0;
	xe += c->l_qseq - ye + bw/2;
	if(xe - xb - c->l_qseq > bw)
		xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;
	{ // glocal
		uint8_t *s, *r, *q, *seq = bam1_seq(b), *bq;
		int *state;
		bq = (uint8_t*)calloc(c->l_qseq + 1, 1);
		memcpy(bq, qual, c->l_qseq);
		s = (uint8_t*)calloc(c->l_qseq, 1);
		for(i = 0; i < c->l_qseq; ++i)
			s[i] = bam_nt16_nt4_table[bam1_seqi(seq, i)];
		r = (uint8_t*)calloc(xe - xb, 1);
		for(i = xb; i < xe; ++i) {
			if(ref[i] == 0) {
				xe = i;
				break;
			}
			r[i-xb] = bam_nt16_nt4_table[bam_nt16_table[(int)ref[i]]];
		}
		state = (int*)calloc(c->l_qseq, sizeof(int));
		q = (uint8_t*)calloc(c->l_qseq, 1);
		kpa_glocal(r, xe-xb, s, c->l_qseq, qual, &conf, state, q);
		if(!extend_baq) { // in this block, bq[] is capped by base quality qual[]
			for(k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for(i = y; i < y + l; ++i) {
						if((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))
							bq[i] = 0;
						else
							bq[i] = bq[i] < q[i] ? bq[i] : q[i];
					}
					x += l;
					y += l;
				}
				else if(op == BAM_CSOFT_CLIP || op == BAM_CINS)
					y += l;
				else if(op == BAM_CDEL)
					x += l;
			}
			for(i = 0; i < c->l_qseq; ++i)
				bq[i] = qual[i] - bq[i] + 64; // finalize BQ
		}
		else { // in this block, bq[] is BAQ that can be larger than qual[] (different from the above!)
			uint8_t *left, *rght;
			left = (uint8_t*)calloc(c->l_qseq, 1);
			rght = (uint8_t*)calloc(c->l_qseq, 1);
			for(k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for(i = y; i < y + l; ++i)
						bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y)) ? 0 : q[i];
					for(left[y] = bq[y], i = y + 1; i < y + l; ++i)
						left[i] = bq[i] > left[i-1] ? bq[i] : left[i-1];
					for(rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
						rght[i] = bq[i] > rght[i+1] ? bq[i] : rght[i+1];
					for(i = y; i < y + l; ++i)
						bq[i] = left[i] < rght[i] ? left[i] : rght[i];
					x += l;
					y += l;
				}
				else if(op == BAM_CSOFT_CLIP || op == BAM_CINS)
					y += l;
				else if(op == BAM_CDEL)
					x += l;
			}
			for(i = 0; i < c->l_qseq; ++i)
				bq[i] = 64 + (qual[i] <= bq[i] ? 0 : qual[i] - bq[i]); // finalize BQ
			free(left);
			free(rght);
		}
		if(apply_baq) {
			for(i = 0; i < c->l_qseq; ++i)
				qual[i] -= bq[i] - 64; // modify qual
			bam_aux_append(b, "ZQ", 'Z', c->l_qseq + 1, bq);
		}
		else
			bam_aux_append(b, "BQ", 'Z', c->l_qseq + 1, bq);
		free(bq);
		free(s);
		free(r);
		free(q);
		free(state);
	}
	return 0;
}
*/

/*------------------------------------------------------------------------------------------------*/
// sam_header.c

/*
struct _HeaderList {
	struct _HeaderList *last;   // Hack: Used and maintained only by list_append_to_end. Maintained in the root node only.
	struct _HeaderList *next;
	void *data;
};
typedef struct _HeaderList list_t;
typedef list_t HeaderDict;

typedef struct {
	char key[2];
	char *value;
}
HeaderTag;

typedef struct {
	char type[2];
	list_t *tags;
}
HeaderLine;

static void list_free(list_t *root) {
	list_t *l = root;
	while(root) {
		l    = root;
		root = root->next;
		free(l);
	}
}

static void sam_header_line_free(HeaderLine *hline) {
	list_t *tags = hline->tags;
	while(tags) {
		HeaderTag *tag = (HeaderTag*)(tags->data);
		free(tag->value);
		free(tag);
		tags = tags->next;
	}
	list_free(hline->tags);
	free(hline);
}

void sam_header_free(void *_header) {
	HeaderDict *header = (HeaderDict*)_header;
	list_t *hlines = header;
	while(hlines) {
		sam_header_line_free((HeaderLine*)(hlines->data));
		hlines = hlines->next;
	}
	list_free(header);
}

KHASH_MAP_INIT_STR(str, const char *)

void sam_tbl_destroy(void *h) {
	khash_t(str) *tbl = (khash_t(str)*)h;
	kh_destroy(str, tbl);
}
*/
























