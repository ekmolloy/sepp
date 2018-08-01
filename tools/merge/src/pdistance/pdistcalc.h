//
// Created by miken on 7/18/2018.
//

#ifndef MISC_C_UTILITIES_PDISTCALC_H
#define MISC_C_UTILITIES_PDISTCALC_H

typedef struct {
    int nPos;
    int nSeq;
    char **names;
    char **seqs;
    int nSaved; /* actual allocated size of names and seqs */
} alignment_t;

typedef struct {
    int pos;
    char *name;
} listentry_t;


alignment_t *ReadAlignment(/*IN*/FILE *fp, bool bQuote);
void *myfree(void *p, size_t sz);
void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy);
void *mymemdup(void *data, size_t sz);
void *mymalloc(size_t sz);



#endif //MISC_C_UTILITIES_PDISTCALC_H
