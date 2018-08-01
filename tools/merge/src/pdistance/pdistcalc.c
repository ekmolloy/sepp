//
// Created by miken on 7/18/2018.
//
/*   ***********************************************
 *   P-Distance Calcuation Utility for SEPP and HIPPI
 *
 *                  Mike Nute
 *                July 19, 2018
 *   ***********************************************
 *
 *   This utility takes an alignment in FASTA form either from stdin or
 *   from a file and calculates the average normalized hamming distance (ANHD)
 *   between the sequences in it. Optionally a list of sequence names can
 *   be given in a stream separated by newlines, and the ANHD can be calculated
 *   just on that subset. This program was able to calculate the ANHD for
 *   an alignment of ~17100 sequences with ~2000 columns in 7 minutes. (That
 *   is from memory).
 *
 *   PLEASE NOTE: some chunks of this code were copied directly from the
 *   source code for FastTree by Morgan Price. The website for that program
 *   is http://www.microbesonline.org/fasttree/. Specifically, the FastTree
 *   struct definition for a multiple sequence alignment was quite nice, es-
 *   pecially combined with a robust function for reading and parsing a
 *   Fasta of Phylip file. So that function is copied almost exactly, and it
 *   uses a few other memory management functions in the FastTree code, so
 *   those have been copied over as well. This file is self-contained.
 *
 * */



#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
// #include "pdistcalc.h"

//
//  Struct Definitions and Prototypes
//
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

// These three functions were copied and pasted from the source code for FastTree
//   by Morgan Price
alignment_t *ReadAlignment(/*IN*/FILE *fp, bool bQuote);
void *myfree(void *p, size_t sz);
void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy);
void *mymemdup(void *data, size_t sz);
void *mymalloc(size_t sz);



#define MAXCODES 20
#define NOCODE 127
/* Note -- sequence lines longer than BUFFER_SIZE are
   allowed, but FASTA header lines must be within this limit */
#define BUFFER_SIZE 5000
#define MAX_SEQ_NAME_SIZE 500
#define MIN(X,Y) ((X) <  (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) >  (Y) ? (X) : (Y))

#ifdef USE_SSE3
#define SSE_STRING "SSE3"
#define ALIGNED __attribute__((aligned(16)))
#define IS_ALIGNED(X) ((((unsigned long) new) & 15L) == 0L)
//#include <xmmintrin.h>

#else

#define ALIGNED
#define IS_ALIGNED(X) 1

#ifndef USE_DOUBLE
#define SSE_STRING "No SSE3"
#endif

#endif /* USE_SSE3 */

long szAllAlloc = 0;
long mymallocUsed = 0;		/* useful allocations by mymalloc */

char *usage =
        "Usage for P-Distance Calculator: (Result is printed as a 8-digit float to    "
        "STDOUT)                                                                    \n"
        "                                                                           \n"
        "basic (alignment & subset via files:                                       \n"
        "    ./pd -aln <aln_path> -sub <sub_path> -n_sub <# seqs>                   \n"
        "                                                                           \n"
        "Required Options:                                                          \n"
        "    -aln        path to FASTA file containing the alignment. If this is    \n"
        "                not given, alignment sequences are assumed to come via     \n"
        "                STDIN. If the alignment comes via stdin, nothing else      \n"
        "                can be given via STDIN.                                    \n"
        "                                                                           \n"
        "    -sub        path to file containing a subset of sequences to be        \n"
        "                used for the computation. File should be one sequence      \n"
        "                name per line, separated by newlines. Specified            \n"
        "                sequences that are not actually in the fasta are           \n"
        "                ignored. If the alignment file is provided as a path,      \n"
        "                this can be omitted and given through STDIN.               \n"
        "                                                                           \n"
        "    -n_sub      Number of sequences in the subset list. If omitted,        \n"
        "                calculates p-dist on the whole file. If this is not        \n"
        "                given, the subset file will not be used at all.            \n"
        "                                                                           \n"
        "Additional Flags:                                                          \n"
        "    -max        outputs the maximum p-distance instead of the average      \n"
        "                pairwise. (default: off)                                   \n"
        "                                                                           \n"
        "    -verbose    if given, program prints several intermediate values to    \n"
        "                stdout. This includes counts of total gapped residues      \n"
        "                and total matching residues, as well as the ``weighted``   \n"
        "                p-distance (rather than averaging all the pairwise PDs).   \n"
        "                                                                           \n"
        "    -debug      Goes into debug mode, prints some stuff and then exits.    \n"
        "                Currently does a small unit test to make sure reading      \n"
        "                the sequences was working. (default: off)                  \n"
        "                                                                           \n"
        "    -h          Prints this help description and then terminates           \n";


//line 3629 in fasttree.c
alignment_t *ReadAlignment(/*IN*/FILE *fp, bool bQuote) {
    /* bQuote supports the -quote option */
    int nSeq = 0;
    int nPos = 0;
    char **names = NULL;
    char **seqs = NULL;
    char buf[BUFFER_SIZE] = "";
    if (fgets(buf,sizeof(buf),fp) == NULL) {
        fprintf(stderr, "Error reading header line\n");
        exit(1);
    }
    int nSaved = 100;
    if (buf[0] == '>') {
        /* FASTA, truncate names at any of these */
        char *nameStop = bQuote ? "'\t\r\n" : "(),: \t\r\n";
        char *seqSkip = " \t\r\n";	/* skip these characters in the sequence */
        seqs = (char**)mymalloc(sizeof(char*) * nSaved);
        names = (char**)mymalloc(sizeof(char*) * nSaved);

        do {
            /* loop over lines */
            if (buf[0] == '>') {
                /* truncate the name */
                char *p, *q;
                for (p = buf+1; *p != '\0'; p++) {
                    for (q = nameStop; *q != '\0'; q++) {
                        if (*p == *q) {
                            *p = '\0';
                            break;
                        }
                    }
                    if (*p == '\0') break;
                }

                /* allocate space for another sequence */
                nSeq++;
                if (nSeq > nSaved) {
                    int nNewSaved = nSaved*2;
                    seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
                    names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
                    nSaved = nNewSaved;
                }
                names[nSeq-1] = (char*)mymemdup(buf+1,strlen(buf));
                seqs[nSeq-1] = NULL;
            } else {
                /* count non-space characters and append to sequence */
                int nKeep = 0;
                char *p, *q;
                for (p=buf; *p != '\0'; p++) {
                    for (q=seqSkip; *q != '\0'; q++) {
                        if (*p == *q)
                            break;
                    }
                    if (*p != *q)
                        nKeep++;
                }
                int nOld = (seqs[nSeq-1] == NULL) ? 0 : strlen(seqs[nSeq-1]);
                seqs[nSeq-1] = (char*)myrealloc(seqs[nSeq-1], nOld, nOld+nKeep+1, /*copy*/false);
                if (nOld+nKeep > nPos)
                    nPos = nOld + nKeep;
                char *out = seqs[nSeq-1] + nOld;
                for (p=buf; *p != '\0'; p++) {
                    for (q=seqSkip; *q != '\0'; q++) {
                        if (*p == *q)
                            break;
                    }
                    if (*p != *q) {
                        *out = *p;
                        out++;
                    }
                }
                assert(out-seqs[nSeq-1] == nKeep + nOld);
                *out = '\0';
            }
        } while(fgets(buf,sizeof(buf),fp) != NULL);

        if (seqs[nSeq-1] == NULL) {
            fprintf(stderr, "No sequence data for last entry %s\n",names[nSeq-1]);
            exit(1);
        }
        names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
        seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
    } else {
        /* PHYLIP interleaved-like format
           Allow arbitrary length names, require spaces between names and sequences
           Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
           or not.
         */
        if (buf[0] == '\n' || buf[0] == '\r') {
            if (fgets(buf,sizeof(buf),fp) == NULL) {
                fprintf(stderr, "Empty header line followed by EOF\n");
                exit(1);
            }
        }
        if (sscanf(buf, "%d%d", &nSeq, &nPos) != 2
            || nSeq < 1 || nPos < 1) {
            fprintf(stderr, "Error parsing header line:%s\n", buf);
            exit(1);
        }
        names = (char **)mymalloc(sizeof(char*) * nSeq);
        seqs = (char **)mymalloc(sizeof(char*) * nSeq);
        nSaved = nSeq;

        int i;
        for (i = 0; i < nSeq; i++) {
            names[i] = NULL;
            seqs[i] = (char *)mymalloc(nPos+1);	/* null-terminate */
            seqs[i][0] = '\0';
        }
        int iSeq = 0;

        while(fgets(buf,sizeof(buf),fp)) {
            if ((buf[0] == '\n' || buf[0] == '\r') && (iSeq == nSeq || iSeq == 0)) {
                iSeq = 0;
            } else {
                int j = 0; /* character just past end of name */
                if (buf[0] == ' ') {
                    if (names[iSeq] == NULL) {
                        fprintf(stderr, "No name in phylip line %s", buf);
                        exit(1);
                    }
                } else {
                    while (buf[j] != '\n' && buf[j] != '\0' && buf[j] != ' ')
                        j++;
                    if (buf[j] != ' ' || j == 0) {
                        fprintf(stderr, "No sequence in phylip line %s", buf);
                        exit(1);
                    }
                    if (iSeq >= nSeq) {
                        fprintf(stderr, "No empty line between sequence blocks (is the sequence count wrong?)\n");
                        exit(1);
                    }
                    if (names[iSeq] == NULL) {
                        /* save the name */
                        names[iSeq] = (char *)mymalloc(j+1);
                        int k;
                        for (k = 0; k < j; k++) names[iSeq][k] = buf[k];
                        names[iSeq][j] = '\0';
                    } else {
                        /* check the name */
                        int k;
                        int match = 1;
                        for (k = 0; k < j; k++) {
                            if (names[iSeq][k] != buf[k]) {
                                match = 0;
                                break;
                            }
                        }
                        if (!match || names[iSeq][j] != '\0') {
                            fprintf(stderr, "Wrong name in phylip line %s\nExpected %s\n", buf, names[iSeq]);
                            exit(1);
                        }
                    }
                }
                int seqlen = strlen(seqs[iSeq]);
                for (; buf[j] != '\n' && buf[j] != '\0'; j++) {
                    if (buf[j] != ' ') {
                        if (seqlen >= nPos) {
                            fprintf(stderr, "Too many characters (expected %d) for sequence named %s\nSo far have:\n%s\n",
                                    nPos, names[iSeq], seqs[iSeq]);
                            exit(1);
                        }
                        seqs[iSeq][seqlen++] = toupper(buf[j]);
                    }
                }
                seqs[iSeq][seqlen] = '\0'; /* null-terminate */
//                if(verbose>10) fprintf(stderr,"Read iSeq %d name %s seqsofar %s\n", iSeq, names[iSeq], seqs[iSeq]);
                iSeq++;
                if (iSeq == nSeq && strlen(seqs[0]) == nPos)
                    break; /* finished alignment */
            } /* end else non-empty phylip line */
        }
        if (iSeq != nSeq && iSeq != 0) {
            fprintf(stderr, "Wrong number of sequences: expected %d\n", nSeq);
            exit(1);
        }
    }
    /* Check lengths of sequences */
    int i;
    for (i = 0; i < nSeq; i++) {
        int seqlen = strlen(seqs[i]);
        if (seqlen != nPos) {
            fprintf(stderr, "Wrong number of characters for %s: expected %d but have %d instead.\n"
                            "This sequence may be truncated, or another sequence may be too long.\n",
                    names[i], nPos, seqlen);
            exit(1);
        }
    }
    /* Replace "." with "-" and warn if we find any */
    /* If nucleotide sequences, replace U with T and N with X */
    bool findDot = false;
    for (i = 0; i < nSeq; i++) {
        char *p;
        for (p = seqs[i]; *p != '\0'; p++) {
            if (*p == '.') {
                findDot = true;
                *p = '-';
            }
//            if (nCodes == 4 && *p == 'U')
//                *p = 'T';
//            if (nCodes == 4 && *p == 'N')
//                *p = 'X';
        }
    }
    if (findDot)
        fprintf(stderr, "Warning! Found \".\" character(s). These are treated as gaps\n");

    if (ferror(fp)) {
        fprintf(stderr, "Error reading input file\n");
        exit(1);
    }

    alignment_t *align = (alignment_t*)mymalloc(sizeof(alignment_t));
    align->nSeq = nSeq;
    align->nPos = nPos;
    align->names = names;
    align->seqs = seqs;
    align->nSaved = nSaved;
    return(align);
}

void *mymalloc(size_t sz) {
    if (sz == 0) return(NULL);
    void *mnew = malloc(sz);
    if (mnew == NULL) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    szAllAlloc += sz;
    mymallocUsed += sz;
#ifdef TRACK_MEMORY
    struct mallinfo mi = mallinfo();
  if (mi.arena+mi.hblkhd > maxmallocHeap)
    maxmallocHeap = mi.arena+mi.hblkhd;
#endif
    /* gcc malloc should always return 16-byte-aligned values... */
    assert(IS_ALIGNED(mnew));
    return (mnew);
}

void *mymemdup(void *data, size_t sz) {
    if(data==NULL) return(NULL);
    void *mnew = mymalloc(sz);
    memcpy(/*to*/mnew, /*from*/data, sz);
    return(mnew);
}

void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy) {
    if (data == NULL && szOld == 0)
        return(mymalloc(szNew));
    if (data == NULL || szOld == 0 || szNew == 0) {
        fprintf(stderr,"Empty myrealloc\n");
        exit(1);
    }
    if (szOld == szNew)
        return(data);
    void *mnew = NULL;
    if (bCopy) {
        /* Try to reduce memory fragmentation by allocating anew and copying
           Seems to help in practice */
        mnew = mymemdup(data, szNew);
        myfree(data, szOld);
    } else {
        mnew = realloc(data,szNew);
        if (mnew == NULL) {
            fprintf(stderr, "Out of memory\n");
            exit(1);
        }
        assert(IS_ALIGNED(mnew));
        szAllAlloc += (szNew-szOld);
        mymallocUsed += (szNew-szOld);
#ifdef TRACK_MEMORY
        struct mallinfo mi = mallinfo();
    if (mi.arena+mi.hblkhd > maxmallocHeap)
      maxmallocHeap = mi.arena+mi.hblkhd;
#endif
    }
    return(mnew);
}

void *myfree(void *p, size_t sz) {
    if(p==NULL) return(NULL);
    free(p);
    mymallocUsed -= sz;
    return(NULL);
}

void zero_2_int_arrs(int* arr1, int* arr2, int len)
{
    int i = 0;
    for (i=0; i<len; i++) {
        arr1[i] = 0;
        arr2[i] = 0;
    }
}

int l_entry_compare(const void *a, const void *b)
{
    int j;
    listentry_t* a_le;
    listentry_t* b_le;
    a_le = *(listentry_t**)(a);
    b_le = *(listentry_t**)(b);
    return strcmp(a_le->name, b_le->name);
}
int l_entry_compare_search(const void *a, const void *b)
{
    char *a_le;
    listentry_t *b_le;
    a_le = (char*) (a);
    b_le = *(listentry_t **) (b);
    return strcmp(a_le, b_le->name);
}

int main(int argc, char *argv[]) {
    int ar;
    bool maxpd = false;
    bool debug = false;
    bool verbose = false;
    int aln_file = -1;
    int sub_file = -1;
    int n_sub = -1;
    int i,j,k, runct, i_ind, j_ind;
    double runtot = 0.;
    double pd;
    runct = 0;

    // Processing some options
    for (ar = 0; ar < argc; ar++)
    {
        if (strcmp(argv[ar],"-h")==0) {
            printf("%s\n",usage);
        } else if (strcmp(argv[ar],"-max")==0) {
            maxpd = true;
        } else if (strcmp(argv[ar],"-aln")==0) {
            aln_file = ar + 1;
        } else if (strcmp(argv[ar],"-sub")==0) {
            sub_file = ar + 1;
        } else if (strcmp(argv[ar],"-n_sub")==0) {
            if (argc>ar+1) {
                n_sub = atoi(argv[ar+1]);
            }
        } else if (strcmp(argv[ar],"-debug")==0) {
            debug = true;
            printf(" **** IN DEBUG MODE **** \n");
            printf("  -max\t%d\n",maxpd);
            printf("  -aln\t%d, %s\n",aln_file, (aln_file!=-1 ? argv[aln_file] : "(stdin)"));
            printf("  -sub\t%d\n", sub_file);
            printf("  -n_sub\t%d\n", n_sub);
        } else if (strcmp(argv[ar],"-verbose")==0) {
            verbose = true;
        }
    }

    // Reading the Input Alignment
    bool bQuote = false;
    alignment_t *aln;
    if (aln_file==-1) {
        aln = ReadAlignment(stdin,bQuote);
    } else {
        FILE *aln_f = fopen(argv[aln_file],"r");
        aln = ReadAlignment(aln_f,bQuote);
    }

    // Making a list of the sequences and noting the order it appears
    listentry_t** larr = (listentry_t**)malloc(aln->nSeq * sizeof(listentry_t*));
    for (i=0; i<aln->nSeq; i++) {
        larr[i] = malloc(sizeof(listentry_t));
        larr[i]->name = aln->names[i];
        larr[i]->pos = i;
//        larr[i]->name = (char*)malloc(200);
    }
    qsort(&larr[0], aln->nSeq, sizeof(listentry_t*),l_entry_compare);


    // Reading the Subset List
    int pos;
    char *myc;
    void *bsres;
    listentry_t *bsres_le;

    char** subset;
    int* subset_seq_inds;
    if (n_sub == -1) {
        if (verbose) printf("No sequence subset count found, computing on whole alignment\n");
        n_sub = aln->nSeq;
        subset_seq_inds = (int*)malloc(n_sub * sizeof(int));
        for (i=0; i<n_sub; i++) {
            subset_seq_inds[i] = i;
        }
        //
        // Make the appropriate dummy index here
        //
    } else {
        subset = (char**)malloc(n_sub * sizeof(char*));
        subset_seq_inds = (int*)malloc(n_sub * sizeof(int));

        FILE *subset_in = sub_file!=-1 ? fopen(argv[sub_file], "r") : stdin;
        for (i=0; i<n_sub; i++) {
            subset[i] = (char*)malloc(MAX_SEQ_NAME_SIZE);
            fgets(subset[i],MAX_SEQ_NAME_SIZE,subset_in);
            myc = strchr(subset[i],10);
            pos = myc-subset[i];
            subset[i][pos]=0;
            bsres = bsearch(subset[i],&larr[0],aln->nSeq, sizeof(listentry_t*),l_entry_compare_search);
            bsres_le = *(listentry_t**)bsres;
            subset_seq_inds[i] = bsres_le->pos;
        }
    }

    if (debug) {
        for (i=0; i<n_sub; i++) {
            printf("%s\t%d\n",subset[i],subset_seq_inds[i]);
            return 0;
        }
    }


    // Compute the average (or max) pairwise p-distance
    int* comp = (int*)malloc(n_sub * sizeof(int)); //(COMP)arable sequence residues (i.e. ungapped)
    int* same = (int*)malloc(n_sub * sizeof(int));
    long long all_comp = 0;
    long long all_same = 0;

    for (i=0; i<n_sub; i++){
        zero_2_int_arrs(comp,same,n_sub);
        i_ind = subset_seq_inds[i];

        for (k=0; k<aln->nPos; k++) {
            if (aln->seqs[i_ind][k]!=45) {
                for (j=0; j<i; j++) {
                    j_ind = subset_seq_inds[j];
                    if (aln->seqs[j_ind][k]!=45) {
                        comp[j]++;
                        if (aln->seqs[j_ind][k]==aln->seqs[i_ind][k]) {
                            same[j]++;
                        }
                    }
                }
            }
        }
        if (maxpd) {
            for (j=0; j<i; j++) {
                if (comp[j]==0) {continue;}
                pd = 1.0 - (double)same[j] / (double)comp[j];
                runtot = ( pd > runtot ? pd : runtot);
            }
        } else {
            for (j=0; j<i; j++) {
                if (comp[j]==0) continue;
                pd = 1.0 - (double)same[j] / (double)comp[j];
                all_comp = all_comp + comp[j];
                all_same = all_same + same[j];
                runtot = runtot + pd;
                runct++;
            }
        }
    }

    if (maxpd) {
        pd = runtot;
    } else {
        pd = runtot / (double) runct;
    }
    printf("%.8f\n",pd);

    if (verbose) {
        printf("   verbose results:\n");
        printf("   Running Total: %e\n",runtot);
        printf("   Running Count: %d\n",runct);
        printf("   Total Ungapped Residues: %ll\n",all_comp);
        printf("   Total Matching Residues: %ll\n", all_same);
        printf("   Weighted ANHD: %.8f\n\n",1.0 - (double)all_same / (double)all_comp);
    }

    free(comp);
    free(same);
    free(subset);
    free(subset_seq_inds);
    return(0);
}
