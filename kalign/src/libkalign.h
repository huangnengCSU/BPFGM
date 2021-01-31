#include "global.h"
#include "msa.h"
#include "parameters.h"
//#include "align_io.h"
#include "alignment_parameters.h"

#include "bisectingKmeans.h"
/* #include "alignment.h" */
#include "weave_alignment.h"

#include "esl_stopwatch.h"
//#include "../src/misc.h"
#include <getopt.h>
#include <string.h>
#include "idata.h"
#include "alphabet.h"

#include "aln_run.h"
#include "tlrng.h"
#include "rwalign.h"
#include "aln_task.h"

typedef struct MSAOutput
{
    char **out_ptr;
    int num_seq;
    int len;
} MSAOut;

MSAOut new_MSAOut(int num_seq, int len);

struct msa *read_input_from_array(char **readArray, int nReads);

MSAOut RunMSA(struct msa *m);

void free_MSAOut(MSAOut *out);
