#include "libkalign.h"

#define OPT_SET 1
#define OPT_ALNPARAM 2
#define OPT_RENAME 3
#define OPT_REFORMAT 4
#define OPT_SHOWW 5
#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8
#define OPT_MATADD 9

#define OPT_DEVTEST 10
#define OPT_CHAOS 11
#define OPT_DUMP_INTERNAL 12
#define OPT_NTHREADS 14

static int check_for_sequences(struct msa *msa);
static int detect_aligned(struct msa *msa);
static int detect_alphabet(struct msa *msa);
static int set_sip_nsip(struct msa *msa);

MSAOut new_MSAOut(int num_seq, int len)
{
    MSAOut out;
    out.num_seq = num_seq;
    out.len = len;
    char **ptr = (char **)malloc(sizeof(char *) * num_seq);
    for (int i = 0; i < num_seq; i++)
    {
        ptr[i] = (char *)malloc(sizeof(char) * len);
    }
    out.out_ptr = ptr;
    return out;
}

void free_MSAOut(MSAOut *out){
    for(int i=0;i<out->num_seq;i++){
        free(out->out_ptr[i]);
    }
    free(out->out_ptr);
}

int check_for_sequences(struct msa *msa)

{
    if (!msa)
    {
        ERROR_MSG("No sequences were found in the input files or standard input.");
    }
    if (msa->numseq < 2)
    {
        if (msa->numseq == 0)
        {
            ERROR_MSG("No sequences were found in the input files or standard input.");
        }
        else if (msa->numseq == 1)
        {
            ERROR_MSG("Only 1 sequence was found in the input files or standard input");
        }
    }
    return OK;
ERROR:
    return FAIL;
}

int detect_alphabet(struct msa *msa)
{
    int i;
    int min, c;
    uint8_t DNA[128];
    uint8_t protein[128];
    int diff[3];
    char DNA_letters[] = "acgtuACGTUnN";
    char protein_letters[] = "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY";

    ASSERT(msa != NULL, "No alignment");

    for (i = 0; i < 128; i++)
    {
        DNA[i] = 0;
        protein[i] = 0;
    }

    for (i = 0; i < (int)strlen(DNA_letters); i++)
    {
        DNA[(int)DNA_letters[i]] = 1;
    }

    for (i = 0; i < (int)strlen(protein_letters); i++)
    {
        protein[(int)protein_letters[i]] = 1;
    }

    diff[0] = 0;
    diff[1] = 0;
    for (i = 0; i < 128; i++)
    {
        if ((msa->letter_freq[i]) && (!DNA[i]))
        {
            diff[0]++;
        }
        if ((msa->letter_freq[i]) && (!protein[i]))
        {
            diff[1]++;
        }
    }
    // for (i = 0; i < 128; i++)
    // {
    //     fprintf(stdout, "%d,", msa->letter_freq[i]);
    // }
    // fprintf(stdout, "\n");
    // fprintf(stdout, "diff0:%d,diff1:%d\n", diff[0], diff[1]);

    if (diff[0] + diff[1] == 0)
    {
        ERROR_MSG("Could not detect any AA or nucleotides.");
    }
    c = -1;
    min = 2147483647;
    for (i = 0; i < 2; i++)
    {
        if (diff[i] < min)
        {
            min = diff[i];
            c = i;
        }
    }
    if (c == 0)
    {
        // LOG_MSG("Detected DNA sequences.");
        msa->L = defDNA;
        RUN(convert_msa_to_internal(msa, defDNA));
    }
    else if (c == 1)
    {
        // LOG_MSG("Detected protein sequences.");
        msa->L = redPROTEIN;
        RUN(convert_msa_to_internal(msa, redPROTEIN));
    }
    else
    {
        ERROR_MSG("Alphabet not recognized.");
    }
    return OK;
ERROR:
    return FAIL;
}

int detect_aligned(struct msa *msa)
{
    int min_len;
    int max_len;
    int l;
    int i;
    int j;
    int n;
    int gaps = 0;
    /* assume that sequences are not aligned */
    msa->aligned = 0;

    /* Improved logic:
           Lets sum up the number of gaps plus sequence length of the first
           X sequences. if min == max it is aligned.
        */
    min_len = INT32_MAX;
    max_len = 0;
    gaps = 0;
    n = MACRO_MIN(50, msa->numseq);
    for (i = 0; i < n; i++)
    {
        l = 0;
        for (j = 0; j <= msa->sequences[i]->len; j++)
        {
            l += msa->sequences[i]->gaps[j];
        }
        gaps += l;
        l += msa->sequences[i]->len;
        min_len = MACRO_MIN(min_len, l);
        max_len = MACRO_MAX(max_len, l);
    }
    if (gaps)
    {
        if (min_len == max_len)
        { /* sequences have gaps and total length is identical - clearly aligned  */
            msa->aligned = ALN_STATUS_ALIGNED;
        }
        else
        { /* odd there are gaps but total length differs - unknown status  */
            msa->aligned = ALN_STATUS_UNKNOWN;
        }
    }
    else
    {
        if (min_len == max_len)
        { /* no gaps and sequences have same length. Can' tell if they are aligned  */
            msa->aligned = ALN_STATUS_UNKNOWN;
        }
        else
        { /* No gaps and sequences have different lengths - unaligned */
            msa->aligned = ALN_STATUS_UNALIGNED;
        }
    }
    //LOG_MSG("Aligned: %d gaps: %d",msa->aligned,gaps);
    return OK;
}

int set_sip_nsip(struct msa *msa)
{
    int i;
    ASSERT(msa != NULL, "No msa");
    if (msa->plen)
    {
        for (i = msa->num_profiles; i--;)
        {
            if (msa->sip[i])
            {
                MFREE(msa->sip[i]);
            }
        }
        if (msa->plen)
        {
            MFREE(msa->plen);
        }
        if (msa->sip)
        {
            MFREE(msa->sip);
        }
        if (msa->nsip)
        {
            MFREE(msa->nsip);
        }
        msa->plen = NULL;
        msa->sip = NULL;
        msa->nsip = NULL;
    }

    msa->num_profiles = (msa->numseq << 1) - 1;

    MMALLOC(msa->sip, sizeof(int *) * msa->num_profiles);
    MMALLOC(msa->nsip, sizeof(int) * msa->num_profiles);
    MMALLOC(msa->plen, sizeof(int) * msa->num_profiles);

    for (i = 0; i < msa->num_profiles; i++)
    {
        msa->sip[i] = NULL;
        msa->nsip[i] = 0;
    }

    for (i = 0; i < msa->numseq; i++)
    {

        MMALLOC(msa->sip[i], sizeof(int));
        msa->nsip[i] = 1;
        msa->sip[i][0] = i;
        msa->plen[i] = 0;
    }
    return OK;
ERROR:
    return FAIL;
}

struct msa *read_input_from_array(char **readArray, int nReads)
{
    struct msa *msa = NULL;
    struct msa_seq *seq_ptr = NULL;
    if (msa == NULL)
    {
        msa = alloc_msa();
    }
    for (int k = 0; k < nReads; k++)
    {
        if (msa->alloc_numseq == msa->numseq)
        {
            resize_msa(msa);
        }
        seq_ptr = msa->sequences[msa->numseq];
        // char name[16] = {};
        // sprintf(name, "%d", k);
        // seq_ptr->name = name;
        msa->numseq++;
        int i = 0;
        while (readArray[k][i] != '\0')
        {
            msa->letter_freq[(int)readArray[k][i]]++;
            if (isalpha((int)readArray[k][i]))
            {
                if (seq_ptr->alloc_len == seq_ptr->len)
                {
                    resize_msa_seq(seq_ptr);
                }
                seq_ptr->seq[seq_ptr->len] = readArray[k][i];
                seq_ptr->len++;
            }
            if (ispunct((int)readArray[k][i]))
            {
                seq_ptr->gaps[seq_ptr->len]++;
            }
            i++;
        }
        msa->letter_freq[(int)readArray[k][i]]++;
    }
    // for (int k = 0; k < nReads; k++)
    // {
    //     fprintf(stdout, "%s\n", msa->sequences[k]->seq);
    // }
    null_terminate_sequences(msa);
    detect_alphabet(msa);
    detect_aligned(msa);
    set_sip_nsip(msa);
    return msa;
}

MSAOut RunMSA(struct msa *m)
{
    struct parameters *param;
    param = init_param();
    param->gpo = 4.0;
    param->gpe = 2.0;
    param->tgpe = 1.0;
    param->nthreads = 2;

    struct aln_param *ap = NULL;
    struct aln_tasks *tasks = NULL;
    int **map = NULL;
    check_for_sequences(m);
    // LOG_MSG("Detected: %d sequences.", m->numseq);
    init_ap(&ap, param, m->numseq, m->L);
    alloc_tasks(&tasks, m->numseq);
    build_tree_kmeans(m, ap, &tasks);
    // LOG_MSG("Aligning");
    create_msa_serial(m,ap, tasks);
    // create_msa_openMP(m, ap, tasks);
    m->aligned = ALN_STATUS_ALIGNED;
    // LOG_MSG("Weaving");
    weave(m, tasks);
    // 输出多序列比对结果
    int i, j, c, f;
    int msa_aligned_len = 0;
    for(j=0;j<m->sequences[0]->len;j++){
        for (c = 0; c < m->sequences[i]->gaps[j]; c++){
            msa_aligned_len++;  //number of gaps
        }
        msa_aligned_len++; // base
    }
    for (c = 0; c < m->sequences[i]->gaps[m->sequences[i]->len]; c++){
        msa_aligned_len++; // end gaps
    }
    msa_aligned_len++; // end symbol \0
    // fprintf(stdout,"len=%d\n",msa_aligned_len);
    MSAOut msa_out = new_MSAOut(m->numseq, msa_aligned_len);
    for (i = 0; i < m->numseq; i++)
    {
        f = 0;
        // fprintf(stdout, ">%s\n", m->sequences[i]->name);
        for (j = 0; j < m->sequences[i]->len; j++)
        {

            for (c = 0; c < m->sequences[i]->gaps[j]; c++)  // read[i][j]前面的gaps
            {
                msa_out.out_ptr[i][f] = '-';
                f++;
                // fprintf(stdout, "-");
            }
            msa_out.out_ptr[i][f] = m->sequences[i]->seq[j];
            f++;
            // fprintf(stdout, "%c", m->sequences[i]->seq[j]); // read[i][j]
        }
        for (c = 0; c < m->sequences[i]->gaps[m->sequences[i]->len]; c++)   // read[i][last j]的gaps
        {
            msa_out.out_ptr[i][f] = '-';
            f++;
            // fprintf(stdout, "-");
        }
        msa_out.out_ptr[i][f] = '\0';
        // fprintf(stdout, "\n");
    }
    // LOG_MSG("free memories.");
    free_tasks(tasks);
    free_ap(ap);
    free_parameters(param);
    return msa_out;
}
