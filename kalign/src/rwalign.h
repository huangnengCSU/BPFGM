#include "msa.h"
#include "config.h"

struct msa* alloc_msa(void);
int resize_msa(struct msa* msa);
void free_msa(struct msa* msa);
int resize_msa_seq(struct msa_seq* seq);
void free_msa_seq(struct msa_seq* seq);
int null_terminate_sequences(struct msa* msa);