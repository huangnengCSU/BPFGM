#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libkalign.h"

int main(){
	char *reads[] = {"TGAGAA","TGAGAA","TGAGAA","TGAGAAA","TGAGA"};
	int num_reads = 0;
	int i=0;
	num_reads = sizeof(reads);
	printf("size of reads: %d\n",num_reads);
	char **p;
	p = (char**)calloc(5,sizeof(char*));
	
	for(int j=0;j<5;j++){
		printf("%s\t, %ld\n",reads[j],strlen(reads[j]));
		char *pc = malloc(strlen(reads[j])*sizeof(char)+1);
		strcpy(pc,reads[j]);
		p[j] = pc;	
	}
	struct msa *msa;
    read_input_from_array(p, 5, msa);
	return 0;
}
