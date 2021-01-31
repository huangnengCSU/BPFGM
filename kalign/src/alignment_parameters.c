/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "alignment_parameters.h"

#include "alphabet.h"

#include "tlrng.h"

int set_subm_gaps(struct aln_param* ap);
int set_subm_gaps_DNA(struct aln_param* ap);
int set_param_number(struct aln_param* ap,int L, int sel);

int new_aln_matrices(struct aln_param* ap);


int init_ap(struct aln_param** aln_param, struct parameters* param, int numseq,int L)
{
        struct aln_param* ap = NULL;
        int i,j;

        if(*aln_param){
                ap = *aln_param;
        }else{

                MMALLOC(ap, sizeof(struct aln_param));

                ap->tree = NULL;
                MMALLOC(ap->tree, sizeof(int) * (numseq*3+1));

                for(i = 0;i < (numseq*3+1);i++){
                        ap->tree[i] = 0;
                }
                ap->subm = NULL;
                ap->chaos = param->chaos;
                RUNP(ap->rng = init_rng(42));
                MMALLOC(ap->subm,sizeof (float*) * 23);

                for (i = 23;i--;){
                        ap->subm[i] = NULL;
                        MMALLOC(ap->subm[i],sizeof(float) * 23);
                        for (j = 23;j--;){
                                ap->subm[i][j] = 0.0f;
                        }
                }
        }
        if(L == defDNA){
                RUN(set_subm_gaps_DNA(ap));
        }else if(L == defPROTEIN){

                RUN(set_subm_gaps(ap));
        }else if(L == redPROTEIN){
                RUN(set_subm_gaps(ap));
        }else if(L == ambigiousPROTEIN){
                RUN(new_aln_matrices(ap));
        }

        if(param->gpo != FLT_MAX){
                ap->gpo = param->gpo;
        }
        if(param->gpe != FLT_MAX){
                ap->gpe = param->gpe;
        }
        if(param->tgpe != FLT_MAX){
                ap->tgpe = param->tgpe;
        }
        if(param->matadd!= 0.0F){
                for(i = 0; i < 23;i++){
                        for(j = 0; j < 23;j++){
                                ap->subm[i][j] += param->matadd;
                        }
                }

        }

        *aln_param = ap;
        return OK;
ERROR:
        free_ap(ap);
        return FAIL;
}



void free_ap(struct aln_param* ap)
{
        int i;
        if(ap){
                if(ap->subm){
                        for (i = 23;i--;){
                                MFREE(ap->subm[i]);
                        }
                        MFREE(ap->subm);
                }
                if(ap->rng){
                        free_rng(ap->rng);
                }
                if(ap->tree){
                        MFREE(ap->tree);
                }
                MFREE(ap);
        }
}

/* These are old parameters from kalign 2 */
int set_subm_gaps_DNA(struct aln_param* ap)
{
        int i,j;
        for(i = 0; i < 5; i++){
                for(j =0; j < 5;j++){
                        ap->subm[i][j] = 283;
                }
        }
        //	A   91 -114  -31 -123    0  -43
        ap->subm[0][0] += 91;
        ap->subm[0][1] += -114;
        ap->subm[0][2] += -31;
        ap->subm[0][3] += -123;

//	C -114  100 -125  -31    0  -43
        ap->subm[1][0] += -114;
        ap->subm[1][1] += 100;
        ap->subm[1][2] += -125;
        ap->subm[1][3] += -31;

//	G  -31 -125  100 -114    0  -43
        ap->subm[2][0] += -31;
        ap->subm[2][1] += -125;
        ap->subm[2][2] += 100;
        ap->subm[2][3] += -114;

//	T -123  -31 -114   91    0  -43
        ap->subm[3][0] += -123;
        ap->subm[3][1] += -31;
        ap->subm[3][2] += -114;
        ap->subm[3][3] += 91;

        ap->gpo = 217;
        ap->gpe = 39.4;
        ap->tgpe =  292.6;
        //param->secret = 28.3;
        //         A    C    G    T    .    N

        /*for(i = 0; i < 5; i++){
                for(j =0; j < 5;j++){
                        if(i == j){
                                ap->subm[i][j] = 5;
                        }else{
                                ap->subm[i][j] = -4;
                        }
                }
        }

        ap->gpo = 7;
        ap->gpe = 4;
        ap->tgpe = 3;*/
        return OK;
}

int set_subm_gaps(struct aln_param* ap)
{
        int i,j;
        int m_pos = 0;
        float *matrix_pointer = 0;
        float balimt[]={
                24.501946,
                5.998169, 115.750240,
                -2.470710, -31.062287, 47.937530,
                0.999786, -29.101076, 28.000000, 36.003891,
                -22.005890, -7.007568, -44.750011, -38.000458, 71.000000,
                6.000000, -19.000000, 2.000000, -7.015625, -51.000000, 66.992218,
                -9.000000, -12.000000, 4.843778, 4.934356, -0.000031, -13.499763, 60.750057,
                -7.000855, -10.015595, -37.000214, -26.249912, 10.985351, -44.001923, -21.030732, 40.753445,
                -3.000214, -27.998062, 6.000000, 12.000229, -32.055085, -10.000061, 6.999969, -20.013794, 32.875029,
                -11.007813, -14.000000, -39.000000, -27.124605, 20.844236, -43.003876, -18.001831, 29.000000, -20.000458, 40.875059,
                -6.015106, -8.986221, -29.128878, -19.062470, 16.875029, -34.029297, -12.000946, 25.503868, -13.000000, 29.000000, 43.938384,
                -2.499519, -17.003632, 22.780331, 10.000000, -30.001923, 4.999786, 12.999542, -27.375036, 9.000000, -31.000000, -21.000000, 38.902403,
                3.999908, -30.249973, -6.060548, -4.000000, -37.003662, -15.000000, -10.029297, -25.246525, -5.001801, -22.015595, -23.124971, -8.500008, 77.000000,
                -1.000214, -23.499855, 9.999786, 17.000473, -25.014832, -9.000092, 12.624781, -18.148531, 15.877928, -15.031189, -9.015595, 7.999786, -1.062470, 27.000473,
                -5.001923, -21.078096, -2.124971, 5.000000, -31.750011, -9.000000, 7.000000, -23.030274, 27.999542, -21.492195, -16.001923, 3.757809, -8.000000, 15.500023, 47.984375,
                11.999054, 1.996338, 5.875120, 3.000000, -27.000000, 4.875029, -1.250919, -17.499977, 2.000000, -20.046876, -13.015564, 9.972198, 4.546899, 2.265614, -1.062013, 22.750027,
                6.993225, -4.031220, 1.000000, -0.499977, -21.000214, -10.000000, -2.062013, -5.000946, 1.985351, -12.999985, -5.000000, 6.000000, 1.562402, -0.500481, -1.000519, 15.960937, 25.986114,
                0.001923, 0.554681, -28.999985, -18.999557, 1.968780, -32.124025, -19.031220, 32.000000, -16.999985, 18.750027, 16.500053, -21.875227, -17.000458, -14.499519, -19.124971, -9.499886, 0.000015, 34.999512,
                -35.249973, -9.000031, -51.062959, -42.996109, 36.996124, -39.048310, -7.503426, -17.015595, -34.124971, -7.984436, -9.063233, -35.187503, -49.496101, -26.000214, -15.000092, -32.265599, -34.937026, -25.499977, 143.000000,
                -21.007782, -4.999985, -27.999985, -26.015595, 51.875029, -39.242649, 22.750027, -6.000458, -20.015595, 0.999969, -1.000000, -13.500008, -30.000000, -16.000458, -17.059052, -18.062470, -18.055146, -10.109377, 41.000107, 78.000961,
                0.750973, 0.621088, 1.000000, 0.750027, 0.999786, 0.937530, 0.937560, 0.984405, 0.999054, 0.991241, 1.000000, 0.871580, 0.999786, 0.031235, 1.000000, 0.265614, 0.097642, 0.969726, 0.999054, 1.000000, 0.999908,
        };
        ap->gpo = 55.918190;
        ap->gpe =  9.335495;
        ap->tgpe =  5.017874;

        matrix_pointer = balimt;
        m_pos = 0;

        for (i = 0;i < 21;i++){
                for (j = 0;j <= i;j++){
                        ap->subm[i][j] = matrix_pointer[m_pos];
                        ap->subm[j][i] = matrix_pointer[m_pos];
                        m_pos++;
                }
        }

        return OK;
}



int new_aln_matrices(struct aln_param* ap)
{
        int i;
        int j;
        //char aacode[20] = "ACDEFGHIKLMNPQRSTVWY";
        //char aa_order[23] = "ARNDCQEGHILKMFPSTWYVBZX";
        //int num_aa = 23;

        int CorBLOSUM66_13plus[23][23] = {
                /*A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X*/
                {5,-1,-1,-2,-2,-1,-1,0,-2,-1,-1,-1,0,-2,-1,1,0,-2,-2,0,-2,-1,0},
                {-1,6,0,-1,-3,1,1,-2,0,-2,-2,3,-1,-3,-1,-1,-1,-1,-1,-2,0,1,-1},
                {-1,0,6,2,-3,1,0,0,0,-3,-3,0,-2,-2,-1,1,0,-2,-1,-2,4,0,-1},
                {-2,-1,2,7,-3,1,2,-1,-1,-3,-3,0,-3,-3,-1,0,-1,-3,-2,-3,5,2,-1},
                {-2,-3,-3,-3,12,-3,-4,-3,-2,-2,-3,-3,-2,-1,-3,-2,-2,-3,-2,-2,-3,-3,-2},
                {-1,1,1,1,-3,5,2,-2,0,-2,-2,1,0,-2,-1,0,0,-1,-1,-2,1,3,0},
                {-1,1,0,2,-4,2,6,-2,-1,-3,-3,1,-2,-3,0,0,-1,-2,-2,-2,1,4,-1},
                {0,-2,0,-1,-3,-2,-2,7,-2,-4,-4,-2,-3,-3,-2,0,-2,-3,-3,-3,-1,-2,-1},
                {-2,0,0,-1,-2,0,-1,-2,10,-3,-3,0,-2,-2,-2,-1,-1,-2,1,-3,0,0,-1},
                {-1,-2,-3,-3,-2,-2,-3,-4,-3,5,2,-2,2,0,-3,-2,-1,-1,-1,3,-3,-2,-1},
                {-1,-2,-3,-3,-3,-2,-3,-4,-3,2,5,-2,3,1,-3,-3,-2,0,-1,1,-3,-2,-1},
                {-1,3,0,0,-3,1,1,-2,0,-2,-2,5,-1,-3,-1,0,0,-2,-2,-2,0,1,-1},
                {0,-1,-2,-3,-2,0,-2,-3,-2,2,3,-1,6,1,-2,-1,-1,0,-1,1,-2,-1,0},
                {-2,-3,-2,-3,-1,-2,-3,-3,-2,0,1,-3,1,7,-3,-2,-2,2,3,0,-3,-3,-1},
                {-1,-1,-1,-1,-3,-1,0,-2,-2,-3,-3,-1,-2,-3,9,0,-1,-2,-2,-2,-1,-1,-1},
                {1,-1,1,0,-2,0,0,0,-1,-2,-3,0,-1,-2,0,4,2,-2,-2,-1,0,0,0},
                {0,-1,0,-1,-2,0,-1,-2,-1,-1,-2,0,-1,-2,-1,2,5,-1,-1,0,0,0,0},
                {-2,-1,-2,-3,-3,-1,-2,-3,-2,-1,0,-2,0,2,-2,-2,-1,13,3,-2,-2,-2,-1},
                {-2,-1,-1,-2,-2,-1,-2,-3,1,-1,-1,-2,-1,3,-2,-2,-1,3,9,-1,-2,-2,-1},
                {0,-2,-2,-3,-2,-2,-2,-3,-3,3,1,-2,1,0,-2,-1,0,-2,-1,4,-3,-2,-1},
                {-2,0,4,5,-3,1,1,-1,0,-3,-3,0,-2,-3,-1,0,0,-2,-2,-3,4,1,-1},
                {-1,1,0,2,-3,3,4,-2,0,-2,-2,1,-1,-3,-1,0,0,-2,-2,-2,1,4,-1},
                {0,-1,-1,-1,-2,0,-1,-1,-1,-1,-1,-1,0,-1,-1,0,0,-1,-1,-1,-1,-1,-1},
        };

        for(i = 0; i < 23;i++){
                for(j = 0; j < 23;j++){
                        ap->subm[i][j] = (float)(CorBLOSUM66_13plus[i][j]);
                }
        }
        ap->gpo = 5.5F;
        ap->gpe = 2.0F;
        ap->tgpe = 1.0F;
        return OK;
}
