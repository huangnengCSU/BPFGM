#include <iostream>
#include <vector>

#ifdef __cplusplus
extern "C"
{
#endif
#include "libkalign.h"
#ifdef __cplusplus
}
#endif

using namespace std;

int main()
{
    vector<string> readBlock;
    readBlock.push_back("TGAGGAA");
    readBlock.push_back("TGAGAA");
    readBlock.push_back("TGAGAAT");
    readBlock.push_back("TGAGTAA");
    readBlock.push_back("TTAGAAA");
    readBlock.push_back("TGAAG");

    char **p;
    p = new char *[readBlock.size()];
    for (int i = 0; i < readBlock.size(); i++)
    {
        char *pc = new char[readBlock[i].size() + 1];
        strcpy(pc, readBlock[i].c_str());
        p[i] = pc;
    }
    struct msa *m;
    // for(int i=0;i<readBlock.size();i++){
    //     cout<<p[i]<<endl;
    // }
    m = read_input_from_array(p, readBlock.size());
    cout << m->numseq << endl;
    MSAOut out;
    out = RunMSA(m);
    cout<<out.num_seq<<","<<out.len<<endl;
    int i,j;
    for(i=0;i<out.num_seq;i++){
        j=0;
        while(out.out_ptr[i][j]!='\0'){
            cout<<out.out_ptr[i][j];
            j++;
        }
        cout<<"\n";
    }
    free_msa(m);
    free_MSAOut(&out);
    return 0;
}
