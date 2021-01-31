/*
 * block multi-task polishing
 */
#include <iostream>
#include <time.h>
#include <assert.h>
#include "./src/consensus.h"
#include <unistd.h>
#include <unistd.h>

struct CmdArgs
{
    char *bamfile = nullptr;
    char *trivialout = nullptr;
    char *complexout = nullptr;
};

CmdArgs getOpt(int argc, char *argv[])
{
    int opt;
    const char *optstring = "b:s:c:";
    CmdArgs args;
    while ((opt = getopt(argc, argv, optstring)) != -1)
    {
        if ((char)opt == 'b')
            args.bamfile = optarg;
        else if ((char)opt == 's')
            args.trivialout = optarg;
        else if ((char)opt == 'c')
        {
            args.complexout = optarg;
        }
    }
    return args;
}

void usage()
{
    fprintf(stderr, "usage: block -b <BAM_FILE> -s <TRIVIAL_BLOCK_FEATURES> -c <COMPLEX_BLOCK_FEATURES>\n");
    fprintf(stderr, "Generate the features for BlockPolish polishing using reads-to-assembly alignments in BAM_FILE.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE is the alignment of reads to the assembly.\n");
    fprintf(stderr, "    TRIVIAL_BLOCK_FEATURES is the output file of the features of the trivial blocks.\n");
    fprintf(stderr, "    COMPLEX_BLOCK_FEATURES is the output file of the features of the trivial blocks.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
    clock_t startTime, endTime, trivialstartTime, trivialendTime, complexstartTime, complexendTime;
    startTime = clock();
    CmdArgs args;
    if (argc != 7)
    {
        usage();
        return 0;
    }
    args = getOpt(argc, argv);
    std::string BamFile = args.bamfile;
    std::string trivial_filename = args.trivialout;
    std::string complex_filename = args.complexout;
    std::string ReadFile = "";
    std::string RefFile = "";

    std::map<int, std::string> refIDMap; // RefID->Contig_name
    std::ofstream trivialfeatureWritter(trivial_filename);
    std::ofstream complexfeatureWritter(complex_filename);
    bool status = CreateConsensusTable(ReadFile, RefFile, BamFile, trivialfeatureWritter, complexfeatureWritter);
    trivialfeatureWritter.close();
    complexfeatureWritter.close();
    std::cout << "CreateConsensusTable finished." << std::endl;
    endTime = clock();
    std::cout << "running time:" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
    return 0;
}
