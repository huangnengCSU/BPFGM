//
// Created by user on 2020/10/31.
//

#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include "src/dataset.h"
#include <unistd.h>

struct CmdArgs
{
    char *ref2assm_bamfile = nullptr;
    char *trivialfeature_filename = nullptr;
    char *trivialtraining_outfile = nullptr;
};

CmdArgs getOpt(int argc, char *argv[])
{
    int opt;
    const char *optstring = "b:i:o:";
    CmdArgs args;
    while ((opt = getopt(argc, argv, optstring)) != -1)
    {
        if ((char)opt == 'b')
            args.ref2assm_bamfile = optarg;
        else if ((char)opt == 'i')
            args.trivialfeature_filename = optarg;
        else if ((char)opt == 'o')
        {
            args.trivialtraining_outfile = optarg;
        }
    }
    return args;
}

void usage()
{
    fprintf(stderr, "usage: trivial_training -b <BAM_FILE> -i <TRIVIAL_BLOCK_FEATURES> -o <TRIVIAL_BLOCK_FEATURES_TRAINING>\n");
    fprintf(stderr, "Generate the features for BlockPolish polishing using reads-to-assembly alignments in BAM_FILE.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE is the alignment of truth assembly to the assembly.\n");
    fprintf(stderr, "    TRIVIAL_BLOCK_FEATURES is the input file of the features of the trivial blocks.\n");
    fprintf(stderr, "    TRIVIAL_BLOCK_FEATURES_TRAINING is the output file of the features of the trivial blocks for training.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
    /*生成训练数据集*/
    CmdArgs args;
    args = getOpt(argc, argv);
    if (argc != 7)
    {
        usage();
        return 0;
    }
    std::string trivialfeatures = args.trivialfeature_filename;
    std::string trivialtrainingdata = args.trivialtraining_outfile;
    // ref2assembly
    std::string Ref2assembly_BamFile = args.ref2assm_bamfile;

    std::string RefFile = "";
    std::string AssemblyFile = "";
    std::map<std::string, int> IDrefMap;
    BamTools::BamReader reader1;
    if (!reader1.Open(Ref2assembly_BamFile))
    {
        std::cerr << "bamtools convert ERROR: could not open input BAM file(s)... Aborting." << std::endl;
        return false;
    }
    BamTools::BamAlignment al;
    BamTools::RefVector references;
    references = reader1.GetReferenceData();
    for (auto refData : references)
    {
        int refid = reader1.GetReferenceID(refData.RefName);
        IDrefMap.insert(std::pair<std::string, int>(refData.RefName, refid));
    }

    ConsensusTables consensusTables;
    std::vector<std::string> samSequence;
    std::vector<Cigar> samCigars;
    std::vector<std::string> samReadNames;
    std::map<int, std::string> refIDMap; // RefID->Contig_name
    bool status = CreateConsensusTableForReference(RefFile, AssemblyFile, Ref2assembly_BamFile, samSequence, samCigars,
                                                   samReadNames, consensusTables, refIDMap);
    std::ifstream ftrivial(trivialfeatures);
    std::ofstream trivialout(trivialtrainingdata);

    GenerateTrainingDataset(ftrivial, consensusTables, samSequence, IDrefMap, trivialout, 40);
    ftrivial.close();
    trivialout.close();
    return 0;
}
