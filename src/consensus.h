//
// Created by lab on 2020/9/22.
//

#ifndef BMP_CONSENSUS_H
#define BMP_CONSENSUS_H

#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include "../bamtools/src/api/BamAux.h"
#include "../bamtools/src/api/BamReader.h"
#include "../bamtools/src/api/BamAlignment.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "../kalign/src/libkalign.h"
#ifdef __cplusplus
}
#endif

typedef struct cigar_summary_of_column {
    std::vector<int> MatchReads;
    std::vector<int> InsertReads;
    std::vector<int> DeleteReads;
} ColumnCigars;

typedef struct consensus_table_elements {
    int number_of_A = 0;
    int number_of_C = 0;
    int number_of_G = 0;
    int number_of_T = 0;
    int number_of_N = 0;
    int number_of_M = 0;
    int number_of_Insertion = 0;
    int number_of_Deletion = 0;
    ColumnCigars CigarVec; // SamSequencesIndex, cigarOp
    std::vector<std::pair<int, int>> InsertionSize;
    std::vector<std::pair<int, int>> SequencePosition;  // SamSequencesIndex, PosOnRead
//    std::vector<int> PositionOnSequence;
//    std::vector<int> SamSequencesIndex;
} TableElem;

typedef std::map<int, std::vector<TableElem>> ConsensusTables;
typedef std::vector<TableElem> SingleRefConsensusTable;
typedef std::pair<int, int> RegionInterval;
typedef std::vector<std::vector<float>> TrivialFeatureMatrix;
typedef std::vector<std::vector<float>> ComplexFeatureMatrix;
typedef std::vector<BamTools::CigarOp> Cigar;

class TrivialFeature {
public:
    TrivialFeatureMatrix mat;
    int refID;
    RegionInterval regionInterval;
};

class ComplexFeature {
public:
    ComplexFeatureMatrix mat;
    std::string RefName;
    RegionInterval regionInterval;
};


bool CreateConsensusTable(const std::string ReadFile, const std::string RefFile, const std::string BamFile,
                          std::ofstream &trivialfeatureWritter, std::ofstream &complexfeatureWritter);

void FilterConsensusTable(SingleRefConsensusTable &consensusTable, std::vector<std::string> &SamSequences,
                          std::vector<Cigar> &SamCigars, std::vector<std::string> &SamReadNames);

void DevideConsensusRegions(SingleRefConsensusTable &consensusTable, std::vector<RegionInterval> &trivial_regions,
                            std::vector<RegionInterval> &complex_regions);

void GetTrivialRegionFeatures(SingleRefConsensusTable &consensusTable, std::string RefName,
                              std::vector<RegionInterval> &trivial_regions, std::ofstream &fout);

void SolveAbnormalComplexRegion(std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                std::vector<std::string> &SamReadNames,
                                RegionInterval &ComplexRegion, std::map<int, int> &leftReadsMap,
                                std::map<int, int> &rightReadsMap, std::vector<std::string> &readBlock);

bool ProcessConsensusTable(SingleRefConsensusTable &table, std::string RefName, std::vector<std::string> &SamSequences,
                           std::vector<Cigar> &SamCigars, std::vector<std::string> &SamReadNames,
                           std::ofstream &trivialfeatureWritter,
                           std::ofstream &complexfeatureWritter);

template<typename W>
void MultiThreadRun(uint32_t thread_size, W work_func);

void GetComplexRegionFeaturesThreads(SingleRefConsensusTable &consensusTable, std::string &RefName,
                                     std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                     std::vector<std::string> &SamReadNames,
                                     std::vector<RegionInterval> &complex_regions, std::ofstream &fout,
                                     uint32_t thread_size);


#endif //BMP_CONSENSUS_H

