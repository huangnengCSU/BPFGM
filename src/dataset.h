//
// Created by user on 2020/11/18.
//

#ifndef BMP_DATASET_H
#define BMP_DATASET_H

#include "consensus.h"
#include <set>

struct FeatureRecord {
    std::string refname, feature_str;
    int region_start, region_end;
};

struct FeatureTargetRecord {
    std::string refname, feature_str, label_str;
    int region_start, region_end;
};

template<typename W>
void MultiThreadRun2(uint32_t thread_size, W work_func);

bool CreateConsensusTableForReference(const std::string ReadFile, const std::string RefFile, const std::string BamFile,
                                      std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                      std::vector<std::string> &SamReadNames,
                                      ConsensusTables &consensusTables,
                                      std::map<int, std::string> &refIDMap);

void GenerateTrainingDataset(std::ifstream &feature_reader, ConsensusTables &consensusTables,
                             std::vector<std::string> &SamSequence, std::map<std::string, int> &IDrefMap,
                             std::ofstream &fout, uint32_t thread_size);

void SplitString(std::string &s, std::vector<std::string> &v, std::string &c);

#endif //BMP_DATASET_H
