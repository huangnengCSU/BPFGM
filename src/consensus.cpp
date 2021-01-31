//
// Created by lab on 2020/9/22.
//

#include "consensus.h"

bool ProcessConsensusTable(SingleRefConsensusTable &table, std::string RefName, std::vector<std::string> &SamSequences,
                           std::vector<Cigar> &SamCigars, std::vector<std::string> &SamReadNames,
                           std::ofstream &trivialfeatureWritter,
                           std::ofstream &complexfeatureWritter)
{
    std::vector<RegionInterval> trivial_regions, complex_regions;
    // filter high coverage
    FilterConsensusTable(table, SamSequences, SamCigars, SamReadNames);
    DevideConsensusRegions(table, trivial_regions, complex_regions);
    GetTrivialRegionFeatures(table, RefName, trivial_regions, trivialfeatureWritter);
    GetComplexRegionFeaturesThreads(table, RefName, SamSequences, SamCigars, SamReadNames, complex_regions,
                                    complexfeatureWritter, 40);
}

bool CreateConsensusTable(const std::string ReadFile, const std::string RefFile, const std::string BamFile,
                          std::ofstream &trivialfeatureWritter,
                          std::ofstream &complexfeatureWritter)
{
    BamTools::BamReader reader1;
    if (!reader1.Open(BamFile))
    {
        std::cerr << "bamtools convert ERROR: could not open input BAM file(s)... Aborting." << std::endl;
        return false;
    }
    BamTools::BamAlignment al;
    BamTools::RefVector references;
    references = reader1.GetReferenceData();
    std::map<int, std::string> refIDMap;
    std::map<int, int> refLenMap;
    for (auto refData : references)
    {
        int refid = reader1.GetReferenceID(refData.RefName);
        refIDMap.insert(std::pair<int, std::string>(refid, refData.RefName));
        refLenMap.insert(std::pair<int, int>(refid, refData.RefLength));
    }
    int SamSequencesIndex = -1;
    std::vector<std::string> SamSequences;
    std::vector<Cigar> SamCigars;
    std::vector<std::string> SamReadNames;
    SingleRefConsensusTable table;
    int cur_refid = -1;
    while (reader1.GetNextAlignment(al))
    {
        std::string Name = al.Name; // read name
        uint32_t AlignmentFlag = al.AlignmentFlag;
        uint16_t MappingQuality = al.MapQuality;
        int RefID = al.RefID;
        std::string Sequence = al.QueryBases;
        int32_t ReadLength = al.Length;
        int32_t StartPositionOnRef = al.Position; // position (0-based) where alignment starts (on the reference)
        std::string ReadQualities = al.Qualities;
        Cigar CigarData = al.CigarData;
        if ((AlignmentFlag & 0x900) != 0)
            // filter reads alignment flagged as secondary or supplementary
            continue;
        if (MappingQuality < 30)
            // mapping quality score not less than 30
            continue;
        if (RefID == -1)
            continue;
        if (Sequence.size() == 0)
            continue;
        if (cur_refid == -1)
        {
            SingleRefConsensusTable newtable(refLenMap.at(RefID));
            table.clear();
            table = newtable;
            cur_refid = RefID;
        }
        if (cur_refid != RefID)
        {
            // 读取到新的染色体的比对信息，把当前染色体的table处理，提取简单区域和复杂区域特征
            ProcessConsensusTable(table, refIDMap.at(cur_refid), SamSequences, SamCigars, SamReadNames,
                                  trivialfeatureWritter, complexfeatureWritter);
            SingleRefConsensusTable newtable(refLenMap.at(RefID));
            table.clear();
            table = newtable;
            cur_refid = RefID;
        }
        SamSequences.push_back(Sequence);
        SamCigars.push_back(CigarData);
        SamReadNames.push_back(Name);
        SamSequencesIndex++;
        if (CigarData[0].Type == 'S')
        {
            // soft clip.
            int clip_size = CigarData[0].Length;
            int PosOnTable = StartPositionOnRef;
            int PosOnRead = clip_size;
            std::vector<TableElem>::iterator te = table.begin() + PosOnTable;
            for (int i = 1; i < CigarData.size(); i++)
            {
                BamTools::CigarOp cop = CigarData[i];
                if (cop.Type == 'H' or cop.Type == 'S')
                    continue; //末尾的clip不计算
                if (cop.Type == 'M')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.MatchReads.push_back(SamSequencesIndex);
                        te->number_of_M++;
                        switch (Sequence[PosOnRead])
                        {
                        case 'A':
                            te->number_of_A++;
                            break;
                        case 'C':
                            te->number_of_C++;
                            break;
                        case 'G':
                            te->number_of_G++;
                            break;
                        case 'T':
                            te->number_of_T++;
                            break;
                        case 'N':
                            te->number_of_N++;
                            break;
                        }
                        //                        te->ReadNames.push_back(Name);
                        //                        te->PositionOnSequence.push_back(PosOnRead);
                        //                        te->SamSequencesIndex.push_back(SamSequencesIndex);
                        te->SequencePosition.push_back(std::pair<int, int>(SamSequencesIndex, PosOnRead));
                        assert(SamSequences[SamSequencesIndex][PosOnRead] == Sequence[PosOnRead]);
                        PosOnRead++;
                        te++;
                    }
                }
                else if (cop.Type == 'I')
                {
                    // (te - 1)->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                    (te - 1)->CigarVec.InsertReads.push_back(SamSequencesIndex);
                    (te - 1)->number_of_Insertion++;
                    (te - 1)->InsertionSize.push_back(std::make_pair(SamSequencesIndex, cop.Length));
                    PosOnRead += cop.Length;
                }
                else if (cop.Type == 'D')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.DeleteReads.push_back(SamSequencesIndex);
                        te->number_of_Deletion++;
                        te++;
                    }
                }
            }
        }
        else if (CigarData[0].Type == 'H')
        {
            // hard clip.
            int clip_size = CigarData[0].Length;
            int PosOnTable = StartPositionOnRef;
            int PosOnRead = 0;
            std::vector<TableElem>::iterator te = table.begin() + PosOnTable;
            for (int i = 1; i < CigarData.size(); i++)
            {
                BamTools::CigarOp cop = CigarData[i];
                if (cop.Type == 'H' or cop.Type == 'S')
                    continue; //末尾的clip不计算
                if (cop.Type == 'M')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.MatchReads.push_back(SamSequencesIndex);
                        te->number_of_M++;
                        switch (Sequence[PosOnRead])
                        {
                        case 'A':
                            te->number_of_A++;
                            break;
                        case 'C':
                            te->number_of_C++;
                            break;
                        case 'G':
                            te->number_of_G++;
                            break;
                        case 'T':
                            te->number_of_T++;
                            break;
                        case 'N':
                            te->number_of_N++;
                            break;
                        }
                        //                        te->ReadNames.push_back(Name);
                        //                        te->PositionOnSequence.push_back(PosOnRead);
                        //                        te->SamSequencesIndex.push_back(SamSequencesIndex);
                        te->SequencePosition.push_back(std::pair<int, int>(SamSequencesIndex, PosOnRead));
                        assert(SamSequences[SamSequencesIndex][PosOnRead] == Sequence[PosOnRead]);
                        PosOnRead++;
                        te++;
                    }
                }
                else if (cop.Type == 'I')
                {
                    // (te - 1)->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                    (te - 1)->CigarVec.InsertReads.push_back(SamSequencesIndex);
                    (te - 1)->number_of_Insertion++;
                    (te - 1)->InsertionSize.push_back(std::make_pair(SamSequencesIndex, cop.Length));
                    PosOnRead += cop.Length;
                }
                else if (cop.Type == 'D')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.DeleteReads.push_back(SamSequencesIndex);
                        te->number_of_Deletion++;
                        te++;
                    }
                }
            }
        }
        else
        {
            // without clip.
            int clip_size = CigarData[0].Length;
            int PosOnTable = StartPositionOnRef;
            int PosOnRead = 0;
            std::vector<TableElem>::iterator te = table.begin() + PosOnTable;
            for (int i = 0; i < CigarData.size(); i++)
            { //cigar起点无clip，则从第0个flag开始计算
                BamTools::CigarOp cop = CigarData[i];
                if (cop.Type == 'H' or cop.Type == 'S')
                    continue; //末尾的clip不计算
                if (cop.Type == 'M')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.MatchReads.push_back(SamSequencesIndex);
                        te->number_of_M++;
                        switch (Sequence[PosOnRead])
                        {
                        case 'A':
                            te->number_of_A++;
                            break;
                        case 'C':
                            te->number_of_C++;
                            break;
                        case 'G':
                            te->number_of_G++;
                            break;
                        case 'T':
                            te->number_of_T++;
                            break;
                        case 'N':
                            te->number_of_N++;
                            break;
                        }
                        //                        te->ReadNames.push_back(Name);
                        //                        te->PositionOnSequence.push_back(PosOnRead);
                        //                        te->SamSequencesIndex.push_back(SamSequencesIndex);
                        te->SequencePosition.push_back(std::pair<int, int>(SamSequencesIndex, PosOnRead));
                        assert(SamSequences[SamSequencesIndex][PosOnRead] == Sequence[PosOnRead]);
                        PosOnRead++;
                        te++;
                    }
                }
                else if (cop.Type == 'I')
                {
                    // (te - 1)->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                    (te - 1)->CigarVec.InsertReads.push_back(SamSequencesIndex);
                    (te - 1)->number_of_Insertion++;
                    (te - 1)->InsertionSize.push_back(std::make_pair(SamSequencesIndex, cop.Length));
                    PosOnRead += cop.Length;
                }
                else if (cop.Type == 'D')
                {
                    for (int j = 0; j < cop.Length; j++)
                    {
                        // te->CigarVec.push_back(std::make_pair(SamSequencesIndex, cop));
                        te->CigarVec.DeleteReads.push_back(SamSequencesIndex);
                        te->number_of_Deletion++;
                        te++;
                    }
                }
            }
        }
    }
    ProcessConsensusTable(table, refIDMap.at(cur_refid), SamSequences, SamCigars, SamReadNames,
                          trivialfeatureWritter, complexfeatureWritter);
    return true;
}

void DevideConsensusRegions(SingleRefConsensusTable &consensusTable, std::vector<RegionInterval> &trivial_regions,
                            std::vector<RegionInterval> &complex_regions)
{
    // TODO:简单区域、复杂区域的划分不是很严格
    SingleRefConsensusTable &table = consensusTable;
    int consecutive_number = 0;
    bool find_trivial = true;
    int watcher = 0;
    for (int i = 0, size = table.size(); i < size; i++)
    {
        if (find_trivial)
        {
            if (table[i].number_of_Insertion < 6 and
                ((float)(table[i].number_of_M) / (table[i].number_of_M + table[i].number_of_Deletion) > 0.8 or
                 (float)(table[i].number_of_Deletion) / (table[i].number_of_M + table[i].number_of_Deletion) >
                     0.8))
            {
                consecutive_number++;
            }
            else
            {
                if (i == 0)
                {
                    //如果第0个位置就是complex的话，直接将标志转为find complex
                    find_trivial = false;
                    consecutive_number = 0;
                    continue;
                }
                trivial_regions.push_back(RegionInterval(watcher, i - 1));
                consecutive_number = 0;
                watcher = i;
                find_trivial = false;
            }
        }
        else
        {
            if (table[i].number_of_Insertion < 6 and
                ((float)(table[i].number_of_M) / (table[i].number_of_M + table[i].number_of_Deletion) > 0.8 or
                 (float)(table[i].number_of_Deletion) / (table[i].number_of_M + table[i].number_of_Deletion) >
                     0.8))
            {
                consecutive_number++;
                if (consecutive_number >= 3)
                {
                    complex_regions.push_back(RegionInterval(watcher, i - 1));
                    consecutive_number = 0;
                    watcher = i;
                    find_trivial = true;
                }
            }
            else
            {
                consecutive_number = 0;
            }
        }
    }
    if (find_trivial)
    {
        trivial_regions.push_back(RegionInterval(watcher, table.size() - 1));
    }
    else
    {
        complex_regions.push_back(RegionInterval(watcher, table.size() - 1));
    }
}

void GetTrivialRegionFeatures(SingleRefConsensusTable &consensusTable, std::string RefName,
                              std::vector<RegionInterval> &trivial_regions, std::ofstream &fout)
{
    std::vector<TrivialFeature> TrivialFeatureVector;
    for (auto interval : trivial_regions)
    {
        TrivialFeature feature;
        feature.regionInterval = interval;
        TrivialFeatureMatrix array(interval.second - interval.first + 1);
        for (int i = interval.first, j = 0; i <= interval.second; i++, j++)
        {
            TableElem elem = consensusTable[i];
            int cov = elem.number_of_A + elem.number_of_C + elem.number_of_G + elem.number_of_T + elem.number_of_N +
                      elem.number_of_Insertion + elem.number_of_Deletion;
            std::vector<float> array1 = {(float)(elem.number_of_A) / cov,
                                         (float)(elem.number_of_C) / cov,
                                         (float)(elem.number_of_G) / cov,
                                         (float)(elem.number_of_T) / cov,
                                         (float)(elem.number_of_N) / cov,
                                         (float)(elem.number_of_Insertion) / cov,
                                         (float)(elem.number_of_Deletion) / cov};
            array[j] = array1;
        }
        feature.mat = array;
        TrivialFeatureVector.push_back(feature);
    }
    for (auto trivialfeature : TrivialFeatureVector)
    {
        fout << '>' << RefName << "\t" << trivialfeature.regionInterval.first << "\t"
             << trivialfeature.regionInterval.second << "\n";
        for (auto trivialarray : trivialfeature.mat)
        {
            for (auto ratio : trivialarray)
            {
                fout << ratio << ",";
            }
        }
        fout << "\n";
    }
}

template <typename W>
void MultiThreadRun(uint32_t thread_size, W work_func)
{
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < thread_size; ++i)
    {
        workers.push_back(std::thread(work_func, i));
    }

    for (uint32_t i = 0; i < workers.size(); ++i)
    {
        workers[i].join();
    }
}

void GetComplexRegionFeaturesThreads(SingleRefConsensusTable &consensusTable, std::string &RefName,
                                     std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                     std::vector<std::string> &SamReadNames,
                                     std::vector<RegionInterval> &complex_regions, std::ofstream &fout,
                                     uint32_t thread_size)
{
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    std::vector<RegionInterval>::iterator RegIntervalVectorIt = complex_regions.begin();

    auto generate_func = [&mutex_gen, &complex_regions, &RegIntervalVectorIt]() {
        std::lock_guard<std::mutex> lock(mutex_gen);
        RegionInterval tmp_itvl;
        if (RegIntervalVectorIt != complex_regions.end())
        {
            tmp_itvl = *RegIntervalVectorIt;
            RegIntervalVectorIt++;
            return std::pair<bool, RegionInterval>(true, tmp_itvl);
        }
        else
        {
            return std::pair<bool, RegionInterval>(false, tmp_itvl);
        }
    };

    auto combine_func = [&mutex_comb](std::vector<ComplexFeature> &ComplexFeatureVector, std::ofstream &fout) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (auto complexMSAOut : ComplexFeatureVector)
        {
            fout << '>' << complexMSAOut.RefName << "\t" << complexMSAOut.regionInterval.first << "\t"
                 << complexMSAOut.regionInterval.second << "\n";
            for (auto array : complexMSAOut.mat)
            {
                for (auto ratio : array)
                {
                    fout << ratio << ",";
                }
            }
            fout << "\n";
        }
    };
    auto work_func = [generate_func, combine_func, &consensusTable, &RefName, &SamSequences, &SamCigars, &SamReadNames, &complex_regions, &fout](
                         size_t) {
        std::vector<ComplexFeature> ComplexFeatureVector;
        while (true)
        {
            std::pair<bool, RegionInterval> product = generate_func();
            if (!product.first)
            {
                break;
            }
            else
            {
                std::vector<std::string> readBlock;
                int contigSize = consensusTable.size();
                ComplexFeature feature;
                RegionInterval interval = product.second;
                feature.RefName = RefName;
                feature.regionInterval = product.second;
                assert(interval.first >= 0);
                assert(interval.second < contigSize);
                int intervalLeftValue = (interval.first >= 0 ? interval.first : 0);
                int intervalRightValue = (interval.second <= contigSize - 1 ? interval.second : contigSize - 1);
                std::map<int, int> leftReadsMap, rightReadsMap;
                for (auto vSeqPos : consensusTable[intervalLeftValue].SequencePosition)
                { //SamSequencesIndex, PosOnRead
                    leftReadsMap.insert(std::pair<int, int>(vSeqPos.first, vSeqPos.second));
                }
                for (auto vSeqPos : consensusTable[intervalRightValue].SequencePosition)
                {
                    rightReadsMap.insert(std::pair<int, int>(vSeqPos.first, vSeqPos.second));
                }
                for (auto mSeqPos : leftReadsMap)
                {
                    if (rightReadsMap.count(mSeqPos.first) == 0)
                    {
                        bool found = false;
                        for (int tidx = intervalRightValue - 1; tidx >= intervalLeftValue + 1; tidx--)
                        {
                            for (auto vSeqPos : consensusTable[tidx].SequencePosition)
                            {
                                if (vSeqPos.first == mSeqPos.first)
                                {
                                    rightReadsMap.insert(std::pair<int, int>(vSeqPos.first, vSeqPos.second));
                                    found = true;
                                    break;
                                }
                            }
                            if (found)
                            {
                                break;
                            }
                        }
                    }
                }
                for (auto mSeqPos : rightReadsMap)
                {
                    if (leftReadsMap.count(mSeqPos.first) == 0)
                    {
                        bool found = false;
                        for (int tidx = intervalLeftValue + 1; tidx <= intervalRightValue - 1; tidx++)
                        {
                            for (auto vSeqPos : consensusTable[tidx].SequencePosition)
                            {
                                if (vSeqPos.first == mSeqPos.first)
                                {
                                    leftReadsMap.insert(std::pair<int, int>(vSeqPos.first, vSeqPos.second));
                                    found = true;
                                    break;
                                }
                            }
                            if (found)
                            {
                                break;
                            }
                        }
                    }
                }
                int segmentLength;
                for (auto mSeqPos : leftReadsMap)
                { //SamSequencesIndex, PosOnRead
                    if (rightReadsMap.count(mSeqPos.first) == 0)
                        continue;
                    segmentLength = rightReadsMap.at(mSeqPos.first) - mSeqPos.second + 1;
                    std::string readSegment = SamSequences[mSeqPos.first].substr(mSeqPos.second, segmentLength);
                    readBlock.push_back(readSegment);
                }
                if (readBlock.size() <= 1)
                {
                    continue;
                }
                char **p;
                p = new char *[readBlock.size()];
                for (int i = 0; i < readBlock.size(); i++)
                {
                    char *pc = new char[readBlock[i].size() + 1];
                    strcpy(pc, readBlock[i].c_str());
                    p[i] = pc;
                }

                //                printf("%s:%d-%d,%ld\n", refIDMap.at(product.second.first).c_str(), product.second.second.first,
                //                       product.second.second.second, readBlock.size());
                // TODO: 复杂区域长度超过1k 或者 复杂区域覆盖度超过1k
                if (readBlock.size() > 1000 or (product.second.second - product.second.first) > 1000)
                {
                    // 复杂区域长度超过1k 或者 复杂区域覆盖度超过1k
                    printf("%s:%d-%d,%ld\n", RefName.c_str(), product.second.first, product.second.second,
                           readBlock.size());
                    SolveAbnormalComplexRegion(SamSequences, SamCigars, SamReadNames, interval, leftReadsMap,
                                               rightReadsMap, readBlock);

                    printf("(filtered)%s:%d-%d,%ld\n", RefName.c_str(), product.second.first, product.second.second,
                           readBlock.size());
                    //                    continue;
                }

                // use kalign to perform multiple seqience alignment
                struct msa *m;
                m = read_input_from_array(p, readBlock.size());
                MSAOut out;
                out = RunMSA(m);
                ComplexFeatureMatrix array(out.len - 1); // MSA比对结果字符串以\0结尾，所以len要大1.
                int i, j;
                int cov;
                std::map<char, int> baseNumber;
                for (i = 0; i < out.len - 1; i++)
                { // out.len-1是因为\0结尾，不处理\0
                    baseNumber.clear();
                    baseNumber.insert(std::pair<char, int>('A', 0));
                    baseNumber.insert(std::pair<char, int>('C', 0));
                    baseNumber.insert(std::pair<char, int>('G', 0));
                    baseNumber.insert(std::pair<char, int>('T', 0));
                    baseNumber.insert(std::pair<char, int>('N', 0));
                    baseNumber.insert(std::pair<char, int>('I', 0));
                    baseNumber.insert(std::pair<char, int>('-', 0));
                    for (j = 0; j < out.num_seq; j++)
                    {
                        if ((int)out.out_ptr[j][i] >= 97 && (int)out.out_ptr[j][i] <= 122)
                            out.out_ptr[j][i] -= 32;
                        baseNumber.at(out.out_ptr[j][i])++;
                    }
                    cov = baseNumber.at('A') + baseNumber.at('C') + baseNumber.at('G') + baseNumber.at('T') +
                          baseNumber.at('N') + baseNumber.at('I') + baseNumber.at('-');
                    std::vector<float> array1 = {(float)(baseNumber.at('A')) / cov,
                                                 (float)(baseNumber.at('C')) / cov,
                                                 (float)(baseNumber.at('G')) / cov,
                                                 (float)(baseNumber.at('T')) / cov,
                                                 (float)(baseNumber.at('N')) / cov,
                                                 (float)(baseNumber.at('I')) / cov,
                                                 (float)(baseNumber.at('-')) / cov};
                    array[i] = array1;
                }
                feature.mat = array;
                ComplexFeatureVector.push_back(feature);
                // 释放内存
                free_msa(m);
                free_MSAOut(&out);
                for (int i = 0; i < readBlock.size(); i++)
                {
                    delete[] p[i];
                }
                delete[] p;
                if (ComplexFeatureVector.size() > 10000)
                {
                    combine_func(ComplexFeatureVector, fout);
                    ComplexFeatureVector.clear();
                }
            }
        }
        if (ComplexFeatureVector.size() > 0)
        {
            combine_func(ComplexFeatureVector, fout);
            ComplexFeatureVector.clear();
        }
    };
    MultiThreadRun(thread_size, work_func);
}

bool cmp(std::pair<int, float> a, std::pair<int, float> b)
{
    /*
     * 按identuty从大到小排序
     */
    return a.second > b.second;
}

void SortCoverCigarIdentity(std::map<int, Cigar> &coverCigars, std::vector<std::string> &SamReadNames,
                            std::vector<int> &sortedIndex)
{
    /*
     * 计算identity时不区分match和mismatch
     */
    std::vector<std::pair<int, float>> identities; //identities: index, identity
    std::map<int, Cigar>::iterator it = coverCigars.begin();
    for (; it != coverCigars.end(); it++)
    {
        int index = it->first;
        Cigar cigar = it->second;
        BamTools::CigarOp cop;
        int overhead_size = 0, M_size = 0, I_size = 0, D_size = 0, read_size = 0;
        for (int i = 0; i < cigar.size(); i++)
        {
            cop = cigar[i];
            if (cop.Type == 'S' or cop.Type == 'H')
            {
                overhead_size += cop.Length;
                read_size += cop.Length;
            }
            else if (cop.Type == 'M')
            {
                M_size += cop.Length;
                read_size += cop.Length;
            }
            else if (cop.Type == 'I')
            {
                I_size += cop.Length;
                read_size += cop.Length;
            }
            else if (cop.Type == 'D')
                D_size += cop.Length;
            else
            {
                printf("Unknown flag:%c\n", cop.Type);
            }
        }
        float overhead_score = 1 - ((float)(overhead_size) / read_size);
        float accuracy_score = (float)(M_size) / (M_size + I_size + D_size);
        float mapping_length_score = (M_size + D_size) * accuracy_score;

        identities.push_back(std::make_pair(index, mapping_length_score));
    }
    std::sort(identities.begin(), identities.end(), cmp);
    for (auto idt : identities)
    {
        //        printf("%d(%s):%f\n", idt.first, SamReadNames[idt.first].c_str(), idt.second);
        sortedIndex.push_back(idt.first);
    }
}

void SolveAbnormalComplexRegion(std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                std::vector<std::string> &SamReadNames,
                                RegionInterval &ComplexRegion, std::map<int, int> &leftReadsMap,
                                std::map<int, int> &rightReadsMap, std::vector<std::string> &readBlock)
{
    readBlock.clear(); // 清除read内容，后面存入过滤后的read
    std::vector<int> indexVec;
    for (auto mSeqPos : leftReadsMap)
    {
        /*SamSequencesIndex, PosOnRead*/
        if (rightReadsMap.count(mSeqPos.first) == 0)
            continue;
        indexVec.push_back(mSeqPos.first);
    }
    /*
     * coverCigars: index, Cigar
     * coverSequences: index, Sequence
     * */
    std::map<int, Cigar> coverCigars;
    std::map<int, std::string> coverSequences;
    std::vector<int> sortedIndex;
    std::set<int> chosenIndex;
    for (int i = 0, vsize = indexVec.size(); i < vsize; i++)
    {
        coverCigars[indexVec[i]] = SamCigars[indexVec[i]];
        coverSequences[indexVec[i]] = SamSequences[indexVec[i]];
    }
    SortCoverCigarIdentity(coverCigars, SamReadNames, sortedIndex);
    int sampleNumber = 30 < sortedIndex.size() ? 30 : sortedIndex.size();
    for (int i = 0; i < sampleNumber; i++)
    {
        chosenIndex.insert(sortedIndex[i]);
    }
    for (auto mSeqPos : leftReadsMap)
    {
        /*SamSequencesIndex, PosOnRead*/
        if (rightReadsMap.count(mSeqPos.first) == 0)
            continue;
        if (chosenIndex.find(mSeqPos.first) != chosenIndex.end())
        {
            int segmentLength = rightReadsMap.at(mSeqPos.first) - mSeqPos.second + 1;
            std::string readSegment = SamSequences[mSeqPos.first].substr(mSeqPos.second, segmentLength);
            readBlock.push_back(readSegment);
        }
    }
}

void ChooseAlignmentIndex(SingleRefConsensusTable &consensusTable, RegionInterval &interval,
                          std::vector<std::string> &SamSequences,
                          std::vector<Cigar> &SamCigars, std::vector<std::string> &SamReadNames,
                          std::set<int> &chosenIndex)
{
    SingleRefConsensusTable &table = consensusTable;
    std::map<int, int> indexMap;
    for (int i = interval.first; i < interval.second; i++)
    {
        for (auto SeqPos : table[i].SequencePosition)
        {
            if (indexMap.find(SeqPos.first) != indexMap.end())
                indexMap[SeqPos.first]++;
            else
                indexMap.insert(std::make_pair(SeqPos.first, 1));
        }
    }
    std::vector<int> indexVec;
    for (auto indexitem : indexMap)
    {
        indexVec.push_back(indexitem.first);
    }
    std::map<int, Cigar> coverCigars;
    std::map<int, std::string> coverSequences;
    std::vector<int> sortedIndex;
    for (int i = 0, vsize = indexVec.size(); i < vsize; i++)
    {
        coverCigars[indexVec[i]] = SamCigars[indexVec[i]];
        coverSequences[indexVec[i]] = SamSequences[indexVec[i]];
    }
    SortCoverCigarIdentity(coverCigars, SamReadNames, sortedIndex);
    int sampleNumber = 30 < sortedIndex.size() ? 30 : sortedIndex.size();
    for (int i = 0; i < sampleNumber; i++)
    {
        chosenIndex.insert(sortedIndex[i]);
    }
}

void ModifyTableItem(SingleRefConsensusTable &consensusTable, std::vector<std::string> &SamSequences,
                     RegionInterval &interval, std::set<int> &chosenIndexSet)
{
    SingleRefConsensusTable &table = consensusTable;
    for (int i = interval.first; i < interval.second; i++)
    {
        int number_of_A = 0;
        int number_of_C = 0;
        int number_of_G = 0;
        int number_of_T = 0;
        int number_of_N = 0;
        int number_of_M = 0;
        int number_of_Insertion = 0;
        int number_of_Deletion = 0;
        ColumnCigars filteredCigarVec;
        std::vector<std::pair<int, int>> filteredInsertionSize;
        std::vector<std::pair<int, int>> filteredSequencePosition;
        for (auto SeqPos : table[i].SequencePosition)
        {
            if (chosenIndexSet.find(SeqPos.first) != chosenIndexSet.end())
            {
                filteredSequencePosition.push_back(std::make_pair(SeqPos.first, SeqPos.second));
                switch (SamSequences[SeqPos.first][SeqPos.second])
                {
                case 'A':
                    number_of_A++;
                    break;
                case 'C':
                    number_of_C++;
                    break;
                case 'G':
                    number_of_G++;
                    break;
                case 'T':
                    number_of_T++;
                    break;
                case 'N':
                    number_of_N++;
                    break;
                }
            }
        }
        for (auto cid : table[i].CigarVec.MatchReads)
        {
            if (chosenIndexSet.find(cid) != chosenIndexSet.end())
            {
                number_of_M++;
                filteredCigarVec.MatchReads.push_back(cid);
            }
        }
        for (auto cid : table[i].CigarVec.InsertReads)
        {
            if (chosenIndexSet.find(cid) != chosenIndexSet.end())
            {
                number_of_Insertion++;
                filteredCigarVec.InsertReads.push_back(cid);
            }
        }
        for (auto cid : table[i].CigarVec.DeleteReads)
        {
            if (chosenIndexSet.find(cid) != chosenIndexSet.end())
            {
                number_of_Deletion++;
                filteredCigarVec.DeleteReads.push_back(cid);
            }
        }
        // for (auto idxCigar : table[i].CigarVec)
        // {
        //     if (chosenIndexSet.find(idxCigar.first) != chosenIndexSet.end())
        //     {
        //         filteredCigarVec.push_back(std::make_pair(idxCigar.first, idxCigar.second));
        //         if (idxCigar.second.Type == 'M')
        //             number_of_M++;
        //         else if (idxCigar.second.Type == 'I')
        //             number_of_Insertion++;
        //         else if (idxCigar.second.Type == 'D')
        //             number_of_Deletion++;
        //         else
        //             continue;
        //     }
        // }
        for (auto idxInsert : table[i].InsertionSize)
        {
            if (chosenIndexSet.find(idxInsert.first) != chosenIndexSet.end())
            {
                filteredInsertionSize.push_back(std::make_pair(idxInsert.first, idxInsert.second));
            }
        }
        table[i].number_of_A = number_of_A;
        table[i].number_of_C = number_of_C;
        table[i].number_of_G = number_of_G;
        table[i].number_of_T = number_of_T;
        table[i].number_of_N = number_of_N;
        table[i].number_of_M = number_of_M;
        table[i].number_of_Insertion = number_of_Insertion;
        table[i].number_of_Deletion = number_of_Deletion;
        table[i].CigarVec = filteredCigarVec;
        table[i].InsertionSize = filteredInsertionSize;
        table[i].SequencePosition = filteredSequencePosition;
    }
}

void FilterConsensusTable(SingleRefConsensusTable &consensusTable, std::vector<std::string> &SamSequences,
                          std::vector<Cigar> &SamCigars, std::vector<std::string> &SamReadNames)
{
    /*
     * 固定滑动窗口大小，对覆盖度超高的窗口，过滤比对结果，从而降低覆盖度
     */
    int shiftWindowSize = 1000;
    SingleRefConsensusTable &table = consensusTable;
    int contigSize = table.size();
    int windowStartPos = 0;
    bool break_flag = false;
    while (true)
    {
        RegionInterval interval;
        if (windowStartPos + 2 * shiftWindowSize - 1 < contigSize)
        {
            interval.first = windowStartPos;
            interval.second = windowStartPos + shiftWindowSize - 1;
        }
        else
        {
            interval.first = windowStartPos;
            interval.second = contigSize - 1;
            break_flag = true;
        }
        int totalCoverage = 0;
        for (int i = interval.first; i < interval.second; i++)
        {
            totalCoverage += table[i].SequencePosition.size();
        }
        float averageCoverage = (float)(totalCoverage) / (interval.second - interval.first + 1);
        if (averageCoverage < 50)
        {
            if (break_flag)
                break;
            windowStartPos += shiftWindowSize;
            continue;
        }
        //            printf("(modify)%s:%d-%d,%f\n", RefIDMap.at(RefID).c_str(), windowStartPos,
        //                   windowStartPos + shiftWindowSize - 1, averageCoverage);
        std::set<int> chosenIndexSet;
        ChooseAlignmentIndex(table, interval, SamSequences, SamCigars, SamReadNames, chosenIndexSet);
        ModifyTableItem(table, SamSequences, interval, chosenIndexSet);
        windowStartPos += shiftWindowSize;
        if (break_flag)
        {
            break;
        }
    }
}
