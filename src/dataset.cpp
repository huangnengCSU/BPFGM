//
// Created by user on 2020/11/18.
//

#include "dataset.h"

void SplitString(std::string &s, std::vector<std::string> &v, std::string &c)
{
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2 - pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

bool CreateConsensusTableForReference(const std::string ReadFile, const std::string RefFile, const std::string BamFile,
                                      std::vector<std::string> &SamSequences, std::vector<Cigar> &SamCigars,
                                      std::vector<std::string> &SamReadNames,
                                      ConsensusTables &consensusTables,
                                      std::map<int, std::string> &refIDMap)
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
    for (auto refData : references)
    {
        int refid = reader1.GetReferenceID(refData.RefName);
        refIDMap.insert(std::pair<int, std::string>(refid, refData.RefName));
    }
    int SamSequencesIndex = -1;
    //    ConsensusTables consensusTables;
    for (BamTools::RefVector::iterator it = references.begin(); it != references.end(); it++)
    {
        std::vector<TableElem> table(it->RefLength);
        consensusTables.insert(std::pair<int, std::vector<TableElem>>(reader1.GetReferenceID(it->RefName), table));
    }
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
        if (RefID == -1)
            continue;
        if (Sequence.size() == 0)
            continue;
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
            std::vector<TableElem>::iterator te = consensusTables.at(RefID).begin() + PosOnTable;
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
            std::vector<TableElem>::iterator te = consensusTables.at(RefID).begin() + PosOnTable;
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
            std::vector<TableElem>::iterator te = consensusTables.at(RefID).begin() + PosOnTable;
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
    return true;
}

template <typename W>
void MultiThreadRun2(uint32_t thread_size, W work_func)
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

void GenerateTrainingDataset(std::ifstream &feature_reader, ConsensusTables &consensusTables,
                             std::vector<std::string> &SamSequence, std::map<std::string, int> &IDrefMap,
                             std::ofstream &fout, uint32_t thread_size)
{
    std::string sep = "\t";

    std::mutex mutex_gen;
    std::mutex mutex_comb;
    if (!feature_reader)
    {
        std::cout << "file not open!\n";
        return;
    }

    auto generate_func = [&mutex_gen, &feature_reader, &sep]() {
        // 读取5000条feature特征
        std::lock_guard<std::mutex> lock(mutex_gen);
        std::string line;
        std::vector<std::string> header_parts;
        if (feature_reader.good() && !feature_reader.eof())
        {
            struct FeatureRecord fr;
            getline(feature_reader, line);
            if (line.size() == 0)
                return std::pair<bool, FeatureRecord>(false, fr);
            SplitString(line, header_parts, sep);
            fr.refname = header_parts[0].substr(1, header_parts[0].size() - 1);
            fr.region_start = std::atoi(header_parts[1].c_str());
            fr.region_end = std::atoi(header_parts[2].c_str());
            getline(feature_reader, line);
            if (line.size() == 0)
                return std::pair<bool, FeatureRecord>(false, fr);
            fr.feature_str = line;
            return std::pair<bool, FeatureRecord>(true, fr);
        }
        else
        {
            struct FeatureRecord fr;
            return std::pair<bool, FeatureRecord>(false, fr);
        }
    };

    auto combine_func = [&mutex_comb](std::vector<FeatureTargetRecord> &finished_pool, std::ofstream &fout) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (auto record : finished_pool)
        {
            fout << ">" << record.refname << "\t" << record.region_start << "\t" << record.region_end << "\t"
                 << record.label_str << std::endl;
            fout << record.feature_str << std::endl;
        }
    };

    auto work_func = [generate_func, combine_func, &consensusTables, &SamSequence, &IDrefMap, &fout](size_t) {
        // 多线程为每条feature特征找对应的target
        std::vector<FeatureTargetRecord> finished_pool;
        while (true)
        {
            std::pair<bool, FeatureRecord> product = generate_func();
            if (!product.first)
                break;
            else
            {
                FeatureRecord record = product.second;
                std::string label_seq;
                if (IDrefMap.find(record.refname) == IDrefMap.end())
                {
                    label_seq = "";
                    continue;
                }
                int refID = IDrefMap.at(record.refname);
                std::vector<TableElem> &table = consensusTables.at(refID);
                for (int target_indx = record.region_start; target_indx <= record.region_end; target_indx++)
                {
                    TableElem &elem = table[target_indx];
                    int delele_and_insert_count = elem.number_of_Insertion + elem.number_of_Deletion;
                    int base_count = elem.number_of_M;
                    if (delele_and_insert_count > 1 or base_count > 1 or base_count + delele_and_insert_count > 2)
                    {
                        label_seq = "";
                        break;
                    }
                }
                // 该位置没有参考序列覆盖到
                if (table[record.region_start].SequencePosition.size() == 0 or
                    table[record.region_end].SequencePosition.size() == 0)
                {
                    label_seq = "";
                    continue;
                }
                int samSequenceIndex_start = table[record.region_start].SequencePosition[0].first;
                int seq_start = table[record.region_start].SequencePosition[0].second;
                int samSequenceIndex_end = table[record.region_end].SequencePosition[0].first;
                int seq_end = table[record.region_end].SequencePosition[0].second;

                if (samSequenceIndex_start != samSequenceIndex_end)
                {
                    // 该区域的左右端点的序列不同，中间断开了
                    printf("%s:%d-%d disconnected!\n", record.refname.c_str(), record.region_start, record.region_end);
                    label_seq = "";
                    continue;
                }
                label_seq = SamSequence[samSequenceIndex_start].substr(seq_start, seq_end - seq_start + 1);
                std::set<char> cters;
                for (auto v : label_seq)
                    cters.insert(toupper(v));
                if (cters.size() > 5)
                {
                    // label sequence has characters besides "ATCGN"
                    printf("%s:%d-%d label_seq:%s\n", record.refname.c_str(), record.region_start, record.region_end, label_seq.c_str());
                    label_seq = "";
                    continue;
                }
                else
                {
                    for (auto v : cters)
                    {
                        if (v != 'A' and v != 'C' and v != 'G' and v != 'T' and v != 'N')
                        {
                            printf("%s:%d-%d label_seq:%s\n", record.refname.c_str(), record.region_start, record.region_end, label_seq.c_str());
                            label_seq = "";
                            continue;
                        }
                    }
                }
                if (label_seq != "")
                {
                    FeatureTargetRecord ftr;
                    ftr.feature_str = record.feature_str;
                    ftr.refname = record.refname;
                    ftr.region_start = record.region_start;
                    ftr.region_end = record.region_end;
                    ftr.label_str = label_seq;
                    finished_pool.push_back(ftr);
                }
                if (finished_pool.size() > 5000)
                {
                    combine_func(finished_pool, fout);
                    finished_pool.clear();
                }
            }
        }
    };

    MultiThreadRun2(thread_size, work_func);
}
