#pragma once

#include <vector>
#include <queue>
#include <list>
#include "patchOrganizerS.h"

namespace PMVS3
{

class CfindMatch;
  
class Cexpand
{
public:
    Cexpand(CfindMatch& findMatch);
    ~Cexpand() {};

    void init(void);
    void run(void);

    float computeRadius(const Patch::Cpatch& patch);

protected:
    int expandSub(const Patch::Ppatch& orgppatch, const int id, const Vec4f& canCoord);

    int updateCounts(const Patch::Cpatch& patch);

    int checkCounts(Patch::Cpatch& patch);

    void findEmptyBlocks(const Patch::Ppatch& ppatch, std::vector<std::vector<Vec4f>>& canCoords);

    std::priority_queue<Patch::Ppatch, std::vector<Patch::Ppatch>, P_compare> m_queue;

    CfindMatch& m_fm;

    void expandThread(void);
    static void* expandThreadTmp(void* arg);

    std::vector<int> m_ecounts;   // Number of trials
    std::vector<int> m_fcounts0;  // Number of failures in the prep
    std::vector<int> m_fcounts1;  // Number of failures in the post processing
    std::vector<int> m_pcounts;   // Number passes
};

};