#pragma once

#include <list>
#include <string>
#include <mutex>

#include "../image/photoSetS.h"
#include "point.h"

namespace Image
{
    class CphotoSetS;
};

namespace PMVS3
{

class CdetectFeatures
{
public:
    CdetectFeatures() = default;
    virtual ~CdetectFeatures();

    void run(const Image::CphotoSetS& pss, const int num, const int csize, const int level, const int CPU = 1);

    std::vector<std::vector<Cpoint>> m_points;

protected:
    const Image::CphotoSetS* m_ppss;
    int m_csize;
    int m_level;

    std::mutex m_rwlock;
    int m_CPU;

    std::list<int> m_jobs;

    void runThread(void);
    static void* runThreadTmp(void*arg);
};

};
