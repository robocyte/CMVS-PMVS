#include <iostream>
#include <fstream>
#include <thread>

#include "../image/image.h"
#include "detectFeatures.h"
#include "harris.h"
#include "dog.h"
#include "point.h"

using namespace PMVS3;
using namespace Image;

void CdetectFeatures::run(const CphotoSetS& pss, const int num, const int csize, const int level, const int CPU)
{
    m_ppss = &pss;
    m_csize = csize;
    m_level = level;
    m_CPU = CPU;

    m_points.clear();
    m_points.resize(num);

    for (int index = 0; index < num; ++index) m_jobs.push_back(index);

    std::vector<std::thread> threads(m_CPU);
    for (auto& t : threads) t = std::thread(&CdetectFeatures::runThread, this);
    for (auto& t : threads) t.join();

    std::cerr << "done" << std::endl;
}

void CdetectFeatures::runThread()
{
    while (1)
    {
        int index = -1;
        m_rwlock.lock();
        if (!m_jobs.empty())
        {
            index = m_jobs.front();
            m_jobs.pop_front();
        }
        m_rwlock.unlock();
        if (index == -1) break;

        const int image = m_ppss->m_images[index];
        std::cerr << image << ' ' << std::flush;

        //?????????????  May need file lock, because targetting images should not overlap among multiple processors.    
        char buffer[1024];
        sprintf(buffer, "%smodels/%08d.affin%d", m_ppss->m_prefix.c_str(), image, m_level);
        std::ifstream ifstr;
        ifstr.open(buffer);
        if (ifstr.is_open())
        {
            ifstr.close();
            continue;
        }
        ifstr.close();

        const float sigma = 4.0f;       // Parameters for harris...
        const float firstScale = 1.0f;  // ... for DoG
        const float lastScale = 3.0f;   // ... for DoG

        // Harris
        {
            Charris harris;
            std::multiset<Cpoint> result;
            harris.run(m_ppss->m_photos[index].getImage(m_level),
                       m_ppss->m_photos[index].Cimage::getMask(m_level),
                       m_ppss->m_photos[index].Cimage::getEdge(m_level),
                       m_ppss->m_photos[index].getWidth(m_level),
                       m_ppss->m_photos[index].getHeight(m_level), m_csize, sigma, result);
      
            auto rbegin = result.rbegin();
            while (rbegin != result.rend())
            {
                m_points[index].push_back(*rbegin);
                rbegin++;
            }
        }

        // DoG
        {
            Cdog dog;
            std::multiset<Cpoint> result;
            dog.run(m_ppss->m_photos[index].getImage(m_level),
                    m_ppss->m_photos[index].Cimage::getMask(m_level),
                    m_ppss->m_photos[index].Cimage::getEdge(m_level),
                    m_ppss->m_photos[index].getWidth(m_level),
                    m_ppss->m_photos[index].getHeight(m_level),
                    m_csize, firstScale, lastScale, result);

            auto rbegin = result.rbegin();
            while (rbegin != result.rend())
            {
                m_points[index].push_back(*rbegin);
                rbegin++;
            }
        }
    }
}
