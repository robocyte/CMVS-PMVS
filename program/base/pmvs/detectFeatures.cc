#include <iostream>
#include <fstream>
#include <numeric>
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
    m_jobs.resize(num);
    std::iota(m_jobs.begin(), m_jobs.end(), 0);

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
        {
            std::lock_guard<std::mutex> lock(m_rwlock);
            if (!m_jobs.empty())
            {
                index = m_jobs.front();
                m_jobs.pop_front();
            }
        }

        if (index == -1) break;

        const int image = m_ppss->m_images[index];
        std::cerr << image << ' ' << std::flush;

        const float sigma      = 4.0f;  // Parameters for harris...
        const float firstScale = 1.0f;  // ... for DoG
        const float lastScale  = 3.0f;  // ... for DoG

        // Harris
        {
            Charris harris;
            std::multiset<Cpoint> result;
            harris.run(m_ppss->m_photos[index].getImage(m_level),
                       m_ppss->m_photos[index].Cimage::getMask(m_level),
                       m_ppss->m_photos[index].Cimage::getEdge(m_level),
                       m_ppss->m_photos[index].getWidth(m_level),
                       m_ppss->m_photos[index].getHeight(m_level), m_csize, sigma, result);
      
            for (const auto& point : result) m_points[index].push_back(point);
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

            for (const auto& point : result) m_points[index].push_back(point);
        }
    }
}
