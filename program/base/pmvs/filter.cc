#include <numeric>
#include <ctime>
#include <time.h>
#include <thread>

#include "Eigen/Dense"

#include "findMatch.h"
#include "filter.h"

using namespace Patch;
using namespace PMVS3;

Cfilter::Cfilter(CfindMatch& findMatch) : m_fm(findMatch)
{
}

void Cfilter::init(void)
{
}

void Cfilter::run(void)
{
    setDepthMapsVGridsVPGridsAddPatchV(0);

    filterOutside();
    setDepthMapsVGridsVPGridsAddPatchV(1);

    filterExact();
    setDepthMapsVGridsVPGridsAddPatchV(1);

    filterNeighbor(1);
    setDepthMapsVGridsVPGridsAddPatchV(1);

    filterSmallGroups();
    setDepthMapsVGridsVPGridsAddPatchV(1);
}

void Cfilter::filterOutside(void)
{
    time_t tv;
    time(&tv); 
    time_t curtime = tv;
    std::cerr << "FilterOutside" << std::endl;
    //??? notice (1) here to avoid removing m_fix=1
    m_fm.m_pos.collectPatches(1);

    const int psize = (int)m_fm.m_pos.m_ppatches.size();  
    m_gains.resize(psize);

    std::cerr << "mainbody: " << std::flush;

    m_fm.m_count = 0;
    std::vector<std::thread> threads(m_fm.m_CPU);
    for (auto& t : threads) t = std::thread(filterOutsideThreadTmp, this);
    for (auto& t : threads) t.join();
    std::cerr << std::endl;

    // Delete patches with positive m_gains
    int count = 0;

    double ave = 0.0f;
    double ave2 = 0.0f;
    int denom = 0;  

    for (int p = 0; p < psize; ++p)
    {
        ave += m_gains[p];
        ave2 += m_gains[p] * m_gains[p];
        ++denom;
    
        if (m_gains[p] < 0.0)
        {
            m_fm.m_pos.removePatch(m_fm.m_pos.m_ppatches[p]);
            count++;
        }
    }

    if (denom == 0) denom = 1;
    ave /= denom;
    ave2 /= denom;
    ave2 = sqrt(std::max(0.0, ave2 - ave * ave));
    std::cerr << "Gain (ave/var): " << ave << ' ' << ave2 << std::endl;

    time(&tv);
    std::cerr << (int)m_fm.m_pos.m_ppatches.size() << " -> "
              << (int)m_fm.m_pos.m_ppatches.size() - count << " ("
              << 100 * ((int)m_fm.m_pos.m_ppatches.size() - count) / (float)m_fm.m_pos.m_ppatches.size()
              << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

float Cfilter::computeGain(const Patch::Cpatch& patch, const int lock)
{
    float gain = patch.score2(m_fm.m_nccThreshold);

    const int size = (int)patch.m_images.size();  
    for (int i = 0; i < size; ++i)
    {
        const int& index = patch.m_images[i];
        if (m_fm.m_tnum <= index) continue;

        const int& ix = patch.m_grids[i][0];
        const int& iy = patch.m_grids[i][1];
        const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;

        float maxpressure = 0.0f;
        if (lock) m_fm.m_imageLocks[index].lock();

        for (int j = 0; j < (int)m_fm.m_pos.m_pgrids[index][index2].size(); ++j)
        {
            if (!m_fm.isNeighbor(patch, *m_fm.m_pos.m_pgrids[index][index2][j], m_fm.m_neighborThreshold1))
            maxpressure = std::max(maxpressure, m_fm.m_pos.m_pgrids[index][index2][j]->m_ncc - m_fm.m_nccThreshold);
        }
        if (lock) m_fm.m_imageLocks[index].unlock();

        gain -= maxpressure;
    }

    const int vsize = (int)patch.m_vimages.size();
    for (int i = 0; i < vsize; ++i)
    {
        const int& index = patch.m_vimages[i];
        if (m_fm.m_tnum <= index) continue;

        const float pdepth = m_fm.m_pss.computeDepth(index, patch.m_coord);    

        const int& ix = patch.m_vgrids[i][0];
        const int& iy = patch.m_vgrids[i][1];
        const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;
        float maxpressure = 0.0f;      

        if (lock) m_fm.m_imageLocks[index].lock();

        for (int j = 0; j < (int)m_fm.m_pos.m_pgrids[index][index2].size(); ++j)
        {
            const float bdepth = m_fm.m_pss.computeDepth(index, m_fm.m_pos.m_pgrids[index][index2][j]->m_coord);
            if (pdepth < bdepth && !m_fm.isNeighbor(patch, *m_fm.m_pos.m_pgrids[index][index2][j], m_fm.m_neighborThreshold1))
            {
                maxpressure = std::max(maxpressure, m_fm.m_pos.m_pgrids[index][index2][j]->m_ncc - m_fm.m_nccThreshold);
            }
        }
        if (lock) m_fm.m_imageLocks[index].unlock();

        gain -= maxpressure;
    }
    return gain;
}

void Cfilter::filterOutsideThread(void)
{
    m_fm.m_lock.lock();
    const int id = m_fm.m_count++;
    m_fm.m_lock.unlock();

    const int size = (int)m_fm.m_pos.m_ppatches.size();  
    const int itmp = (int)ceil(size / (float)m_fm.m_CPU);
    const int begin = id * itmp;
    const int end = std::min(size, (id + 1) * itmp);

    for (int p = begin; p < end; ++p)
    {
        Ppatch& ppatch = m_fm.m_pos.m_ppatches[p];
        m_gains[p] = ppatch->score2(m_fm.m_nccThreshold);

        const int size = (int)ppatch->m_images.size();  
        for (int i = 0; i < size; ++i)
        {
            const int& index = ppatch->m_images[i];
            if (m_fm.m_tnum <= index) continue;

            const int& ix = ppatch->m_grids[i][0];
            const int& iy = ppatch->m_grids[i][1];
            const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;

            float maxpressure = 0.0f;
            for (int j = 0; j < (int)m_fm.m_pos.m_pgrids[index][index2].size(); ++j)
            {
                if (!m_fm.isNeighbor(*ppatch, *m_fm.m_pos.m_pgrids[index][index2][j], m_fm.m_neighborThreshold1))
                    maxpressure = std::max(maxpressure, m_fm.m_pos.m_pgrids[index][index2][j]->m_ncc - m_fm.m_nccThreshold);
            }

            m_gains[p] -= maxpressure;
        }

        const int vsize = (int)ppatch->m_vimages.size();
        for (int i = 0; i < vsize; ++i)
        {
            const int& index = ppatch->m_vimages[i];
            if (m_fm.m_tnum <= index) continue;

            const float pdepth = m_fm.m_pss.computeDepth(index, ppatch->m_coord);    

            const int& ix = ppatch->m_vgrids[i][0];
            const int& iy = ppatch->m_vgrids[i][1];
            const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;
            float maxpressure = 0.0f;      

            for (int j = 0; j < (int)m_fm.m_pos.m_pgrids[index][index2].size(); ++j)
            {
                const float bdepth = m_fm.m_pss.computeDepth(index, m_fm.m_pos.m_pgrids[index][index2][j]->m_coord);
                if (pdepth < bdepth && !m_fm.isNeighbor(*ppatch, *m_fm.m_pos.m_pgrids[index][index2][j], m_fm.m_neighborThreshold1))
                {
                    maxpressure = std::max(maxpressure, m_fm.m_pos.m_pgrids[index][index2][j]->m_ncc - m_fm.m_nccThreshold);
                }
            }
            m_gains[p] -= maxpressure;
        }
    }
}

void* Cfilter::filterOutsideThreadTmp(void* arg)
{
    ((Cfilter*)arg)->filterOutsideThread();
    return NULL;
}

void Cfilter::filterExact(void)
{
    time_t tv;
    time(&tv); 
    time_t curtime = tv;
    std::cerr << "Filter Exact: " << std::flush;

    //??? cannot use (1) because we use patch.m_id to set newimages,....
    m_fm.m_pos.collectPatches();
    const int psize = (int)m_fm.m_pos.m_ppatches.size();

    // dis associate images
    m_newimages.clear();            m_newgrids.clear();
    m_removeimages.clear();         m_removegrids.clear();
    m_newimages.resize(psize);      m_newgrids.resize(psize);
    m_removeimages.resize(psize);   m_removegrids.resize(psize);

    m_fm.m_count = 0;
    std::vector<std::thread> threads(m_fm.m_CPU);
    for (auto& t : threads) t = std::thread(filterExactThreadTmp, this);
    for (auto& t : threads) t.join();

    std::cerr << std::endl;

    for (int p = 0; p < psize; ++p)
    {
        if (m_fm.m_pos.m_ppatches[p]->m_fix) continue;

        for (int i = 0; i < (int)m_removeimages[p].size(); ++i)
        {
            const int index = m_removeimages[p][i];
            if (m_fm.m_tnum <= index)
            {
                std::cerr << "MUST NOT COME HERE" << std::endl;
                exit (1);
            }

            const int ix = m_removegrids[p][i][0];      const int iy = m_removegrids[p][i][1];
            const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;

            m_fm.m_pos.m_pgrids[index][index2].
            erase(remove(m_fm.m_pos.m_pgrids[index][index2].begin(),
            m_fm.m_pos.m_pgrids[index][index2].end(),
            m_fm.m_pos.m_ppatches[p]),
            m_fm.m_pos.m_pgrids[index][index2].end());
        }
    }

    m_fm.m_debug = 1;
  
    int count = 0;
    for (int p = 0; p < psize; ++p)
    {
        if (m_fm.m_pos.m_ppatches[p]->m_fix) continue;
    
        Cpatch& patch = *m_fm.m_pos.m_ppatches[p];

        // This should be images in targetting images. Has to come before the next for-loop.
        patch.m_timages = (int)m_newimages[p].size();
    
        for (int i = 0; i < (int)patch.m_images.size(); ++i)
        {
            const int& index = patch.m_images[i];
            if (m_fm.m_tnum <= index)
            {
                m_newimages[p].push_back(patch.m_images[i]);
                m_newgrids[p].push_back(patch.m_grids[i]);
            }
        }

        patch.m_images.swap(m_newimages[p]);
        patch.m_grids.swap(m_newgrids[p]);

        if (m_fm.m_minImageNumThreshold <= (int)patch.m_images.size())
        {
            m_fm.m_optim.setRefImage(patch, 0);
            m_fm.m_pos.setGrids(patch);
        }

        if ((int)patch.m_images.size() < m_fm.m_minImageNumThreshold)
        {
            m_fm.m_pos.removePatch(m_fm.m_pos.m_ppatches[p]);
            count++;
        }
    }

    time(&tv); 
    std::cerr << (int)m_fm.m_pos.m_ppatches.size() << " -> "
              << (int)m_fm.m_pos.m_ppatches.size() - count << " ("
              << 100 * ((int)m_fm.m_pos.m_ppatches.size() - count) / (float)m_fm.m_pos.m_ppatches.size()
              << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

void Cfilter::filterExactThread(void)
{
    const int psize = (int)m_fm.m_pos.m_ppatches.size();
    std::vector<std::vector<int> > newimages, removeimages;
    std::vector<std::vector<TVec2<int> > > newgrids, removegrids;
    newimages.resize(psize);  removeimages.resize(psize);
    newgrids.resize(psize);   removegrids.resize(psize);

    while (1)
    {
        m_fm.m_lock.lock();
        const int image = m_fm.m_count++;
        m_fm.m_lock.unlock();

        if (m_fm.m_tnum <= image) break;

        std::cerr << '*' << std::flush;

        const int& w = m_fm.m_pos.m_gwidths[image];
        const int& h = m_fm.m_pos.m_gheights[image];
        int index = -1;
        for (int y = 0; y < h; ++y)
        {
            for (int x = 0; x < w; ++x)
            {
                ++index;
                for (int i = 0; i < (int)m_fm.m_pos.m_pgrids[image][index].size(); ++i)
                {
                    const Cpatch& patch = *m_fm.m_pos.m_pgrids[image][index][i];
                    if (patch.m_fix) continue;
          
                    int safe = 0;

                    if (m_fm.m_pos.isVisible(patch, image, x, y, m_fm.m_neighborThreshold1, 0))                         safe = 1;
                    // use 4 neighbors?
                    else if (0 < x && m_fm.m_pos.isVisible(patch, image, x - 1, y, m_fm.m_neighborThreshold1, 0))       safe = 1;
                    else if (x < w - 1 && m_fm.m_pos.isVisible(patch, image, x + 1, y, m_fm.m_neighborThreshold1, 0))   safe = 1;
                    else if (0 < y && m_fm.m_pos.isVisible(patch, image, x, y - 1, m_fm.m_neighborThreshold1, 0))       safe = 1;
                    else if (y < h - 1 && m_fm.m_pos.isVisible(patch, image, x, y + 1, m_fm.m_neighborThreshold1, 0))   safe = 1;

                    if (safe)
                    {
                        newimages[patch.m_id].push_back(image);
                        newgrids[patch.m_id].push_back(TVec2<int>(x, y));
                    } else
                    {
                        removeimages[patch.m_id].push_back(image);
                        removegrids[patch.m_id].push_back(TVec2<int>(x, y));
                    }
                }
            }
        }
    }

    m_fm.m_lock.lock();
    for (int p = 0; p < psize; ++p)
    {
        m_newimages[p].insert(m_newimages[p].end(), newimages[p].begin(), newimages[p].end());
        m_newgrids[p].insert(m_newgrids[p].end(), newgrids[p].begin(), newgrids[p].end());
        m_removeimages[p].insert(m_removeimages[p].end(), removeimages[p].begin(), removeimages[p].end());
        m_removegrids[p].insert(m_removegrids[p].end(), removegrids[p].begin(), removegrids[p].end());
    }
    m_fm.m_lock.unlock();
}

void* Cfilter::filterExactThreadTmp(void* arg)
{
    ((Cfilter*)arg)->filterExactThread();
    return NULL;
}

void Cfilter::filterNeighborThread(void)
{
    const int size = (int)m_fm.m_pos.m_ppatches.size();  
    while (1)
    {
        int jtmp = -1;
        m_fm.m_lock.lock();
        if (!m_fm.m_jobs.empty())
        {
            jtmp = m_fm.m_jobs.front();
            m_fm.m_jobs.pop_front();
        }
        m_fm.m_lock.unlock();
        if (jtmp == -1) break;

        const int begin = m_fm.m_junit * jtmp;
        const int end = std::min(size, m_fm.m_junit * (jtmp + 1));

        for (int p = begin; p < end; ++p)
        {
            Ppatch& ppatch = m_fm.m_pos.m_ppatches[p];
            if (m_rejects[p]) continue;

            std::vector<Ppatch> neighbors;
            m_fm.m_pos.findNeighbors(*ppatch, neighbors, 0, 4, 2, 1);

            if ((int)neighbors.size() < 6) m_rejects[p] = m_time + 1;
            else
            {
                // Fit a quadratic surface
                if (filterQuad(*ppatch, neighbors)) m_rejects[p] = m_time + 1;
            }
        }
    }
}

int Cfilter::filterQuad(const Patch::Cpatch& patch, const std::vector<Ppatch>& neighbors) const
{
    const int nsize = (int)neighbors.size();

    Eigen::MatrixXf A(nsize, 5);
    Eigen::VectorXf b(nsize);
    Eigen::VectorXf x(nsize);

    Vec4f xdir, ydir;
    ortho(patch.m_normal, xdir, ydir);

    float h = 0.0f;
    for (int n = 0; n < nsize; ++n) h += norm(neighbors[n]->m_coord - patch.m_coord);
    h /= nsize;

    std::vector<float> fxs, fys, fzs;
    fxs.resize(nsize);
    fys.resize(nsize);
    fzs.resize(nsize);

    for (int n = 0; n < nsize; ++n)
    {
        Vec4f diff = neighbors[n]->m_coord - patch.m_coord;
        fxs[n] = diff * xdir / h;
        fys[n] = diff * ydir / h;
        fzs[n] = diff * patch.m_normal;

        A(n, 0) = fxs[n] * fxs[n];
        A(n, 1) = fys[n] * fys[n];
        A(n, 2) = fxs[n] * fys[n];
        A(n, 3) = fxs[n];
        A(n, 4) = fys[n];
        b[n]    = fzs[n];
    }

    x = A.colPivHouseholderQr().solve(b);

    // Compute residual divided by m_dscale
    const int inum = std::min(m_fm.m_tau, (int)patch.m_images.size());
    float unit = 0.0;
    for (int i = 0; i < inum; ++i) unit += m_fm.m_optim.getUnit(patch.m_images[i], patch.m_coord);
    unit /= inum;

    float residual = 0.0f;
    for (int n = 0; n < nsize; ++n)
    {
        const float res =   x[0] * (fxs[n] * fxs[n]) +
                            x[1] * (fys[n] * fys[n]) +
                            x[2] * (fxs[n] * fys[n]) +
                            x[3] * fxs[n]  +
                            x[4] * fys[n]  - fzs[n];
        residual += fabs(res) / unit;
    }

    residual /= (nsize - 5);

    if (residual < m_fm.m_quadThreshold) return 0;
    else                                 return 1;
}

void* Cfilter::filterNeighborThreadTmp(void* arg)
{
    ((Cfilter*)arg)->filterNeighborThread();
    return NULL;
}

void Cfilter::filterNeighbor(const int times)
{
    time_t tv;
    time(&tv); 
    time_t curtime = tv;
    std::cerr << "FilterNeighbor:\t" << std::flush;

    //??? notice (1) to avoid removing m_fix=1
    m_fm.m_pos.collectPatches(1);
    if (m_fm.m_pos.m_ppatches.empty())
    return;

    m_rejects.resize((int)m_fm.m_pos.m_ppatches.size());
    fill(m_rejects.begin(), m_rejects.end(), 0);

    // Lapack is not thread-safe? Sometimes, the code gets stuck here.
    int count = 0;
    for (m_time = 0; m_time < times; ++m_time)
    {
        m_fm.m_count = 0;

        m_fm.m_jobs.clear();
        const int jtmp = (int)ceil(m_fm.m_pos.m_ppatches.size() / (float)m_fm.m_junit);
        for (int j = 0; j < jtmp; ++j) m_fm.m_jobs.push_back(j);

        std::vector<std::thread> threads(m_fm.m_CPU);
        for (auto& t : threads) t = std::thread(filterNeighborThreadTmp, this);
        for (auto& t : threads) t.join();

        auto bpatch  = m_fm.m_pos.m_ppatches.begin();
        auto epatch  = m_fm.m_pos.m_ppatches.end();
        auto breject = m_rejects.begin();

        while (bpatch != epatch)
        {
            if ((*breject) == m_time + 1)
            {
                count++;
                m_fm.m_pos.removePatch(*bpatch);;
            }

            ++bpatch;
            ++breject;
        }
    }

    time(&tv);
    std::cerr << (int)m_fm.m_pos.m_ppatches.size() << " -> "
              << (int)m_fm.m_pos.m_ppatches.size() - count << " ("
              << 100 * ((int)m_fm.m_pos.m_ppatches.size() - count) / (float)m_fm.m_pos.m_ppatches.size()
              << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

void Cfilter::filterSmallGroups(void)
{
    time_t tv;
    time(&tv); 
    time_t curtime = tv;
    std::cerr << "FilterGroups:\t" << std::flush;
    m_fm.m_pos.collectPatches();
    if (m_fm.m_pos.m_ppatches.empty()) return;

    const int psize = (int)m_fm.m_pos.m_ppatches.size();
    std::vector<int> label;
    label.resize(psize);
    std::fill(label.begin(), label.end(), -1);

    std::list<int> untouch;
    auto bpatch = m_fm.m_pos.m_ppatches.begin();
    for (int p = 0; p < psize; ++p, ++bpatch)
    {
        untouch.push_back(p);
        (*bpatch)->m_flag = p;
    }

    int id = -1;
    while (!untouch.empty())
    {
        const int pid = untouch.front();
        untouch.pop_front();

        if (label[pid] != -1) continue;

        label[pid] = ++id;
        std::list<int> ltmp;
        ltmp.push_back(pid);

        while (!ltmp.empty())
        {
            const int ptmp = ltmp.front();
            ltmp.pop_front();

            filterSmallGroupsSub(ptmp, id, label, ltmp);
        }
    }
    id++;

    std::vector<int> size;
    size.resize(id);
    auto bite = label.begin();
    auto eite = label.end();
    while (bite != eite)
    {
        ++size[*bite];
        ++bite;
    }

    const int threshold = std::max(20, psize / 10000);
    std::cerr << threshold << std::endl;

    bite = size.begin();
    eite = size.end();
    while (bite != eite)
    {
        if (*bite < threshold)  *bite = 0;
        else                    *bite = 1;
        ++bite;
    }

    int count = 0;

    bite = label.begin();
    eite = label.end();
    bpatch = m_fm.m_pos.m_ppatches.begin();
    while (bite != eite)
    {
        if ((*bpatch)->m_fix)
        {
            ++bite;
            ++bpatch;
            continue;
        }

        if (size[*bite] == 0)
        {
            m_fm.m_pos.removePatch(*bpatch);
            count++;
        }
        ++bite;
        ++bpatch;
    }

    time(&tv);
    std::cerr << (int)m_fm.m_pos.m_ppatches.size() << " -> "
              << (int)m_fm.m_pos.m_ppatches.size() - count << " ("
              << 100 * ((int)m_fm.m_pos.m_ppatches.size() - count) / (float)m_fm.m_pos.m_ppatches.size()
              << "%)\t" << (tv - curtime)/CLOCKS_PER_SEC << " secs" << std::endl;
}

void Cfilter::filterSmallGroupsSub(const int pid, const int id, std::vector<int>& label, std::list<int>& ltmp) const
{
    // Find neighbors of ptmp and set their ids
    const Cpatch& patch = *m_fm.m_pos.m_ppatches[pid];

    const int index = patch.m_images[0];
    const int ix = patch.m_grids[0][0];
    const int iy = patch.m_grids[0][1];
    const int gwidth = m_fm.m_pos.m_gwidths[index];
    const int gheight = m_fm.m_pos.m_gheights[index];

    for (int y = -1; y <= 1; ++y)
    {
        const int iytmp = iy + y;
        if (iytmp < 0 || gheight <= iytmp) continue;
        for (int x = -1; x <= 1; ++x)
        {
            const int ixtmp = ix + x;
            if (ixtmp < 0 || gwidth <= ixtmp)
            continue;

            const int index2 = iytmp * gwidth + ixtmp;
            auto bgrid = m_fm.m_pos.m_pgrids[index][index2].begin();
            auto egrid = m_fm.m_pos.m_pgrids[index][index2].end();
            while (bgrid != egrid)
            {
                const int itmp = (*bgrid)->m_flag;
                if (label[itmp] != -1)
                {
                    ++bgrid;
                    continue;
                }

                if (m_fm.isNeighbor(patch, **bgrid, m_fm.m_neighborThreshold2))
                {
                    label[itmp] = id;
                    ltmp.push_back(itmp);
                }
                ++bgrid;
            }

            bgrid = m_fm.m_pos.m_vpgrids[index][index2].begin();
            egrid = m_fm.m_pos.m_vpgrids[index][index2].end();
            while (bgrid != egrid)
            {
                const int itmp = (*bgrid)->m_flag;
                if (label[itmp] != -1)
                {
                    ++bgrid;
                    continue;
                }

                if (m_fm.isNeighbor(patch, **bgrid, m_fm.m_neighborThreshold2))
                {
                    label[itmp] = id;
                    ltmp.push_back(itmp);
                }
                ++bgrid;
            }
        }
    }
}

void Cfilter::setDepthMaps(void)
{
    for (int index = 0; index < m_fm.m_tnum; ++index)
    {
        std::fill(m_fm.m_pos.m_dpgrids[index].begin(), m_fm.m_pos.m_dpgrids[index].end(), m_fm.m_pos.m_MAXDEPTH);
    }

    m_fm.m_count = 0;
    std::vector<std::thread> threads(m_fm.m_CPU);
    for (auto& t : threads) t = std::thread(setDepthMapsThreadTmp, this);
    for (auto& t : threads) t.join();
}

void* Cfilter::setDepthMapsThreadTmp(void* arg)
{
    ((Cfilter*)arg)->setDepthMapsThread();
    return NULL;
}

void Cfilter::setDepthMapsThread(void)
{
    while (1)
    {
        m_fm.m_lock.lock();
        const int index = m_fm.m_count++;
        m_fm.m_lock.unlock();

        if (m_fm.m_tnum <= index) break;

        const int gwidth  = m_fm.m_pos.m_gwidths[index];
        const int gheight = m_fm.m_pos.m_gheights[index];

        auto bpatch = m_fm.m_pos.m_ppatches.begin();
        auto epatch = m_fm.m_pos.m_ppatches.end();

        while (bpatch != epatch)
        {
            Ppatch& ppatch = *bpatch;
            const Vec3f icoord = m_fm.m_pss.project(index, ppatch->m_coord, m_fm.m_level);

            const float fx  = icoord[0] / m_fm.m_csize;
            const int xs[2] = {(int)floor(fx), (int)ceil(fx)};
            const float fy  = icoord[1] / m_fm.m_csize;
            const int ys[2] = {(int)floor(fy), (int)ceil(fy)};
    
            const float depth = m_fm.m_pss.m_photos[index].m_oaxis * ppatch->m_coord;

            for (int j = 0; j < 2; ++j)
            {
                for (int i = 0; i < 2; ++i)
                {
                    if (xs[i] < 0 || gwidth <= xs[i] || ys[j] < 0 || gheight <= ys[j]) continue;
                    const int index2 = ys[j] * gwidth + xs[i];

                    if (m_fm.m_pos.m_dpgrids[index][index2] == m_fm.m_pos.m_MAXDEPTH)
                    {
                        m_fm.m_pos.m_dpgrids[index][index2] = ppatch;
                    } else
                    {
                        const float dtmp = m_fm.m_pss.m_photos[index].m_oaxis * m_fm.m_pos.m_dpgrids[index][index2]->m_coord;

                        if (depth < dtmp) m_fm.m_pos.m_dpgrids[index][index2] = ppatch;
                    }
                }
            }
            ++bpatch;
        }
    }
}

void Cfilter::setDepthMapsVGridsVPGridsAddPatchV(const int additive)
{
    m_fm.m_pos.collectPatches();
    setDepthMaps();

    // clear m_vpgrids
    for (int index = 0; index < m_fm.m_tnum; ++index)
    {
        auto bvvp = m_fm.m_pos.m_vpgrids[index].begin();
        auto evvp = m_fm.m_pos.m_vpgrids[index].end();
        while (bvvp != evvp)
        {
            (*bvvp).clear();
            ++bvvp;
        }
    }

    if (additive == 0)
    {
        auto bpatch = m_fm.m_pos.m_ppatches.begin();
        auto epatch = m_fm.m_pos.m_ppatches.end();
        while (bpatch != epatch)
        {
            (*bpatch)->m_vimages.clear();
            (*bpatch)->m_vgrids.clear();
            ++bpatch;
        }
    }

    m_fm.m_count = 0;

    std::vector<std::thread> threads0(m_fm.m_CPU);
    for (auto& t : threads0) t = std::thread(setVGridsVPGridsThreadTmp, this);
    for (auto& t : threads0) t.join();

    std::vector<std::thread> threads1(m_fm.m_CPU);
    for (auto& t : threads1) t = std::thread(addPatchVThreadTmp, this);
    for (auto& t : threads1) t.join();
}

void* Cfilter::setVGridsVPGridsThreadTmp(void* arg)
{
    ((Cfilter*)arg)->setVGridsVPGridsThread();
    return NULL;
}

void* Cfilter::addPatchVThreadTmp(void* arg)
{
    ((Cfilter*)arg)->addPatchVThread();
    return NULL;
}

void Cfilter::setVGridsVPGridsThread(void)
{
    const int noj = 1000;
    const int size = (int)m_fm.m_pos.m_ppatches.size();    
    const int job = std::max(1, size / (noj - 1));

    while (1)
    {
        m_fm.m_lock.lock();
        const int id = m_fm.m_count++;
        m_fm.m_lock.unlock();

        const int begin = id * job;
        const int end   = std::min(size, (id + 1) * job);

        if (size <= begin)
        break;

        // add patches to m_vpgrids
        for (int p = begin; p < end; ++p)
        {
            Ppatch& ppatch = m_fm.m_pos.m_ppatches[p];
            m_fm.m_pos.setVImagesVGrids(ppatch);
        }
    }
}

void Cfilter::addPatchVThread(void)
{
    while (1)
    {
        m_fm.m_lock.lock();
        const int index = m_fm.m_count++;
        m_fm.m_lock.unlock();

        if (m_fm.m_tnum <= index) break;

        auto bpatch = m_fm.m_pos.m_ppatches.begin();
        auto epatch = m_fm.m_pos.m_ppatches.end();
        while (bpatch != epatch)
        {
            Ppatch& ppatch = *bpatch;
            auto bimage = ppatch->m_vimages.begin();
            auto eimage = ppatch->m_vimages.end();
            auto bgrid = ppatch->m_vgrids.begin();
            while (bimage != eimage)
            {
                if (*bimage == index)
                {
                    const int& ix = (*bgrid)[0];
                    const int& iy = (*bgrid)[1];
                    const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;
                    m_fm.m_pos.m_vpgrids[index][index2].push_back(ppatch);
                    break;
                }
                ++bimage;
                ++bgrid;
            }
            ++bpatch;
        }
    }
}
