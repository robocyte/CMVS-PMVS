#include <map>
#include <ctime>

#include "findMatch.h"
#include "detectFeatures.h"

using namespace PMVS3;
using namespace Patch;

CfindMatch::CfindMatch(void)
    : m_pos(*this), m_seed(*this), m_expand(*this), m_filter(*this), m_optim(*this)
{
    m_debug = 0;
}

CfindMatch::~CfindMatch()
{
}

void CfindMatch::updateThreshold(void)
{
    m_nccThreshold -= 0.05f;
    m_nccThresholdBefore -= 0.05f;

    m_countThreshold1 = 2;
}

void CfindMatch::init(const Soption& option)
{
    m_timages = option.m_timages;
    m_oimages = option.m_oimages;
    m_images.clear();
    m_images.insert(m_images.end(), m_timages.begin(), m_timages.end());
    m_images.insert(m_images.end(), m_oimages.begin(), m_oimages.end());

    m_tnum = (int)m_timages.size();
    m_num  = (int)m_images.size();

    m_prefix                = option.m_prefix;
    m_level                 = option.m_level;
    m_csize                 = option.m_csize;
    m_nccThreshold          = option.m_threshold;
    m_wsize                 = option.m_wsize;
    m_minImageNumThreshold  = option.m_minImageNum;
    m_CPU                   = option.m_CPU;
    m_setEdge               = option.m_setEdge;
    m_sequenceThreshold     = option.m_sequence;

    m_junit = 100;
    m_visibleThreshold = 0.0f;
    m_visibleThresholdLoose = 0.0f;

    m_tau = std::min(option.m_minImageNum * 2, m_num);
    m_depth = 0;

    // Set target images and other images
    m_bindexes = option.m_bindexes;
    m_visdata  = option.m_visdata;
    m_visdata2 = option.m_visdata2;

    m_imageLocks.resize(m_num);
    m_countLocks.resize(m_num);

    // We set m_level + 3, to use multi-resolutional texture grabbing
    m_pss.init(m_images, m_prefix, m_level + 3, m_wsize, 1);

    if (m_setEdge != 0.0f) m_pss.setEdge(m_setEdge);
    m_pss.setDistances();

    // Detect features if not yet done
    CdetectFeatures df;
    const int fcsize = 16;
    df.run(m_pss, m_num, fcsize, m_level, m_CPU);  

    // Initialize each core member. m_pos should be first
    m_pos.init();
    m_seed.init(df.m_points);
    m_expand.init();
    m_filter.init();
    m_optim.init();

    // Init thresholds
    m_angleThreshold0 = 60.0f * (float)M_PI / 180.0f;
    m_angleThreshold1 = 60.0f * (float)M_PI / 180.0f;

    m_countThreshold0 = 2;
    m_countThreshold1 = 4;
    m_countThreshold2 = 2;

    m_neighborThreshold = 0.5f;
    m_neighborThreshold1 = 1.0f;
    m_neighborThreshold2 = 1.0f;

    m_maxAngleThreshold  = option.m_maxAngleThreshold;
    m_nccThresholdBefore = m_nccThreshold - 0.3f;
    m_quadThreshold      = option.m_quadThreshold;
    m_epThreshold        = 2.0f;
}

int CfindMatch::insideBimages(const Vec4f& coord) const
{
    for (int i = 0; i < (int)m_bindexes.size(); ++i)
    {
        const int index = m_bindexes[i];
        const Vec3f icoord = m_pss.project(index, coord, m_level);
        if (icoord[0] < 0.0 || m_pss.getWidth(index, m_level) - 1 < icoord[0] || icoord[1] < 0.0 || m_pss.getHeight(index, m_level) - 1 < icoord[1]) return 0;
    }

    return 1;
}

int CfindMatch::isNeighbor(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float neighborThreshold) const
{
    const float hunit = (m_optim.getUnit(lhs.m_images[0], lhs.m_coord) + m_optim.getUnit(rhs.m_images[0], rhs.m_coord)) / 2.0f * m_csize;
    return isNeighbor(lhs, rhs, hunit, neighborThreshold);
}

int CfindMatch::isNeighbor(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float hunit, const float neighborThreshold) const
{
    if (lhs.m_normal * rhs.m_normal < cos(120.0 * M_PI / 180.0)) return 0;
    const Vec4f diff = rhs.m_coord - lhs.m_coord;

    const float vunit = lhs.m_dscale + rhs.m_dscale;
    const float f0 = lhs.m_normal * diff;
    const float f1 = rhs.m_normal * diff;   
    float ftmp = (fabs(f0) + fabs(f1)) / 2.0f;
    ftmp /= vunit;

    // This may loosen the isneighbor testing. Need to tighten (decrease) threshold?
    const float hsize = norm(2.0f * diff - lhs.m_normal * f0 - rhs.m_normal * f1) / 2.0f / hunit;
    if (1.0f < hsize) ftmp /= std::min(2.0f, hsize);

    if (ftmp < neighborThreshold)   return 1;
    else                            return 0;
}

int CfindMatch::isNeighborRadius(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float hunit, const float neighborThreshold, const float radius) const
{
    if (lhs.m_normal * rhs.m_normal < cos(120.0 * M_PI / 180.0)) return 0;
    const Vec4f diff = rhs.m_coord - lhs.m_coord;

    const float vunit = lhs.m_dscale + rhs.m_dscale;

    const float f0 = lhs.m_normal * diff;
    const float f1 = rhs.m_normal * diff;   
    float ftmp = (fabs(f0) + fabs(f1)) / 2.0f;
    ftmp /= vunit;

    // This may loosen the isneighbor testing. Need to tighten (decrease) threshold?  
    const float hsize = norm(2.0f * diff - lhs.m_normal * f0 - rhs.m_normal * f1) / 2.0f / hunit;

    // Rradius check
    if (radius / hunit < hsize) return 0;

    if (1.0 < hsize) ftmp /= std::min(2.0f, hsize);

    if (ftmp < neighborThreshold)   return 1;
    else                            return 0;
}

void CfindMatch::run(void)
{
    time_t tv;
    time(&tv); 
    time_t curtime = tv;

    // Seed generation
    m_seed.run();
    m_seed.clear();

    ++m_depth;
    m_pos.collectPatches();

    // Expansion
    const int TIME = 3;
    for (int t = 0; t < TIME; ++t)
    {
        m_expand.run();
        m_filter.run();
        
        updateThreshold();

        std::cout << "STATUS: ";
        for (int i = 0; i < (int)m_optim.m_status.size(); ++i)
        {
            std::cout << m_optim.m_status[i] << ' ';
            if (i % 10 == 9) std::cout << std::endl;
        }
        std::cout << std::endl;

        ++m_depth;
    }
    time(&tv);
    std::cerr << "---- Total: " << (tv - curtime)/CLOCKS_PER_SEC << " secs ----" << std::endl;
}

void CfindMatch::write(const std::string prefix)
{
    m_pos.writePatches2(prefix);
}
