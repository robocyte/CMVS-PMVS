#pragma once

#include "patch.h"
#include <queue>

namespace PMVS3
{

class CfindMatch;

class P_compare
{
public:
    bool operator()(const Patch::Ppatch& lhs, const Patch::Ppatch& rhs) const { return lhs->m_tmp < rhs->m_tmp; }
};

class CpatchOrganizerS
{
public:
    CpatchOrganizerS(CfindMatch& findMatch);

    void init(void);
    void collectPatches(const int target = 0);
    void collectPatches(std::priority_queue<Patch::Ppatch, std::vector<Patch::Ppatch>, P_compare>& pqpatches);
    
    void writePatches2(const std::string prefix);
    void writePLY(const std::vector<Patch::Ppatch>& patches, const std::string filename);

    void clearCounts(void);
    void clearFlags(void);

    void setGridsImages(Patch::Cpatch& patch, const std::vector<int>& images) const;
    void addPatch(Patch::Ppatch& ppatch);
    void removePatch(const Patch::Ppatch& ppatch);
    void setGrids(Patch::Ppatch& ppatch) const;
    void setGrids(Patch::Cpatch& patch) const;
    void setVImagesVGrids(Patch::Ppatch& ppatch);
    void setVImagesVGrids(Patch::Cpatch& patch);
    void updateDepthMaps(Patch::Ppatch& ppatch);

    int isVisible(const Patch::Cpatch& patch, const int image, const int& ix, const int& iy, const float strict, const int lock);
    int isVisible0(const Patch::Cpatch& patch, const int image, int& ix, int& iy, const float strict, const int lock);

    void findNeighbors(const Patch::Cpatch& patch, std::vector<Patch::Ppatch>& neighbors, const int lock, const float scale = 1.0f, const int margin = 1, const int skipvis = 0);

    void setScales(Patch::Cpatch& patch) const;

    // Change the contents of m_images from images to indexes
    void image2index(Patch::Cpatch& patch);
    void index2image(Patch::Cpatch& patch);

    std::vector<int>                                     m_gwidths;
    std::vector<int>                                     m_gheights;
    std::vector<std::vector<std::vector<Patch::Ppatch>>> m_pgrids;
    std::vector<std::vector<std::vector<Patch::Ppatch>>> m_vpgrids;
    std::vector<std::vector<Patch::Ppatch>>              m_dpgrids;     // Closest patch
    std::vector<Patch::Ppatch>                           m_ppatches;    // All the patches in the current level of m_pgrids 
    std::vector<std::vector<unsigned char>>              m_counts;       // Check how many times patch optimization was performed for expansion

    static Patch::Ppatch m_MAXDEPTH;
    static Patch::Ppatch m_BACKGROUND;

protected:
    CfindMatch& m_fm;
};

};
