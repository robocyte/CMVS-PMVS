#pragma once

#include <deque>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include <queue>
#include <string>
#include <vector>

#include "patch.h"
#include "../image/photoSetS.h"
#include "patchOrganizerS.h"
#include "seed.h"
#include "expand.h"
#include "filter.h"
#include "optim.h"
#include "option.h"

namespace PMVS3
{

class CfindMatch
{
public:
    CfindMatch();
    virtual ~CfindMatch() {}

    void init(const PMVS3::Soption& option);
    void run(void);
    void write(const std::string prefix);

    int insideBimages(const Vec4f& coord) const;

    int isNeighborRadius(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float hunit, const float neighborThreshold, const float radius) const;
    int isNeighbor(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float hunit, const float neighborThreshold) const;
    int isNeighbor(const Patch::Cpatch& lhs, const Patch::Cpatch& rhs, const float neighborThreshold) const;

    int m_CPU;
    int m_tnum;                                 // num of target images
    int m_num;                                  // num of total images

    std::vector<int> m_timages;                 // target images
    std::vector<int> m_oimages;                 // other images where patches are not computed
    std::vector<int> m_images;                  // total images

    std::string m_prefix;
    int m_level;
    int m_csize;                                // cellsize
    int m_wsize;                                // windows size
    float m_setEdge;                            // use edge detection or not

    std::vector<int>                m_bindexes; // bounding images
    std::vector<std::vector<int>>   m_visdata;  // visdata from SfM. m_num x m_num matrix
    std::vector<std::vector<int>>   m_visdata2; // an array of relavant images

    int m_tau;                                  // Maximum number of images used in the optimization
    int m_depth = 0;                            // If patches are dense or not, that is, if we use check(patch) after patch optimization

    float m_quadThreshold;                      // Threshold on filterQuad
    float m_nccThreshold;
    int   m_minImageNumThreshold;
    int   m_sequenceThreshold;
    float m_angleThreshold0;                    // For first feature matching. Images within this angle are used in matching.
    float m_angleThreshold1;                    // tigher angle
    int   m_countThreshold0 = 2;                // Number of success generation from each seed point
    int   m_countThreshold1 = 4;                // Number of counts, expansion can be tried
    int   m_countThreshold2 = 2;                // Number of trials for each cell in seed
    float m_neighborThreshold = 0.5f;           // Parameter for isNeighbor in findemptyblocks
    float m_neighborThreshold1 = 1.0f;          // Parameter for isNeighbor in filterOutside
    float m_neighborThreshold2 = 1.0f;          // Parameter for filterNeighbor
    float m_nccThresholdBefore;                 // ncc threshold before optim
    float m_maxAngleThreshold;                  // Maximum angle of images must be at least as large as this
    float m_visibleThreshold = 0.0f;
    float m_visibleThresholdLoose = 0.0f;
    float m_epThreshold = 2.0f;                 // Maximum angle of images must be at least as large as this

    std::mutex                      m_lock;     // General lock
    std::deque<std::mutex>          m_imageLocks;
    std::deque<std::mutex>          m_countLocks;

    int                             m_count;
    int                             m_junit = 100;
    std::list<int>                  m_jobs;

    Image::CphotoSetS               m_pss;
    CpatchOrganizerS                m_pos;
    Cseed                           m_seed;
    Cexpand                         m_expand;

public:
    Cfilter                         m_filter;
    Coptim                          m_optim;

    int m_debug = 0;

protected:
    void init(void);
    void initTargets(void);
    void updateThreshold(void);
    void initImages(void);
};

};