#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "../numeric/vec4.h"

namespace Patch
{

class Cpatch
{
public:
    Cpatch(void) = default;

    // Associated image ids. first image id is the reference one. images can be non-targetting image.
    std::vector<int>        m_images;
    std::vector<TVec2<int>> m_grids;

    // Visible images. m_vimages must be targetting images.
    std::vector<int>        m_vimages;
    std::vector<TVec2<int>> m_vgrids;
  
    inline float score(const float threshold) const  { return std::max(0.0f, m_ncc - threshold) * (int)m_images.size(); }
    inline float score2(const float threshold) const { return std::max(0.0f, m_ncc - threshold) * m_timages; }

    Vec4f m_coord;              // 3D coordinates of the center of the patch
    Vec4f m_normal;             // Patch outward normal vector

    int m_timages = 0;          // Number of targetting images in m_images
    int m_flag;                 // Flat for expansion - 0: not yet tested, 1: done
    int m_id;                   // Id number in m_ppatches

    unsigned char m_dflag = 0;  // For directional flag

    char m_fix = 0;             // Fixed patch or not

    float m_dscale;             // Scaling factor corresponding to one pixel difference
    float m_ascale;
    float m_tmp;
    float m_ncc = -1.0;         // Average ncc
};

typedef std::shared_ptr<Cpatch> Ppatch;

struct Spatchcmp
{
    bool operator()(const Ppatch& lhs, const Ppatch& rhs)
    {
        if (lhs.get() < rhs.get())  return true;
        else                        return false;
    }
};

std::istream& operator >>(std::istream& istr, Patch::Cpatch& rhs);
std::ostream& operator <<(std::ostream& ostr, const Patch::Cpatch& rhs);

};