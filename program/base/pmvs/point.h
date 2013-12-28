#pragma once

#include "../numeric/vec4.h"
#include "../numeric/mat3.h"

namespace PMVS3
{

class Cpoint
{
public:
    Cpoint() = default;
    virtual ~Cpoint();

    Vec3f m_icoord;
    Vec4f m_coord;      // 3D coordinate
    float m_response = -1.0f;

    int m_type = 1;     // 0: Harris, 1: DoG
    int m_itmp;         // Temporary variable, used to store original imageid in initial match

    bool operator < (const Cpoint& rhs) const { return m_response < rhs.m_response; }

    friend std::istream& operator >>(std::istream& istr, Cpoint& rhs);
    friend std::ostream& operator <<(std::ostream& ostr, const Cpoint& rhs);
};

bool SortCpoint(const Cpoint& a, const Cpoint& b);

std::istream& operator >>(std::istream& istr, Cpoint& rhs);
std::ostream& operator <<(std::ostream& ostr, const Cpoint& rhs);

};