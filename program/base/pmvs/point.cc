#include <iostream>

#include "point.h"

using namespace PMVS3;

Cpoint::Cpoint(void)
{
    m_response = -1.0;
    m_type = -1;
}

Cpoint::~Cpoint()
{
}

std::istream& PMVS3::operator >>(std::istream& istr, Cpoint& rhs)
{
    std::string header;
    char str[1024];
    istr >> str;
    header = std::string(str);
    istr >> rhs.m_icoord[0] >> rhs.m_icoord[1] >> rhs.m_response >> rhs.m_type;
    rhs.m_icoord[2] = 1.0f;

    return istr;
}

std::ostream& PMVS3::operator <<(std::ostream& ostr, const Cpoint& rhs)
{
    ostr << "POINT0" << std::endl << rhs.m_icoord[0] << ' ' << rhs.m_icoord[1] << ' ' << rhs.m_response << ' ' << rhs.m_type;

    return ostr;
}

bool SortCpoint(const Cpoint& a, const Cpoint& b)
{
    return a.m_response < b.m_response;
}
