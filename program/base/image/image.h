#pragma once

#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>

#include "../numeric/vec3.h"

namespace Image
{

class Cimage
{
public:
    Cimage(void);
    virtual ~Cimage();

    virtual void init(const std::string name, const std::string mname, const int maxLevel = 1);
    virtual void init(const std::string name, const std::string mname, const std::string ename, const int maxLevel = 1);

    void setEdge(const float threshold);

    // access to image/masks
    inline Vec3f getColor(const float fx, const float fy, const int level) const;
    inline Vec3f getColor(const int ix, const int iy, const int level) const;

    inline void setColor(const int ix, const int iy, const int level, const Vec3f& rgb);
  
    inline int getMask(const float fx, const float fy, const int level) const;
    inline int getMask(const int ix, const int iy, const int level) const;

    inline int getEdge(const float fx, const float fy, const int level) const;
    inline int getEdge(const int ix, const int iy, const int level) const;  

    inline int getWidth(const int level = 0) const;
    inline int getHeight(const int level = 0) const;

    inline const std::vector<unsigned char>& getImage(const int level) const;
    inline const std::vector<unsigned char>& getMask(const int level) const;
    inline const std::vector<unsigned char>& getEdge(const int level) const;
    inline std::vector<unsigned char>& getImage(const int level);
    inline std::vector<unsigned char>& getMask(const int level);
    inline std::vector<unsigned char>& getEdge(const int level);

    inline int isSafe(const Vec3f& icoord, const int level) const;
    inline int isMask(void) const;        // Check if a mask image exists
    inline int isEdge(void) const;        // Check if an edge image exists

    void alloc(const int fast = 0, const int filter = 0); // This function is also called when you call getColor/getMask when the memory is not allocated
    void free(void);
    void free(const int freeLevel);   // free memory below the specified level

    static int readPBMImage(const std::string file, std::vector<unsigned char>& image, int& width, int& height, const int fast);
    static int writePBMImage(const std::string file, std::vector<unsigned char>& image, int& width, int& height, const int fast);

    static int readPGMImage(const std::string file, std::vector<unsigned char>& image, int& width, int& height, const int fast);
    static int writePGMImage(const std::string file, const std::vector<unsigned char>& image, const int width, const int height);

    static int readPPMImage(const std::string file, std::vector<unsigned char>& image, int& width, int& height, const int fast);
    static int writePPMImage(const std::string file, const std::vector<unsigned char>& image, const int width, const int height);

    static int readJpegImage(const std::string file, std::vector<unsigned char>& image, int& width, int& height, const int fast);

    static void writeJpegImage(const std::string filename, const std::vector<unsigned char>& buffer, const int width, const int height, const int flip = 0);

    static void rgb2hsv(const float r, const float g, const float b, float& hr, float& sr, float& vr);
    static void rgb2hsv(const Vec3f& rgb, Vec3f& hsv);
    static void rgb2hsv(const Vec3f& rgb, float& hr, float& sr, float& vr);
    static void rgb2hsv(const float r, const float g, const float b, Vec3f& hsv);  

    static void rgb2hs(const float r, const float g, const float b, float& h, float& s);
    static void rgb2hs(const Vec3f& rgb, float& h, float& s);
    static void rgb2hs(const Vec3f& rgb, Vec2f& hs);
    static void rgb2hs(const float r, const float g, const float b, Vec2f& hs);

    static void gray2rgb(const float gray, float& r, float& g, float& b);

    // Create a 1d gaussian filter based on sigma
    static void createFilter(const float sigma, std::vector<float>& filter);

    static void filterG(const std::vector<float>& filter, std::vector<std::vector<float> >& data);
    static void filterG(const std::vector<float>& filter, std::vector<std::vector<float> >& data, std::vector<std::vector<float> >& buffer);
    static void filterG(const std::vector<float>& filter, const int width, const int height, std::vector<float>& data);
    static void filterG(const std::vector<float>& filter, const int width, const int height, std::vector<float>& data, std::vector<float>& buffer);                      

    // non maximum surpression
    static void nms(std::vector<std::vector<float> >& data);
    static void nms(std::vector<std::vector<float> >& data, std::vector<std::vector<float> >& buffer);

protected:
    static void completeName(const std::string& lhs, std::string& rhs, const int color);    // complete the name of an image file

    void buildImageMaskEdge(const int filter);
    void buildImage(const int filter);
    void buildMask(void);
    void buildEdge(void);

    //----------------------------------------------------------------------
    // Variables updated at every alloc/free
    //----------------------------------------------------------------------  
    int m_alloc;                                        // 0: nothing allocated; 1: width/height allocated; 2: memory allocated
    std::vector<std::vector<unsigned char>> m_images;   // a pyramid of images
    std::vector<std::vector<unsigned char>> m_masks;    // a pyramid of masks
    std::vector<std::vector<unsigned char>> m_edges;    // a pyramid of images specifying regions with edges(texture)

    std::vector<int> m_widths;      // width of an image in each level
    std::vector<int> m_heights;     // height of an image in each level

    //----------------------------------------------------------------------
    // Variables keep fixed
    //----------------------------------------------------------------------  
    std::string m_name;       // a name of an image
    std::string m_mname;      // a name of a mask image
    std::string m_ename;      // a name of an image specifying regions with edges(texture)
    int m_maxLevel;           // number of levels
};

inline int Cimage::isSafe(const Vec3f& icoord, const int level) const
{
    if (icoord[0] < 0.0 || m_widths[level] - 2 < icoord[0] || icoord[1] < 0.0 || m_heights[level] - 2 < icoord[1]) return 0;
    else                                                                                                           return 1;
};

inline int Cimage::isMask(void) const
{
    if (m_masks[0].empty()) return 0;
    else                    return 1;
};

inline int Cimage::isEdge(void) const
{
    if (m_edges[0].empty()) return 0;
    else                    return 1;
};

inline const std::vector<unsigned char>& Cimage::getImage(const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_images[level];
};

inline const std::vector<unsigned char>& Cimage::getMask(const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_masks[level];
};

inline const std::vector<unsigned char>& Cimage::getEdge(const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_edges[level];
};

inline std::vector<unsigned char>& Cimage::getImage(const int level)
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_images[level];
};

inline std::vector<unsigned char>& Cimage::getMask(const int level)
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_masks[level];
};

inline std::vector<unsigned char>& Cimage::getEdge(const int level)
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    return m_edges[level];
};

int Cimage::getWidth(const int level) const
{
    if (m_alloc == 0)
    {
        std::cerr << "First allocate (getWidth)" << std::endl;
        exit (1);
    }

    return m_widths[level];
};

int Cimage::getHeight(const int level) const
{
    if (m_alloc == 0)
    {
        std::cerr << "First allocate (getHeight)" << std::endl;
        exit (1);
    }

    return m_heights[level];
};

Vec3f Cimage::getColor(const float x, const float y, const int level) const
{
    // Bilinear case
    const int lx = (int)floor(x);
    const int ly = (int)floor(y);
    const int index = 3 * (ly * m_widths[level] + lx);

    const float dx1 = x - lx;  const float dx0 = 1.0f - dx1;
    const float dy1 = y - ly;  const float dy0 = 1.0f - dy1;

    const float f00 = dx0 * dy0;  const float f01 = dx0 * dy1;
    const float f10 = dx1 * dy0;  const float f11 = dx1 * dy1;
    const int index2 = index + 3 * m_widths[level];

    const unsigned char* ucp0 = &m_images[level][index]  - 1;
    const unsigned char* ucp1 = &m_images[level][index2] - 1;
    float r = 0.0f;  float g = 0.0f;  float b = 0.0f;
    r += *(++ucp0) * f00 + *(++ucp1) * f01;
    g += *(++ucp0) * f00 + *(++ucp1) * f01;
    b += *(++ucp0) * f00 + *(++ucp1) * f01;
    r += *(++ucp0) * f10 + *(++ucp1) * f11;
    g += *(++ucp0) * f10 + *(++ucp1) * f11;
    b += *(++ucp0) * f10 + *(++ucp1) * f11;

    return Vec3f(r, g, b);
};

void Cimage::setColor(const int ix, const int iy, const int level, const Vec3f& rgb)
{
    const int index = (iy * m_widths[level] + ix) * 3;

    m_images[level][index]   = (unsigned char)floor(rgb[0] + 0.5f);
    m_images[level][index+1] = (unsigned char)floor(rgb[1] + 0.5f);
    m_images[level][index+2] = (unsigned char)floor(rgb[2] + 0.5f);
};

Vec3f Cimage::getColor(const int ix, const int iy, const int level) const
{
    const int index = (iy * m_widths[level] + ix) * 3;

    return Vec3f(m_images[level][index], m_images[level][index+1], m_images[level][index+2]);
};

int Cimage::getMask(const float fx, const float fy, const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }    
  
    if (m_masks[level].empty()) return 1;
  
    const int ix = (int)floor(fx + 0.5f);
    const int iy = (int)floor(fy + 0.5f);
    return getMask(ix, iy, level);
};

int Cimage::getMask(const int ix, const int iy, const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    if (m_masks[level].empty()) return 1;

    if (ix < 0 || m_widths[level] <= ix || iy < 0 || m_heights[level] <= iy) return 1;

    const int index = iy * m_widths[level] + ix;
    return m_masks[level][index];
};

int Cimage::getEdge(const float fx, const float fy, const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    if (m_edges[level].empty()) return 1;

    const int ix = (int)floor(fx + 0.5f);
    const int iy = (int)floor(fy + 0.5f);

    return getEdge(ix, iy, level);
};

int Cimage::getEdge(const int ix, const int iy, const int level) const
{
    if (m_alloc != 2)
    {
        std::cerr << "First allocate" << std::endl;
        exit (1);
    }

    if (m_edges[level].empty()) return 1;

    if (ix < 0 || m_widths[level] <= ix || iy < 0 || m_heights[level] <= iy) return 1;

    const int index = iy * m_widths[level] + ix;
    return m_edges[level][index];
};

};