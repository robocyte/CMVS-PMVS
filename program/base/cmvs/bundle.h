#pragma once

#include <vector>
#include <string>
#include <list>
#include <mutex>

#include <boost/pending/disjoint_sets.hpp>

#include "../stann/sfcnn.hpp"
#include "../image/photoSetS.h"

namespace CMVS
{

struct Sadd
{
    Sadd(const int image, const float gain)
        : m_image(image)
        , m_gain(gain)
    {
    };

    int   m_image;
    float m_gain;
};

struct Ssfm2
{
    Ssfm2() = default;

    int m_cluster = -1;                 // which cluster it belongs to currently
    float m_score = -2.0f;
    float m_scoreThreshold = -1.0f;

    // If SFM is satisfied or not.
    // 1: satisfied, 0: not satisfied
    //
    // In adding images,
    // 2: currently not satisfied, 1: satisfied,
    // 0: not satisfied from the begining and no hope
    char m_satisfied = 1;
  
    // For SfM point that has not bee satisfied, compute several number of images that can be added to gain more info
    std::vector<Sadd> m_adds;

    std::vector<int> m_uimages;         // best images
};

class Cbundle
{
public:
    Cbundle() = default;
    virtual ~Cbundle() {}

    void run(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU);

    std::string m_prefix;                       // Root dir

    int m_cnum;                                 // # of cameras
    int m_pnum;                                 // # of points

    // Point params
    std::vector<Vec4f> m_coords;
    std::vector<std::vector<int>> m_visibles;

    std::vector<Vec3f> m_colors;
    std::vector<std::vector<int>> m_vpoints;    // A set of point ids visible in each camera

    std::vector<int> m_pweights;

    std::vector<std::vector<int>> m_neighbors;  // A list of connected images, for each camera.

    // Width and height of depth map
    std::vector<int> m_widths;
    std::vector<int> m_heights;

    std::vector<int> m_levels; // scale

    std::vector<std::vector<int> > m_timages;   // m_timages need to be sorted for addImages.
    std::vector<std::vector<int> > m_oimages;

protected:
    void prep(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU);
    void prep2(void);

    void readBundle(const std::string file);
    void setWidthsHeightsLevels(void);
    void setNeighbors(void);

    int totalNum(void) const;

    void setScoreThresholds(void);

    void resetVisibles(void);

    void setNewImages(const int pid, const int rimage, std::vector<int>& newimages); // set new images without image while taking into account m_removed
    void sRemoveImages(void);
    void checkImage(const int image);

    void setCluster(const int p);
    void setClusters(void); // For unsatisfied sfm points, update cluster
    void setScoresClusters(void);

    void slimNeighborsSetLinks(void);

    float computeLink(const int image0, const int image1);

    void addImagesP(void);
    int addImages(void);
    int addImagesSub(const std::vector<std::map<int, float> >& cands);

    static float angleScore(const Vec4f& ray0, const Vec4f& ray1);

    void mergeSfMP(void);
    void mergeSfMPThread(void);

    std::vector<char> m_merged;

    void findPNeighbors(sfcnn<const float*, 3, float>& tree, const int pid, std::vector<int>& pneighbors);

    void resetPoints(void);

    static void mymerge(const std::vector<int>& lhs, const std::vector<int>& rhs, std::vector<int>& output);

    static int my_isIntersect(const std::vector<int>& lhs, const std::vector<int>& rhs);
  
    // Cluster images
    void setTimages(void);
    void divideImages(const std::vector<int>& lhs, std::vector<std::vector<int> >& rhs);

    float computeScore2(const Vec4f& coord, const std::vector<int>& images) const;
    float computeScore2(const Vec4f& coord, const std::vector<int>& images, std::vector<int>& uimages) const;
    // Enforce the specified image to be inside
    float computeScore2(const Vec4f& coord, const std::vector<int>& images, const int index) const;

    void writeCameraCenters(void);
    void writeVis(void);
    void writeGroups(void);

    std::vector<std::vector<float>> m_links; // Link info

    std::vector<int> m_removed; // Removed or not for images

    // Number of SFM points above threshold for each image. We can remove images until m_allows is non-negative for all the images.
    std::vector<int> m_allows; // Used in removing images
    std::vector<int> m_lacks;  // Used in adding images

    // For an image, how sfm point (m_vpoints) changes if the image is removed.
    //  0: unsatisfy
    //  1: satisfy->satisfy
    //  2: satisfy->unsatisfy
    std::vector<char> m_statsT;

    int m_imageT; // image under consideration
    int m_lacksT; // The value of lacks

    std::vector<Ssfm2> m_sfms2; // sfm information used in addimages

    int m_tau; // Number of images used in computeScore2
    std::vector<std::vector<std::vector<int>>> m_ufsT; // union find operations to be executed
    std::vector<float> m_minScales;

    std::vector<int> m_addnums;

    float m_dscale;                 // scaling factor for depth
    float m_dscale2;                // scaling for kdtree version

    Image::CphotoSetS m_pss;

    int m_dlevel;                   // depth level
    int m_maxLevel;                 // maxLevel in m_pss

    int m_imageThreshold;
    int m_pnumThreshold;            // Num of points for images to be connected
    float m_linkThreshold;          // link threshold for neighbor
    float m_scoreRatioThreshold;    // Score ratio threshold. Optimal score using all the visible images times this threshold is the mimimum possible score to be satisfied.
    float m_coverageThreshold;      // How much SFM must be satisfied in each image.

    boost::disjoint_sets_with_storage<>* m_puf = nullptr; // union find for sfm points

    sfcnn<const float*, 3, float>* m_ptree = nullptr;

    // Threads
    int             m_CPU = 8;
    std::mutex      m_lock;
    std::list<int>  m_jobs;
    int             m_junit = 100;
    int             m_thread;
    int             m_count;

    int m_debug = 0;

    void startTimer(void);
    time_t curTimer(void);

    time_t m_tv; //PM
    time_t m_curtime;
};

};