#pragma once

#include <vector>

#include "patch.h"

#include "unsupported/Eigen/NonLinearOptimization"

namespace PMVS3
{

class CfindMatch;

class Coptim
{
public:
    Coptim(CfindMatch& findMatch);

    void init(void);

    void collectImages(const int index, std::vector<int>& indexes) const;
    void addImages(Patch::Cpatch& patch) const;
    void removeImagesEdge(Patch::Cpatch& patch) const;

    float getUnit(const int index, const Vec4f& coord) const;

    void computeUnits(const Patch::Cpatch& patch, std::vector<int>& indexes, std::vector<float>& fineness, std::vector<Vec4f>& rays) const;
    void computeUnits(const Patch::Cpatch& patch, std::vector<float>& fineness) const;

    int preProcess(Patch::Cpatch& patch, const int id, const int seed);
    int postProcess(Patch::Cpatch& patch, const int id, const int seed);

    bool refinePatchBFGS(Patch::Cpatch& patch, const int id);

    void setRefImage(Patch::Cpatch& patch, const int id);

    int check(Patch::Cpatch& patch);

    std::vector<int> m_status;

protected:
    void filterImagesByAngle(Patch::Cpatch& patch);

    void sortImages(Patch::Cpatch& patch) const;
    void constraintImages(Patch::Cpatch& patch, const float nccThreshold, const int id);

    void setINCCs(const Patch::Cpatch& patch, std::vector<float> & nccs, const std::vector<int>& indexes, const int id, const int robust);
    void setINCCs(const Patch::Cpatch& patch, std::vector<std::vector<float>>& nccs, const std::vector<int>& indexes, const int id, const int robust);

    int grabTex(const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, const int index, const int size, std::vector<float>& tex) const;
    int grabSafe(const int index, const int size, const Vec3f& center, const Vec3f& dx, const Vec3f& dy, const int level) const;

    double computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const Vec4f& pxaxis, const Vec4f& pyaxis, const int id, const int robust);

public:
    static void normalize(std::vector<float>& tex);
    static void normalize(std::vector<std::vector<float>>& texs, const int size);

    float dot(const std::vector<float>& tex0, const std::vector<float>& tex1) const;

    //BFGS
    static void my_f_lm(const Eigen::VectorXd &par, Eigen::VectorXd &fvec, int id);
    static void my_f_lm(const Eigen::VectorXf &par, Eigen::VectorXf &fvec, int id);

    void encode(const Vec4f& coord, double* const vect, const int id) const;
    void encode(const Vec4f& coord, const Vec4f& normal, double* const vect, const int id) const;
    void decode(Vec4f& coord, Vec4f& normal, const double* const vect, const int id) const;
    void decode(Vec4f& coord, const double* const vect, const int id) const;

    void setWeightsT(const Patch::Cpatch& patch, const int id);

    void getPAxes(const int index, const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis) const;

    double computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const int id, const int robust);
    static inline float robustincc(const float rhs) { return rhs / (1 + 3 * rhs); }
    static inline float unrobustincc(const float rhs) { return rhs / (1 - 3 * rhs); }

protected:
    void setAxesScales(void);

    static Coptim* m_one;
    CfindMatch& m_fm;

    std::vector<Vec3f> m_xaxes;
    std::vector<Vec3f> m_yaxes;
    std::vector<Vec3f> m_zaxes;
    std::vector<float> m_ipscales;

    std::vector<float> m_vect0T;  
    std::vector<Vec4f> m_centersT;
    std::vector<Vec4f> m_raysT;
    std::vector<std::vector<int> > m_indexesT;
    std::vector<float> m_dscalesT;
    std::vector<float> m_ascalesT;

    std::vector<Vec3f>                              m_paramsT;      // stores current parameters for derivative computation
    std::vector<std::vector<std::vector<float>>>    m_texsT;        // Grabbed texture, last is 7x7x3 patch
    std::vector<std::vector<float>>                 m_weightsT;     // weights for refineDepthOrientationWeighed
    std::vector<std::vector<double>>                m_worksT;       // Working array for levmar
};

};