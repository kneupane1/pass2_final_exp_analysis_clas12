#ifndef EFF_CORR_H_GUARD
#define EFF_CORR_H_GUARD
#include "constants.hpp"
class EffCorr
{
protected:
    static const int CD_SEC = 3;
    static const int FD_SEC = 6;

    const float cd_phi_bin_ranges[CD_SEC + 1] = {0, 120, 240, 360};
    const float fd_phi_bin_ranges[FD_SEC + 1] = {0, 60, 120, 180, 240, 300, 360};

    //////////// theta maps
    const float ProtThMap[3] = {-0.00006, 0.00405, 0.95412};
    const float PipThMap[4] = {0.0000006, -0.0000924, 0.0041458, 0.9670324};
    // const float PipThMap[4] = {0, 0, 1, 0};

    //////////// phi maps
    const float ProtPhiMap[2][3] = {{2.24862356e-07, -1.11751535e-04, 1.01446988e+00}, {0, 0, 1.0}};
    const float PipPhiMap[2][3] = {{4.30687779e-07, -1.90534394e-04, 1.02032149e+00}, {0, 0, 1.0}};
    /// // these are from the mapping of mes/miss momentum, and fitting them
    //// a, b ,c for ax2 + bx+c
    const float CDProtMap[3][3] = {{0.02395, -0.04266, 0.99411}, {0.02328, -0.03952, 1.00249}, {0.01712, -0.00029, 0.95352}};

    const float FDProtMap[6][3] = {{-0.00321, 0.02041, 0.96814},
                                   {-0.00145, 0.01266, 0.97548},
                                   {-0.00269, 0.01893, 0.96422},
                                   {-0.00061, 0.00929, 0.97292},
                                   {-0.00129, 0.01430, 0.96118},
                                   {-0.00264, 0.02051, 0.95482}};

    /// // coefficients are from the fit of efficiencies..
    const float CDProtCoef[3][3] = {{0.20052, -0.79964, 1.38699},
                                    {0.16842, -0.64970, 1.31246},
                                    {0.18845, -0.75924, 1.41677}};

    const float FDProtCoef[6][3] = {{-0.00136, 0.01822, 0.92778},
                                    {-0.01176, 0.07598, 0.87691},
                                    {-0.00186, 0.03132, 0.92515},
                                    {-0.00447, 0.03279, 0.93003},
                                    {0.00483, -0.02757, 1.02468},
                                    {-0.00658, 0.04342, 0.93986}};

    const float CDPipMap[3][3] = {{-0.03214, 0.08964, 0.92059},
                                  {-0.04959, 0.13683, 0.89443},
                                  {-0.04007, 0.13766, 0.87928}};

    const float FDPipMap[6][3] = {{-0.01193, 0.05671, 0.94107},
                                  {-0.00850, 0.04254, 0.95453},
                                  {-0.00633, 0.03247, 0.96096},
                                  {-0.00663, 0.03284, 0.95924},
                                  {-0.00958, 0.04797, 0.93805},
                                  {-0.01003, 0.05244, 0.93133}};

    const float CDPipCoef[3][3] = {{-0.14107, 0.19469, 0.66807},
                                   {-0.19271, 0.32662, 0.67771},
                                   {-0.19560, 0.39853, 0.61896}};

    const float FDPipCoef[6][4] = {{-0.06346, 0.30654, 0.63456},
                                   {-0.09727, 0.40524, 0.58949},
                                   {0.06756, -0.00928, 0.79671},
                                   {0.08739, -0.10334, 0.91378},
                                   {-0.08930, 0.35192, 0.64703},
                                   {-0.01981, 0.24775, 0.61119}};

    const float FDPipCoef1[6][2] = {{-0.02222, 1.03444},
                                    {-0.04444, 1.09889},
                                    {0.01111, 0.99778},
                                    {-0.03333, 1.10667},
                                    {-0.01111, 1.00222},
                                    {-0.02222, 1.05444}};

public:
    EffCorr() {};
    ~EffCorr();

    float PROT_EFF_CORR_FACT(float mom_, float theta_, float phi_);
    float PIP_EFF_CORR_FACT(float mom_, float theta_, float phi_);
    float PIM_EFF_CORR_FACT(float mom_, float theta_, float phi_);

    // float EFF_CORR_FACT(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip, float mom_pim, float theta_pim, float phi_pim);
    float EFF_CORR_FACT1(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip);
};

#endif