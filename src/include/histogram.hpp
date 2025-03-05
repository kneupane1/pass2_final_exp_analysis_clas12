
#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"
#include <mutex>

using namespace std;

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram
{
protected:
    std::shared_ptr<TFile> RootOutputFile;
    std::shared_ptr<TCanvas> def;

    int bins = 500;
    double p_min = 0.0;
    double p_max = 6.0;
    double Dt_max = 10.0;
    double Dt_min = -Dt_max;
    double q2_min = 1.0;
    double q2_max = 10.0;

    double w_max = 2.5;
    double w_min = 1.0;

    double zero = 0.0;

    static const short particle_num = 3; // 0-e 1-Pi 2-P 3-K
    std::string particle_name[particle_num] = {"e", "pi", "P"};
    static const short charge_num = 2; // 0-pos 1-neg
    std::string charge_name[charge_num] = {"positive", "negative"};
    static const short with_id_num = 3; // 0-without 1-with 2-anti
    std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

    static const short num_sectors = 6;
    std::string sec_name[num_sectors] = {"1", "2", "3", "4", "5", "6"};

    // std::string cut_name[2] = {"bc", "ac"};

    static const short CUTS = 4;
    enum cuts
    {
        before_any_cuts,
        with_one_cut,
        outside_one_cut,
        after_all_cuts
    };
    std::mutex mutex;

    static const short w_range_num = 3;
    std::string w_range_name[w_range_num] = {
        "all_W_range", " W < 2.5 ",
        " W > 2.5 "}; //{" W<2.0 ", " 2.0<W<2.5 ", " 2.5<W<3.0 ", " 3.0<W<3.5 "};
    static const short q2_range_num = 3;
    std::string q2_range_name[q2_range_num] = /*{" Q2<1.0 ",     " 1.0<Q2<2.0 ",
                                               " 2.0<Q2<3.0 ", " 3.0<Q2<4.0 ",
                                               " 4.0<Q2<5.0 ", " 5.0<Q2<6.0 ",
                                               " 6.0<Q2<7.0 ", " 7.0<Q2<8.0 ",
                                               " 8.0<Q2<9.0 ", " Q2>9.0 "};*/
        {"all_Q2_range ", " Q2 < 4.5 ", " Q2 > 4.5 "};

    static const short inv_Ppip_range_num = 2;
    std::string inv_Ppip_range_name[inv_Ppip_range_num] = {" M[Ppip]<1.5 GeV ",
                                                           " M[Ppip]>1.5 GeV "};

    static const short inv_pip_pim_range_num = 2;
    std::string inv_pip_pim_range_name[inv_pip_pim_range_num] = {
        " M[pi+pi-]<0.9 GeV ", " M[pi+pi-]>0.9 GeV "};

    static const short theta_pim_range_num = 2;
    std::string theta_pim_range_name[theta_pim_range_num] = {
        " theta[pi-]<90 deg ", " theta[pi-]>90 deg "};

    static const short phi_pim_range_num = 2;
    std::string phi_pim_range_name[phi_pim_range_num] = {" phi[pi-]<90 deg ",
                                                         " phi[pi-]>90 deg "};
    static const short alpha_pim_range_num = 2;
    std::string alpha_pim_range_name[alpha_pim_range_num] = {
        " alpha[pi-]<180 deg ", " alpha[pi-]>180 deg "};

    static const short NUM_CONDITIONS = 1;
    std::string NUM_CONDITIONS_NAME[NUM_CONDITIONS] = {
        /*"twoPi_event" ,*/ "missingPim events"};
    // Kinematics

    static const short th_bin_size = 3;
    float th_low_values[3] = {5, 25, 37};
    float th_up_values[3] = {25, 37, 60};

    static const short mom_bin_size = 9;
    float mom_low_values[3][9] = {{0, 0.4, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0}, {0, 0.4, 0.8, 1.1, 1.3, 1.5, 1.75, 2.0, 2.4}, {0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.2, 2.5}};
    float mom_up_values[3][9] = {{0.4, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0}, {0.4, 0.8, 1.1, 1.3, 1.5, 1.75, 2.0, 2.4, 4.0}, {0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.2, 2.5, 3.0}};

    static const short phi_bin_size = 3;
    float phi_low_values[3] = {0, 12, 240};
    float phi_up_values[3] = {120, 240, 360};

    // //////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////

    static const short q2_bin = 11;
    float q2_low_values[10] = {1.0, 2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0};
    float q2_up_values[10] = {2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0, 9.0};
    int q2_bin_size = 10;
    int w_lower_bin = 8;
    int w_higher_bin = 24;

    // // // ////////////// backgraound multiplication factors obtained from AND logic in the background fitting using exclusive topology data /////////
    // float background_fact[2][16][9] = {{{3.50, 3.50, 3.50, 3.32, 2.58, 2.53, 2.60, 1.00, 2.00},
    //                                     {3.50, 3.50, 3.22, 3.19, 2.78, 1.78, 1.33, 1.33, 2.00},
    //                                     {3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 2.00},
    //                                     {3.50, 3.50, 3.50, 3.50, 2.93, 3.50, 3.43, 2.00, 2.00},
    //                                     {3.50, 3.50, 3.50, 3.50, 2.94, 2.62, 2.46, 3.00, 2.75},
    //                                     {3.50, 3.50, 3.50, 3.50, 3.18, 3.37, 3.50, 1.56, 2.87},
    //                                     {3.04, 3.50, 3.20, 2.86, 2.25, 2.34, 3.13, 2.00, 2.00},
    //                                     {3.38, 3.50, 3.13, 3.10, 2.86, 3.50, 1.80, 3.50, 2.00},
    //                                     {3.50, 3.50, 3.50, 3.02, 2.97, 3.01, 2.78, 2.89, 2.33},
    //                                     {3.32, 3.30, 3.14, 3.06, 2.77, 2.76, 3.29, 2.74, 2.50},
    //                                     {3.42, 3.20, 3.40, 3.07, 3.27, 3.14, 2.42, 2.17, 2.58},
    //                                     {3.50, 3.21, 3.10, 3.20, 3.43, 2.97, 3.50, 3.06, 2.00},
    //                                     {3.50, 3.36, 3.40, 3.07, 2.88, 2.77, 3.46, 3.50, 2.00},
    //                                     {3.50, 3.50, 3.29, 3.09, 3.43, 3.21, 3.21, 2.78, 2.00},
    //                                     {3.47, 3.32, 3.17, 3.25, 3.25, 3.50, 3.50, 3.45, 2.00}},
    //                                    {{2.69, 3.14, 3.45, 3.42, 3.50, 3.50, 3.50, 2.50, 2.50},
    //                                     {3.50, 3.50, 3.40, 3.50, 3.47, 3.50, 3.50, 3.50, 2.50},
    //                                     {3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 2.50},
    //                                     {3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 2.50},
    //                                     {3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 2.50},
    //                                     {2.92, 2.62, 2.61, 2.44, 2.41, 2.61, 2.69, 2.50, 2.50},
    //                                     {2.78, 2.69, 2.61, 2.55, 2.46, 2.62, 2.78, 2.91, 2.50},
    //                                     {2.44, 2.26, 2.22, 2.18, 2.05, 2.22, 2.26, 2.50, 2.50},
    //                                     {2.30, 2.53, 2.17, 2.26, 2.19, 2.22, 2.32, 2.50, 2.50},
    //                                     {2.55, 2.24, 2.13, 2.25, 2.16, 2.27, 2.31, 2.51, 2.50},
    //                                     {2.49, 2.45, 2.36, 2.26, 2.34, 2.40, 2.52, 2.65, 2.50},
    //                                     {2.39, 2.53, 2.38, 2.46, 2.37, 2.49, 2.60, 2.80, 2.50},
    //                                     {2.45, 2.62, 2.51, 2.49, 2.49, 2.59, 2.74, 2.83, 2.50},
    //                                     {2.80, 2.57, 2.61, 2.55, 2.53, 2.63, 2.76, 2.98, 2.50},
    //                                     {2.50, 2.70, 2.49, 2.53, 2.57, 2.74, 2.87, 3.08, 2.50}}};

    // // ////////////// backgraound multiplication factors obtained from AND logic in the background fitting using exclusive topology data /////////
    float background_fact[15][9] = {{7.06, 5.34, 4.89, 3.27, 2.44, 1.96, 1.70, 1.36, 0.75},
                                    {6.74, 5.02, 3.74, 4.03, 2.98, 2.14, 1.78, 1.80, 1.00},
                                    {6.77, 5.90, 4.39, 3.86, 2.12, 2.17, 1.55, 1.30, 0.50},
                                    {6.42, 5.44, 5.02, 3.73, 2.73, 1.81, 1.84, 1.00, 1.50},
                                    {5.80, 5.49, 4.35, 3.91, 2.64, 2.19, 1.87, 1.38, 1.57},
                                    {4.48, 4.20, 3.97, 3.35, 3.27, 3.57, 1.99, 2.10, 1.21},
                                    {4.39, 3.95, 3.79, 3.54, 3.33, 2.84, 2.86, 2.13, 1.36},
                                    {4.25, 3.95, 3.81, 3.68, 3.31, 3.47, 2.99, 2.79, 1.48},
                                    {3.96, 3.76, 3.62, 3.43, 3.27, 3.32, 3.26, 2.78, 1.58},
                                    {3.84, 3.68, 3.45, 3.42, 3.42, 3.27, 3.23, 3.12, 2.00},
                                    {3.89, 3.65, 3.48, 3.33, 3.44, 3.31, 3.38, 2.99, 1.75},
                                    {3.93, 3.61, 3.53, 3.35, 3.37, 3.35, 3.55, 3.28, 1.57},
                                    {3.89, 3.70, 3.64, 3.52, 3.27, 3.35, 3.64, 3.43, 1.25},
                                    {3.93, 3.80, 3.58, 3.56, 3.57, 3.43, 3.75, 3.18, 1.00},
                                    {3.85, 3.79, 3.52, 3.58, 3.34, 3.58, 3.83, 3.45, 1.00}};

    float background_fact_loose[15][9] = {{8.50, 4.71, 4.20, 3.12, 2.38, 2.21, 1.67, 1.00, 1.00},
                                          {6.32, 5.26, 3.11, 3.63, 2.91, 2.15, 1.07, 1.71, 1.00},
                                          {6.01, 6.91, 3.98, 4.79, 1.83, 2.01, 1.35, 1.20, 0.33},
                                          {5.62, 5.22, 4.53, 3.96, 2.58, 1.50, 2.14, 1.00, 1.50},
                                          {7.03, 5.43, 4.56, 4.03, 2.54, 2.30, 1.32, 1.60, 2.00},
                                          {4.92, 5.27, 4.41, 3.41, 3.34, 2.93, 1.39, 1.62, 1.50},
                                          {4.09, 4.32, 3.61, 3.31, 2.88, 2.53, 2.14, 2.12, 1.00},
                                          {3.93, 3.59, 3.42, 3.35, 3.17, 3.61, 2.64, 2.96, 1.42},
                                          {3.81, 3.54, 3.61, 3.09, 3.00, 3.34, 3.08, 2.89, 1.08},
                                          {3.54, 3.49, 3.39, 3.21, 2.99, 2.95, 2.95, 2.96, 1.67},
                                          {3.67, 3.47, 3.35, 3.20, 3.12, 2.92, 2.95, 2.64, 1.78},
                                          {3.77, 3.49, 3.32, 3.30, 3.31, 3.06, 3.41, 3.26, 1.50},
                                          {3.67, 3.56, 3.50, 3.31, 3.01, 3.13, 3.68, 3.24, 1.00},
                                          {3.89, 3.67, 3.39, 3.35, 3.52, 3.46, 3.49, 2.97, 1.00},
                                          {3.70, 3.65, 3.34, 3.42, 3.21, 3.71, 3.75, 3.49, 1.00}};

    float background_fact_tight[15][9] = {{6.23, 5.24, 4.85, 3.63, 2.61, 2.05, 1.94, 1.35, 0.75},
                                          {7.02, 5.39, 3.97, 3.97, 3.19, 2.21, 1.89, 1.92, 2.25},
                                          {6.43, 5.94, 4.65, 3.76, 2.66, 2.27, 1.75, 1.33, 0.92},
                                          {5.87, 5.24, 4.92, 3.43, 2.71, 2.34, 1.96, 1.66, 1.50},
                                          {5.09, 5.06, 4.10, 3.92, 3.24, 2.70, 2.26, 1.72, 1.75},
                                          {4.66, 4.70, 3.99, 4.09, 3.66, 3.79, 2.78, 2.00, 1.48},
                                          {4.59, 4.28, 4.14, 3.87, 3.56, 3.31, 3.07, 2.76, 1.95},
                                          {4.60, 4.11, 3.98, 3.87, 3.54, 3.80, 3.41, 3.11, 1.76},
                                          {4.19, 4.06, 3.91, 3.64, 3.66, 3.74, 3.62, 3.26, 1.72},
                                          {4.14, 3.94, 3.76, 3.66, 3.69, 3.51, 3.34, 3.24, 2.10},
                                          {4.09, 3.93, 3.72, 3.57, 3.66, 3.61, 3.63, 3.27, 1.83},
                                          {4.17, 3.85, 3.80, 3.61, 3.54, 3.53, 3.77, 3.45, 1.44},
                                          {4.15, 3.91, 3.84, 3.72, 3.52, 3.56, 3.68, 3.68, 0.89},
                                          {4.16, 4.00, 3.80, 3.77, 3.75, 3.64, 3.98, 3.51, 1.00},
                                          {4.01, 3.96, 3.75, 3.75, 3.62, 3.78, 3.93, 3.77, 1.00}};

    ////////////  exp: low: {excl, mPip, mProt}, high:{excl, mPip, mProt},) ; sim: low: {excl, mPip, mProt}, high: {excl, mPip, mProt}
    // float mmsq_low_values_for_bkg[2][2][3] = {{{-0.004, -0.028, 0.763}, {0.002, 0.071, 1.003}}, {{-0.004, -0.024, 0.79}, {0.002, 0.079, 1.025}}};

    int bin_val = -1;

    // Define q2 bin ranges
    std::vector<float> q2_bins = {2.0, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0, 9.0};

    // Find the bin based on q2

    int q2_bining(float q2)
    {

        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                bin_val = i;
                break;
            }
        }
        return bin_val;
    }

    // /// MMSQ cuts
    int q2_bin_val = -1;
    // /////////////////////////// exp, sim data mmsq cuts [Q2][up/down][a,b,c]   /////////////////  updated nov 2024 with 0.65 sm * w dependent sm
    // double mmsq_cuts[2][9][2][3] = {
    //     {{{-0.1330, 0.5196, -0.4026},
    //       {0.0719, -0.3014, 0.2581}},
    //      {{-0.1176, 0.4645, -0.3538},
    //       {0.0410, -0.1844, 0.1470}},
    //      {{-0.1461, 0.5680, -0.4461},
    //       {0.0863, -0.3422, 0.2793}},
    //      {{-0.1028, 0.4131, -0.3115},
    //       {0.0642, -0.2602, 0.2022}},
    //      {{-0.0591, 0.2539, -0.1706},
    //       {0.0615, -0.2450, 0.1804}},
    //      {{-0.0358, 0.1777, -0.1085},
    //       {0.0511, -0.2008, 0.1339}},
    //      {{-0.0256, 0.1207, -0.0358},
    //       {0.0627, -0.2115, 0.1135}},
    //      {{-0.0500, 0.2273, -0.1382},
    //       {-0.0289, 0.1211, -0.1872}},
    //      {{-0.0822, 0.3223, -0.2004},
    //       {0.1153, -0.3838, 0.2522}}},
    //     {{{-0.1648, 0.6717, -0.5615},
    //       {0.0845, -0.3438, 0.2921}},
    //      {{-0.1332, 0.5590, -0.4635},
    //       {0.0880, -0.3502, 0.2887}},
    //      {{-0.1135, 0.4879, -0.3995},
    //       {0.0647, -0.2682, 0.2166}},
    //      {{-0.1316, 0.5516, -0.4541},
    //       {0.0614, -0.2510, 0.1951}},
    //      {{-0.1120, 0.4788, -0.3874},
    //       {0.0496, -0.1989, 0.1396}},
    //      {{-0.0901, 0.3937, -0.3071},
    //       {0.0307, -0.1256, 0.0684}},
    //      {{-0.0512, 0.2545, -0.1847},
    //       {0.0150, -0.0601, 0.0004}},
    //      {{-0.0414, 0.2194, -0.1526},
    //       {0.0160, -0.0525, -0.0182}},
    //      {{-0.0117, 0.1203, -0.0703},
    //       {-0.0040, 0.0208, -0.0896}}}};

    /// MMSQ 3 SIGMA mid cuts //////////////////////////////////////
    double mmsq_cuts[2][9][2][3] =
        {{{{-0.1036, 0.4002, -0.2877}, {0.0207, -0.1156, 0.0831}},
          {{-0.1261, 0.4775, -0.3517}, {0.0693, -0.2855, 0.2266}},
          {{-0.1368, 0.5127, -0.3829}, {0.0784, -0.3141, 0.2477}},
          {{-0.0655, 0.2737, -0.1872}, {0.0569, -0.2435, 0.1881}},
          {{-0.0669, 0.2756, -0.1876}, {0.0286, -0.1445, 0.1002}},
          {{-0.0785, 0.3229, -0.2386}, {0.1088, -0.4230, 0.3420}},
          {{-0.0461, 0.1819, -0.0921}, {0.0514, -0.1975, 0.1264}},
          {{-0.0584, 0.2435, -0.1596}, {0.0227, -0.0900, 0.0309}},
          {{-0.0998, 0.3329, -0.2000}, {0.0623, -0.1873, 0.0918}}},

         {{{-0.1424, 0.5468, -0.4204}, {0.0517, -0.2097, 0.1605}},
          {{-0.1435, 0.5510, -0.4240}, {0.0660, -0.2556, 0.1940}},
          {{-0.1352, 0.5214, -0.3979}, {0.0577, -0.2257, 0.1666}},
          {{-0.1148, 0.4494, -0.3357}, {0.0458, -0.1825, 0.1278}},
          {{-0.0913, 0.3686, -0.2679}, {0.0264, -0.1129, 0.0667}},
          {{-0.0799, 0.3270, -0.2319}, {0.0152, -0.0692, 0.0245}},
          {{-0.0662, 0.2791, -0.1911}, {0.0043, -0.0263, -0.0179}},
          {{-0.0652, 0.2780, -0.1918}, {-0.0011, -0.0024, -0.0438}},
          {{-0.0800, 0.3320, -0.2403}, {-0.0762, 0.2432, -0.2462}}}};

    ///////// tight 2.5 sigma ////////////
    double mmsq_cuts_tight[2][9][2][3] =
        {
            {{{-0.0810, 0.3161, -0.2230}, {0.0122, -0.0786, 0.0571}},
             {{-0.1021, 0.3882, -0.2822}, {0.0566, -0.2341, 0.1885}},
             {{-0.1037, 0.3932, -0.2885}, {0.0578, -0.2361, 0.1876}},
             {{-0.0613, 0.2503, -0.1720}, {0.0569, -0.2342, 0.1844}},
             {{-0.0801, 0.3102, -0.2204}, {0.0521, -0.2131, 0.1607}},
             {{-0.0701, 0.2846, -0.2101}, {0.0937, -0.3625, 0.2950}},
             {{-0.0301, 0.1236, -0.0517}, {0.0499, -0.1883, 0.1269}},
             {{-0.0379, 0.1693, -0.1050}, {0.0040, -0.0219, -0.0187}},
             {{-0.0695, 0.2339, -0.1302}, {0.0460, -0.1347, 0.0599}}},

            {{{-0.1167, 0.4519, -0.3456}, {0.0416, -0.1668, 0.1288}},
             {{-0.1137, 0.4421, -0.3379}, {0.0546, -0.2087, 0.1594}},
             {{-0.1100, 0.4283, -0.3251}, {0.0556, -0.2107, 0.1590}},
             {{-0.0940, 0.3718, -0.2762}, {0.0503, -0.1904, 0.1399}},
             {{-0.0833, 0.3348, -0.2454}, {0.0489, -0.1826, 0.1311}},
             {{-0.0702, 0.2878, -0.2054}, {0.0388, -0.1433, 0.0934}},
             {{-0.0623, 0.2604, -0.1823}, {0.0372, -0.1330, 0.0790}},
             {{-0.0649, 0.2718, -0.1940}, {0.0332, -0.1145, 0.0584}},
             {{-0.0795, 0.3238, -0.2401}, {-0.0766, 0.2509, -0.2460}}}};

    /////////  loose 3.5 sigma cuts /////////////////////
    double mmsq_cuts_loose[2][9][2][3] = {
        {{{-0.0996, 0.3951, -0.2790}, {0.0308, -0.1576, 0.1132}},
         {{-0.1338, 0.5126, -0.3764}, {0.0884, -0.3586, 0.2827}},
         {{-0.1360, 0.5191, -0.3837}, {0.0901, -0.3619, 0.2828}},
         {{-0.0849, 0.3472, -0.2433}, {0.0805, -0.3311, 0.2557}},
         {{-0.1066, 0.4148, -0.2966}, {0.0785, -0.3178, 0.2369}},
         {{-0.1028, 0.4141, -0.3111}, {0.1264, -0.4919, 0.3960}},
         {{-0.0461, 0.1860, -0.0874}, {0.0659, -0.2507, 0.1627}},
         {{-0.0463, 0.2076, -0.1222}, {0.0123, -0.0601, -0.0015}},
         {{-0.0926, 0.3076, -0.1682}, {0.0691, -0.2084, 0.0979}}},

        {{{-0.1483, 0.5756, -0.4404}, {0.0732, -0.2905, 0.2237}},
         {{-0.1473, 0.5723, -0.4373}, {0.0882, -0.3388, 0.2588}},
         {{-0.1431, 0.5560, -0.4220}, {0.0888, -0.3385, 0.2558}},
         {{-0.1229, 0.4842, -0.3594}, {0.0792, -0.3028, 0.2231}},
         {{-0.1098, 0.4383, -0.3207}, {0.0753, -0.2861, 0.2064}},
         {{-0.0919, 0.3741, -0.2651}, {0.0606, -0.2295, 0.1532}},
         {{-0.0822, 0.3390, -0.2345}, {0.0571, -0.2116, 0.1313}},
         {{-0.0845, 0.3491, -0.2445}, {0.0528, -0.1918, 0.1088}},
         {{-0.0800, 0.3383, -0.2390}, {-0.0760, 0.2364, -0.2471}}}};
    //////////////////////////////////////////////////
    int inv_mass_binning(float inv_mass, float inv_pPip_llim, float bin_size_inv)
    {
        for (int i = 0; i < 11; ++i)
        {
            float lower_limit = inv_pPip_llim + i * bin_size_inv;
            float upper_limit = inv_pPip_llim + (i + 1) * bin_size_inv;

            if (lower_limit <= inv_mass && inv_mass < upper_limit)
            {
                return i;
            }
        }
        return -1; // Return -1 if the mass is outside of the bins
    }
    /////////////////////////////////////////////////

    int inv_binning(float w, float inv_mass, bool hp)
    {
        // Find the index of the bin that w falls into
        int w_bin_index = static_cast<int>((w - 1.4) / 0.05);

        // Calculate the midpoint of the bin
        float w_mid = 1.4 + (w_bin_index + 0.5) * 0.05;
        float inv_ulim;
        float inv_llim;
        // Double_t xmin_5D[ndims_5D] = {((0.938272 + 0.13957) - 2 * Bin_size_pPip), (0.13957 + 0.13957) - 2 * Bin_size_pipPim, 0., 0., 0.};
        // Double_t xmax_5D[ndims_5D] = {((1.0 + 0.05 * w + 0.025 - MASS_PIM) + 2 * Bin_size_pPip), ((1.0 + 0.05 * w + 0.025 - MASS_P) + 2 * Bin_size_pipPim), 180, 360, 360};

        if (hp)
        {
            inv_ulim = w_mid - MASS_PIM;
            inv_llim = MASS_P + MASS_PIP;
            // std::cout << " w  " << w << "  w_bin_index  " << w_bin_index << "   w mid  " << w_mid << " inv_llim  " << inv_llim << "   inv_ulim  " << inv_ulim << std::endl;
        }
        else
        {
            inv_ulim = w_mid - MASS_P;
            inv_llim = MASS_PIP + MASS_PIM;
        }
        int inv_bin_val = -1;

        for (int j = 0; j < 7; ++j)
        {
            float lim_value = inv_llim + (j + 1) * (inv_ulim - inv_llim) / 7.0;
            if (inv_mass <= lim_value)
            {
                inv_bin_val = j;
                break;
            }
        }
        return inv_bin_val;
    }

    ////////////////////////////////////////////
    int alpha_bin_val = -1;

    // Define q2 bin ranges
    std::vector<float> alpha_bins = {36.0, 72.0, 108.0, 144.0, 180.0, 216.0, 252.0, 288.0, 324.0, 360.0};

    // Find the bin based on q2

    int alpha_bining(float alpha)
    {

        for (int i = 0; i < alpha_bins.size(); ++i)
        {
            if (alpha < alpha_bins[i])
            {
                alpha_bin_val = i;
                break;
            }
        }
        return alpha_bin_val;
    }

    // double dCosTh(float theta)
    // {
    //     double theta_res = 18.0;
    //     int theta_bin_index = static_cast<int>((theta) / theta_res);
    //     double lower_edge = theta_res * theta_bin_index;
    //     double upper_edge = lower_edge + theta_res;
    //     // Calculate the cosine values of the upper and lower bin edges
    //     double cos_lower = TMath::Cos(lower_edge * TMath::Pi() / 180.0);
    //     double cos_upper = TMath::Cos(upper_edge * TMath::Pi() / 180.0);

    //     double difference = std::abs(cos_upper - cos_lower);

    //     return difference;
    // }
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    static const short w_bin = 16;
    THnSparse *threeDHist[q2_bin][w_bin];
    THnSparse *sevenDHist_pim[q2_bin][w_bin];
    THnSparse *sevenDHist_pip[q2_bin][w_bin];
    THnSparse *sevenDHist_prot[q2_bin][w_bin];

    THnSparse *h_5dim_prot_evt[q2_bin][w_bin];
    THnSparse *h_5dim_pip_evt[q2_bin][w_bin];
    THnSparse *h_5dim_pim_evt[q2_bin][w_bin];

    THnSparse *sevenDHist_pim_tight[q2_bin][w_bin];
    THnSparse *sevenDHist_pip_tight[q2_bin][w_bin];
    THnSparse *sevenDHist_prot_tight[q2_bin][w_bin];
    THnSparse *h_5dim_prot_evt_tight[q2_bin][w_bin];
    THnSparse *h_5dim_pip_evt_tight[q2_bin][w_bin];
    THnSparse *h_5dim_pim_evt_tight[q2_bin][w_bin];

    THnSparse *sevenDHist_pim_loose[q2_bin][w_bin];
    THnSparse *sevenDHist_pip_loose[q2_bin][w_bin];
    THnSparse *sevenDHist_prot_loose[q2_bin][w_bin];
    THnSparse *h_5dim_prot_evt_loose[q2_bin][w_bin];
    THnSparse *h_5dim_pip_evt_loose[q2_bin][w_bin];
    THnSparse *h_5dim_pim_evt_loose[q2_bin][w_bin];

    TH1D_ptr w_gen_hist[q2_bin][w_bin];
    TH1D_ptr q2_gen_hist[q2_bin][w_bin];

    TH1D_ptr w_gen_hist_inv_pPip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pPip[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_inv_pPim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pPim[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_inv_pipPim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pipPim[q2_bin][w_bin][11];

    TH1D_ptr inv_pPip_hist[q2_bin][w_bin][11];
    TH1D_ptr inv_pPim_hist[q2_bin][w_bin][11];
    TH1D_ptr inv_pipPim_hist[q2_bin][w_bin][11];
    // TH1D *histogram = new TH1D("histogram", "Title", 100, 0, 100);

    TH1D_ptr w_gen_hist_th_prot[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_prot[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_th_pip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_pip[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_th_pim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_pim[q2_bin][w_bin][11];

    TH1D_ptr prot_theta_hist[q2_bin][w_bin][11];
    TH1D_ptr pip_theta_hist[q2_bin][w_bin][11];
    TH1D_ptr pim_theta_hist[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_al_prot[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_prot[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_al_pip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_pip[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_al_pim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_pim[q2_bin][w_bin][11];

    TH1D_ptr prot_alpha_hist[q2_bin][w_bin][11];
    TH1D_ptr pip_alpha_hist[q2_bin][w_bin][11];
    TH1D_ptr pim_alpha_hist[q2_bin][w_bin][11];

    static const short W_BIN_CHECK_NUM = 11;

    // std::string W_BIN_CHECK_NAME[W_BIN_CHECK_NUM] = {" All_W_range "," <1.30W<1.35 ",     " 1.35<W<1.40 ",
    //                                            " 1.40<W<1.45 ", " 1.45<W<1.50 ",
    //                                            " 1.50<W<1.55 ", " 1.55<W<1.60 ",
    //                                            " 1.60<W<1.65 ", " 1.65<W<1.70 ",
    //                                            " 1.70<W<1.75 ", " 1.75<W<1.80 "};

    // for inclusive
    std::string W_BIN_CHECK_NAME[W_BIN_CHECK_NUM] = {" ALL Q2 ", " Q2= 2.37 GeV2 ", " Q2= 2.77 GeV2 ", " Q2= 3.24 GeV2 ",
                                                     " Q2= 3.79 GeV2 ", " Q2= 4.44 GeV2 ",
                                                     " Q2= 5.19 GeV2 ", " Q2= 6.07 GeV2 ",
                                                     " Q2= 7.09 GeV2 ", " Q2= 8.29 GeV2 ",
                                                     " Q2= 9.70 GeV2 "};
    TH1D_ptr mm2_mPim_hist_check_th[11];
    TH1D_ptr W_hist_check_th[11];
    TH2D_ptr W_VS_Q2_hist_check_th[11];

    TH1D_ptr mm2_mPim_hist_check[11];
    TH1D_ptr W_hist_check[11];
    TH2D_ptr W_VS_Q2_hist_check[11];

    static const short NUM_CUT = 2;

    static const short theta_bin_NUM = 18;
    TH1D_ptr theta_pim_measured_3_sigma[theta_bin_NUM];
    std::string theta_bin_NAME[theta_bin_NUM] = {
        "0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180"};

    //  TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500,
    //  0, 0.5);
    //
    TH1D_ptr missing_Energy_hist = std::make_shared<TH1D>("e#pi^{+}#pi^{-}pX", "missing Energy", 500, -0.6, 0.6);

    TH1D_ptr diff_E_P_x_mu_hist_ = std::make_shared<TH1D>(
        "diff_E_P_x_mu ", "diff_E_P_x_mu ", 500, -1.0, 1.0);
    TH1D_ptr P_x_mu = std::make_shared<TH1D>("Mom P ", "Mom P ", 500, -2.0, 10.0);

    TH2D_ptr theta_vs_mom_elec[num_sectors];
    TH2D_ptr theta_vs_mom_prot[num_sectors];
    TH2D_ptr theta_vs_mom_pip[num_sectors];
    TH2D_ptr theta_vs_mom_pim[num_sectors];

    TH2D_ptr pim_phi_vs_theta_rec_FD_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_rec_FD_after_exclusive_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_measured_FD_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_measured_FD_after_exclusive_sec[num_sectors];

    // missingPiP
    static const short HADRON_NUM = 3;
    std::string HADRON_NAME[HADRON_NUM] = {
        "mProt", "mPip", "mPim"};
    static const short EFF_CONDITIONS_NUM_ALL = 3;
    std::string EFF_CONDITIONS_NAME_ALL[EFF_CONDITIONS_NUM_ALL] = {
        "missing", "exclusive", "after_MMSQ_exclusive_cuts"};
    TH2D_ptr pip_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr prot_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr pip_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr prot_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];

    TH1D_ptr weight_hist;
    TH1D_ptr pid_at_zero;
    TH1D_ptr mc_pid_at_zero;
    TH1D_ptr inv_mass_pPip;
    TH1D_ptr inv_mass_pPim;
    TH1D_ptr inv_mass_pipPim;
    TH1D_ptr theta_Prot_cm;
    TH1D_ptr theta_Pip_cm;
    TH1D_ptr theta_Pim_cm;
    TH1D_ptr phi_Prot_cm;
    TH1D_ptr phi_Pip_cm;
    TH1D_ptr phi_Pim_cm;
    TH1D_ptr alpha_Prot_cm;
    TH1D_ptr alpha_Pip_cm;
    TH1D_ptr alpha_Pim_cm;

    TH1D_ptr dp_prot_cdfd_hist;
    TH1D_ptr dp_pip_cdfd_hist;
    TH1D_ptr dth_prot_cdfd_hist;
    TH1D_ptr dth_pip_cdfd_hist;
    TH1D_ptr dphi_prot_cdfd_hist;
    TH1D_ptr dphi_pip_cdfd_hist;

    TH1D_ptr inv_mass_pPip_swapped;
    TH1D_ptr inv_mass_pPim_swapped;
    TH1D_ptr inv_mass_pipPim_swapped;
    TH1D_ptr theta_Prot_cm_swapped;
    TH1D_ptr theta_Pip_cm_swapped;
    TH1D_ptr theta_Pim_cm_swapped;
    TH1D_ptr phi_Prot_cm_swapped;
    TH1D_ptr phi_Pip_cm_swapped;
    TH1D_ptr phi_Pim_cm_swapped;
    TH1D_ptr alpha_Prot_cm_swapped;
    TH1D_ptr alpha_Pip_cm_swapped;
    TH1D_ptr alpha_Pim_cm_swapped;

    TH1D_ptr dp_prot_hist;
    TH1D_ptr dp_pip_hist;
    TH1D_ptr dp_ambi_prot_all_hist;
    TH1D_ptr dp_ambi_pip_all_hist;
    TH1D_ptr dp_sum_hist;
    TH1D_ptr dp_sum_hist_twoPi;
    TH1D_ptr dp_prot_for_pip_hist;
    TH1D_ptr dp_pip_for_prot_hist;

    TH1D_ptr entries_in_each_event;
    TH1D_ptr entries_prot;
    TH1D_ptr entries_pip;
    TH1D_ptr MM2_mPim_all_comb;
    TH1D_ptr MM2_mPim_1_comb;
    TH1D_ptr MM2_mPim_2_comb;
    TH1D_ptr MM2_mPim_3_comb;
    TH1D_ptr MM2_mPim_4_or_more_comb;

    TH1D_ptr p_gen_prot_hist;
    TH1D_ptr p_gen_pip_hist;
    TH1D_ptr p_gen_ambi_prot_hist;
    TH1D_ptr p_gen_ambi_pip_hist;
    TH1D_ptr p_gen_prot_for_pip_hist;
    TH1D_ptr p_gen_pip_for_prot_hist;

    TH1D_ptr p_rec_prot_hist;
    TH1D_ptr p_rec_pip_hist;
    TH1D_ptr p_rec_ambi_prot_hist;
    TH1D_ptr p_rec_ambi_pip_hist;
    TH1D_ptr p_rec_prot_for_pip_hist;
    TH1D_ptr p_rec_pip_for_prot_hist;

    TH1D_ptr scalar_triple_product_hist;
    TH1D_ptr scalar_triple_product_hist_with_mmsq_cuts;
    TH2D_ptr W_vs_sclar_product;
    TH2D_ptr W_vs_scalar_product_after_mmsq_cuts;

    TH1D_ptr W_hist;
    TH1D_ptr Q2_hist;
    TH2D_ptr W_vs_q2;
    TH1D_ptr W_P2pi_hist;

    TH1D_ptr W_thrown;
    TH2D_ptr W_vs_Q2_thrown;
    TH1D_ptr Q2_thrown;

    TH1D_ptr vz_position[CUTS];
    TH2D_ptr pcal_sec[CUTS];
    TH2D_ptr pcal_sec_ineff_cuts[CUTS];

    TH2D_ptr pcal_hx_hy_sec[CUTS];
    TH2D_ptr dcr1_sec[CUTS];
    TH2D_ptr dcr2_sec[CUTS];
    TH2D_ptr dcr3_sec[CUTS];

    // EC Sampling Fraction
    TH2D_ptr EC_sampling_fraction[CUTS];
    TH2D_ptr ECin_sf_vs_PCAL_sf[CUTS];
    TH1D_ptr momentum[CUTS];

    //// Hadron pid
    TH1D_ptr prot_Delta_vz_cut_fd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_fd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pim_Delta_vz_cut[CUTS];
    TH1D_ptr pim_Chi2pid_cut[CUTS];

    TH1D_ptr prot_Delta_vz_cut_cd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_cd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_cd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_cd[CUTS];

    TH2D_ptr dcr1_sec_prot[CUTS];
    TH2D_ptr dcr2_sec_prot[CUTS];
    TH2D_ptr dcr3_sec_prot[CUTS];
    TH2D_ptr dcr1_sec_pip[CUTS];
    TH2D_ptr dcr2_sec_pip[CUTS];
    TH2D_ptr dcr3_sec_pip[CUTS];
    TH2D_ptr dcr1_sec_pim[CUTS];
    TH2D_ptr dcr2_sec_pim[CUTS];
    TH2D_ptr dcr3_sec_pim[CUTS];
    TH2D_ptr phi_vs_mom_prot_fd[CUTS];
    TH2D_ptr phi_vs_mom_pip_fd[CUTS];
    TH2D_ptr phi_vs_mom_pim_fd[CUTS];
    TH1D_ptr theta_prot_fd[CUTS];
    TH1D_ptr theta_pip_fd[CUTS];
    TH1D_ptr theta_pim_fd[CUTS];
    TH2D_ptr Theta_prot_lab_vs_mom_prot_fd[CUTS];
    TH2D_ptr Theta_pip_lab_vs_mom_pip_fd[CUTS];
    TH2D_ptr Theta_pim_lab_vs_mom_pim_fd[CUTS];

    TH2D_ptr phi_vs_momT_prot_cd[CUTS];
    TH2D_ptr phi_vs_momT_pip_cd[CUTS];
    TH2D_ptr phi_vs_momT_pim_cd[CUTS];

    TH1D_ptr theta_prot_cd[CUTS];
    TH1D_ptr theta_pip_cd[CUTS];
    TH1D_ptr theta_pim_cd[CUTS];
    TH2D_ptr Theta_prot_lab_vs_mom_prot_cd[CUTS];
    TH2D_ptr Theta_pip_lab_vs_mom_pip_cd[CUTS];
    TH2D_ptr Theta_pim_lab_vs_mom_pim_cd[CUTS];

    TH2D_ptr Theta_fd_prot_lab_vs_mom_prot[num_sectors];
    TH2D_ptr Theta_fd_pip_lab_vs_mom_pip[num_sectors];
    TH2D_ptr Theta_fd_elec_lab_vs_mom_elec[num_sectors];

    TH2D_ptr Theta_prot_cm_vs_mom_prot;
    TH2D_ptr Theta_pip_cm_vs_mom_pip;
    TH2D_ptr Theta_pim_cm_vs_mom_pim;

    TH2D_ptr Theta_prot_thrown_cm_vs_mom_prot;
    TH2D_ptr Theta_pip_thrown_cm_vs_mom_pip;
    TH2D_ptr Theta_pim_thrown_cm_vs_mom_pim;

    TH2D_ptr Theta_prot_thrown_lab_vs_mom_prot;
    TH2D_ptr Theta_pip_thrown_lab_vs_mom_pip;
    TH2D_ptr Theta_pim_thrown_lab_vs_mom_pim;

    TH1D_ptr Phi_gamma;
    TH1D_ptr Phi_prot;
    TH1D_ptr Phi_pip;
    TH1D_ptr Phi_pim;

    TH1D_ptr alpha_pim;
    TH1D_ptr alpha_pip;
    TH1D_ptr alpha_prot;

    TH1D_ptr theta_prot_mc;
    TH1D_ptr theta_pip_mc;
    TH1D_ptr theta_pim_mc;

    TH1D_ptr Phi_gamma_mc;
    TH1D_ptr Phi_prot_mc;
    TH1D_ptr Phi_pip_mc;
    TH1D_ptr Phi_pim_mc;

    TH1D_ptr alpha_pim_mc;
    TH1D_ptr alpha_pip_mc;
    TH1D_ptr alpha_prot_mc;

    TH1D_ptr theta_prot_thrown;
    TH1D_ptr theta_pip_thrown;
    TH1D_ptr theta_pim_thrown;

    TH1D_ptr Phi_gamma_thrown;
    TH1D_ptr Phi_prot_thrown;
    TH1D_ptr Phi_pip_thrown;
    TH1D_ptr Phi_pim_thrown;

    TH1D_ptr alpha_pim_thrown;
    TH1D_ptr alpha_pip_thrown;
    TH1D_ptr alpha_prot_thrown;

    TH2D_ptr W_vs_q2_sec[num_sectors];
    TH1D_ptr W_sec[num_sectors];

    TH1D_ptr W_det[3];
    TH2D_ptr WQ2_det[3];

    TH1D_ptr W_hist_singlePi;
    TH1D_ptr Q2_hist_singlePi;
    TH2D_ptr W_vs_q2_singlePi;
    TH2D_ptr W_vs_q2_singlePi_sec[num_sectors];
    TH1D_ptr W_singlePi_sec[num_sectors];
    TH2D_ptr W_vs_MM_singlePi[num_sectors];

    TH2D_ptr W_vs_q2_Npip_sec[num_sectors];
    TH1D_ptr W_Npip_sec[num_sectors];
    TH1D_ptr MM_Npip_sec[num_sectors];

    TH1D_ptr MM_neutron;
    TH1D_ptr MM_neutron_sec[num_sectors];

    TH2D_ptr W_vs_MM;
    TH2D_ptr W_vs_MM2;

    //////////////////////////   electron pid cuts ///////////////
    TH1D_ptr htcc_nphe_sec[CUTS][num_sectors];
    TH1D_ptr elec_Chi2pid_sec[CUTS][num_sectors];
    TH1D_ptr vz_sec[CUTS][num_sectors];
    TH2D_ptr ECAL_VS_PCAL[CUTS][num_sectors][9];
    TH2D_ptr SF_VS_MOM[CUTS][num_sectors];
    std::string ECIN_ECOUT_MOM_NAME[9] = {" p<2 GeV ", " 2<p<=3 GeV ", " 3<p<=4 GeV ",
                                          " 4<p<=5 GeV ", " 5<=p<6 GeV ",
                                          " 6<p<=7 GeV ", " 7<p<=8 GeV ",
                                          " 8<p<=9 GeV ", " p>9 GeV "};

    TH2D_ptr W_vs_q2_twoPi_sec[num_sectors];
    TH1D_ptr MM_twoPi;
    TH1D_ptr MM_mPim_twoPi_sec[num_sectors];
    TH1D_ptr MM2_twoPi_excl;
    TH1D_ptr MM_twoPi_excl;

    TH1D_ptr MM2_twoPi_mPim;
    TH1D_ptr MM_twoPi_mPim;
    TH1D_ptr MM2_mPim_twoPi_sec[num_sectors];
    TH1D_ptr MM2_twoPi_missingPip;
    TH1D_ptr MM2_twoPi_missingPip_sec[num_sectors];
    TH1D_ptr MM2_twoPi_missingProt;
    TH1D_ptr MM2_twoPi_missingProt_sec[num_sectors];
    TH1D_ptr W_hist_twoPi;
    TH1D_ptr Q2_hist_twoPi;
    TH2D_ptr W_vs_q2_twoPi;
    TH2D_ptr W_vs_q2_twoPi_thrown;

    TH1D_ptr MMSQ_mPim_hist[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_with_cut[q2_bin][15];

    TH1D_ptr MMSQ_mPim_hist_3D[3][9][3];
    TH1D_ptr MMSQ_mPim_hist_with_cut_3D[3][9][3];

    TH1D_ptr MMSQ_mPim_hist_1_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_2_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_3_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_4_or_more_comb[q2_bin][15];

    TH1D_ptr Inv_mass_pPip[q2_bin][15];
    TH1D_ptr Inv_mass_pPim[q2_bin][15];
    TH1D_ptr Inv_mass_pipPim[q2_bin][15];
    TH1D_ptr Alpha_Prot_cm[q2_bin][15];
    TH1D_ptr Alpha_Pip_cm[q2_bin][15];
    TH1D_ptr Alpha_Pim_cm[q2_bin][15];

    static const short FDmomArraySize = 24;
    static const short CDmomArraySize = 15;

    TH1D_ptr dt_prot_fd_hist[FDmomArraySize];
    TH1D_ptr dt_prot_cd_hist[CDmomArraySize];
    TH1D_ptr dt_pip_fd_hist[FDmomArraySize];
    TH1D_ptr dt_pip_cd_hist[CDmomArraySize];
    TH1D_ptr dt_pim_fd_hist[FDmomArraySize];
    TH1D_ptr dt_pim_cd_hist[CDmomArraySize];

    TH2D_ptr W_vs_q2_twoPi_sec_thrown[num_sectors];

    // Mom vs Beta
    TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
    TH2D_ptr momvsbeta_hist_prot[2][3];
    TH2D_ptr momvsbeta_hist_pip[2][3];
    TH2D_ptr momvsbeta_hist_pim[2][3];

    // Mom vs Beta

    // Delta T
    TH2D_ptr delta_t_hist[3][2][3];

    // Delta T

public:
    Histogram(const std::string &output_file);
    ~Histogram();

    bool MM_cut(float w, float q2, float mm2)
    {
        int is_mc = 0;
        if (_mc)
        {
            is_mc = 1;
        }
        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                q2_bin_val = i;
                break;
            }
        }

        if ((mm2 < (mmsq_cuts[is_mc][q2_bin_val - 1][0][0] * pow(w, 2) + mmsq_cuts[is_mc][q2_bin_val - 1][0][1] * pow(w, 1) + mmsq_cuts[is_mc][q2_bin_val - 1][0][2])) &&
            (mm2 > (mmsq_cuts[is_mc][q2_bin_val - 1][1][0] * pow(w, 2) + mmsq_cuts[is_mc][q2_bin_val - 1][1][1] * pow(w, 1) + mmsq_cuts[is_mc][q2_bin_val - 1][1][2])))

        {
            return true;
        }
        else
            return false;
    }

    bool MM_cut_tight(float w, float q2, float mm2)
    {
        int is_mc = 0;
        if (_mc)
        {
            is_mc = 1;
        }
        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                q2_bin_val = i;
                break;
            }
        }

        if ((mm2 < (mmsq_cuts_tight[is_mc][q2_bin_val - 1][0][0] * pow(w, 2) + mmsq_cuts_tight[is_mc][q2_bin_val - 1][0][1] * pow(w, 1) + mmsq_cuts_tight[is_mc][q2_bin_val - 1][0][2])) &&
            (mm2 > (mmsq_cuts_tight[is_mc][q2_bin_val - 1][1][0] * pow(w, 2) + mmsq_cuts_tight[is_mc][q2_bin_val - 1][1][1] * pow(w, 1) + mmsq_cuts_tight[is_mc][q2_bin_val - 1][1][2])))

        {
            return true;
        }
        else
            return false;
    }

    bool MM_cut_loose(float w, float q2, float mm2)
    {
        int is_mc = 0;
        if (_mc)
        {
            is_mc = 1;
        }
        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                q2_bin_val = i;
                break;
            }
        }

        if ((mm2 < (mmsq_cuts_loose[is_mc][q2_bin_val - 1][0][0] * pow(w, 2) + mmsq_cuts_loose[is_mc][q2_bin_val - 1][0][1] * pow(w, 1) + mmsq_cuts_loose[is_mc][q2_bin_val - 1][0][2])) &&
            (mm2 > (mmsq_cuts_loose[is_mc][q2_bin_val - 1][1][0] * pow(w, 2) + mmsq_cuts_loose[is_mc][q2_bin_val - 1][1][1] * pow(w, 1) + mmsq_cuts_loose[is_mc][q2_bin_val - 1][1][2])))

        {
            return true;
        }
        else
            return false;
    }

    void populate_eff_check_mPim(const std::shared_ptr<Reaction> &_e, double min_w, double max_w, double min_theta, double max_theta, double min_phi, double max_phi, short index_w, short index_theta, short index_phi);
    // void populate_eff_check_exclusive(const std::shared_ptr<Reaction> &_e, double min, double max, short index_w);
    void Fill_eff_ckeck_mPim(const std::shared_ptr<Reaction> &_e);
    // void Fill_eff_ckeck_exclusive(const std::shared_ptr<Reaction> &_e);
    void writeHist_eff_check_mPim();
    // void writeHist_eff_check_mPim_after_MMSQ_cuts();
    // void writeHist_eff_check_exclusive();
    // void writeHist_eff_check_exclusive_MMSQ_cuts();

    void Fill_histthreeD(const std::shared_ptr<Reaction> &_e);
    void writeHists3D();
    void Fill_hists4D_background(const std::shared_ptr<Reaction> &_e);
    void writehists4D_background();

    void Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e);

    void Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e);

    void Fill_histSevenD_prot_evt(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pim_evt(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip_evt(const std::shared_ptr<Reaction> &_e);

    void Fill_histSevenD_prot_tight(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pim_tight(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip_tight(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_prot_evt_tight(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pim_evt_tight(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip_evt_tight(const std::shared_ptr<Reaction> &_e);

    void Fill_histSevenD_prot_loose(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pim_loose(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip_loose(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_prot_evt_loose(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pim_evt_loose(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_pip_evt_loose(const std::shared_ptr<Reaction> &_e);

    void writeHists7D_prot();
    void writeHists7D_thrown_prot();
    void writeHists7D_pim();
    void writeHists7D_thrown_pim();
    void writeHists7D_pip();
    void writeHists7D_thrown_pip();
    void writeHists7D_prot_evt();
    void writeHists7D_thrown_prot_evt();
    void writeHists7D_pip_evt();
    void writeHists7D_thrown_pip_evt();
    void writeHists7D_pim_evt();
    void writeHists7D_thrown_pim_evt();

    void writeHists7D_prot_tight();
    void writeHists7D_pim_tight();
    void writeHists7D_pip_tight();
    void writeHists7D_prot_evt_tight();
    void writeHists7D_pip_evt_tight();
    void writeHists7D_pim_evt_tight();

    void writeHists7D_prot_loose();
    void writeHists7D_pim_loose();
    void writeHists7D_pip_loose();
    void writeHists7D_prot_evt_loose();
    void writeHists7D_pip_evt_loose();
    void writeHists7D_pim_evt_loose();

    //// bin centering corrections
    void Fill_hist1D_thrown_w_q2(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_inv_mass(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_theta(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_alpha(const std::shared_ptr<MCReaction> &_e);

    void writeHists1D_thrown_w_gen();
    void writeHists1D_thrown_q2_gen();

    void writeHists1D_thrown_w_gen_inv_pPip();
    void writeHists1D_thrown_q2_gen_inv_pPip();

    void writeHists1D_thrown_w_gen_inv_pPim();
    void writeHists1D_thrown_q2_gen_inv_pPim();

    void writeHists1D_thrown_w_gen_inv_pipPim();
    void writeHists1D_thrown_q2_gen_inv_pipPim();

    void writeHists1D_thrown_w_gen_th_prot();
    void writeHists1D_thrown_q2_gen_th_prot();

    void writeHists1D_thrown_w_gen_th_pip();
    void writeHists1D_thrown_q2_gen_th_pip();

    void writeHists1D_thrown_w_gen_th_pim();
    void writeHists1D_thrown_q2_gen_th_pim();

    void writeHists1D_thrown_w_gen_al_prot();
    void writeHists1D_thrown_q2_gen_al_prot();

    void writeHists1D_thrown_w_gen_al_pip();
    void writeHists1D_thrown_q2_gen_al_pip();

    void writeHists1D_thrown_w_gen_al_pim();
    void writeHists1D_thrown_q2_gen_al_pim();

    void writeHists1D_thrown_protPip();
    void writeHists1D_thrown_protPim();
    void writeHists1D_thrown_pipPim();
    void writeHists1D_thrown_th_prot();
    void writeHists1D_thrown_th_pip();
    void writeHists1D_thrown_th_pim();
    void writeHists1D_thrown_alpha_prot();
    void writeHists1D_thrown_alpha_pip();
    void writeHists1D_thrown_alpha_pim();
    double CosTheta(int theta_bin_);

    // W and Q^2
    void makeHists_sector();
    void Fill_WvsQ2(const std::shared_ptr<Reaction> &_e);
    // void Fill_WvsQ2(const std::shared_ptr<MCReaction> &_e);
    void Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction> &_e);
    void Fill_WvsQ2_Npip(const std::shared_ptr<Reaction> &_e);
    void Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction> &_e);
    // void Fill_WvsQ2_twoPi(const std::shared_ptr<MCReaction>& _e);
    void Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<MCReaction> &_e);
    void Write_WvsQ2();

    void Fill_W_vs_Q2_all_sec();
    void Fill_W_vs_Q2_thrown();
    void Fill_inv_mass_hist();

    // void makeHistMMSQ_mPim_3D();
    // void Fill_MMSQ_mPim_3D(const std::shared_ptr<Reaction> &_e);
    // void writeMMSQ_mPim_3D();

    void makeHistMMSQ_mPim();
    void Fill_MMSQ_mPim(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim();

    void Fill_MMSQ_mPim_1_comb(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_1_comb();

    void Fill_MMSQ_mPim_2_comb(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_2_comb();

    void Fill_MMSQ_mPim_3_comb(float dv2, const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_3_comb();

    void Fill_MMSQ_mPim_4_or_more_comb(float dv2, const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_4_or_more_comb();

    void write_Inv_Mass_hist();
    void Fill_cdfd_prot(float dp, float dth, float dphi, const std::shared_ptr<Reaction> &_e);
    void Fill_cdfd_pip(float dp, float dth, float dphi, const std::shared_ptr<Reaction> &_e);

    // P and E
    // ecectron cuts
    void makeHists_electron_cuts();
    void FillHists_electron_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);
    void FillHists_electron_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);

    void FillHists_prot_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_prot_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &prot);

    void FillHists_pip_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pip_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &pip);
    void FillHists_pim_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pim_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);

    void Write_Electron_cuts();
    void Write_Hadrons_cuts();
    void write_hist_cd_fid();
    void Fill_hist_cd_fid(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Reaction> &_e, int i);

    void Fill_pi0(const std::shared_ptr<Reaction> &_e);

    void makeHists_MomVsBeta();
    void Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part, const std::shared_ptr<Reaction> &_e);
    void Write_MomVsBeta();

    // Delta T
    void makeHists_deltat();
    void Fill_deltat_pip(const std::shared_ptr<Branches12> &data,
                         const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pim(const std::shared_ptr<Branches12> &data,
                         const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_prot(const std::shared_ptr<Branches12> &data,
                          const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_before_cut(const std::shared_ptr<Branches12> &data,
                                const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_prot_after_cut(const std::shared_ptr<Branches12> &data,
                                    const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pip_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pim_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Write_deltat();

    //////////////
    void Fill_deltaP_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_pip(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_prot_for_pip(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_pip_for_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_ambi_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_ambi_pip(const std::shared_ptr<Reaction> &_e, double dp);

    void Fill_deltaP_sum_twoPi(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_sum(const std::shared_ptr<Reaction> &_e, double dp);
    void Write_deltaP();

    void Fill_Entries(int num_entries);
    void Fill_Entries_prot(int num_entries);
    void Fill_Entries_pip(int num_entries);
    void Fill_all_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_1_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_2_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_3_Combi(float dv2, const std::shared_ptr<Reaction> &_e);
    void Fill_4_or_more_Combi(float dv2, const std::shared_ptr<Reaction> &_e);

    ///////////////

    void makeHistTheta_pim_measured();
    void populate_theta_pim_measured(const std::shared_ptr<Reaction> &_e, double min, double max, short index_theta_pim);
    void Fill_theta_pim_measured(const std::shared_ptr<Reaction> &_e);
    void write_hist_theta_pim_measured();

    void Write();

    // Function to get the momentum range index based on the value of p

    int getMomRange(double p)
    {
        const double boundaries[] = {2, 3, 4, 5, 6, 7, 8, 9};
        const int numBoundaries = sizeof(boundaries) / sizeof(boundaries[0]);

        for (int i = 0; i < numBoundaries; ++i)
        {
            if (p < boundaries[i])
            {
                return i;
            }
        }
        return numBoundaries; // For p >= 9
    }
};

#endif
