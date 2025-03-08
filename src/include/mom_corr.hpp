#ifndef MOM_CORR_H_GUARD
#define MOM_CORR_H_GUARD
#include <cstdlib>
#include <ctime> // add this line to include <ctime> header
#include <iostream>
#include "TROOT.h"
#include "constants.hpp"
// using namespace std;
class mom_corr
{
private:
  // float alpha_CD[3][3] = {{1, 1, 1}, {0.9, 0.9, 0.9}, {0.9, 0.9, 0.9}};
  // float alpha_FD[3] = {1.0, 0.9, 0.9};

  float alpha_CD[3][3] = {{0.9, 0.9, 0.9}, {0.8, 0.8, 0.8}, {0.8, 0.8, 0.8}};
  float alpha_FD[3] = {0.9, 0.8, 0.8};

  double sys_fac = 0.9;

public:
  mom_corr() {};
  ~mom_corr();

  bool is_FD(int prot_status);
  bool is_CD(int prot_status);
  bool is_lower_band(float mom_, float theta_DCr1_, int status_);

  float CD_prot_Emom_corr(float mom_, float theta_);
  float FD_prot_Emom_corr_lower(float mom_, float theta_);
  float FD_prot_Emom_corr_upper(float mom_, float theta_);

  float CD_prot_Eth_corr(float mom_, float theta_);
  float FD_prot_Eth_corr_lower(float mom_, float theta_);
  float FD_prot_Eth_corr_upper(float mom_, float theta_);

  float CD_prot_Eph_corr(float mom_, float theta_, float phi_);
  float FD_prot_Eph_corr_lower(float mom_, float theta_, float phi_);
  float FD_prot_Eph_corr_upper(float mom_, float theta_, float phi_);

  // // Calcuating Energy loss corr parameters
  float A_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
  float B_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec);

  // float A_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
  // float B_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
  // float C_th(float mom_, float theta_, int dc_sec);

  // float A_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
  // float B_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
  // float C_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);

  // pip energy loss correction functions
  float CD_pip_Emom_corr(float mom_, float theta_);
  float FD_pip_Emom_corr_lower(float mom_, float theta_);
  float FD_pip_Emom_corr_upper(float mom_, float theta_);

  float CD_pip_Eth_corr(float mom_, float theta_);
  float FD_pip_Eth_corr_lower(float mom_, float theta_);
  float FD_pip_Eth_corr_upper(float mom_, float theta_);

  float CD_pip_Eph_corr(float mom_, float theta_, float phi_);
  float FD_pip_Eph_corr_lower(float mom_, float theta_, float phi_);
  float FD_pip_Eph_corr_upper(float mom_, float theta_, float phi_);

  // pim energy loss correction functions
  float CD_pim_Emom_corr(float mom_, float theta_);
  float FD_pim_Emom_corr_lower(float mom_, float theta_);
  float FD_pim_Emom_corr_upper(float mom_, float theta_);

  float CD_pim_Eth_corr(float mom_, float theta_);
  float FD_pim_Eth_corr_lower(float mom_, float theta_);
  float FD_pim_Eth_corr_upper(float mom_, float theta_);

  float CD_pim_Eph_corr(float mom_, float theta_, float phi_);
  float FD_pim_Eph_corr_lower(float mom_, float theta_, float phi_);
  float FD_pim_Eph_corr_upper(float mom_, float theta_, float phi_);

  // hadron mom corrections
  double dppC(float Px, float Py, float Pz, int sec, int ivec);
  double elossPipFD(double pion_p, double pip_theta);
  double elossPipCD(double pion_p, double pip_theta);

  // 4-vector method dp corrections:
  float CD_prot_Hmom_corr(float mom_, float phi_);
  float FD_prot_Hmom_corr(float mom_, float dc_sec);

  float CD_pip_Hmom_corr(float mom_, float phi_);
  float FD_pip_Hmom_corr(float mom_, float dc_sec);

  float CD_pim_Hmom_corr(float mom_, float phi_);
  float FD_pim_Hmom_corr(float mom_, float dc_sec);

  //// mes - missing dp corrections

  // float CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot);
  // float FD_prot_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_prot);
  // float FD_prot_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_prot);
  // float FD_prot_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_prot);
  // float FD_prot_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_prot);

  // float CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip);
  // float FD_pip_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pip);
  // float FD_pip_Hmom_corr_lower_Except_All_FD(float mom_, float , float alpha_pipdc_sec);
  // float FD_pip_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pip);
  // float FD_pip_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pip);

  // float CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim);
  // float FD_pim_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pim);
  // float FD_pim_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pim);
  // float FD_pim_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pim);
  // float FD_pim_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pim);

  float CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot_CD[3]);
  float FD_prot_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_prot_FD);
  float FD_prot_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_prot_FD);
  float FD_prot_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_prot_FD);
  float FD_prot_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_prot_FD);

  float CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip_CD[3]);
  float FD_pip_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pip_FD);
  float FD_pip_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pip_FD);
  float FD_pip_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pip_FD);
  float FD_pip_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pip_FD);

  float CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim_CD[3]);
  float FD_pim_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pim_FD);
  float FD_pim_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pim_FD);
  float FD_pim_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pim_FD);
  float FD_pim_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pim_FD);

  // float alpha_prot_mom_corr_FD[4] = {0.5, 0.6, 0.5, 0.5};

  // void random_no_gen() {
  //   std::srand(std::time(nullptr));  // seed the random number generator
  //                                    // rest of your code here
  //   for (int i = 0; i < 4; i++) {
  //     float variation =
  //         -0.1 + (static_cast<float>(std::rand()) / RAND_MAX) * 0.2;  // generate random percentage between -10%
  //         and +10%
  //     alpha_prot_mom_corr_FD[i] *= (1.0 + variation);  // apply the random percentage to the current element

  //     // std::cout << "New value in element " << i << " of array: " << alpha_prot_mom_corr_FD[i] << std::endl;
  //   }
  // }
};

#endif
