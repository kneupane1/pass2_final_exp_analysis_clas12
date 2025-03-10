#include "cuts.hpp"
#include <iostream>
#include "TFile.h"
#include "histogram.hpp"
#include "reaction.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12> &data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : _data(data), _dt(dt) {}
// size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists,
//            int thread_id);

Cuts::~Cuts() {}

// auto RxnClas = std::make_shared<Reaction>(_data, beam_energy);

// // // // /////////////////////////// exp data dt cuts ////////////// prot, pip, pim /////////////////
// double dt_cut_fd_up[3][6] = {{-0.0055, 0.10236, -0.7266, 2.447, -3.926, 2.893},
//                              {-0.003983, 0.0719, -0.4958, 1.634, -2.598, 2.059},
//                              {-0.002695, 0.0391, -0.2146, 0.56, -0.7153, 0.743}};
// double dt_cut_fd_down[3][6] = {{0.005783, -0.1055, 0.7295, -2.377, 3.672, -2.697},
//                                {0.005386, -0.0987, 0.6816, -2.203, 3.326, -2.277},
//                                {0.00419, -0.06027, 0.3281, -0.8403, 1.025, -0.8545}};
// double dt_cut_cd_up[3][3] = {{0.06, -0.3303, 0.7656},
//                              {0.01909, 0.00434, 0.428},
//                              {-0.005173, 0.09735, 0.3801}};
// double dt_cut_cd_down[3][3] = {{-0.04974, 0.286, -0.73},
//                                {-0.03662, 0.1555, -0.4675},
//                                {-0.02998, 0.1083, -0.4429}};

// // /////////////////////////// both data dt cuts ////////////// {exp,sim}-> prot, pip, pim /////////////////

double dt_cut_fd_up[2][3][6] = {{{-0.0055, 0.10236, -0.7266, 2.447, -3.926, 2.893},
                                 {-0.003983, 0.0719, -0.4958, 1.634, -2.598, 2.059},
                                 {-0.002695, 0.0391, -0.2146, 0.56, -0.7153, 0.743}},
                                {{-0.003056, 0.05838, -0.429, 1.513, -2.592, 2.297},
                                 {-0.0013075, 0.02531, -0.1893, 0.6865, -1.227, 1.377},
                                 {-0.0009727, 0.0131, -0.0642, 0.1498, -0.2065, 0.6343}}};

double dt_cut_fd_down[2][3][6] = {{{0.005783, -0.1055, 0.7295, -2.377, 3.672, -2.697},
                                   {0.005386, -0.0987, 0.6816, -2.203, 3.326, -2.277},
                                   {0.00419, -0.06027, 0.3281, -0.8403, 1.025, -0.8545}},
                                  {{0.004017, -0.07574, 0.5454, -1.866, 3.047, -2.426},
                                   {0.002022, -0.0358, 0.2452, -0.823, 1.369, -1.379},
                                   {0.002928, -0.04208, 0.2229, -0.5396, 0.6313, -0.7734}}};

double dt_cut_cd_up[2][3][3] = {{{0.06, -0.3303, 0.7656},
                                 {0.01909, 0.00434, 0.428},
                                 {-0.005173, 0.09735, 0.3801}},
                                {{0.05585, -0.2876, 0.714},
                                 {0.0458, -0.1715, 0.5796},
                                 {0.014305, -0.04828, 0.5063}}};

double dt_cut_cd_down[2][3][3] = {{{-0.04974, 0.286, -0.73},
                                   {-0.03662, 0.1555, -0.4675},
                                   {-0.02998, 0.1083, -0.4429}},
                                  {{-0.02893, 0.1836, -0.649},
                                   {-0.0731, 0.303, -0.58745},
                                   {-0.03882, 0.1859, -0.519}}};

//////////////////////////////////////////////////////////////////////////////////////////////

// double dt_cut_fd_up[3][6] = {{-0.003056, 0.05838, -0.429, 1.513, -2.592, 2.297},
//                              {-0.0013075, 0.02531, -0.1893, 0.6865, -1.227, 1.377},
//                              {-0.0009727, 0.0131, -0.0642, 0.1498, -0.2065, 0.6343}};
// double dt_cut_fd_down[3][6] = {{0.004017, -0.07574, 0.5454, -1.866, 3.047, -2.426},
//                                {0.002022, -0.0358, 0.2452, -0.823, 1.369, -1.379},
//                                {0.002928, -0.04208, 0.2229, -0.5396, 0.6313, -0.7734}};
// double dt_cut_cd_up[3][3] = {{0.05585, -0.2876, 0.714},
//                              {0.0458, -0.1715, 0.5796},
//                              {0.014305, -0.04828, 0.5063}};
// double dt_cut_cd_down[3][3] = {{-0.02893, 0.1836, -0.649},
//                                {-0.0731, 0.303, -0.58745},
//                                {-0.03882, 0.1859, -0.519}};

bool Pass2_Cuts::IsPip(int i, std::string condition)
{
        float dt_shift = -0.06;
        int is_mc = 0;
        if (_mc)
        {
                is_mc = 1;
        }
        if (_data->gpart() <= i)
                return false;
        bool _pip = true;
        _pip &= (_data->charge(i) == POSITIVE);
        _pip &= (_data->pid(i) == PIP);
        // _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
        _pip &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

        // // // // // min/max mom cuts
        if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
        // {
        // if (!(_dt->isCtof()))
        {
                // if (!std::isnan(_dt->dt_Pi(i)))
                // std::cout << "  particle status at cut level :  " << _data->status(i) << "   dt pi at ftof  " << _dt->dt_Pi(i) << std::endl;

                _pip &= (_data->p(i) > 0.5);
                // _pip &= (_data->p(i) < 4.6);
                _pip &= (_dt->dt_Pi(i) < ((dt_cut_fd_up[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][1][1] * pow(_data->p(i), 4) +
                                           dt_cut_fd_up[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][1][3] * pow(_data->p(i), 2) +
                                           dt_cut_fd_up[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][1][5]) +
                                          dt_shift));
                _pip &= (_dt->dt_Pi(i) > ((dt_cut_fd_down[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][1][1] * pow(_data->p(i), 4) +
                                           dt_cut_fd_down[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][1][3] * pow(_data->p(i), 2) +
                                           dt_cut_fd_down[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][1][5]) -
                                          dt_shift));

                _pip &= DC_fiducial_cut_XY_PIP(i, 2, condition);
                _pip &= DC_Ineff_cut_X_Y(i, 2, condition);
        }
        // }
        else if (abs(_data->status(i)) >= 4000)
        // else if ((_dt->isCtof()))
        {
                // std::cout << "   dt pi at ctof  " << _dt->dt_Pi(i) << std::endl;

                _pip &= (_data->p(i) > 0.2);
                // _pip &= (_data->p(i) < 1.7);
                _pip &= (_dt->dt_Pi(i) < ((dt_cut_cd_up[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][1][1] * _data->p(i) + dt_cut_cd_up[is_mc][1][2]) + dt_shift));
                _pip &= (_dt->dt_Pi(i) > ((dt_cut_cd_down[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][1][1] * _data->p(i) + dt_cut_cd_down[is_mc][1][2]) - dt_shift));
                _pip &= CD_fiducial_had(i, condition);
        }
        // _pip &= (_data->p(i) > 0.2);
        _pip &= Hadron_Delta_vz_cut(i, condition);
        _pip &= Hadron_Chi2pid_cut(i, condition);
        return _pip;
}
bool Pass2_Cuts::IsProton(int i, std::string condition)
{
        float dt_shift = -0.06;
        int is_mc = 0;
        if (_mc)
        {
                is_mc = 1;
        }
        if (_data->gpart() <= i)
                return false;
        bool _proton = true;
        _proton &= (_data->charge(i) == POSITIVE);
        _proton &= (_data->pid(i) == PROTON);
        // _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
        // // // _proton &= !(abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
        _proton &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
        // // // // min/max mom cuts
        if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
        // {
        // if (!(_dt->isCtof()))
        {
                // std::cout << "  particle status at cut level :  " << _data->status(i) << "   dt prot at ftof  " << _dt->dt_P(i) << std::endl;

                _proton &= (_data->p(i) > 0.4);
                // _proton &= (_data->p(i) < 4.5);
                _proton &= (_dt->dt_P(i) < ((dt_cut_fd_up[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][0][1] * pow(_data->p(i), 4) +
                                             dt_cut_fd_up[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][0][3] * pow(_data->p(i), 2) +
                                             dt_cut_fd_up[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][0][5]) +
                                            dt_shift));

                _proton &= (_dt->dt_P(i) > ((dt_cut_fd_down[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][0][1] * pow(_data->p(i), 4) +
                                             dt_cut_fd_down[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][0][3] * pow(_data->p(i), 2) +
                                             dt_cut_fd_down[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][0][5]) -
                                            dt_shift));

                _proton &= DC_fiducial_cut_XY_PROT(i, 1, condition);
                _proton &= DC_Ineff_cut_X_Y(i, 1, condition);
        }
        // }
        else if (abs(_data->status(i)) >= 4000)
        // else if ((_dt->isCtof()))
        {
                // std::cout << "   dt prot at ctof  " << _dt->dt_P(i) << std::endl;

                _proton &= (_data->p(i) > 0.2); /// this 0.4 look harse when we do missing Pim channel
                // _proton &= (_data->p(i) < 2.0);
                _proton &= (_dt->dt_P(i) < ((dt_cut_cd_up[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][0][1] * _data->p(i) + dt_cut_cd_up[is_mc][0][2]) + dt_shift));
                _proton &= (_dt->dt_P(i) > ((dt_cut_cd_down[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][0][1] * _data->p(i) + dt_cut_cd_down[is_mc][0][2]) - dt_shift));
                _proton &= CD_fiducial_had(i, condition);
        }
        // _proton &= (_data->p(i) > 0.2);
        // _proton &= (abs(_data->chi2d(i)) < 0.5);
        _proton &= Hadron_Delta_vz_cut(i, condition);
        _proton &= Hadron_Chi2pid_cut(i, condition);
        return _proton;
}

bool Pass2_Cuts::IsPim(int i)
{
        int is_mc = 0;
        if (_mc)
        {
                is_mc = 1;
        }
        if (_data->gpart() <= i)
                return false;
        bool _pim = true;
        _pim &= (_data->charge(i) == NEGATIVE);
        _pim &= (_data->pid(i) == PIM);

        _pim &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
        // min / max mom cuts
        if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
        {
                _pim &= (_data->p(i) > 0.4);
                // _pim &= (_data->p(i) < 4.5);
                _pim &= (_dt->dt_Pi(i) < (dt_cut_fd_up[is_mc][2][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][2][1] * pow(_data->p(i), 4) +
                                          dt_cut_fd_up[is_mc][2][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][2][3] * pow(_data->p(i), 2) +
                                          dt_cut_fd_up[is_mc][2][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][2][5]));

                _pim &= (_dt->dt_Pi(i) > (dt_cut_fd_down[is_mc][2][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][2][1] * pow(_data->p(i), 4) +
                                          dt_cut_fd_down[is_mc][2][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][2][3] * pow(_data->p(i), 2) +
                                          dt_cut_fd_down[is_mc][2][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][2][5]));
        }
        else if (abs(_data->status(i)) >= 4000)
        {

                _pim &= (_data->p(i) > 0.2);
                // _pim &= (_data->p(i) < 1.9);
                _pim &= (_dt->dt_Pi(i) < (dt_cut_cd_up[is_mc][2][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][2][1] * _data->p(i) + dt_cut_cd_up[is_mc][2][2]));
                _pim &= (_dt->dt_Pi(i) > (dt_cut_cd_down[is_mc][2][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][2][1] * _data->p(i) + dt_cut_cd_down[is_mc][2][2]));
        }
        // _pim &= (_data->p(i) > 0.2);

        _pim &= Hadron_Delta_vz_cut(i, "mid"); /// this because of the fact that hadron cuts are removed for pim
        _pim &= Hadron_Chi2pid_cut(i, "mid");  /// this because of the fact that hadron cuts are removed for pim
        // _pim &= DC_fiducial_cut_XY(0, 0); /// this because of the fact that hadron cuts are removed for pim, same dc cuts as electrons

        return _pim;
}
// /////////////////////// Pass2_Cuts ///////////////////////
bool Pass2_Cuts::ElectronCuts(std::string condition)
{
        bool cut = true;
        if (!cut)
                return false;
        cut &= (_data->gpart() > 0);
        cut &= (_data->gpart() < 20);
        // // //
        cut &= (_data->pid(0) == ELECTRON);
        cut &= (_data->charge(0) == NEGATIVE);
        cut &= DC_z_vertex_cut(condition);
        cut &= (_data->p(0) > 1.50);
        // cut &= (_data->vz(0) > -(2.78 + 3 * 2.16) && _data->vz(0) < (-2.78 + 3 * 2.16)); // 3 sigma cut
        cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
        // // cut &= (abs(_data->chi2pid(0)) < 3); ////////////// check it....... along with simulations
        // // cut &= CC_nphe_cut();
        cut &= DC_fiducial_cut_XY_E(condition);
        cut &= EC_sampling_fraction_cut(condition);
        cut &= PCAL_minimum_energy();
        cut &= PCAL_fiducial_cut_X_Y(condition);
        // // //cut &= EC_inner_vs_EC_outer();
        cut &= EC_hit_position_fiducial_cut_homogeneous(condition);
        cut &= PCAL_Ineff_cut_X_Y();

        return cut;
}
// bool Pass2_Cuts::HadronsCuts(int i)
// {
//         bool cut = true;
//         // if (_data->pid(i) == PROTON || _data->pid(i) == PIP)
//         if (_data->charge(i) == POSITIVE)
//         {
//                 if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
//                         cut &= DC_fiducial_cut_XY(i);
//                 // cut &= DC_fiducial_cut_theta_phi(i);
//                 else if (4000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000)
//                         cut &= CD_fiducial_had(i);
//         }

//         cut &= Hadron_Delta_vz_cut(i);
//         cut &= Hadron_Chi2pid_cut(i);

//         return cut;
// }

bool Pass2_Cuts::CC_nphe_cut()
{
        float nphe_min = 2;
        return (_data->cc_nphe_tot(0) > nphe_min);
}
bool Pass2_Cuts::PCAL_minimum_energy()
{
        double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
        return (_data->ec_pcal_energy(0) > edep_medium);
}

bool Pass2_Cuts::EC_sampling_fraction_cut(std::string condition)
{
        int isec = (_data->ec_pcal_sec(0) - 1);
        double upper_lim_total = 0;
        double lower_lim_total = 0;

        double mean_minus_3_5_sigma[6][3] = {{-0.001124, 0.01212, 0.1616}, {-0.001773, 0.02165, 0.1334}, {-0.002983, 0.03586, 0.0957}, {-0.001408, 0.01909, 0.1344}, {-0.0007095, 0.00882, 0.1687}, {-0.002068, 0.02394, 0.1298}};
        double mean_plus_3_5_sigma[6][3] = {{-0.001306, 0.012856, 0.2598}, {-0.00082, 0.00801, 0.271}, {-0.00118, 0.01078, 0.269}, {-0.0009036, 0.01098, 0.2607}, {-0.0001322, 0.002972, 0.2744}, {-0.001183, 0.01233, 0.257}};
        // double mean_minus_3_sigma[6][3] = {{-0.001138, 0.01218, 0.1686}, {-0.001705, 0.02068, 0.1433}, {-0.002853, 0.03406, 0.1081}, {-0.001371, 0.01851, 0.1434}, {-0.0006685, 0.00841, 0.1763}, {-0.002005, 0.02312, 0.1389}};
        // double mean_plus_3_sigma[6][3] = {{-0.001293, 0.0128, 0.2527}, {-0.0008883, 0.00899, 0.261}, {-0.001308, 0.01257, 0.2566}, {-0.0009394, 0.01156, 0.2517}, {-0.0001734, 0.00339, 0.2668}, {-0.0012455, 0.01316, 0.2479}};
        // double mean_minus_4_sigma[6][3] = {{-0.001112, 0.01207, 0.1545}, {-0.001841, 0.02263, 0.1236}, {-0.003113, 0.03763, 0.0833}, {-0.001444, 0.01967, 0.1255}, {-0.000751, 0.00924, 0.1613}, {-0.00213, 0.02478, 0.12067}};
        // double mean_plus_4_sigma[6][3] = {{-0.001319, 0.01291, 0.2668}, {-0.0007524, 0.00704, 0.2808}, {-0.001051, 0.008995, 0.2815}, {-0.0008674, 0.0104, 0.2698}, {-9.09e-05, 0.002552, 0.282}, {-0.00112, 0.0115, 0.266}};

        for (Int_t k = 0; k < 6; k++)
        {
                if (isec == k)
                {
                        // if (condition == "mid")
                        // {
                        upper_lim_total = mean_plus_3_5_sigma[k][0] * pow(_data->p(0), 2) + (mean_plus_3_5_sigma[k][1]) * _data->p(0) +
                                          mean_plus_3_5_sigma[k][2];

                        lower_lim_total = mean_minus_3_5_sigma[k][0] * pow(_data->p(0), 2) + (mean_minus_3_5_sigma[k][1]) * _data->p(0) +
                                          mean_minus_3_5_sigma[k][2];
                        // }

                        // if (condition == "tight")
                        // {
                        //         upper_lim_total = mean_plus_3_sigma[k][0] * pow(_data->p(0), 2) + (mean_plus_3_sigma[k][1]) * _data->p(0) +
                        //                           mean_plus_3_sigma[k][2];

                        //         lower_lim_total = mean_minus_3_sigma[k][0] * pow(_data->p(0), 2) + (mean_minus_3_sigma[k][1]) * _data->p(0) +
                        //                           mean_minus_3_sigma[k][2];
                        // }

                        // if (condition == "loose")
                        // {
                        //         upper_lim_total = mean_plus_4_sigma[k][0] * pow(_data->p(0), 2) + (mean_plus_4_sigma[k][1]) * _data->p(0) +
                        //                           mean_plus_4_sigma[k][2];

                        //         lower_lim_total = mean_minus_4_sigma[k][0] * pow(_data->p(0), 2) + (mean_minus_4_sigma[k][1]) * _data->p(0) +
                        //                           mean_minus_4_sigma[k][2];
                        // }
                }
        }
        bool pass_band = _data->ec_tot_energy(0) / _data->p(0) <= upper_lim_total &&
                         _data->ec_tot_energy(0) / _data->p(0) >= lower_lim_total;
        // bool pass_band = true;
        bool pass_triangle = false;

        if (_data->p(0) < 4.5)
        {
                pass_triangle = true;
        }
        else
        {
                // pass_triangle = (_data->ec_ecin_energy(0) / _data->p(0)) > (0.2 - _data->ec_pcal_energy(0) / _data->p(0));
                pass_triangle = true;
        }

        if (pass_band && pass_triangle)
                return true;
        else
                return false;
}
/////////////////////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_inner_vs_EC_outer()
{

        short isector = (_data->ec_pcal_sec(0) - 1);

        double param_a_exp[9][6] = {
            {0.201232, 0.197423, 0.197381, 0.18777, 0.197711, 0.198306},  // <2 GeV
            {0.207553, 0.20437, 0.209949, 0.198804, 0.208762, 0.209544},  // 2−3 GeV
            {0.213829, 0.215049, 0.215857, 0.209816, 0.218481, 0.221003}, // 3−4 GeV
            {0.217145, 0.218486, 0.21988, 0.215048, 0.221863, 0.22733},   // 4−5 GeV
            {0.220458, 0.218645, 0.220504, 0.218401, 0.222802, 0.230317}, // 5−6 GeV
            {0.22359, 0.219856, 0.226815, 0.221764, 0.222474, 0.23842},   // 6−7 GeV
            {0.226479, 0.219286, 0.22881, 0.225656, 0.22175, 0.242428},   // 7−8 GeV
            {0.226668, 0.219087, 0.226924, 0.228833, 0.21954, 0.24472},   // 8−9 GeV
            {0.225488, 0.221795, 0.21997, 0.228161, 0.216218, 0.246185}   // >9 GeV
        };

        double param_b_exp[9][6] = {
            {-1.00947, -0.975077, -0.988142, -0.892892, -0.97406, -0.981174}, // <2 GeV
            {-1.07716, -1.01246, -1.10285, -1.00481, -1.07982, -1.08578},     // 2−3 GeV
            {-1.10165, -1.08631, -1.1188, -1.0906, -1.15156, -1.17655},       // 3−4 GeV
            {-1.12186, -1.1006, -1.14603, -1.11891, -1.16233, -1.22691},      // 4−5 GeV
            {-1.16741, -1.10554, -1.15965, -1.13694, -1.1562, -1.26216},      // 5−6 GeV
            {-1.22825, -1.13945, -1.26374, -1.16584, -1.14528, -1.38061},     // 6−7 GeV
            {-1.30402, -1.16559, -1.33568, -1.2161, -1.13725, -1.47588},      // 7−8 GeV
            {-1.37011, -1.21059, -1.38956, -1.28444, -1.11712, -1.57429},     // 8−9 GeV
            {-1.45688, -1.30088, -1.4082, -1.34519, -1.09426, -1.67625}       // >9 GeV
        };

        int momRangeIndex = getMomRangeIndex(_data->p(0));
        double a = param_a_exp[momRangeIndex][isector];
        double b = param_b_exp[momRangeIndex][isector];
        double sf_ecin = _data->ec_ecin_energy(0) / _data->p(0);
        double sf_pcal = _data->ec_pcal_energy(0) / _data->p(0);
        // std::cout << "ec sec " << isector << " a  " << a << " b  " << b << '\n';

        return sf_pcal > a + b * sf_ecin;
}
/////////////////////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_hit_position_fiducial_cut_homogeneous(std::string condition)
{
        // Cut using the natural directions of the scintillator bars/ fibers:
        ///////////////////////////////////////////////////////////////////
        /// inbending:

        if (condition == "tight")
        {
                return (_data->ec_pcal_lv(0) > 14 && _data->ec_pcal_lv(0) < 400 && _data->ec_pcal_lw(0) > 14 &&
                        _data->ec_pcal_lw(0) < 400);
        }

        else if (condition == "loose")
        {
                return (_data->ec_pcal_lv(0) > 9 && _data->ec_pcal_lv(0) < 400 && _data->ec_pcal_lw(0) > 9 &&
                        _data->ec_pcal_lw(0) < 400);
        }
        else
        {
                return (_data->ec_pcal_lv(0) > 11.5 && _data->ec_pcal_lv(0) < 400 && _data->ec_pcal_lw(0) > 11.5 &&
                        _data->ec_pcal_lw(0) < 400);
        }
}

bool Pass2_Cuts::PCAL_Ineff_cut_X_Y()
{
        bool pcal_ineff_cuts = true;

        short pcal_sector = (_data->dc_sec(0));

        double X = _data->ec_pcal_x(0);
        double Y = _data->ec_pcal_y(0);

        if (pcal_sector == 1)
        {

                pcal_ineff_cuts &= (Y > 0.56575 * X - 92 + 0.25) || (Y < 0.56575 * X - 94.4 - 0.25 - 2);
                pcal_ineff_cuts &= (Y > 0.56575 * X - 101.1 - 0.25 + 2) || (Y < 0.56575 * X - 103.5 - 0.25);
                pcal_ineff_cuts &= (Y > 0.56575 * X - 219 + 0.25) || (Y < 0.56575 * X - 221.4 - 0.25);
                pcal_ineff_cuts &= (Y > 0.56575 * X - 227 + 0.25) || (Y < 0.56575 * X - 229.4 - 0.25);
        }
        else if (pcal_sector == 2)
        {

                pcal_ineff_cuts &= (Y > 0.5897 * X + 120 + 0.25 + 3.0) || (Y < 0.5913 * X + 114.3872 - 0.25);
        }
        else if (pcal_sector == 4)
        {

                pcal_ineff_cuts &= (Y > (-0.568) * X - 232.8 + 0.25 + 2) || (Y < (-0.568) * X - 236.3 - 0.25 - 2);
        }

        else if (pcal_sector == 6)
        {
                // std::cout << "   6  " << std::endl;

                pcal_ineff_cuts &= (Y > (-0.591377) * X - 185 + 0.25 + 1) || (Y < (-0.591377) * X - 187 - 0.25 - 1);
                pcal_ineff_cuts &= (Y > (-0.591377) * X - 193.3 + 0.25 + 1) || (Y < (-0.591377) * X - 195.5 - 0.25 - 0.25);
        }
        return pcal_ineff_cuts;
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::PCAL_fiducial_cut_X_Y(std::string condition)
{
        double minparams_pcal_in[6][2];
        double maxparams_pcal_in[6][2];
        double min_radious;
        double max_radious;

        if (condition == "mid")
        {
                // double minparams_pcal_in_m[6][2] = {{-0.53233, 21.92333}, {-0.53876, 22.35121}, {-0.50470, 21.58152}, {-0.52421, 24.02394}, {-0.52700, 25.79000}, {-0.52324, 21.14879}};
                // double maxparams_pcal_in_m[6][2] = {{0.51997, -20.34515}, {0.50942, -22.93788}, {0.52876, -25.21121}, {0.53930, -22.61848}, {0.52876, -25.21121}, {0.54148, -22.88758}};

                double minparams_pcal_in_m[6][2] = {{-0.55148, 22.58758}, {-0.55427, 22.51364}, {-0.52912, 23.30939}, {-0.54088, 23.79061}, {-0.54403, 25.97485}, {-0.54233, 21.62333}};
                double maxparams_pcal_in_m[6][2] = {{0.53815, -20.69424}, {0.52906, -23.67970}, {0.53876, -24.91121}, {0.55912, -23.52939}, {0.54227, -24.87364}, {0.55912, -23.52939}};

                min_radious = 72;
                max_radious = 270;
                std::copy(&minparams_pcal_in_m[0][0], &minparams_pcal_in_m[0][0] + 6 * 2, &minparams_pcal_in[0][0]);
                std::copy(&maxparams_pcal_in_m[0][0], &maxparams_pcal_in_m[0][0] + 6 * 2, &maxparams_pcal_in[0][0]);
        }
        else if (condition == "tight")
        {
                // double minparams_pcal_in_t[6][2] = {{-0.50197, 21.60515}, {-0.50633, 21.74333}, {-0.47652, 21.53242}, {-0.48470, 22.18152}, {-0.49209, 24.28455}, {-0.49221, 20.74394}};
                // double maxparams_pcal_in_t[6][2] = {{0.48403, -19.53485}, {0.48258, -23.14212}, {0.49688, -24.31061}, {0.50961, -22.54697}, {0.49221, -23.46394}, {0.51106, -22.37970}};

                double minparams_pcal_in_t[6][2] = {{-0.53233, 21.92333}, {-0.53876, 22.35121}, {-0.50470, 21.58152}, {-0.52421, 24.02394}, {-0.52700, 25.79000}, {-0.52324, 21.14879}};
                double maxparams_pcal_in_t[6][2] = {{0.51997, -20.34515}, {0.50942, -22.93788}, {0.52876, -25.21121}, {0.53930, -22.61848}, {0.52876, -25.21121}, {0.54148, -22.88758}};

                min_radious = 75;
                max_radious = 260;
                std::copy(&minparams_pcal_in_t[0][0], &minparams_pcal_in_t[0][0] + 6 * 2, &minparams_pcal_in[0][0]);
                std::copy(&maxparams_pcal_in_t[0][0], &maxparams_pcal_in_t[0][0] + 6 * 2, &maxparams_pcal_in[0][0]);
        }
        else if (condition == "loose")
        {
                double minparams_pcal_in_l[6][2] = {{-0.57476, 22.79121}, {-0.57403, 22.51485}, {-0.54324, 22.54879}, {-0.55803, 23.79485}, {-0.57082, 27.02091}, {-0.56615, 22.25424}};
                double maxparams_pcal_in_l[6][2] = {{0.55433, -20.38333}, {0.54524, -23.36879}, {0.56615, -25.85424}, {0.58221, -23.72394}, {0.56258, -24.98212}, {0.58494, -24.26030}};
                min_radious = 69;
                max_radious = 280;
                std::copy(&minparams_pcal_in_l[0][0], &minparams_pcal_in_l[0][0] + 6 * 2, &minparams_pcal_in[0][0]);
                std::copy(&maxparams_pcal_in_l[0][0], &maxparams_pcal_in_l[0][0] + 6 * 2, &maxparams_pcal_in[0][0]);
        }
        else
        {
                std::cerr << "Invalid condition: " << condition << ". Choose 'loose', 'mid', or 'tight'.\n";
                return false;
        }

        short pcal_sector = (_data->ec_pcal_sec(0) - 1);
        double X = _data->ec_pcal_x(0);
        double Y = _data->ec_pcal_y(0);

        float X_new = X * cos(DEG2RAD * (-60 * (pcal_sector))) - Y * sin(DEG2RAD * (-60 * (pcal_sector)));
        Y = X * sin(DEG2RAD * (-60 * (pcal_sector))) + Y * cos(DEG2RAD * (-60 * (pcal_sector)));
        X = X_new;

        double radious = sqrt(X * X + Y * Y);
        double calc_min = minparams_pcal_in[pcal_sector][0] * X + minparams_pcal_in[pcal_sector][1];
        double calc_max = maxparams_pcal_in[pcal_sector][0] * X + maxparams_pcal_in[pcal_sector][1];

        return ((Y > calc_min) && (Y < calc_max) && (radious > min_radious) && (radious < max_radious));
}

bool Pass2_Cuts::DC_Ineff_cut_X_Y(int i, int pid, std::string condition)
{
        bool _dc_ineff_cut = true;

        short dc_sector = (_data->dc_sec(i));
        // region 3
        double X3 = _data->dc_r3_x(i);
        double Y3 = _data->dc_r3_y(i);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector - 1))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector - 1)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector - 1))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector - 1)));
        X3 = X3_new;
        int region_3 = 3;

        // double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        // double Min_r3 = min_r[pid][region_3 - 1][dc_sector];
        // double Max_r3 = max_r[pid][region_3 - 1][dc_sector];
        float del_intercept = 0;

        if (condition == "tight")
        {
                del_intercept = 5;
        }
        else if (condition == "loose")
        {
                del_intercept = -5;
        }

        if (pid == 1) ///// for proton
        {
                _dc_ineff_cut &= X3 < 330;
                if (dc_sector == 1 || dc_sector == 4)
                {
                        _dc_ineff_cut &= ((Y3 > (1.36 * X3 - 515.8 + del_intercept)) && (Y3 < (-1.36 * X3 + 515.8 - del_intercept)));
                }
                else if (dc_sector == 2 || dc_sector == 5)
                {
                        _dc_ineff_cut &= ((Y3 > (1.11 * X3 - 425.7 + del_intercept)) && (Y3 < (-1.11 * X3 + 425.7 - del_intercept)));
                }
        }
        if (pid == 2) ////// for pip
        {
                _dc_ineff_cut &= X3 < 330;

                if (dc_sector == 1 || dc_sector == 4)
                {
                        _dc_ineff_cut &= ((Y3 > (1.36 * X3 - 505.8 + del_intercept)) && (Y3 < (-1.36 * X3 + 505.8 - del_intercept)));
                }
                else if (dc_sector == 2 || dc_sector == 5)
                {
                        _dc_ineff_cut &= ((Y3 > (1.11 * X3 - 425.7 + del_intercept)) && (Y3 < (-1.11 * X3 + 425.7 - del_intercept)));
                }
        }
        return _dc_ineff_cut;
}

bool Pass2_Cuts::DC_fiducial_cut_XY_E(std::string condition)

{
        bool _dc_fid_cut = true;

        /// supergaus old mid cuts, new tight
        double minparams_in_t[6][3][2] = {

            {{-0.62170, 15.22182}, {-0.62585, 22.49044}, {-0.56401, 22.36676}},
            {{-0.61648, 14.61909}, {-0.61152, 20.72934}, {-0.56335, 21.28434}},
            {{-0.58788, 14.42273}, {-0.60119, 22.33538}, {-0.53533, 21.37390}},
            {{-0.59939, 14.92364}, {-0.60844, 22.53813}, {-0.54571, 22.24121}},
            {{-0.61624, 16.14455}, {-0.61723, 24.04363}, {-0.55016, 23.68214}},
            {{-0.61624, 14.84455}, {-0.61829, 21.46714}, {-0.54973, 20.49643}}

        };

        double maxparams_in_t[6][3][2] = {
            {{0.60679, -14.11727}, {0.61275, -20.82582}, {0.56231, -21.30000}},
            {{0.60255, -15.79273}, {0.61130, -23.74527}, {0.54742, -23.90027}},
            {{0.66073, -18.50364}, {0.61798, -24.23231}, {0.56500, -25.25962}},
            {{0.62242, -15.04545}, {0.62954, -22.36560}, {0.56780, -21.79451}},
            {{0.64133, -17.42000}, {0.62492, -24.55736}, {0.56440, -24.99945}},
            {{0.61624, -14.84455}, {0.62725, -22.55989}, {0.56951, -22.53049}}};
        ////////// new mid cuts supergaus 0.8
        double minparams_in_m[6][3][2] = {
            {{-0.63636, 15.21182}, {-0.6298, 21.76077}, {-0.5783, 22.24478}},
            {{-0.65442, 15.88545}, {-0.62505, 20.63352}, {-0.57495, 20.93352}},
            {{-0.62739, 15.91364}, {-0.60237, 21.50648}, {-0.54385, 21.5}},
            {{-0.62764, 15.64818}, {-0.61288, 21.85912}, {-0.56181, 22.43049}},
            {{-0.6423, 17.05818}, {-0.61758, 23.27527}, {-0.56005, 23.53379}},
            {{-0.62933, 14.61}, {-0.61196, 20.12604}, {-0.56379, 20.73929}},
        };

        double maxparams_in_m[6][3][2] = {
            {{0.62388, -14.23273}, {0.60374, -19.16484}, {0.56929, -20.57225}},
            {{0.6423, -17.05818}, {0.60347, -22.08396}, {0.56544, -24.38379}},
            {{0.6423, -17.05818}, {0.65499, -26.27758}, {0.57165, -24.6522}},
            {{0.65442, -15.88545}, {0.64022, -22.33407}, {0.58049, -21.91951}},
            {{0.64824, -17.48455}, {0.63811, -24.87275}, {0.58137, -25.49863}},
            {{0.65988, -16.26273}, {0.63244, -21.9267}, {0.58978, -23.51099}},
        };

        // // ########### tight cuts 0.5 ##########
        // double maxparams_in_t[6][3][2] = {
        //     {{0.57079, -13.46727}, {0.58879, -20.73407}, {0.52813, -20.14341}},
        //     {{0.56073, -14.83364}, {0.58448, -23.36066}, {0.52522, -23.29478}},
        //     {{0.62473, -17.85364}, {0.59886, -24.36143}, {0.53742, -24.79643}},
        //     {{0.58242, -14.04545}, {0.59934, -21.59780}, {0.55747, -23.14945}},
        //     {{0.58776, -16.14545}, {0.59525, -23.87989}, {0.53313, -24.02225}},
        //     {{0.59309, -15.06545}, {0.60242, -22.46044}, {0.55495, -23.27198}}};
        // double minparams_in_t[6][3][2] = {
        //     {{-0.58933, 15.03000}, {-0.59903, 22.10582}, {-0.54181, 21.96126}},
        //     {{-0.58242, 14.04545}, {-0.59240, 20.85846}, {-0.55192, 22.29423}},
        //     {{-0.55927, 14.26636}, {-0.57437, 21.95077}, {-0.51302, 21.37775}},
        //     {{-0.56533, 14.35000}, {-0.58022, 22.08407}, {-0.52533, 22.67005}},
        //     {{-0.56412, 15.01727}, {-0.58703, 23.30440}, {-0.53396, 24.24066}},
        //     {{-0.57624, 13.84455}, {-0.59547, 21.72110}, {-0.53071, 20.85082}}};

        // ########### loose cuts 0.9 ##########
        double maxparams_in_l[6][3][2] = {
            {{0.63055, -13.52273}, {0.63187, -20.28242}, {0.58192, -20.20577}},
            {{0.63939, -15.92364}, {0.63354, -23.33275}, {0.57808, -24.01731}},
            {{0.68824, -18.48455}, {0.65174, -24.94198}, {0.58791, -24.84670}},
            {{0.69600, -17.11000}, {0.64369, -21.33945}, {0.59495, -21.86429}},
            {{0.68824, -18.48455}, {0.65174, -24.94198}, {0.59764, -25.69313}},
            {{0.67636, -16.21182}, {0.64655, -21.63231}, {0.60016, -22.64753}}};
        double minparams_in_l[6][3][2] = {
            {{-0.67636, 16.21182}, {-0.64941, 22.01088}, {-0.59093, 21.87830}},
            {{-0.69442, 16.88545}, {-0.63956, 20.79615}, {-0.59115, 21.10577}},
            {{-0.64776, 16.01545}, {-0.63086, 22.59857}, {-0.55813, 21.37802}},
            {{-0.66582, 16.68909}, {-0.63613, 22.65901}, {-0.57835, 22.67473}},
            {{-0.68230, 18.05818}, {-0.64479, 24.20264}, {-0.57643, 23.77280}},
            {{-0.64158, 14.47455}, {-0.63886, 21.11857}, {-0.57665, 20.33874}}};

        double min_r[3][6] =
            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {54.504, 54.487, 54.405, 54.644, 54.858, 54.992},
             {63.980, 64.513, 63.854, 63.939, 63.953, 64.713}

            };
        double max_r[3][6] =
            {{101, 101, 101, 101, 101, 101},
             {155, 155, 155, 155, 155, 155},
             {210, 210, 210, 210, 210, 210}};

        double minparams_in[6][3][2]; // Array to hold selected min cut
        double maxparams_in[6][3][2]; // Array to hold selected max cut

        if (condition == "tight")
        {
                std::copy(&minparams_in_t[0][0][0], &minparams_in_t[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_t[0][0][0], &maxparams_in_t[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }
        else if (condition == "loose")
        {
                std::copy(&minparams_in_l[0][0][0], &minparams_in_l[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_l[0][0][0], &maxparams_in_l[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }
        else
        {
                // Default to medium
                std::copy(&minparams_in_m[0][0][0], &minparams_in_m[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_m[0][0][0], &maxparams_in_m[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }

        short dc_sector = (_data->dc_sec(0) - 1);

        // region 1
        double X1 = _data->dc_r1_x(0);
        double Y1 = _data->dc_r1_y(0);
        float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

        X1 = X1_new;
        int region_1 = 1;

        double calc_min1 = minparams_in[dc_sector][region_1 - 1][0] * X1 + minparams_in[dc_sector][region_1 - 1][1];
        double calc_max1 = maxparams_in[dc_sector][region_1 - 1][0] * X1 + maxparams_in[dc_sector][region_1 - 1][1];
        double DC_r1 = sqrt(X1 * X1 + Y1 * Y1);
        double Min_r1 = min_r[region_1 - 1][dc_sector];
        double Max_r1 = max_r[region_1 - 1][dc_sector];

        // region 2
        double X2 = _data->dc_r2_x(0);
        double Y2 = _data->dc_r2_y(0);
        float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));
        X2 = X2_new;
        int region_2 = 2;
        double calc_min2 = minparams_in[dc_sector][region_2 - 1][0] * X2 + minparams_in[dc_sector][region_2 - 1][1];
        double calc_max2 = maxparams_in[dc_sector][region_2 - 1][0] * X2 + maxparams_in[dc_sector][region_2 - 1][1];
        double DC_r2 = sqrt(X2 * X2 + Y2 * Y2);
        double Min_r2 = min_r[region_2 - 1][dc_sector];
        double Max_r2 = max_r[region_2 - 1][dc_sector];

        // region 3
        double X3 = _data->dc_r3_x(0);
        double Y3 = _data->dc_r3_y(0);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));
        X3 = X3_new;
        int region_3 = 3;
        double calc_min3 = minparams_in[dc_sector][region_3 - 1][0] * X3 + minparams_in[dc_sector][region_3 - 1][1];
        double calc_max3 = maxparams_in[dc_sector][region_3 - 1][0] * X3 + maxparams_in[dc_sector][region_3 - 1][1];
        double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        double Min_r3 = min_r[region_3 - 1][dc_sector];
        double Max_r3 = max_r[region_3 - 1][dc_sector];

        // double Max_line = max_line[region_3 - 1][dc_sector];

        // /// data should be above
        // y = 1.36x−505.8
        //     /// data should be below
        //     y =−1.36x + 505.8

        return ((Y1 > calc_min1) && (Y1 < calc_max1) && (Y2 > calc_min2) && (Y2 < calc_max2) && (Y3 > calc_min3) &&
                (Y3 < calc_max3) && (DC_r1 > Min_r1) && (DC_r2 > Min_r2) && (DC_r3 > Min_r3) // && (DC_r1 < Max_r1) && (DC_r2 < Max_r2)
                && (DC_r3 < Max_r3));
        //////    //&& (Y3 > (1.36 * X3 - 505.8)) && (Y3 < (-1.36 * X3 + 505.8)));
}

bool Pass2_Cuts::DC_fiducial_cut_XY_PIP(int i, int pid, std::string condition)
{
        bool _dc_fid_cut = true;
        // bool isinbending = true;
        // new cut parameters for the linear cut based on x and y coordinates (inbending field):
        // replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
        // (optimized for electrons, do not use it for hadrons)

        // maxparams_in[1][6][3][2] -> [pid][sec][regions][a,b]->a*x+b

        /// supergaus ld mid cuts but new tight cuts
        double minparams_in_pip_t[6][3][2] = {
            {{-0.48371, 2.71544}, {-0.56539, 14.80074}, {-0.56644, 28.90231}},
            {{-0.49800, 3.45000}, {-0.56931, 14.84191}, {-0.56402, 27.89462}},
            {{-0.49165, 3.19853}, {-0.55265, 13.67279}, {-0.56789, 30.39264}},
            {{-0.48529, 2.72206}, {-0.56480, 14.85221}, {-0.57158, 30.11088}},
            {{-0.46906, 1.67941}, {-0.55137, 13.52941}, {-0.54323, 24.85703}},
            {{-0.51318, 4.11324}, {-0.56931, 14.84191}, {-0.56697, 28.84165}}

        };

        double maxparams_in_pip_t[6][3][2] = {
            {{0.49447, -3.13529}, {0.55402, -12.47426}, {0.57062, -26.22209}},
            {{0.49800, -3.45000}, {0.56431, -13.96176}, {0.57756, -28.62341}},
            {{0.49235, -2.92647}, {0.55402, -12.61544}, {0.56486, -26.74143}},
            {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56947, -26.46780}},
            {{0.48971, -2.89044}, {0.55402, -12.61544}, {0.55198, -24.45110}},
            {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56565, -25.70758}}};

        ///// new loose cuts supergaus 0.7
        double minparams_in_pip_l[6][3][2] = {
            {{-0.499, 2.4125}, {-0.56167, 12.55809}, {-0.57233, 26.73923}},
            {{-0.51329, 3.14706}, {-0.58863, 15.01471}, {-0.57629, 27.02714}},
            {{-0.50694, 2.89559}, {-0.56863, 13.47059}, {-0.57268, 27.64165}},
            {{-0.50112, 2.62132}, {-0.57157, 13.90735}, {-0.57431, 26.9689}},
            {{-0.48324, 1.40515}, {-0.54333, 10.50735}, {-0.55136, 22.63615}},
            {{-0.53182, 4.19926}, {-0.58167, 14.10221}, {-0.57035, 25.99527}},
        };

        double maxparams_in_pip_l[6][3][2] = {
            {{0.51171, -3.14044}, {0.57235, -12.46618}, {0.58156, -25.07769}},
            {{0.51382, -3.34926}, {0.57814, -13.3875}, {0.5758, -24.56846}},
            {{0.50865, -2.86103}, {0.57814, -13.3875}, {0.57778, -25.31242}},
            {{0.51612, -3.30882}, {0.57235, -12.46618}, {0.58174, -25.60033}},
            {{0.50035, -2.30147}, {0.56824, -12.07353}, {0.56051, -21.95615}},
            {{0.50512, -2.54632}, {0.57167, -12.51838}, {0.57145, -23.52604}},
        };

        double minparams_in_pip_m[6][3][2] = {
            {{-0.48929, 2.29706}, {-0.55794, 13.11544}, {-0.57725, 29.25879}},
            {{-0.51312, 3.52132}, {-0.56824, 13.95588}, {-0.57756, 28.62341}},
            {{-0.49588, 2.89118}, {-0.56137, 13.56029}, {-0.5747, 29.75198}},
            {{-0.48982, 2.22426}, {-0.56039, 13.46176}, {-0.56569, 27.24538}},
            {{-0.48765, 2.09853}, {-0.55402, 12.61544}, {-0.55198, 24.4511}},
            {{-0.51529, 3.52206}, {-0.58, 15.10294}, {-0.58033, 29.70495}}};

        double maxparams_in_pip_m[6][3][2] = {
            {{0.50465, -3.11103}, {0.57686, -13.75}, {0.57831, -25.9089}},
            {{0.50771, -3.39044}, {0.5751, -14.375}, {0.57114, -25.81857}},
            {{0.50788, -3.19118}, {0.5499, -11.41103}, {0.57114, -25.81857}},
            {{0.51053, -3.27721}, {0.57686, -13.75}, {0.58574, -27.74033}},
            {{0.48671, -2.05294}, {0.55775, -12.05809}, {0.54938, -21.8922}},
            {{0.50788, -3.19118}, {0.57686, -13.75}, {0.5533, -21.87088}}};

        double min_r[3][6] =

            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {66, 66, 66, 66, 66, 66},
             {135, 135, 135, 135, 135, 135}

            };

        double max_r[3][6] =
            {{140, 145, 145, 145, 145, 140},
             {210, 220, 220, 220, 220, 210},
             {325, 340, 340, 340, 340, 325}};

        double minparams_pip_in[6][3][2]; // Array to hold selected min cut
        double maxparams_pip_in[6][3][2]; // Array to hold selected max cut

        if (condition == "tight")
        {
                std::copy(&minparams_in_pip_t[0][0][0], &minparams_in_pip_t[0][0][0] + 6 * 3 * 2, &minparams_pip_in[0][0][0]);
                std::copy(&maxparams_in_pip_t[0][0][0], &maxparams_in_pip_t[0][0][0] + 6 * 3 * 2, &maxparams_pip_in[0][0][0]);
        }
        else if (condition == "loose")
        {
                std::copy(&minparams_in_pip_l[0][0][0], &minparams_in_pip_l[0][0][0] + 6 * 3 * 2, &minparams_pip_in[0][0][0]);
                std::copy(&maxparams_in_pip_l[0][0][0], &maxparams_in_pip_l[0][0][0] + 6 * 3 * 2, &maxparams_pip_in[0][0][0]);
        }
        else
        {
                // Default to medium
                std::copy(&minparams_in_pip_m[0][0][0], &minparams_in_pip_m[0][0][0] + 6 * 3 * 2, &minparams_pip_in[0][0][0]);
                std::copy(&maxparams_in_pip_m[0][0][0], &maxparams_in_pip_m[0][0][0] + 6 * 3 * 2, &maxparams_pip_in[0][0][0]);
        }

        short dc_sector = (_data->dc_sec(i) - 1);

        // region 1
        double X1 = _data->dc_r1_x(i);
        double Y1 = _data->dc_r1_y(i);
        float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

        X1 = X1_new;
        int region_1 = 1;

        double calc_min1 = minparams_pip_in[dc_sector][region_1 - 1][0] * X1 + minparams_pip_in[dc_sector][region_1 - 1][1];
        double calc_max1 = maxparams_pip_in[dc_sector][region_1 - 1][0] * X1 + maxparams_pip_in[dc_sector][region_1 - 1][1];
        double DC_r1 = sqrt(X1 * X1 + Y1 * Y1);
        double Min_r1 = min_r[region_1 - 1][dc_sector];
        double Max_r1 = max_r[region_1 - 1][dc_sector];

        // region 2
        double X2 = _data->dc_r2_x(i);
        double Y2 = _data->dc_r2_y(i);
        float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));
        X2 = X2_new;
        int region_2 = 2;
        double calc_min2 = minparams_pip_in[dc_sector][region_2 - 1][0] * X2 + minparams_pip_in[dc_sector][region_2 - 1][1];
        double calc_max2 = maxparams_pip_in[dc_sector][region_2 - 1][0] * X2 + maxparams_pip_in[dc_sector][region_2 - 1][1];
        double DC_r2 = sqrt(X2 * X2 + Y2 * Y2);
        double Min_r2 = min_r[region_2 - 1][dc_sector];
        double Max_r2 = max_r[region_2 - 1][dc_sector];

        // region 3
        double X3 = _data->dc_r3_x(i);
        double Y3 = _data->dc_r3_y(i);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));
        X3 = X3_new;
        int region_3 = 3;
        double calc_min3 = minparams_pip_in[dc_sector][region_3 - 1][0] * X3 + minparams_pip_in[dc_sector][region_3 - 1][1];
        double calc_max3 = maxparams_pip_in[dc_sector][region_3 - 1][0] * X3 + maxparams_pip_in[dc_sector][region_3 - 1][1];
        double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        double Min_r3 = min_r[region_3 - 1][dc_sector];
        double Max_r3 = max_r[region_3 - 1][dc_sector];
        // double Max_line = max_line[region_3 - 1][dc_sector];

        // /// data should be above
        // y = 1.36x−505.8
        //     /// data should be below
        //     y =−1.36x + 505.8

        return ((Y1 > calc_min1) && (Y1 < calc_max1) && (Y2 > calc_min2) && (Y2 < calc_max2) && (Y3 > calc_min3) &&
                (Y3 < calc_max3) && (DC_r1 > Min_r1) && (DC_r2 > Min_r2) && (DC_r3 > Min_r3) // && (DC_r1 < Max_r1) && (DC_r2 < Max_r2)
                && (DC_r3 < Max_r3));
        //////    //&& (Y3 > (1.36 * X3 - 505.8)) && (Y3 < (-1.36 * X3 + 505.8)));

        // // if(inbending == true) pid = 0; // use only for electrons in inbending case
        // double calc_min = minparams[dc_sector - 1][region - 1] + minparams[dc_sector - 1][region - 1][1] *
        // X; double calc_max = maxparams[dc_sector - 1][region - 1] + maxparams[dc_sector - 1][region -
        // 1][1] * X; return (Y > calc_min) && (Y < calc_max);
}

bool Pass2_Cuts::DC_fiducial_cut_XY_PROT(int i, int pid, std::string condition)
{
        bool _dc_fid_cut = true;
        // bool isinbending = true;
        // new cut parameters for the linear cut based on x and y coordinates (inbending field):
        // replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
        // (optimized for electrons, do not use it for hadrons)

        // maxparams_in[1][6][3][2] -> [sec][regions][a,b]->a*x+b

        /// supergaus
        double minparams_in_prot_t[6][3][2] = {

            {{-0.51812, 4.80882}, {-0.60660, 19.32861}, {-0.60392, 34.35024}},
            {{-0.52912, 5.39632}, {-0.61995, 20.57325}, {-0.62008, 37.14595}},
            {{-0.50012, 4.10882}, {-0.59236, 18.29823}, {-0.58771, 32.53500}},
            {{-0.52647, 5.16029}, {-0.60426, 18.88900}, {-0.60434, 34.83024}},
            {{-0.50600, 4.00000}, {-0.58717, 17.29909}, {-0.57177, 28.34119}},
            {{-0.52388, 4.91618}, {-0.60842, 19.32450}, {-0.61696, 36.42214}}

        };

        double maxparams_in_prot_t[6][3][2] = {

            {{0.53865, -5.86103}, {0.61408, -19.00965}, {0.64226, -38.73976}},
            {{0.52876, -5.36985}, {0.61309, -19.36563}, {0.61519, -34.45357}},
            {{0.52388, -4.91618}, {0.60462, -18.49771}, {0.61353, -34.83833}},
            {{0.53606, -5.61691}, {0.61106, -18.60558}, {0.62413, -35.58310}},
            {{0.51347, -4.49779}, {0.60016, -18.02074}, {0.59738, -31.99500}},
            {{0.52918, -5.16324}, {0.60592, -18.18320}, {0.61468, -33.58214}}};

        double minparams_in_prot_m[6][3][2] = {
            {{-0.52224, 4.54265}, {-0.61242, 19.00117}, {-0.60912, 34.15024}},
            {{-0.53994, 5.65809}, {-0.62951, 20.80061}, {-0.63005, 36.81833}},
            {{-0.51629, 4.38456}, {-0.59927, 18.18736}, {-0.60413, 34.00929}},
            {{-0.52135, 4.41397}, {-0.6119, 19.04316}, {-0.61312, 35.025}},
            {{-0.50665, 3.48603}, {-0.59262, 17.04866}, {-0.59332, 31.3769}},
            {{-0.53229, 4.78456}, {-0.61808, 19.38156}, {-0.61987, 35.31548}},
        };

        double maxparams_in_prot_m[6][3][2] = {
            {{0.53441, -5.04338}, {0.62712, -19.26043}, {0.63847, -36.1431}},
            {{0.53765, -5.52353}, {0.61662, -18.66104}, {0.61862, -33.53262}},
            {{0.53271, -5.07794}, {0.61496, -18.45255}, {0.62236, -34.3669}},
            {{0.53235, -4.97647}, {0.61366, -18.15753}, {0.62465, -34.37833}},
            {{0.52753, -4.78971}, {0.60468, -17.54113}, {0.61862, -33.53262}},
            {{0.53235, -4.97647}, {0.6121, -17.94065}, {0.60418, -30.45738}},
        };

        double minparams_in_prot_l[6][3][2] = {
            {{-0.52806, 4.16691}, {-0.61506, 18.32987}, {-0.60979, 31.99452}},
            {{-0.54824, 5.39265}, {-0.63055, 19.91662}, {-0.6259, 34.885}},
            {{-0.51006, 3.46691}, {-0.60047, 17.34792}, {-0.60839, 33.02929}},
            {{-0.536, 4.5}, {-0.61532, 18.2803}, {-0.62636, 35.24167}},
            {{-0.52124, 3.65515}, {-0.60244, 17.05502}, {-0.58501, 27.34833}},
            {{-0.54182, 4.77426}, {-0.62618, 19.33602}, {-0.62382, 34.62786}},
        };

        double maxparams_in_prot_l[6][3][2] = {
            {{0.54341, -4.75588}, {0.62384, -18.18212}, {0.64081, -35.13595}},
            {{0.55029, -5.60956}, {0.62223, -18.57896}, {0.6461, -36.40357}},
            {{0.53388, -4.29118}, {0.62088, -18.33576}, {0.6233, -33.05167}},
            {{0.547, -5.1375}, {0.62275, -18.00364}, {0.64008, -34.99595}},
            {{0.52329, -3.72206}, {0.60743, -16.84238}, {0.61312, -31.29167}},
            {{0.54129, -4.74706}, {0.62275, -18.00364}, {0.63312, -33.59881}},
        };

        double min_r[3][6] =

            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {56, 56, 56, 56, 56, 56},
             {135, 135, 135, 135, 135, 135}

            };

        double max_r[3][6] =
            {{145, 145, 145, 145, 145, 145},
             {220, 220, 220, 220, 220, 220},
             {350, 350, 350, 350, 350, 350}};

        double minparams_in[6][3][2]; // Array to hold selected min cut
        double maxparams_in[6][3][2]; // Array to hold selected max cut

        if (condition == "tight")
        {
                std::copy(&minparams_in_prot_t[0][0][0], &minparams_in_prot_t[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_prot_t[0][0][0], &maxparams_in_prot_t[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }
        else if (condition == "loose")
        {
                std::copy(&minparams_in_prot_l[0][0][0], &minparams_in_prot_l[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_prot_l[0][0][0], &maxparams_in_prot_l[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }
        else
        {
                // Default to medium
                std::copy(&minparams_in_prot_m[0][0][0], &minparams_in_prot_m[0][0][0] + 6 * 3 * 2, &minparams_in[0][0][0]);
                std::copy(&maxparams_in_prot_m[0][0][0], &maxparams_in_prot_m[0][0][0] + 6 * 3 * 2, &maxparams_in[0][0][0]);
        }
        short dc_sector = (_data->dc_sec(i) - 1);

        // region 1
        double X1 = _data->dc_r1_x(i);
        double Y1 = _data->dc_r1_y(i);
        float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

        X1 = X1_new;
        int region_1 = 1;

        double calc_min1 = minparams_in[dc_sector][region_1 - 1][0] * X1 + minparams_in[dc_sector][region_1 - 1][1];
        double calc_max1 = maxparams_in[dc_sector][region_1 - 1][0] * X1 + maxparams_in[dc_sector][region_1 - 1][1];
        double DC_r1 = sqrt(X1 * X1 + Y1 * Y1);
        double Min_r1 = min_r[region_1 - 1][dc_sector];
        double Max_r1 = max_r[region_1 - 1][dc_sector];

        // region 2
        double X2 = _data->dc_r2_x(i);
        double Y2 = _data->dc_r2_y(i);
        float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));
        X2 = X2_new;
        int region_2 = 2;
        double calc_min2 = minparams_in[dc_sector][region_2 - 1][0] * X2 + minparams_in[dc_sector][region_2 - 1][1];
        double calc_max2 = maxparams_in[dc_sector][region_2 - 1][0] * X2 + maxparams_in[dc_sector][region_2 - 1][1];
        double DC_r2 = sqrt(X2 * X2 + Y2 * Y2);
        double Min_r2 = min_r[region_2 - 1][dc_sector];
        double Max_r2 = max_r[region_2 - 1][dc_sector];

        // region 3
        double X3 = _data->dc_r3_x(i);
        double Y3 = _data->dc_r3_y(i);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));
        X3 = X3_new;
        int region_3 = 3;
        double calc_min3 = minparams_in[dc_sector][region_3 - 1][0] * X3 + minparams_in[dc_sector][region_3 - 1][1];
        double calc_max3 = maxparams_in[dc_sector][region_3 - 1][0] * X3 + maxparams_in[dc_sector][region_3 - 1][1];
        double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        double Min_r3 = min_r[region_3 - 1][dc_sector];
        double Max_r3 = max_r[region_3 - 1][dc_sector];
        // double Max_line = max_line[region_3 - 1][dc_sector];

        // /// data should be above
        // y = 1.36x−505.8
        //     /// data should be below
        //     y =−1.36x + 505.8

        return ((Y1 > calc_min1) && (Y1 < calc_max1) && (Y2 > calc_min2) && (Y2 < calc_max2) && (Y3 > calc_min3) &&
                (Y3 < calc_max3) && (DC_r1 > Min_r1) && (DC_r2 > Min_r2) && (DC_r3 > Min_r3) // && (DC_r1 < Max_r1) && (DC_r2 < Max_r2)
                && (DC_r3 < Max_r3));
        //////    //&& (Y3 > (1.36 * X3 - 505.8)) && (Y3 < (-1.36 * X3 + 505.8)));

        // // if(inbending == true) pid = 0; // use only for electrons in inbending case
        // double calc_min = minparams[pid][dc_sector - 1][region - 1][pid] + minparams[pid][dc_sector - 1][region - 1][1] *
        // X; double calc_max = maxparams[pid][dc_sector - 1][region - 1][pid] + maxparams[pid][dc_sector - 1][region -
        // 1][1] * X; return (Y > calc_min) && (Y < calc_max);
}

bool Pass2_Cuts::DC_z_vertex_cut(std::string condition)
{
        float partvz = _data->vz(0);

        if (condition == "tight")
        {
                return partvz > -10 && partvz < 5;
        }

        else if (condition == "loose")
        {
                return partvz > -12 && partvz < 7;
        }
        else
        {
                return partvz > -11 && partvz < 6;
        }
}

// public class HadronCuts {
//
// /**
//  * DC fiducial cut for hadrons
//  * @param dc_sector sector of hits in DC
//  * @param region specify fiducial Pass2_Cuts for which region to use
//  * @param trajx x for region 1 or 2 or 3 from REC::Traj
//  * @param trajy y for region 1 or 2 or 3 from REC::Traj
//  * @param trajz z for region 1 or 2 or 3 from REC::Traj
//  * @param partpid pid assigned to particle candidate
//  * @param isinbending True if magnetic field is inbending
//  */

/** Delta VZ cut for hadrons
 * @param pid hadron PID code
 * @param dvz difference between Vz of hadron candidate and electron
 */
bool Pass2_Cuts::Hadron_Delta_vz_cut(int i, std::string condition)
{
        int pid = _data->pid(i);
        float dvz = (_data->vz(i) - _data->vz(0));

        return dvz > -20 && dvz < 20;

        if (condition == "tight")
        {
                return dvz > -20 && dvz < 20;
        }

        else if (condition == "loose")
        {
                return dvz > -22 && dvz < 22;
        }
        else
        {
                return dvz > -21 && dvz < 21;
        }
}

/** chi2pid cut for hadrons
 * @param chi2pid chi2pid value
 * @param pid hadron PID code
 */
bool Pass2_Cuts::Hadron_Chi2pid_cut(int i, std::string condition)
{
        bool isstrict = false;
        float chi2pid = _data->chi2pid(i);
        float p = _data->p(i);
        int pid = _data->pid(i);
        int status = abs(_data->status(i));
        if (condition == "loose")
        {
                if (status < 4000)
                        return abs(chi2pid) < 6.0;
                else
                {
                        return abs(chi2pid) < 8.0;
                }
        }
        if (condition == "tight")
        {
                if (status < 4000)
                        return abs(chi2pid) < 5.0;
                else
                {
                        return abs(chi2pid) < 7.0;
                }
        }
        else
        {
                if (status < 4000)
                        return abs(chi2pid) < 5.5;
                else
                {
                        return abs(chi2pid) < 7.5;
                }
        }
        // double coef;
        // if (pid == 211)
        //         coef = 0.88;
        // else if (pid == -211)
        //         coef = 0.93;

        // else if (pid == 2212)
        // {
        //         if (status < 4000)
        //                 return abs(chi2pid) < 3.0; /// please confirm this first 2.64 is given for rga fall 2018
        //         else
        //         {
        //                 return abs(chi2pid) < 6.0;
        //         }
        // }

        // else
        //         return false;

        // bool chi2cut = false;
        // if (status < 4000)
        // {
        //         if (isstrict)
        //         {
        //                 if (p < 2.44)
        //                         chi2cut = chi2pid < 3 * coef;
        //                 else if (p < 4.6)
        //                         chi2cut = chi2pid < coef * (0.00869 + 14.98587 * exp(-p / 1.18236) + 1.81751 * exp(-p / 4.86394));
        //                 else
        //                         chi2cut = chi2pid < coef * (-1.14099 + 24.14992 * exp(-p / 1.36554) + 2.66876 * exp(-p / 6.80522));
        //         }
        //         else
        //         {
        //                 if (p < 2.44)
        //                         chi2cut = chi2pid < 3 * coef;
        //                 else
        //                         chi2cut = chi2pid < coef * (0.00869 + 14.98587 * exp(-p / 1.18236) + 1.81751 * exp(-p / 4.86394));
        //         }

        //         return chi2cut && chi2pid > coef * -3;
        // }
        // else
        // {
        //         return abs(chi2pid) < 6.0;
        // }
}
// bool Pass2_Cuts::CD_fiducial_Prot(double phi, double theta, double momT)
bool Pass2_Cuts::CD_fiducial_had(int i, std::string condition)
{
        bool pass_fiducial = true;
        // int pid = _data->pid(i);
        //        if (pid == 2212)
        {
                double momT = sqrt(_data->px(i) * _data->px(i) + _data->py(i) * _data->py(i));
                double theta = atan2(momT, _data->pz(i)) * 180 / PI;
                double phi = atan2(_data->py(i), _data->px(i)) * 180 / PI;

                double fiducial_phi_width; // = 3; // 3 is used by andrew
                if (condition == "tight")
                {
                        fiducial_phi_width = 4;
                }
                if (condition == "loose")
                {
                        fiducial_phi_width = 2;
                }
                else
                {
                        fiducial_phi_width = 3;
                }

                double fiducial_phi_shift = 0;
                double fiducial_momT_start = 0.15;
                double fiducial_phi_central = (-asin(fiducial_momT_start / momT) - (PI / 2)) * 180 / PI;

                if ((fabs(phi - fiducial_phi_central - fiducial_phi_shift) < fiducial_phi_width) ||
                    (fabs(phi - fiducial_phi_central - fiducial_phi_shift - 120) < fiducial_phi_width) ||
                    (fabs(phi - fiducial_phi_central - fiducial_phi_shift - 240) < fiducial_phi_width)) // || (theta < 40) || (theta > 125))
                {
                        pass_fiducial = false;
                        // std::cout << "phi   " << phi << "   momT   " << momT << std::endl;
                        // std::cout << " phi - fiducial_phi_central   :   " << abs(phi - fiducial_phi_central) << std::endl;
                        // std::cout << " phi - fiducial_phi_cent - 120:   " << abs(phi - fiducial_phi_central - 120) << std::endl;
                        // std::cout << " phi - fiducial_phi_cent -240 :   " << abs(phi - fiducial_phi_central - 240) << std::endl;
                }
        }
        return pass_fiducial;
}

///////////////////// Pass2_Cuts ///////////////////////
