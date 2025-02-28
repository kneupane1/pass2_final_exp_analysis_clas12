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

// // // /////////////////////////// exp data dt cuts ////////////// prot, pip, pim /////////////////
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

bool Pass2_Cuts::IsPip(int i)
{
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
                _pip &= (_dt->dt_Pi(i) < (dt_cut_fd_up[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][1][1] * pow(_data->p(i), 4) +
                                          dt_cut_fd_up[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][1][3] * pow(_data->p(i), 2) +
                                          dt_cut_fd_up[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][1][5]));
                _pip &= (_dt->dt_Pi(i) > (dt_cut_fd_down[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][1][1] * pow(_data->p(i), 4) +
                                          dt_cut_fd_down[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][1][3] * pow(_data->p(i), 2) +
                                          dt_cut_fd_down[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][1][5]));

                _pip &= DC_fiducial_cut_XY(i, 2);
                _pip &= DC_Ineff_cut_X_Y(i, 2);
        }
        // }
        else if (abs(_data->status(i)) >= 4000)
        // else if ((_dt->isCtof()))
        {
                // std::cout << "   dt pi at ctof  " << _dt->dt_Pi(i) << std::endl;

                _pip &= (_data->p(i) > 0.2);
                // _pip &= (_data->p(i) < 1.7);
                _pip &= (_dt->dt_Pi(i) < (dt_cut_cd_up[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][1][1] * _data->p(i) + dt_cut_cd_up[is_mc][1][2]));
                _pip &= (_dt->dt_Pi(i) > (dt_cut_cd_down[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][1][1] * _data->p(i) + dt_cut_cd_down[is_mc][1][2]));
                _pip &= CD_fiducial_had(i);
        }
        // _pip &= (_data->p(i) > 0.2);
        _pip &= Hadron_Delta_vz_cut(i);
        _pip &= Hadron_Chi2pid_cut(i);
        return _pip;
}
bool Pass2_Cuts::IsProton(int i)
{
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
                _proton &= (_dt->dt_P(i) < (dt_cut_fd_up[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][0][1] * pow(_data->p(i), 4) +
                                            dt_cut_fd_up[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][0][3] * pow(_data->p(i), 2) +
                                            dt_cut_fd_up[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][0][5]));

                _proton &= (_dt->dt_P(i) > (dt_cut_fd_down[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][0][1] * pow(_data->p(i), 4) +
                                            dt_cut_fd_down[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][0][3] * pow(_data->p(i), 2) +
                                            dt_cut_fd_down[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][0][5]));

                _proton &= DC_fiducial_cut_XY(i, 1);
                _proton &= DC_Ineff_cut_X_Y(i, 1);
        }
        // }
        else if (abs(_data->status(i)) >= 4000)
        // else if ((_dt->isCtof()))
        {
                // std::cout << "   dt prot at ctof  " << _dt->dt_P(i) << std::endl;

                _proton &= (_data->p(i) > 0.2); /// this 0.4 look harse when we do missing Pim channel
                // _proton &= (_data->p(i) < 2.0);
                _proton &= (_dt->dt_P(i) < (dt_cut_cd_up[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][0][1] * _data->p(i) + dt_cut_cd_up[is_mc][0][2]));
                _proton &= (_dt->dt_P(i) > (dt_cut_cd_down[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][0][1] * _data->p(i) + dt_cut_cd_down[is_mc][0][2]));
                _proton &= CD_fiducial_had(i);
        }
        // _proton &= (_data->p(i) > 0.2);
        // _proton &= (abs(_data->chi2d(i)) < 0.5);
        _proton &= Hadron_Delta_vz_cut(i);
        _proton &= Hadron_Chi2pid_cut(i);
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

        _pim &= Hadron_Delta_vz_cut(i); /// this because of the fact that hadron cuts are removed for pim
        _pim &= Hadron_Chi2pid_cut(i);  /// this because of the fact that hadron cuts are removed for pim
        // _pim &= DC_fiducial_cut_XY(0, 0); /// this because of the fact that hadron cuts are removed for pim, same dc cuts as electrons

        return _pim;
}
// /////////////////////// Pass2_Cuts ///////////////////////
bool Pass2_Cuts::ElectronCuts()
{
        bool cut = true;
        if (!cut)
                return false;
        cut &= (_data->gpart() > 0);
        cut &= (_data->gpart() < 20);
        // // //
        cut &= (_data->pid(0) == ELECTRON);
        cut &= (_data->charge(0) == NEGATIVE);
        cut &= DC_z_vertex_cut();
        cut &= (_data->p(0) > 1.50);
        // // cut &= (_data->vz(0) > -(2.78 + 3 * 2.16) && _data->vz(0) < (-2.78 + 3 * 2.16)); // 3 sigma cut
        cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
        // // cut &= (abs(_data->chi2pid(0)) < 3); ////////////// check it....... along with simulations
        // // cut &= CC_nphe_cut();
        cut &= DC_fiducial_cut_XY(0, 0);
        cut &= EC_sampling_fraction_cut();
        cut &= PCAL_minimum_energy();
        cut &= PCAL_fiducial_cut_X_Y();
        // // //cut &= EC_inner_vs_EC_outer();
        cut &= EC_hit_position_fiducial_cut_homogeneous();
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

bool Pass2_Cuts::EC_sampling_fraction_cut()
{
        int isec = (_data->ec_pcal_sec(0) - 1);
        double upper_lim_total = 0;
        double lower_lim_total = 0;
        int is_mc = 0;
        if (_mc)
        {
                is_mc = 1;
        }

        ///// Both  ////////
        double mean_minus_3_5_sigma[2][6][3] = {{{-0.001124, 0.01212, 0.1616}, {-0.001773, 0.02165, 0.1334}, {-0.002983, 0.03586, 0.0957}, {-0.001408, 0.01909, 0.1344}, {-0.0007095, 0.00882, 0.1687}, {-0.002068, 0.02394, 0.1298}},
                                                {{-0.00058, 0.00687, 0.19312}, {-0.00088, 0.01022, 0.18360}, {-0.00089, 0.00941, 0.18832}, {-0.00066, 0.00888, 0.18466}, {-0.00066, 0.00798, 0.18884}, {-0.00055, 0.00685, 0.19319}}};
        double mean_plus_3_5_sigma[2][6][3] = {{{-0.001306, 0.012856, 0.2598}, {-0.00082, 0.00801, 0.271}, {-0.00118, 0.01078, 0.269}, {-0.0009036, 0.01098, 0.2607}, {-0.0001322, 0.002972, 0.2744}, {-0.001183, 0.01233, 0.257}},
                                               {{-0.00002, -0.00078, 0.29991}, {0.00023, -0.00396, 0.31026}, {0.00010, -0.00156, 0.30077}, {0.00017, -0.00400, 0.31052}, {0.00018, -0.00342, 0.30823}, {0.00012, -0.00297, 0.30706}}};

        for (Int_t k = 0; k < 6; k++)
        {
                if (isec == k)
                {
                        upper_lim_total = mean_plus_3_5_sigma[is_mc][k][0] * pow(_data->p(0), 2) + (mean_plus_3_5_sigma[is_mc][k][1]) * _data->p(0) +
                                          mean_plus_3_5_sigma[is_mc][k][2];

                        lower_lim_total = mean_minus_3_5_sigma[is_mc][k][0] * pow(_data->p(0), 2) + (mean_minus_3_5_sigma[is_mc][k][1]) * _data->p(0) +
                                          mean_minus_3_5_sigma[is_mc][k][2];
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

        // double param_a_sim[9][6] = {
        //     {0.20233, 0.202647, 0.200986, 0.201774, 0.201084, 0.201648},  // <2 GeV
        //     {0.212437, 0.21288, 0.211826, 0.213158, 0.212129, 0.212346},  // 2−3 GeV
        //     {0.219554, 0.220162, 0.219478, 0.220034, 0.219012, 0.219218}, // 3−4 GeV
        //     {0.224078, 0.224914, 0.223814, 0.224357, 0.224371, 0.224192}, // 4−5 GeV
        //     {0.22785, 0.228319, 0.227713, 0.228105, 0.22709, 0.227346},   // 5−6 GeV
        //     {0.230326, 0.230837, 0.230451, 0.230969, 0.229402, 0.229665}, // 6−7 GeV
        //     {0.23258, 0.233087, 0.2319, 0.233145, 0.232148, 0.232072},    // 7−8 GeV
        //     {0.232341, 0.233287, 0.231747, 0.232848, 0.231751, 0.232098}, // 8−9 GeV
        //     {0.22484, 0.233669, 0.236727, 0.233346, 0.233825, 0.23524}    // >9 GeV
        // };

        // double param_b_sim[9][6] = {
        //     {-0.949695, -0.956382, -0.934262, -0.940217, -0.933617, -0.936363}, // <2 GeV
        //     {-1.04219, -1.04446, -1.035, -1.04832, -1.03597, -1.03963},         // 2−3 GeV
        //     {-1.08307, -1.09327, -1.08529, -1.08885, -1.07962, -1.08145},       // 3−4 GeV
        //     {-1.11223, -1.11772, -1.10889, -1.11142, -1.11528, -1.11228},       // 4−5 GeV
        //     {-1.13177, -1.13684, -1.13048, -1.13521, -1.12421, -1.12983},       // 5−6 GeV
        //     {-1.14596, -1.14755, -1.14677, -1.15218, -1.13488, -1.14008},       // 6−7 GeV
        //     {-1.16058, -1.16113, -1.1519, -1.1639, -1.15457, -1.1564},          // 7−8 GeV
        //     {-1.15733, -1.16306, -1.14965, -1.15888, -1.15131, -1.15527},       // 8−9 GeV
        //     {-1.08707, -1.16674, -1.21746, -1.18771, -1.17542, -1.19943}        // >9 GeV
        // };
}
/////////////////////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_hit_position_fiducial_cut_homogeneous()
{
        // Cut using the natural directions of the scintillator bars/ fibers:
        ///////////////////////////////////////////////////////////////////
        /// inbending:
        //
        double min_v_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double min_v_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double min_v_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        //
        double max_v_tight_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_v_med_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_v_loose_inb[6] = {400, 400, 400, 400, 400, 400};
        //
        double min_w_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double min_w_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double min_w_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        //
        double max_w_tight_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_w_med_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_w_loose_inb[6] = {400, 400, 400, 400, 400, 400};
        //////////////////////////////////////////////////////////////
        int isec = (_data->ec_pcal_sec(0) - 1);
        double min_v = min_v_med_inb[isec];
        double max_v = max_v_med_inb[isec];
        double min_w = min_w_med_inb[isec];
        double max_w = max_w_med_inb[isec];
        return (_data->ec_pcal_lv(0) > min_v && _data->ec_pcal_lv(0) < max_v && _data->ec_pcal_lw(0) > min_w &&
                _data->ec_pcal_lw(0) < max_w);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::PCAL_fiducial_cut_X_Y()
{

        double minparams_pcal_in[6][2] = {{-0.53233, 21.92333}, {-0.53876, 22.35121}, {-0.50470, 21.58152}, {-0.52421, 24.02394}, {-0.52700, 25.79000}, {-0.52324, 21.14879}};
        double maxparams_pcal_in[6][2] = {{0.51997, -20.34515}, {0.50942, -22.93788}, {0.52876, -25.21121}, {0.53930, -22.61848}, {0.52876, -25.21121}, {0.54148, -22.88758}};

        // double min_radious[6] = {71.245, 71.587, 72.142, 73.101, 72.025, 72.921};
        double min_radious = 75;
        double max_radious = 260;

        short pcal_sector = (_data->ec_pcal_sec(0) - 1);

        double X = _data->ec_pcal_x(0);
        double Y = _data->ec_pcal_y(0);

        float X_new = X * cos(DEG2RAD * (-60 * (pcal_sector))) - Y * sin(DEG2RAD * (-60 * (pcal_sector)));
        Y = X * sin(DEG2RAD * (-60 * (pcal_sector))) + Y * cos(DEG2RAD * (-60 * (pcal_sector)));

        X = X_new;

        double radious = sqrt(X * X + Y * Y);
        double Min_radious = min_radious; //[pcal_sector];
        double Max_radious = max_radious;

        double calc_min = minparams_pcal_in[pcal_sector][0] * X + minparams_pcal_in[pcal_sector][1];
        double calc_max = maxparams_pcal_in[pcal_sector][0] * X + maxparams_pcal_in[pcal_sector][1];

        return ((Y > calc_min) && (Y < calc_max) && (radious > Min_radious) && (radious < Max_radious));
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

bool Pass2_Cuts::DC_Ineff_cut_X_Y(int i, int pid)
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

        if (pid == 1)
        {
                _dc_ineff_cut &= X3 < 330;
                if (dc_sector == 1 || dc_sector == 4)
                {
                        _dc_ineff_cut &= ((Y3 > (1.36 * X3 - 515.8)) && (Y3 < (-1.36 * X3 + 515.8)));
                }
                else if (dc_sector == 2 || dc_sector == 5)
                {
                        _dc_ineff_cut &= ((Y3 > (1.11 * X3 - 425.7)) && (Y3 < (-1.11 * X3 + 425.7)));
                }
        }
        if (pid == 2)
        {
                _dc_ineff_cut &= X3 < 330;

                if (dc_sector == 1 || dc_sector == 4)
                {
                        _dc_ineff_cut &= ((Y3 > (1.36 * X3 - 505.8)) && (Y3 < (-1.36 * X3 + 505.8)));
                }
                else if (dc_sector == 2 || dc_sector == 5)
                {
                        _dc_ineff_cut &= ((Y3 > (1.11 * X3 - 425.7)) && (Y3 < (-1.11 * X3 + 425.7)));
                }
        }
        return _dc_ineff_cut;
}

bool Pass2_Cuts::DC_fiducial_cut_XY(int i, int pid)
{
        bool _dc_fid_cut = true;
        // bool isinbending = true;
        // new cut parameters for the linear cut based on x and y coordinates (inbending field):
        // replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
        // (optimized for electrons, do not use it for hadrons)

        // maxparams_in[1][6][3][2] -> [pid][sec][regions][a,b]->a*x+b

        /// supergaus
        double minparams_in[3][6][3][2] = {
            //     {{{-0.58933, 15.03000}, {-0.60189, 22.39868}, {-0.54181, 21.96126}},
            //      {{-0.58242, 14.04545}, {-0.59240, 20.85846}, {-0.55192, 22.29423}},
            //      {{-0.55927, 14.26636}, {-0.57437, 21.95077}, {-0.51302, 21.37775}},
            //      {{-0.56533, 14.35000}, {-0.57253, 21.57033}, {-0.52352, 22.35879}},
            //      {{-0.58582, 16.10909}, {-0.58703, 23.30440}, {-0.52440, 23.26868}},
            //      {{-0.57624, 13.84455}, {-0.59262, 21.42824}, {-0.53071, 20.85082}}},

            {{{-0.62170, 15.22182}, {-0.62585, 22.49044}, {-0.56401, 22.36676}},
             {{-0.61648, 14.61909}, {-0.61152, 20.72934}, {-0.56335, 21.28434}},
             {{-0.58788, 14.42273}, {-0.60119, 22.33538}, {-0.53533, 21.37390}},
             {{-0.59939, 14.92364}, {-0.60844, 22.53813}, {-0.54571, 22.24121}},
             {{-0.61624, 16.14455}, {-0.61723, 24.04363}, {-0.55016, 23.68214}},
             {{-0.61624, 14.84455}, {-0.61829, 21.46714}, {-0.54973, 20.49643}}},

            {{{-0.51812, 4.80882}, {-0.60660, 19.32861}, {-0.60392, 34.35024}},
             {{-0.52912, 5.39632}, {-0.61995, 20.57325}, {-0.62008, 37.14595}},
             {{-0.50012, 4.10882}, {-0.59236, 18.29823}, {-0.58771, 32.53500}},
             {{-0.52647, 5.16029}, {-0.60426, 18.88900}, {-0.60434, 34.83024}},
             {{-0.50600, 4.00000}, {-0.58717, 17.29909}, {-0.57177, 28.34119}},
             {{-0.52388, 4.91618}, {-0.60842, 19.32450}, {-0.61696, 36.42214}}},

            {{{-0.48371, 2.71544}, {-0.56539, 14.80074}, {-0.56644, 28.90231}},
             {{-0.49800, 3.45000}, {-0.56931, 14.84191}, {-0.56402, 27.89462}},
             {{-0.49165, 3.19853}, {-0.55265, 13.67279}, {-0.56789, 30.39264}},
             {{-0.48529, 2.72206}, {-0.56480, 14.85221}, {-0.57158, 30.11088}},
             {{-0.46906, 1.67941}, {-0.55137, 13.52941}, {-0.54323, 24.85703}},
             {{-0.51318, 4.11324}, {-0.56931, 14.84191}, {-0.56697, 28.84165}}}

        };

        double maxparams_in[3][6][3][2] = {
            //     {{{0.57079, -13.46727}, {0.57108, -19.36121}, {0.52978, -20.38791}},
            //      {{0.56255, -14.79273}, {0.58448, -23.36066}, {0.52522, -23.29478}},
            //      {{0.62473, -17.85364}, {0.59042, -23.65901}, {0.53396, -24.24066}},
            //      {{0.58242, -14.04545}, {0.59934, -21.59780}, {0.55005, -22.42995}},
            //      {{0.57624, -15.64455}, {0.59525, -23.87989}, {0.53313, -24.02225}},
            //      {{0.59309, -15.06545}, {0.60044, -22.17527}, {0.54753, -22.55247}}},

            {{{0.60679, -14.11727}, {0.61275, -20.82582}, {0.56231, -21.30000}},
             {{0.60255, -15.79273}, {0.61130, -23.74527}, {0.54742, -23.90027}},
             {{0.66073, -18.50364}, {0.61798, -24.23231}, {0.56500, -25.25962}},
             {{0.62242, -15.04545}, {0.62954, -22.36560}, {0.56780, -21.79451}},
             {{0.64133, -17.42000}, {0.62492, -24.55736}, {0.56440, -24.99945}},
             {{0.61624, -14.84455}, {0.62725, -22.55989}, {0.56951, -22.53049}}},

            {{{0.53865, -5.86103}, {0.61408, -19.00965}, {0.64226, -38.73976}},
             {{0.52876, -5.36985}, {0.61309, -19.36563}, {0.61519, -34.45357}},
             {{0.52388, -4.91618}, {0.60462, -18.49771}, {0.61353, -34.83833}},
             {{0.53606, -5.61691}, {0.61106, -18.60558}, {0.62413, -35.58310}},
             {{0.51347, -4.49779}, {0.60016, -18.02074}, {0.59738, -31.99500}},
             {{0.52918, -5.16324}, {0.60592, -18.18320}, {0.61468, -33.58214}}},

            {{{0.49447, -3.13529}, {0.55402, -12.47426}, {0.57062, -26.22209}},
             {{0.49800, -3.45000}, {0.56431, -13.96176}, {0.57756, -28.62341}},
             {{0.49235, -2.92647}, {0.55402, -12.61544}, {0.56486, -26.74143}},
             {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56947, -26.46780}},
             {{0.48971, -2.89044}, {0.55402, -12.61544}, {0.55198, -24.45110}},
             {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56565, -25.70758}}}};

        // // simu
        // double maxparams_in_sim[3][6][3][2] = {
        //     {{0.56012, -12.48727}, {0.57116, -18.59769}, {0.53214, -18.98324}},
        //     {{0.55394, -12.28636}, {0.56703, -18.46868}, {0.51978, -18.26868}},
        //     {{0.53236, -10.51182}, {0.55455, -16.73088}, {0.51912, -16.84011}},
        //     {{0.55079, -11.41727}, {0.56356, -17.47758}, {0.52621, -17.83379}},
        //     {{0.53236, -10.01182}, {0.55960, -16.57868}, {0.51181, -15.30357}},
        //     {{0.55467, -12.07000}, {0.57042, -18.40901}, {0.53214, -18.98324}}};

        // double minparams_in_sim[3][6][3][2] = {
        //     {{-0.55079, 11.41727}, {-0.57226, 17.96088}, {-0.52852, 17.96071}},
        //     {{-0.55079, 11.41727}, {-0.56356, 17.47758}, {-0.52621, 17.83379}},
        //     {{-0.53394, 10.23636}, {-0.55960, 16.57868}, {-0.51462, 15.83846}},
        //     {{-0.55988, 12.71273}, {-0.57815, 19.81956}, {-0.51929, 18.59148}},
        //     {{-0.53830, 10.43818}, {-0.55178, 16.16022}, {-0.51462, 15.83846}},
        //     {{-0.55079, 11.41727}, {-0.57424, 18.24604}, {-0.52852, 17.96071}}};

        double min_r[3][3][6] = {
            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {54.504, 54.487, 54.405, 54.644, 54.858, 54.992},
             {63.980, 64.513, 63.854, 63.939, 63.953, 64.713}},

            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {56, 56, 56, 56, 56, 56},
             {135, 135, 135, 135, 135, 135}},

            {{36.754, 36.386, 36.528, 36.946, 36.680, 36.903},
             {56, 56, 56, 56, 56, 56},
             {135, 135, 135, 135, 135, 135}},

        };

        double max_r[3][3][6] = {
            {{101, 101, 101, 101, 101, 101},
             {155, 155, 155, 155, 155, 155},
             {210, 210, 210, 210, 210, 210}},

            {{145, 145, 145, 145, 145, 145},
             {220, 220, 220, 220, 220, 220},
             {350, 350, 350, 350, 350, 350}},

            {{140, 145, 145, 145, 145, 140},
             {210, 220, 220, 220, 220, 210},
             {325, 340, 340, 340, 340, 325}}};

        short dc_sector = (_data->dc_sec(i) - 1);

        // region 1
        double X1 = _data->dc_r1_x(i);
        double Y1 = _data->dc_r1_y(i);
        float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

        X1 = X1_new;
        int region_1 = 1;

        double calc_min1 = minparams_in[pid][dc_sector][region_1 - 1][0] * X1 + minparams_in[pid][dc_sector][region_1 - 1][1];
        double calc_max1 = maxparams_in[pid][dc_sector][region_1 - 1][0] * X1 + maxparams_in[pid][dc_sector][region_1 - 1][1];
        double DC_r1 = sqrt(X1 * X1 + Y1 * Y1);
        double Min_r1 = min_r[pid][region_1 - 1][dc_sector];
        double Max_r1 = max_r[pid][region_1 - 1][dc_sector];

        // region 2
        double X2 = _data->dc_r2_x(i);
        double Y2 = _data->dc_r2_y(i);
        float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));
        X2 = X2_new;
        int region_2 = 2;
        double calc_min2 = minparams_in[pid][dc_sector][region_2 - 1][0] * X2 + minparams_in[pid][dc_sector][region_2 - 1][1];
        double calc_max2 = maxparams_in[pid][dc_sector][region_2 - 1][0] * X2 + maxparams_in[pid][dc_sector][region_2 - 1][1];
        double DC_r2 = sqrt(X2 * X2 + Y2 * Y2);
        double Min_r2 = min_r[pid][region_2 - 1][dc_sector];
        double Max_r2 = max_r[pid][region_2 - 1][dc_sector];

        // region 3
        double X3 = _data->dc_r3_x(i);
        double Y3 = _data->dc_r3_y(i);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));
        X3 = X3_new;
        int region_3 = 3;
        double calc_min3 = minparams_in[pid][dc_sector][region_3 - 1][0] * X3 + minparams_in[pid][dc_sector][region_3 - 1][1];
        double calc_max3 = maxparams_in[pid][dc_sector][region_3 - 1][0] * X3 + maxparams_in[pid][dc_sector][region_3 - 1][1];
        double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        double Min_r3 = min_r[pid][region_3 - 1][dc_sector];
        double Max_r3 = max_r[pid][region_3 - 1][dc_sector];
        // double Max_line = max_line[pid][region_3 - 1][dc_sector];

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

bool Pass2_Cuts::DC_z_vertex_cut()
{
        int pcal_sector = _data->ec_pcal_sec(0);
        float partvz = _data->vz(0);

        double vz_min_sect_inb[] = {-10, -10, -10, -10, -10, -10};
        double vz_max_sect_inb[] = {5, 5, 5, 5, 5, 5};

        float vz_min_sect[6];
        float vz_max_sect[6];

        for (int i = 0; i < 6; i++)
        {
                vz_min_sect[i] = vz_min_sect_inb[i];
                vz_max_sect[i] = vz_max_sect_inb[i];
        }

        int isec = pcal_sector - 1;
        float vz_min = vz_min_sect[isec];
        float vz_max = vz_max_sect[isec];

        return partvz > vz_min && partvz < vz_max;
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
bool Pass2_Cuts::Hadron_Delta_vz_cut(int i)
{
        int pid = _data->pid(i);
        // if(pid==PROTON){
        float dvz = (_data->vz(i) - _data->vz(0));

        // std::cout<<"dvz  "<<dvz<<std::endl;

        return dvz > -20 && dvz < 20;
        // switch (pid)
        // {
        // case 2212:
        //         return dvz > -20 && dvz < 20;
        // case 22:
        //         return dvz > -20 && dvz < 20;
        // case 2112:
        //         return dvz > -20 && dvz < 20;
        // case 211:
        //         return dvz > -20 && dvz < 20;
        // case -211:
        //         return dvz > -20 && dvz < 20;
        // case 321:
        //         return dvz > -20 && dvz < 20;
        // case -321:
        //         return dvz > -20 && dvz < 20;
        // }
        // return false;
}

/** chi2pid cut for hadrons
 * @param chi2pid chi2pid value
 * @param pid hadron PID code
 */
bool Pass2_Cuts::Hadron_Chi2pid_cut(int i)
{
        bool isstrict = false;
        float chi2pid = _data->chi2pid(i);
        float p = _data->p(i);
        int pid = _data->pid(i);
        int status = abs(_data->status(i));

        if (status < 4000)
                return abs(chi2pid) < 5.0; /// trying very loose cuts
        else
        {
                return abs(chi2pid) < 7.0;
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
bool Pass2_Cuts::CD_fiducial_had(int i)
{
        bool pass_fiducial = true;
        // int pid = _data->pid(i);
        //        if (pid == 2212)
        {
                double momT = sqrt(_data->px(i) * _data->px(i) + _data->py(i) * _data->py(i));
                double theta = atan2(momT, _data->pz(i)) * 180 / PI;
                double phi = atan2(_data->py(i), _data->px(i)) * 180 / PI;

                double fiducial_phi_width = 3; // 3 is used by andrew
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
