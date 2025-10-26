
#include "histogram.hpp"

Histogram::Histogram(const std::string &output_file)
{
        RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
        def = std::make_shared<TCanvas>("def");
        if (getenv("BEAM_E") != NULL)
        {
                if (atof(getenv("BEAM_E")) < 8)
                {
                        q2_max = 4.0;
                        w_max = 4.0;
                        p_max = 4.0;
                }
                else if (atof(getenv("BEAM_E")) < 9)
                {
                        q2_max = 7.0;
                        w_max = 7.0;
                        p_max = 7.0;
                }
        }
        ////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////
        ////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////
        ////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////

        // for (short q2 = 3; q2 < 4; q2++)
        for (short q2 = 1; q2 < q2_bin_size; q2++)
        {
                const Int_t ndims_5D = 5;

                // Int_t bins_5D[ndims_5D] = {22, 22, 10, 6, 10};
                // in our expected range
                // Int_t bins_5D_original[ndims_5D] = {7, 7, 10, 6, 10};
                Int_t bins_5D[ndims_5D] = {15, 15, 10, 6, 10};

                // Mlower = mh1 + mh2
                //  Double_t xmin_5D_original[ndims_5D] = {(0.938272 + 0.13957), (0.13957 + 0.13957), 0., 0., 0.};

                float q2_lower_lim = q2_low_values[q2];
                float q2_upper_lim = q2_up_values[q2];

                // for (short w = w_lower_bin; w < w_higher_bin; w++)
                for (short w = w_lower_bin; w < w_higher_bin; w++)
                {
                        // Mupper(W) = W − mh3
                        // //50 MeV w bin

                        // // // // //adding extra bins in each end of invariant mass hist
                        // Double_t Bin_size_pPip0 = ((1.0 + 0.05 * w + 0.025 - MASS_PIM) - (0.938272 + 0.13957)) / 14.0;
                        // Double_t Bin_size_pipPim0 = ((1.0 + 0.05 * w + 0.025 - MASS_P) - (0.13957 + 0.13957)) / 14.0;

                        // // //adding extra bins in each end of invariant mass hist
                        Double_t Bin_size_pPip0 = ((1.0 + 0.05 * w + 0.025 - MASS_PIM) - (0.938272 + 0.13957)) / 7.0;
                        Double_t Bin_size_pipPim0 = ((1.0 + 0.05 * w + 0.025 - MASS_P) - (0.13957 + 0.13957)) / 7.0;

                        // Double_t xmin_5D[ndims_5D] = {(0.938272 + 0.13957), (0.13957 + 0.13957), 0., 0., 0.};
                        // Double_t xmax_5D[ndims_5D] = {(1.0 + 0.05 * w + 0.025 - MASS_PIM), (1.0 + 0.05 * w + 0.025 - MASS_P), 180, 360, 360};

                        Double_t xmin_5D[ndims_5D] = {((0.938272 + 0.13957) - 4 * Bin_size_pPip0), (0.13957 + 0.13957) - 4 * Bin_size_pipPim0, 0., 0.0, 0.};
                        Double_t xmax_5D[ndims_5D] = {((1.0 + 0.05 * w + 0.025 - MASS_PIM) + 4 * Bin_size_pPip0), ((1.0 + 0.05 * w + 0.025 - MASS_P) + 4 * Bin_size_pipPim0), 180, 360, 360};

                        // Double_t xmin_5D[ndims_5D] = {(0.938272 + 0.13957) - 3 * Bin_size_pPip, (0.13957 + 0.13957) - 3 * Bin_size_pipPim, 0., 0., 0.};
                        // Double_t xmax_5D[ndims_5D] = {(1.0 + 0.05 * w - MASS_PIM) + 4 * Bin_size_pPip, (1.0 + 0.05 * w - MASS_P) + 4 * Bin_size_pipPim, 180, 360, 360};

                        // 25 MeV w bin
                        //  Double_t xmax_5D[ndims_5D] = {(1.0 + 0.025 * w - MASS_PIM), (1.0 + 0.025 * w - MASS_P), 180, 360, 360};

                        // Double_t xmax_5D[ndims_5D] = {(1.0 + 0.05 * w + 0.05), (1.0 + 0.05 * w + 0.05) / 2.0 + 0.1, 180, 360, 360};

                        // auto name = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (1.0 + 1.0 * q2), (1.0 + 1.0 * q2 + 1.0), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05));
                        // auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (1.0 + 1.0 * q2), (1.0 + 1.0 * q2 + 1.0), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05));

                        // //50 MeV w bin

                        auto name = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05));
                        auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05));

                        // std::cout << " name is :" << name << std::endl;

                        // 25 MeV w bin

                        // auto name = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.0 + 0.025 * w), (1.0 + 0.025 * w + 0.025));
                        // auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.0 + 0.025 * w), (1.0 + 0.025 * w + 0.025));

                        sevenDHist_prot[q2][w] = new THnSparseD(name, name,
                                                                ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_prot[q2][w]->Sumw2();

                        sevenDHist_pip[q2][w] = new THnSparseD(name, name,
                                                               ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenDHist_pip[q2][w]->Sumw2();

                        sevenDHist_pim[q2][w] = new THnSparseD(name, name,
                                                               ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenDHist_pim[q2][w]->Sumw2();

                        h_5dim_prot_evt[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pip_evt[q2][w] = new THnSparseD(name_evt, name_evt,
                                                               ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pim_evt[q2][w] = new THnSparseD(name_evt, name_evt,
                                                               ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        /////////////////////  tight cuts /////////////////////////

                        sevenDHist_prot_tight[q2][w] = new THnSparseD(name, name,
                                                                      ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_prot_tight[q2][w]->Sumw2();

                        sevenDHist_pip_tight[q2][w] = new THnSparseD(name, name,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_pip_tight[q2][w]->Sumw2();

                        sevenDHist_pim_tight[q2][w] = new THnSparseD(name, name,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_pim_tight[q2][w]->Sumw2();

                        h_5dim_prot_evt_tight[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                      ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pip_evt_tight[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pim_evt_tight[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        /////////////////////  loose cuts /////////////////////////

                        sevenDHist_prot_loose[q2][w] = new THnSparseD(name, name,
                                                                      ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_prot_loose[q2][w]->Sumw2();

                        sevenDHist_pip_loose[q2][w] = new THnSparseD(name, name,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_pip_loose[q2][w]->Sumw2();

                        sevenDHist_pim_loose[q2][w] = new THnSparseD(name, name,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_pim_loose[q2][w]->Sumw2();

                        h_5dim_prot_evt_loose[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                      ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pip_evt_loose[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        h_5dim_pim_evt_loose[q2][w] = new THnSparseD(name_evt, name_evt,
                                                                     ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        for (size_t xi = 0; xi < 14; xi++)
                        {
                                float xmin_pPip = xmin_5D[0] + Bin_size_pPip0 * xi;
                                float xmax_pPip = xmin_5D[0] + Bin_size_pPip0 * (xi + 1);

                                float xmin_pipPim = xmin_5D[1] + Bin_size_pipPim0 * xi;
                                float xmax_pipPim = xmin_5D[1] + Bin_size_pipPim0 * (xi + 1);

                                auto name_mm2_inv_pPip = Form("h_mm2_mPim_inv_pPip_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV_%.3f<=M_pPip<=%.3f GeV",
                                                              (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05), xmin_pPip, xmax_pPip);

                                auto name_mm2_inv_pPim = Form("h_mm2_mPim_inv_pPim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV_%.3f<=M_pPim<=%.3f GeV",
                                                              (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05), xmin_pPip, xmax_pPip);

                                auto name_mm2_inv_pipPim = Form("h_mm2_mPim_inv_pipPim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV_%.3f<=M_pipPi<=%.3f GeV",
                                                                (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05), xmin_pipPim, xmax_pipPim);
                                mm2_mPim_hist_inv_pPip[q2][w][xi] = std::make_shared<TH1D>(name_mm2_inv_pPip, name_mm2_inv_pPip, 100, -0.4, 0.4);
                                mm2_mPim_hist_inv_pPim[q2][w][xi] = std::make_shared<TH1D>(name_mm2_inv_pPim, name_mm2_inv_pPim, 100, -0.4, 0.4);
                                mm2_mPim_hist_inv_pipPim[q2][w][xi] = std::make_shared<TH1D>(name_mm2_inv_pipPim, name_mm2_inv_pipPim, 100, -0.4, 0.4);
                        }
                }
        }

        ////////////////////////////// Done THNSPARSE and related ////////////////////////////////
        ////////////////////////////// Done THNSPARSE and related //////////////////////////////////////
        //////////////////////// DONE THNSPARSE and related ////////////////////////////////
        mc_pid_at_zero = std::make_shared<TH1D>("mc_pid_at_zero", "mc pid at zero", 1000, -2500, 2500);
        pid_at_zero = std::make_shared<TH1D>("pid_at_zero", "pid at zero", 1000, -2500, 2500);
        weight_hist = std::make_shared<TH1D>("weight", "weight", bins, 0.9, 2.2);

        dp_prot_cdfd_hist = std::make_shared<TH1D>("P_fd-P_cd_Prot", "P_fd-P_cd_Prot", 200, -1, 1);
        dp_pip_cdfd_hist = std::make_shared<TH1D>("P_fd-P_cd_Pip", "P_fd-P_cd_Pip", 200, -1, 1);

        dth_prot_cdfd_hist = std::make_shared<TH1D>("th_fd-th_cd_Prot", "th_fd-th_cd_Prot", 200, -100, 25);
        dth_pip_cdfd_hist = std::make_shared<TH1D>("th_fd-th_cd_Pip", "th_fd-Pthcd_Pip", 200, -100, 25);

        dphi_prot_cdfd_hist = std::make_shared<TH1D>("phi_fd-P_cd_Prot", "phi_fd-phi_cd_Prot", 200, -100, 100);
        dphi_pip_cdfd_hist = std::make_shared<TH1D>("phi_fd-P_cd_Pip", "phi_fd-phi_cd_Pip", 200, -100, 100);

        inv_mass_pPip = std::make_shared<TH1D>("pPip_mass", "Prot-Pip mass", bins, 1.0, 2.25);
        inv_mass_pPim = std::make_shared<TH1D>("pPim_mass", "Prot-Pim mass", bins, 0.75, 2.25);
        inv_mass_pipPim = std::make_shared<TH1D>("pipPim_mass", "Pip-Pim mass", bins, 0.0, 1.5);

        theta_Prot_cm = std::make_shared<TH1D>("theta_prot_cm", "theta prot cm", bins, 0.0, 180.0);
        theta_Pip_cm = std::make_shared<TH1D>("theta_pip_cm", "theta pip cm", bins, 0.0, 180.0);
        theta_Pim_cm = std::make_shared<TH1D>("theta_pim_cm", "theta pim cm", bins, 0.0, 180.0);

        phi_Prot_cm = std::make_shared<TH1D>("phi_prot_cm", "phi prot cm", bins, 0.0, 360.0);
        phi_Pip_cm = std::make_shared<TH1D>("phi_pip_cm", "phi pip cm", bins, 0.0, 360.0);
        phi_Pim_cm = std::make_shared<TH1D>("phi_pim_cm", "phi pim cm", bins, 0.0, 360.0);

        alpha_Prot_cm = std::make_shared<TH1D>("alpha_prot_cm", "alpha prot cm", bins, 0.0, 360.0);
        alpha_Pip_cm = std::make_shared<TH1D>("alpha_pip_cm", "alpha pip cm", bins, 0.0, 360.0);
        alpha_Pim_cm = std::make_shared<TH1D>("alpha_pim_cm", "alpha pim cm", bins, 0.0, 360.0);

        inv_mass_pPip_swapped = std::make_shared<TH1D>("pPip_mass_swapped", "Prot-Pip mass swapped", bins, 1.0, 2.25);
        inv_mass_pPim_swapped = std::make_shared<TH1D>("pPim_mass_swapped", "Prot-Pim mass swapped", bins, 0.0, 1.5);
        inv_mass_pipPim_swapped = std::make_shared<TH1D>("pipPim_mass_swapped", "Pip-Pim mass swapped", bins, 0.75, 2.25);

        theta_Prot_cm_swapped = std::make_shared<TH1D>("theta_prot_cm_swapped", "theta prot cm swapped", bins, 0.0, 180.0);
        theta_Pip_cm_swapped = std::make_shared<TH1D>("theta_pip_cm_swapped", "theta pip cm swapped", bins, 0.0, 180.0);
        theta_Pim_cm_swapped = std::make_shared<TH1D>("theta_pim_cm_swapped", "theta pim cm swapped", bins, 0.0, 180.0);

        phi_Prot_cm_swapped = std::make_shared<TH1D>("phi_prot_cm_swapped", "phi prot cm swapped", bins, 0.0, 360.0);
        phi_Pip_cm_swapped = std::make_shared<TH1D>("phi_pip_cm_swapped", "phi pip cm swapped", bins, 0.0, 360.0);
        phi_Pim_cm_swapped = std::make_shared<TH1D>("phi_pim_cm_swapped", "phi pim cm swapped", bins, 0.0, 360.0);

        alpha_Prot_cm_swapped = std::make_shared<TH1D>("alpha_prot_cm_swapped", "alpha prot cm swapped", bins, 0.0, 360.0);
        alpha_Pip_cm_swapped = std::make_shared<TH1D>("alpha_pip_cm_swapped", "alpha pip cm swapped", bins, 0.0, 360.0);
        alpha_Pim_cm_swapped = std::make_shared<TH1D>("alpha_pim_cm_swapped", "alpha pim cm swapped", bins, 0.0, 360.0);

        dp_prot_hist = std::make_shared<TH1D>("P_gen-P_rec_Prot", "Prot (Gen - Rec) Mom mPim events", 200, -0.005, 0.01);
        dp_pip_hist = std::make_shared<TH1D>("P_gen-P_rec_Pip", "Pip (Gen - Rec) Mom mPim events", 200, -0.005, 0.01);

        dp_prot_for_pip_hist = std::make_shared<TH1D>("P_gen-P_rec_Prot_for_pip", "Prot (Gen - Rec) Mom (for pip)", 200, -0.005, .01);
        dp_pip_for_prot_hist = std::make_shared<TH1D>("P_gen-P_rec_Pip_for_Proton", "Pip (Gen - Rec) Mom (for Prot)", 200, -0.005, .01);

        dp_ambi_prot_all_hist = std::make_shared<TH1D>("P_gen-P_rec_ambi_Prot", " without event selection Prot (Gen - Rec) Mom", 200, -0.005, .01);
        dp_ambi_pip_all_hist = std::make_shared<TH1D>("P_gen-P_rec_ambi_Pip", " without event selection Pip (Gen - Rec) Mom", 200, -0.005, .01);

        dp_sum_hist = std::make_shared<TH1D>("P_gen-P_rec_sum", "sum (Gen - Rec) Mom", 100, -0.005, 0.01);
        dp_sum_hist_twoPi = std::make_shared<TH1D>("P_gen-P_rec_sum_twoPi", "sum (Gen - Rec) Mom Selected Events", 100, -0.005, 0.01);

        entries_in_each_event = std::make_shared<TH1D>("Entries_per_event", "No of entries per event", 20, -1, 10);
        entries_prot = std::make_shared<TH1D>("protons_per_event", "No of prot per event", 20, -1, 10);
        entries_pip = std::make_shared<TH1D>("pips_per_event", "No of pip per event", 20, -1, 10);
        MM2_mPim_all_comb = std::make_shared<TH1D>("MMSQ_all_Combination", "MMSQ all Combination", bins, -0.4, 0.4);
        MM2_mPim_1_comb = std::make_shared<TH1D>("MMSQ_1_Combination", "MMSQ 1 Combination", bins, -0.4, 0.4);
        MM2_mPim_2_comb = std::make_shared<TH1D>("MMSQ_mPim_Swapped_P_pip", "MMSQ mPim Swapped Proton Pip", bins, -0.4, 0.4);
        // MM2_mPim_3_comb = std::make_shared<TH1D>("MMSQ_diff_swapped_unswapped", "Diff MMSQ mPim Unswapped - Swapped", bins, -0.4, 0.4);
        MM2_mPim_3_comb = std::make_shared<TH1D>("dv2_original-swapped_pip", "dv2 original - swapped Pip", bins, -0.4, 2.4);
        MM2_mPim_4_or_more_comb = std::make_shared<TH1D>("dv2_original-swapped_prot", "dv2 original - swapped Proton", bins, -0.4, 2.4);

        // p_gen_prot_hist = std::make_shared<TH1D>("P_gen_Prot", "Prot (Gen) Mom", bins, 0, 5);
        // p_gen_pip_hist = std::make_shared<TH1D>("P_gen_Pip", "Pip (Gen) Mom", bins, 0, 5);
        // p_gen_prot_for_pip_hist = std::make_shared<TH1D>("P_gen_Prot_for_pip", "Prot (Gen ) Mom (for pip)", bins, 0, 5);
        // p_gen_pip_for_prot_hist = std::make_shared<TH1D>("P_gen_Pip_for_Proton", "Pip (Gen ) Mom (for Prot)", bins, 0, 5);
        // p_gen_ambi_prot_hist = std::make_shared<TH1D>("P_gen_ambi_Prot", "Ambi Prot (Gen) Mom", bins, 0, 5);
        // p_gen_ambi_pip_hist = std::make_shared<TH1D>("P_gen_ambi_Pip", "Ambi Pip (Gen) Mom", bins, 0, 5);

        // p_rec_prot_hist = std::make_shared<TH1D>("P_rec_Prot", "Prot (rec) Mom", bins, 0, 5);
        // p_rec_pip_hist = std::make_shared<TH1D>("P_rec_Pip", "Pip (rec) Mom", bins, 0, 5);
        // p_rec_prot_for_pip_hist = std::make_shared<TH1D>("P_rec_Prot_for_pip", "Prot (rec ) Mom (for pip)", bins, 0, 5);
        // p_rec_pip_for_prot_hist = std::make_shared<TH1D>("P_rec_Pip_for_Proton", "Pip (rec ) Mom (for Prot)", bins, 0, 5);
        // p_rec_ambi_prot_hist = std::make_shared<TH1D>("P_rec_ambi_Prot", "Ambi Prot (rec) Mom", bins, 0, 5);
        // p_rec_ambi_pip_hist = std::make_shared<TH1D>("P_rec_ambi_Pip", "Ambi Pip (rec) Mom", bins, 0, 5);

        W_hist = std::make_shared<TH1D>("W", "W", bins, w_min, w_max);

        // Theta_vs_mom_x_mu = std::make_shared<TH2D>("theta_vs_mom_x_mu_all_sec",
        // "theta_vs_mom_x_mu_all_sec", bins, zero,
        //                                            p_max, bins, zero, 120);

        W_P2pi_hist = std::make_shared<TH1D>("W_P2pi", "W_P2pi", bins, zero, w_max);

        Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, zero, q2_max);
        W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2", bins, w_min, w_max,
                                         bins, q2_min, q2_max);

        W_thrown = std::make_shared<TH1D>("W_thrown", "W_thrown", bins, w_min, w_max);
        Q2_thrown = std::make_shared<TH1D>("Q2_thrown", "Q2_thrown", bins, q2_min, q2_max);

        W_vs_Q2_thrown =
            std::make_shared<TH2D>("W_vs_q2_thrown", "W_vs_q2_thrown", bins, w_min,
                                   w_max, bins, q2_min, q2_max);
        W_vs_MM =
            std::make_shared<TH2D>("W_vs_MM", "W vs MM mPim", bins, w_min,
                                   w_max, bins, -1, 2.5);
        W_vs_MM2 =
            std::make_shared<TH2D>("W_vs_MMSQ", "W vs MMSQ mPim", bins, w_min,
                                   w_max, bins, -1, 2.5);

        Phi_gamma = std::make_shared<TH1D>("Phi_gamma", "Phi_gamma", bins, 0, 360);
        Phi_prot = std::make_shared<TH1D>("Phi_prot", "Phi_prot", bins, 0, 360);
        Phi_pip = std::make_shared<TH1D>("Phi_pip", "Phi_pip", bins, 0, 360);
        Phi_pim = std::make_shared<TH1D>("Phi_pim", "Phi_pim", bins, 0, 360);

        alpha_prot =
            std::make_shared<TH1D>("alpha_prot", "alpha_prot", bins, 0.5, 360);
        alpha_pip = std::make_shared<TH1D>("alpha_pip", "alpha_pip", bins, 0.5, 360);
        alpha_pim = std::make_shared<TH1D>("alpha_pim", "alpha_pim", bins, 0.5, 360);

        //
        // theta_prot_mc = std::make_shared<TH1D>("theta_prot_mc", "theta_prot_mc",
        // bins, 0, 180); theta_pip_mc = std::make_shared<TH1D>("theta_pip_mc",
        // "theta_pip_mc", bins, 0, 180); theta_pim_mc =
        // std::make_shared<TH1D>("theta_pim_mc", "theta_pim_mc", bins, 0, 180);
        //
        // Phi_gamma_mc = std::make_shared<TH1D>("Phi_gamma_mc", "Phi_gamma_mc", bins,
        // 0, 360); Phi_prot_mc = std::make_shared<TH1D>("Phi_prot_mc", "Phi_prot_mc",
        // bins, 0, 360); Phi_pip_mc = std::make_shared<TH1D>("Phi_pip_mc",
        // "Phi_pip_mc", bins, 0, 360); Phi_pim_mc =
        // std::make_shared<TH1D>("Phi_pim_mc", "Phi_pim_mc", bins, 0, 360);
        //
        // alpha_prot_mc = std::make_shared<TH1D>("alpha_prot_mc", "alpha_prot_mc",
        // bins, 0, 360); alpha_pip_mc = std::make_shared<TH1D>("alpha_pip_mc",
        // "alpha_pip_mc", bins, 0, 360); alpha_pim_mc =
        // std::make_shared<TH1D>("alpha_pim_mc", "alpha_pim_mc", bins, 0, 360);

        theta_prot_thrown = std::make_shared<TH1D>("theta_prot_thrown",
                                                   "theta_prot_thrown", bins, 0, 180);
        theta_pip_thrown = std::make_shared<TH1D>("theta_pip_thrown",
                                                  "theta_pip_thrown", bins, 0, 180);
        theta_pim_thrown = std::make_shared<TH1D>("theta_pim_thrown",
                                                  "theta_pim_thrown", bins, 0, 180);

        Phi_gamma_thrown = std::make_shared<TH1D>("Phi_gamma_thrown",
                                                  "Phi_gamma_thrown", bins, 0, 360);
        Phi_prot_thrown = std::make_shared<TH1D>("Phi_prot_thrown", "Phi_prot_thrown",
                                                 bins, 0, 360);
        Phi_pip_thrown =
            std::make_shared<TH1D>("Phi_pip_thrown", "Phi_pip_thrown", bins, 0, 360);
        Phi_pim_thrown =
            std::make_shared<TH1D>("Phi_pim_thrown", "Phi_pim_thrown", bins, 0, 360);

        alpha_prot_thrown = std::make_shared<TH1D>("alpha_prot_thrown",
                                                   "alpha_prot_thrown", bins, 0, 360);
        alpha_pip_thrown = std::make_shared<TH1D>("alpha_pip_thrown",
                                                  "alpha_pip_thrown", bins, 0, 360);
        alpha_pim_thrown = std::make_shared<TH1D>("alpha_pim_thrown",
                                                  "alpha_pim_thrown", bins, 0, 360);

        // MM_neutron = std::make_shared<TH1D>("missMass", "missMass", bins,
        // zero, 4.0);
        MM2_twoPi_excl = std::make_shared<TH1D>("MMSQ_excl", "MMSQ excl", bins,
                                                -0.1, 0.1);
        MM_twoPi_excl = std::make_shared<TH1D>("MM_excl", "MM excl", bins,
                                               -0.25, 0.25);
        MM_twoPi_mPim = std::make_shared<TH1D>("MM_e#pi^{+}pX", "MM: expecting #pi^{-}", bins,
                                               -0.4, 0.4);
        MM2_twoPi_mPim = std::make_shared<TH1D>("MMSQ_e#pi^{+}pX", "MMSQ: expecting #pi^{-}", bins,
                                                -0.4, 0.4);
        MM2_twoPi_missingPip = std::make_shared<TH1D>(
            "e#pi^{-}pX", "MMSQ: expecting #pi^{+}", bins, -0.4, 0.6);

        MM2_twoPi_missingProt = std::make_shared<TH1D>(
            "e#pi^{-}#pi^{+}X", "MMSQ: expecting proton", bins, 0.6, 1.3);
        W_hist_twoPi =
            std::make_shared<TH1D>("W_twoPi", "W_twoPi", bins, 1.0, 3.0);
        Q2_hist_twoPi =
            std::make_shared<TH1D>("Q2_twoPi", "Q2_twoPi", bins, zero, q2_max);
        W_vs_q2_twoPi = std::make_shared<TH2D>("W_vs_q2_twoPi", "W_vs_q2_twoPi", bins,
                                               w_min, w_max, bins, q2_min, q2_max);
        W_vs_q2_twoPi_thrown = std::make_shared<TH2D>("W_vs_q2_twoPi_thrown", "W_vs_q2_twoPi_thrown", bins,
                                                      w_min, w_max, bins, q2_min, q2_max);
        Theta_prot_cm_vs_mom_prot = std::make_shared<TH2D>(
            "Theta_prot_cm_vs_mom_prot", "Theta_prot_cm_vs_mom_prot", bins, zero, 8.5,
            bins, zero, 180);
        Theta_pip_cm_vs_mom_pip = std::make_shared<TH2D>(
            "Theta_pip_cm_vs_mom_pip", "Theta_pip_cm_vs_mom_pip", bins, zero, 8.5,
            bins, zero, 180);
        Theta_pim_cm_vs_mom_pim = std::make_shared<TH2D>(
            "Theta_pim_cm_vs_mom_pim", "Theta_pim_cm_vs_mom_pim", bins, zero, 8.5,
            bins, zero, 180);

        Theta_prot_thrown_cm_vs_mom_prot = std::make_shared<TH2D>(
            "Theta_prot_thrown_cm_vs_mom_prot", "Theta_prot_thrown_cm_vs_mom_prot",
            bins, zero, 8.5, bins, zero, 180);
        Theta_pip_thrown_cm_vs_mom_pip = std::make_shared<TH2D>(
            "Theta_pip_thrown_cm_vs_mom_pip", "Theta_pip_thrown_cm_vs_mom_pip", bins,
            zero, 8.5, bins, zero, 180);
        Theta_pim_thrown_cm_vs_mom_pim = std::make_shared<TH2D>(
            "Theta_pim_thrown_cm_vs_mom_pim", "Theta_pim_thrown_cm_vs_mom_pim", bins,
            zero, 8.5, bins, zero, 180);

        Theta_prot_thrown_lab_vs_mom_prot = std::make_shared<TH2D>(
            "Theta_prot_thrown_lab_vs_mom_prot", "Theta_prot_thrown_cm_vs_mom_prot",
            bins, zero, 8.5, bins, zero, 80);
        Theta_pip_thrown_lab_vs_mom_pip = std::make_shared<TH2D>(
            "Theta_pip_thrown_lab_vs_mom_pip", "Theta_pip_thrown_cm_vs_mom_pip", bins,
            zero, 8.5, bins, zero, 130);
        Theta_pim_thrown_lab_vs_mom_pim = std::make_shared<TH2D>(
            "Theta_pim_thrown_lab_vs_mom_pim", "Theta_pim_thrown_cm_vs_mom_pim", bins,
            zero, 8.5, bins, zero, 130);

        // W_hist_singlePi = std::make_shared<TH1D>("W_singlePi", "W_singlePi", bins,
        // zero, w_max); Q2_hist_singlePi = std::make_shared<TH1D>("Q2_singlePi",
        // "Q2_singlePi", bins, zero, q2_max); W_vs_q2_singlePi =
        //     std::make_shared<TH2D>("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins,
        //     zero, w_max, bins, zero, q2_max);

        makeHists_sector();
        makeHists_deltat();
        makeHists_MomVsBeta();
        makeHists_electron_cuts();
        makeHistMMSQ_mPim();
        // makeHists_x_mu();
        // makeHistTheta_pim_measured();
}

Histogram::~Histogram()
{
        this->Write();
}

void Histogram::Write()
{
        std::cout << GREEN << "Writting" << DEF << std::endl;
        // THnSparse   7D HIST ///////
        // THnSparse   7D HIST ///////
        // THnSparse   7D HIST ///////
        // THnSparse   7D HIST ///////

        std::cerr << BOLDBLUE << " Hists_7D()" << DEF << std::endl;
        TDirectory *THnSparse_7D_prot_folder =
            RootOutputFile->mkdir("THnSparse_7D_prot");
        THnSparse_7D_prot_folder->cd();
        writeHists7D_prot();

        TDirectory *THnSparse_7D_prot_evt_folder =
            RootOutputFile->mkdir("THnSparse_7D_prot_evt");
        THnSparse_7D_prot_evt_folder->cd();
        writeHists7D_prot_evt();

        TDirectory *THnSparse_7D_pim_folder =
            RootOutputFile->mkdir("THnSparse_7D_pim");
        THnSparse_7D_pim_folder->cd();
        writeHists7D_pim();

        TDirectory *THnSparse_7D_pim_evt_folder =
            RootOutputFile->mkdir("THnSparse_7D_pim_evt");
        THnSparse_7D_pim_evt_folder->cd();
        writeHists7D_pim_evt();

        TDirectory *THnSparse_7D_pip_folder =
            RootOutputFile->mkdir("THnSparse_7D_pip");
        THnSparse_7D_pip_folder->cd();
        writeHists7D_pip();

        TDirectory *THnSparse_7D_pip_evt_folder =
            RootOutputFile->mkdir("THnSparse_7D_pip_evt");
        THnSparse_7D_pip_evt_folder->cd();
        writeHists7D_pip_evt();

        // // ///////////////////  tight cuts /////////////////

        // std::cerr << BOLDBLUE << " Hists_7D() tight " << DEF << std::endl;
        // TDirectory *THnSparse_7D_prot_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_prot_tight");
        // THnSparse_7D_prot_folder_tight->cd();
        // writeHists7D_prot_tight();

        // TDirectory *THnSparse_7D_prot_evt_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_prot_evt_tight");
        // THnSparse_7D_prot_evt_folder_tight->cd();
        // writeHists7D_prot_evt_tight();

        // TDirectory *THnSparse_7D_pim_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_pim_tight");
        // THnSparse_7D_pim_folder_tight->cd();
        // writeHists7D_pim_tight();

        // TDirectory *THnSparse_7D_pim_evt_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_pim_evt_tight");
        // THnSparse_7D_pim_evt_folder_tight->cd();
        // writeHists7D_pim_evt_tight();

        // TDirectory *THnSparse_7D_pip_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_pip_tight");
        // THnSparse_7D_pip_folder_tight->cd();
        // writeHists7D_pip_tight();

        // TDirectory *THnSparse_7D_pip_evt_folder_tight =
        //     RootOutputFile->mkdir("THnSparse_7D_pip_evt_tight");
        // THnSparse_7D_pip_evt_folder_tight->cd();
        // writeHists7D_pip_evt_tight();

        // ///////////////////  loose cuts /////////////////

        // std::cerr << BOLDBLUE << " Hists_7D() loose " << DEF << std::endl;
        // TDirectory *THnSparse_7D_prot_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_prot_loose");
        // THnSparse_7D_prot_folder_loose->cd();
        // writeHists7D_prot_loose();

        // TDirectory *THnSparse_7D_prot_evt_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_prot_evt_loose");
        // THnSparse_7D_prot_evt_folder_loose->cd();
        // writeHists7D_prot_evt_loose();

        // TDirectory *THnSparse_7D_pim_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_pim_loose");
        // THnSparse_7D_pim_folder_loose->cd();
        // writeHists7D_pim_loose();

        // TDirectory *THnSparse_7D_pim_evt_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_pim_evt_loose");
        // THnSparse_7D_pim_evt_folder_loose->cd();
        // writeHists7D_pim_evt_loose();

        // TDirectory *THnSparse_7D_pip_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_pip_loose");
        // THnSparse_7D_pip_folder_loose->cd();
        // writeHists7D_pip_loose();

        // TDirectory *THnSparse_7D_pip_evt_folder_loose =
        //     RootOutputFile->mkdir("THnSparse_7D_pip_evt_loose");
        // THnSparse_7D_pip_evt_folder_loose->cd();
        // writeHists7D_pip_evt_loose();

        /////////////////// PID CHECKS //////////////////////
        /////////////////// PID CHECKS //////////////////////
        /////////////////// PID CHECKS //////////////////////
        /////////////////// PID CHECKS //////////////////////

        TDirectory *Hists1D_mm2_mPim_inv_pPip_folder =
            RootOutputFile->mkdir("Hists1D_mm2_mPim_inv_pPip_folder");
        Hists1D_mm2_mPim_inv_pPip_folder->cd();
        writeHists1D_mm2_mPim_inv_pPip();

        TDirectory *Hists1D_mm2_mPim_inv_pPim_folder =
            RootOutputFile->mkdir("Hists1D_mm2_mPim_inv_pPim_folder");
        Hists1D_mm2_mPim_inv_pPim_folder->cd();
        writeHists1D_mm2_mPim_inv_pPim();

        TDirectory *Hists1D_mm2_mPim_inv_pipPim_folder =
            RootOutputFile->mkdir("Hists1D_mm2_mPim_inv_pipPim_folder");
        Hists1D_mm2_mPim_inv_pipPim_folder->cd();
        writeHists1D_mm2_mPim_inv_pipPim();

        std::cerr
            << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
        TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
        WvsQ2_folder->cd();
        Write_WvsQ2();

        // std::cerr << BOLDBLUE << "write_hist_x_mu()" << DEF << std::endl;
        // TDirectory *hists_x_mu = RootOutputFile->mkdir("hists_x_mu");
        // hists_x_mu->cd();
        // write_hist_x_mu();

        // std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
        // TDirectory *Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
        // Write_MomVsBeta_folder->cd();
        // Write_MomVsBeta();

        std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
        TDirectory *Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
        Electron_Cuts->cd();
        Write_Electron_cuts();

        std::cerr << BOLDBLUE << "Write_Hadrons_cuts()" << DEF << std::endl;
        TDirectory *Hadrons_Cuts = RootOutputFile->mkdir("Hadrons_Cuts");
        Hadrons_Cuts->cd();
        Write_Hadrons_cuts();

        std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
        TDirectory *Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
        Write_deltat_folder->cd();
        Write_deltat();

        // std::cerr << BOLDBLUE << "Write_MMSQ_mPim_3D()" << DEF << std::endl;
        // TDirectory *MMSQ_mPim_folder_3D = RootOutputFile->mkdir("MMSQ_mPim_3D");
        // MMSQ_mPim_folder_3D->cd();
        // writeMMSQ_mPim_3D();

        std::cerr << BOLDBLUE << "Write_MMSQ_mPim()" << DEF << std::endl;
        TDirectory *MMSQ_mPim_folder = RootOutputFile->mkdir("MMSQ_mPim");
        MMSQ_mPim_folder->cd();
        writeMMSQ_mPim();

        // std::cerr << BOLDBLUE << "Write_MMSQ_mPim_1_combi()" << DEF << std::endl;
        // TDirectory *MMSQ_mPim_folder_1_combi = RootOutputFile->mkdir("MMSQ_mPim_1_combi");
        // MMSQ_mPim_folder_1_combi->cd();
        // writeMMSQ_mPim_1_comb();

        // std::cerr << BOLDBLUE << "Write_MMSQ_mPim_2_combi()" << DEF << std::endl;
        // TDirectory *MMSQ_mPim_folder_2_combi = RootOutputFile->mkdir("MMSQ_mPim_2_combi");
        // MMSQ_mPim_folder_2_combi->cd();
        // writeMMSQ_mPim_2_comb();

        // std::cerr << BOLDBLUE << "Write_MMSQ_mPim_3_combi()" << DEF << std::endl;
        // TDirectory *MMSQ_mPim_folder_3_combi = RootOutputFile->mkdir("MMSQ_mPim_3_combi");
        // MMSQ_mPim_folder_3_combi->cd();
        // writeMMSQ_mPim_3_comb();

        // std::cerr << BOLDBLUE << "Write_MMSQ_mPim_4_or_more_combi()" << DEF << std::endl;
        // TDirectory *MMSQ_mPim_folder_4_or_more_combi = RootOutputFile->mkdir("MMSQ_mPim_4_or_more_combi");
        // MMSQ_mPim_folder_4_or_more_combi->cd();
        // writeMMSQ_mPim_4_or_more_comb();

        // std::cerr << BOLDBLUE << "Inv_Mass_and_Alpha_cm()" << DEF << std::endl;
        // TDirectory *Inv_Mass_and_Alpha_cm = RootOutputFile->mkdir("Inv_Mass_and_Alpha_cm");
        // Inv_Mass_and_Alpha_cm->cd();
        // write_Inv_Mass_hist();

        // std::cerr << BOLDBLUE << "Write_deltaP()" << DEF << std::endl;
        // TDirectory *Write_deltaP_folder = RootOutputFile->mkdir("DelatP");
        // Write_deltaP_folder->cd();
        // Write_deltaP();
        // // // //
        // // // // std::cerr << BOLDBLUE << "write_hist_theta_pim_measured()" << DEF << std::endl;
        // // // // TDirectory* theta_pim_measured = RootOutputFile->mkdir("theta_pim_measured");
        // // // // theta_pim_measured->cd();
        // // // // write_hist_theta_pim_measured();

        std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

//////////////////////////////////////////////////////////
// void Histogram::makeHistMMSQ_mPim_3D()
// {
//         for (short th = 0; th < th_bin_size; th++)
//         {
//                 float th_lower_lim = th_low_values[th];
//                 float th_upper_lim = th_up_values[th];

//                 for (short mom = 0; mom < mom_bin_size; mom++)
//                 {
//                         float mom_lower_lim = mom_low_values[th][mom];
//                         float mom_upper_lim = mom_up_values[th][mom];

//                         for (short phi = 0; phi < phi_bin_size; phi++)
//                         {
//                                 float phi_lower_lim = phi_low_values[phi];
//                                 float phi_upper_lim = phi_up_values[phi];

//                                 auto name_mmsq = Form(
//                                     "MMSQ_mPim_%.0f<=th<=%.0f_deg_%.1f<=mom<=%.1f_GeV_%.0f<=phi<=%.0f_deg",
//                                     th_lower_lim, th_upper_lim, mom_lower_lim, mom_upper_lim, phi_lower_lim, phi_upper_lim);

//                                 auto name_mmsq_cut = Form(
//                                     "MMSQ_mPim_th_cut_%.0f<=th<=%.0f_deg_%.1f<=mom<=%.1f_GeV_%.0f<=phi<=%.0f_deg",
//                                     th_lower_lim, th_upper_lim, mom_lower_lim, mom_upper_lim, phi_lower_lim, phi_upper_lim);

//                                 MMSQ_mPim_hist_3D[mom][th][phi] = std::make_shared<TH1D>(name_mmsq, name_mmsq, 200, -0.3, 0.3);
//                                 MMSQ_mPim_hist_with_cut_3D[mom][th][phi] = std::make_shared<TH1D>(name_mmsq_cut, name_mmsq_cut, 200, -0.3, 0.3);
//                         }
//                 }
//         }
// }

// void Histogram::Fill_MMSQ_mPim_1_comb(const std::shared_ptr<Reaction> &_e)
// {
//         if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
//         {
//                 MMSQ_mPim_hist_1_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());
//         }
// }

// void Histogram::makeHistMMSQ_mPim_3D()
// {
//         for (short th = 0; th < th_bin_size; th++)
//         {

//                 for (size_t mom = 1; mom < mom_bin_size; mom++)
//                 {

//                         for (short phi = 0; phi < phi_bin_size; phi++)
//                         {
//                                 MMSQ_mPim_hist[mom][th][phi]->SetXTitle("MMSQ(GeV2/c2)");
//                                 if (MMSQ_mPim_hist[mom][th][phi]->GetEntries())
//                                         MMSQ_mPim_hist[mom][th][phi]->Write();
//                         }
//                 }
//         }
// }

/// MMSQ Mpim //////////////////////////////
void Histogram::makeHistMMSQ_mPim()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                float q2_lower_lim = q2_low_values[q2];
                float q2_upper_lim = q2_up_values[q2];

                for (short w = 0; w < 15; w++)
                {
                        // //50 MeV w bin

                        auto name_mmsq = Form("MMSQ_mPim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_mmsq_cut = Form("MMSQ_mPim_w_cut_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_mmsq1 = Form("MMSQ_mPim1_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_mmsq2 = Form("MMSQ_mPim2_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_mmsq3 = Form("MMSQ_mPim3_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_mmsq4 = Form("MMSQ_mPim4_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_inv_pPip = Form("inv_pPip_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_inv_pPim = Form("inv_pPim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_inv_pipPim = Form("inv_pipPim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_alpha_Prot = Form("alpha_Prot_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_alpha_Pip = Form("alpha_Pip_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));
                        auto name_alpha_Pim = Form("alpha_Pim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.4 + 0.05 * w), (1.4 + 0.05 * w + 0.05));

                        MMSQ_mPim_hist[q2][w] = std::make_shared<TH1D>(name_mmsq, name_mmsq, 200, -0.3, 0.3);
                        // std::cout << "   bins  =  " << bins << std::endl;
                        // MMSQ_mPim_hist_with_cut[q2][w] = std::make_shared<TH1D>(name_mmsq_cut, name_mmsq_cut, bins, -0.2, 0.2);

                        MMSQ_mPim_hist_1_comb[q2][w] = std::make_shared<TH1D>(name_mmsq1, name_mmsq1, 200, -0.3, 0.3);

                        MMSQ_mPim_hist_2_comb[q2][w] = std::make_shared<TH1D>(name_mmsq2, name_mmsq2, 200, -0.3, 0.3);

                        MMSQ_mPim_hist_3_comb[q2][w] = std::make_shared<TH1D>(name_mmsq3, name_mmsq3, 200, -0.3, 0.3);

                        MMSQ_mPim_hist_4_or_more_comb[q2][w] = std::make_shared<TH1D>(name_mmsq4, name_mmsq4, 200, -0.3, 1.3);

                        Inv_mass_pPip[q2][w] = std::make_shared<TH1D>(name_inv_pPip, name_inv_pPip, 200, 1.0, 2.25);
                        Inv_mass_pPim[q2][w] = std::make_shared<TH1D>(name_inv_pPim, name_inv_pPim, 200, 0.5, 2.25);
                        Inv_mass_pipPim[q2][w] = std::make_shared<TH1D>(name_inv_pipPim, name_inv_pipPim, 200, 0.0, 1.5);
                        Alpha_Prot_cm[q2][w] = std::make_shared<TH1D>(name_alpha_Prot, name_alpha_Prot, 200, 0.0, 360);
                        Alpha_Pip_cm[q2][w] = std::make_shared<TH1D>(name_alpha_Pip, name_alpha_Pip, 200, 0.0, 360);
                        Alpha_Pim_cm[q2][w] = std::make_shared<TH1D>(name_alpha_Pim, name_alpha_Pim, 200, 0.0, 360);
                }
        }
}
void Histogram::Fill_MMSQ_mPim(const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                // MMSQ_mPim_hist[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());
                MMSQ_mPim_hist[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim());

                // if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                //         MMSQ_mPim_hist_with_cut[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());

                // Inv_mass_pPip[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->inv_Ppip(), _e->weight());
                // Inv_mass_pPim[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->inv_Ppim(), _e->weight());
                // Inv_mass_pipPim[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->inv_pip_pim(), _e->weight());
                // Alpha_Prot_cm[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->alpha_pippim_pipf(), _e->weight());
                // Alpha_Pip_cm[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->alpha_ppim_pipip(), _e->weight());
                // Alpha_Pim_cm[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->alpha_ppip_pipim(), _e->weight());
        }
}

void Histogram::writeMMSQ_mPim()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = 0; w < 15; w++)
                {
                        MMSQ_mPim_hist[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        if (MMSQ_mPim_hist[q2][w]->GetEntries())
                                MMSQ_mPim_hist[q2][w]->Write();

                        // MMSQ_mPim_hist_with_cut[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        // // if (MMSQ_mPim_hist_with_cut[q2][w]->GetEntries())
                        // MMSQ_mPim_hist_with_cut[q2][w]->Write();
                }
        }
}
void Histogram::write_Inv_Mass_hist()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = 0; w < 15; w++)
                {
                        Inv_mass_pPip[q2][w]->Write();
                        Inv_mass_pPim[q2][w]->Write();
                        Inv_mass_pipPim[q2][w]->Write();

                        Alpha_Prot_cm[q2][w]->Write();
                        Alpha_Pip_cm[q2][w]->Write();
                        Alpha_Pim_cm[q2][w]->Write();
                }
        }
}
void Histogram::Fill_MMSQ_mPim_1_comb(const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                MMSQ_mPim_hist_1_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());
        }
}

void Histogram::writeMMSQ_mPim_1_comb()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = 0; w < 15; w++)
                {
                        MMSQ_mPim_hist_1_comb[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        if (MMSQ_mPim_hist_1_comb[q2][w]->GetEntries())
                                MMSQ_mPim_hist_1_comb[q2][w]->Write();
                }
        }
}

void Histogram::Fill_MMSQ_mPim_2_comb(const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                MMSQ_mPim_hist_2_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());
                // MMSQ_mPim_hist_2_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim_swapped(), _e->weight());
        }
}

void Histogram::writeMMSQ_mPim_2_comb()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = 0; w < 15; w++)
                {
                        MMSQ_mPim_hist_2_comb[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        if (MMSQ_mPim_hist_2_comb[q2][w]->GetEntries())
                                MMSQ_mPim_hist_2_comb[q2][w]->Write();
                }
        }
}

void Histogram::Fill_MMSQ_mPim_3_comb(float dv2, const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                // MMSQ_mPim_hist_3_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim() - _e->MM2_mPim_swapped(), _e->weight());
                MMSQ_mPim_hist_3_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(dv2, _e->weight());
        }
}

void Histogram::writeMMSQ_mPim_3_comb()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = 0; w < 15; w++)
                {
                        MMSQ_mPim_hist_3_comb[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        if (MMSQ_mPim_hist_3_comb[q2][w]->GetEntries())
                                MMSQ_mPim_hist_3_comb[q2][w]->Write();
                }
        }
}

void Histogram::Fill_MMSQ_mPim_4_or_more_comb(float dv2, const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                // MMSQ_mPim_hist_4_or_more_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(_e->MM2_mPim(), _e->weight());
                MMSQ_mPim_hist_4_or_more_comb[q2_bining(_e->Q2())][int((_e->W() - 1.4) / 0.05)]->Fill(dv2, _e->weight());
        }
}

void Histogram::writeMMSQ_mPim_4_or_more_comb()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = 0; w < 15; w++)
                {
                        MMSQ_mPim_hist_4_or_more_comb[q2][w]->SetXTitle("MMSQ(GeV2/c2)");
                        if (MMSQ_mPim_hist_4_or_more_comb[q2][w]->GetEntries())
                                MMSQ_mPim_hist_4_or_more_comb[q2][w]->Write();
                }
        }
}

//////////////////////////////////////////  mmsq 3-d part   //////////////////////////////////
//////////////////////////////////////////  mmsq 3-d part   //////////////////////////////////
//////////////////////////////////////////  mmsq 3-d part   //////////////////////////////////
//////////////////////////////////////////  mmsq 3-d part   //////////////////////////////////

void Histogram::Fill_hist1D_mm2_mPim_inv_mass(const std::shared_ptr<Reaction> &_e)
{
        // inv_mass_pPim->Fill(_e->inv_Ppim(), _e->weight());

        if (_e->W() <= 2.2 && _e->W() >= 1.4)
        {

                if (_e->Q2() >= 2.0 && _e->Q2() <= 9.0)
                {
                        // inv_mass_pipPim->Fill(_e->inv_pip_pim(), _e->weight());

                        int inv_pPip_bin_val = inv_binning(_e->W(), _e->inv_Ppip(), 1);
                        if (inv_pPip_bin_val != -1)
                        {
                                mm2_mPim_hist_inv_pPip[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)][inv_pPip_bin_val]->Fill(_e->MM2_mPim(), _e->weight());
                        }

                        // int inv_pPim_bin_val = inv_binning(_e->W(), _e->inv_Ppim(), 1);
                        // if (inv_pPim_bin_val != -1)
                        // {
                        //         mm2_mPim_hist_inv_pPim[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)][inv_pPim_bin_val]->Fill(_e->MM2_mPim(), _e->weight());
                        // }
                        // int inv_pipPim_bin_val = inv_binning(_e->W(), _e->inv_pip_pim(), 0);
                        // if (inv_pipPim_bin_val != -1)
                        // {
                        //         mm2_mPim_hist_inv_pipPim[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)][inv_pipPim_bin_val]->Fill(_e->MM2_mPim(), _e->weight());
                        // }
                }
        }
}

void Histogram::writeHists1D_mm2_mPim_inv_pPip()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        for (size_t xi = 0; xi < 14; xi++)
                        {
                                mm2_mPim_hist_inv_pPip[q2][w][xi]->Write();
                        }
                }
        }
}
void Histogram::writeHists1D_mm2_mPim_inv_pPim()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        for (size_t xi = 0; xi < 14; xi++)
                        {
                                mm2_mPim_hist_inv_pPim[q2][w][xi]->Write();
                        }
                }
        }
}
void Histogram::writeHists1D_mm2_mPim_inv_pipPim()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        for (size_t xi = 0; xi < 14; xi++)
                        {
                                mm2_mPim_hist_inv_pipPim[q2][w][xi]->Write();
                        }
                }
        }
}

////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////
////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////
////// //////////////////////// This section is for THNSPARSE and related ////////////////////////////////

void Histogram::Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        // std::cout << "q2 value = " << _e->Q2() << "  q2 bin = " << q2_bining(_e->Q2()) << " inv pPip is outside ...   = " << _e->inv_Ppip() << std::endl;

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ///////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_prot[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_prot[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_prot[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_prot[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}

void Histogram::writeHists7D_prot()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        // std::cout << "q2 in write  " << q2 << " w is " << w << std::endl;

                        sevenDHist_prot[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_prot_evt(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {

                        ///////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_prot_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_prot_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_prot[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_prot_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_prot_evt()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_prot_evt[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        // float mmsq_low_values_for_bkg[2][2][3] = {{{-0.004, -0.028, 0.763}, {0.002, 0.071, 1.003}}};

                        ////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pip[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pip[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_pip[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pip[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pip[q2][w]->Write();
                }
        }
}
void Histogram::Fill_histSevenD_pip_evt(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////   if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pip_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pip_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_pip[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pip_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip_evt()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pip_evt[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////      if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pim[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pim[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_pim[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pim[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pim[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim_evt(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pim_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pim_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_pim[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pim_evt[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim_evt()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pim_evt[q2][w]->Write();
                }
        }
}

/////////////////////////////////////////////  tight cut  //////////////////////////////////////////////////
/////////////////////////////////////////////  tight cut //////////////////////////////////////////////////
/////////////////////////////////////////////  tight cut  //////////////////////////////////////////////////
/////////////////////////////////////////////  tight cut //////////////////////////////////////////////////

void Histogram::Fill_histSevenD_prot_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        // std::cout << "q2 value = " << _e->Q2() << "  q2 bin = " << q2_bining(_e->Q2()) << " inv pPip is outside ...   = " << _e->inv_Ppip() << std::endl;

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////// if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_prot_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_prot_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_prot_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_prot_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}

void Histogram::writeHists7D_prot_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        // std::cout << "q2 in write  " << q2 << " w is " << w << std::endl;

                        sevenDHist_prot_tight[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_prot_evt_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {

                        //////////    if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_prot_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_prot_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_prot_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_prot_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_prot_evt_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_prot_evt_tight[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pip_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        //////////    if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pip_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pip_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pi_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pip_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pip_tight[q2][w]->Write();
                }
        }
}
void Histogram::Fill_histSevenD_pip_evt_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ///////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pip_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pip_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pip_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pip_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip_evt_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pip_evt[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        /////////////// if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pim_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pim_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pim_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pim_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pim_tight[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim_evt_tight(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pim_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pim_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_tight[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pim_tight[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pim_evt_tight[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim_evt_tight()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pim_evt_tight[q2][w]->Write();
                }
        }
}

/////////////////////////////////////////////  loose cut  //////////////////////////////////////////////////
/////////////////////////////////////////////  loose cut //////////////////////////////////////////////////
/////////////////////////////////////////////  loose cut  //////////////////////////////////////////////////
/////////////////////////////////////////////  loose cut //////////////////////////////////////////////////

void Histogram::Fill_histSevenD_prot_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        // std::cout << "q2 value = " << _e->Q2() << "  q2 bin = " << q2_bining(_e->Q2()) << " inv pPip is outside ...   = " << _e->inv_Ppip() << std::endl;

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_prot_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_prot_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_prot_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_prot_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}

void Histogram::writeHists7D_prot_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {

                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        // std::cout << "q2 in write  " << q2 << " w is " << w << std::endl;

                        sevenDHist_prot_loose[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_prot_evt_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {

                        ////////////// if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_prot_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_prot_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // sevenDHist_prot_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_prot_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_prot_evt_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_prot_evt_loose[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pip_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pip_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pip_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pi_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pip_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pip_loose[q2][w]->Write();
                }
        }
}
void Histogram::Fill_histSevenD_pip_evt_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ///////////////  if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pip_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pip_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pip_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pip_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pip_evt_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pip_evt[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();

        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        /////////////////   if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        sevenDHist_pim_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight());
                                else
                                        sevenDHist_pim_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pim_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                sevenDHist_pim_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        sevenDHist_pim_loose[q2][w]->Write();
                }
        }
}

void Histogram::Fill_histSevenD_pim_evt_loose(const std::shared_ptr<Reaction> &_e)
{
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] = _e->alpha_ppip_pipim();
        if (_e->W() <= 2.2 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 9.0)
        {
                if (MM_cut(_e->W(), _e->Q2(), _e->MM2_mPim()))
                {
                        ///////////// if (((_e->MM2_exclusive() < -0.004) || (_e->MM2_exclusive() > 0.002)) && ((_e->MM2_mpip() < -0.028) || (_e->MM2_mpip() > 0.071)) && ((_e->MM2_mprot() < 0.763) || (_e->MM2_mprot() > 1.003)))
                        {
                                TThread::Lock();
                                if (!_bkg)
                                        h_5dim_pim_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, _e->weight() * _e->weight());
                                else
                                        h_5dim_pim_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->Fill(x, pow(_e->weight(), 2) * background_fact_loose[int((_e->W() - 1.0) / 0.05) - 8][q2_bining(_e->Q2()) - 1]);

                                // // sevenDHist_pim_loose[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                                TThread::UnLock();
                                h_5dim_pim_evt_loose[q2_bining(_e->Q2())][int((_e->W() - 1.0) / 0.05)]->GetNbins();
                        }
                }
        }
}
void Histogram::writeHists7D_pim_evt_loose()
{
        for (size_t q2 = 1; q2 < q2_bin_size; q2++)
        {
                for (size_t w = w_lower_bin; w < w_higher_bin; w++)
                {
                        h_5dim_pim_evt_loose[q2][w]->Write();
                }
        }
}

///////////////////////////////////////////  w-q2 and fundamental part //////////////////////////////////////////////////
///////////////////////////////////////////  w-q2 and fundamental part //////////////////////////////////////////////////
///////////////////////////////////////////  w-q2 and fundamental part //////////////////////////////////////////////////
///////////////////////////////////////////  w-q2 and fundamental part //////////////////////////////////////////////////
void Histogram::Fill_cdfd_prot(float dp, float dth, float dphi, const std::shared_ptr<Reaction> &_e)
{
        dp_prot_cdfd_hist->Fill(dp);
        dth_prot_cdfd_hist->Fill(dth);
        dphi_prot_cdfd_hist->Fill(dphi);
}

void Histogram::Fill_cdfd_pip(float dp, float dth, float dphi, const std::shared_ptr<Reaction> &_e)
{
        dp_pip_cdfd_hist->Fill(dp);
        dth_pip_cdfd_hist->Fill(dth);
        dphi_pip_cdfd_hist->Fill(dphi);
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction> &_e)
{
        short sec = _e->sec();
        TThread::Lock(); // Lock the thread to ensure exclusive access to the histograms
        // if (_e->MM2_mPim() > -0.06 && _e->MM2_mPim() < 0.08)
        {
                W_vs_q2->Fill(_e->W(), _e->Q2());
                W_hist->Fill(_e->W());
                Q2_hist->Fill(_e->Q2());
                // weight_hist->Fill(_e->weight());
                weight_hist->Fill(_e->weight());

                inv_mass_pPip->Fill(_e->inv_Ppip(), _e->weight());
                inv_mass_pPim->Fill(_e->inv_Ppim(), _e->weight());
                inv_mass_pipPim->Fill(_e->inv_pip_pim(), _e->weight());

                theta_Prot_cm->Fill(_e->prot_theta(), _e->weight());
                theta_Pip_cm->Fill(_e->pip_theta(), _e->weight());
                theta_Pim_cm->Fill(_e->pim_theta(), _e->weight());

                phi_Prot_cm->Fill(_e->prot_Phi(), _e->weight());
                phi_Pip_cm->Fill(_e->pip_Phi(), _e->weight());
                phi_Pim_cm->Fill(_e->pim_Phi(), _e->weight());

                alpha_Prot_cm->Fill(_e->alpha_pippim_pipf(), _e->weight());
                alpha_Pip_cm->Fill(_e->alpha_ppim_pipip(), _e->weight());
                alpha_Pim_cm->Fill(_e->alpha_ppip_pipim(), _e->weight());

                // inv_mass_pPip_swapped->Fill(_e->inv_Ppip_swapped(), _e->weight());
                // inv_mass_pPim_swapped->Fill(_e->inv_Ppim_swapped(), _e->weight());
                // inv_mass_pipPim_swapped->Fill(_e->inv_pip_pim_swapped(), _e->weight());

                // theta_Prot_cm_swapped->Fill(_e->prot_theta_swapped(), _e->weight());
                // theta_Pip_cm_swapped->Fill(_e->pip_theta_swapped(), _e->weight());
                // theta_Pim_cm_swapped->Fill(_e->pim_theta_swapped(), _e->weight());

                // // phi_Prot_cm_swapped->Fill(_e->prot_Phi_swapped(), _e->weight());
                // // phi_Pip_cm_swapped->Fill(_e->pip_Phi_swapped(), _e->weight());
                // // phi_Pim_cm_swapped->Fill(_e->pim_Phi_swapped(), _e->weight());

                // alpha_Prot_cm_swapped->Fill(_e->alpha_pippim_pipf_swapped(), _e->weight());
                // alpha_Pip_cm_swapped->Fill(_e->alpha_ppim_pipip_swapped(), _e->weight());
                // alpha_Pim_cm_swapped->Fill(_e->alpha_ppip_pipim_swapped(), _e->weight());

                // // if (_e->TwoPion_missingPim())
                //// {

                // if(_e->cut1(_e->Energy_excl()) &&  _e->cut3(_e->MM2_mpip()) &&  _e->cut4(_e->MM2_mprot())) {
                //         MM_twoPi->Fill(_e->MM_mPim(), _e->weight());
                MM_twoPi_mPim->Fill(_e->MM_mPim(), _e->weight());
                MM2_twoPi_mPim->Fill(_e->MM2_mPim());
                // MM2_twoPi_mPim->Fill(_e->MM2_mPim(), _e->weight());

                W_vs_MM->Fill(_e->W(), _e->MM_mPim());
                W_vs_MM2->Fill(_e->W(), _e->MM2_mPim());
                TThread::UnLock(); // Unlock after the operation
        }
        // }
        /*
        if (_e->TwoPion_exclusive())
        {
                // MM_twoPi->Fill(_e->MM2_mPim(), _e->weight());

                // if(_e->cut2(_e->MM2_mPim()) &&  _e->cut3(_e->MM2_mpip()) &&  _e->cut4(_e->MM2_mprot())) {
                missing_Energy_hist->Fill(_e->Energy_excl(), _e->weight());
                MM2_twoPi_excl->Fill(_e->MM2_exclusive(), _e->weight());
                MM_twoPi_excl->Fill(_e->MM_exclusive(), _e->weight());
                // }

                // if (_e->TwoPion_missingPip())
                // {
                // if(_e->cut1(_e->Energy_excl()) &&  _e->cut2(_e->MM2_mPim()) &&  _e->cut4(_e->MM2_mprot())) {
                MM2_twoPi_missingPip->Fill(_e->MM2_mpip(), _e->weight());
                //}
                // }
                // if (_e->TwoPion_missingProt())
                // {
                // if(_e->cut1(_e->Energy_excl()) &&  _e->cut2(_e->MM2_mPim()) &&  _e->cut3(_e->MM2_mpip())) {
                MM2_twoPi_missingProt->Fill(_e->MM2_mprot(), _e->weight());
                //}
        }
        */
        /*
        if (sec > 0 && sec <= 6)
        {
                W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
                W_sec[sec - 1]->Fill(_e->W(), _e->weight());

                if (_e->TwoPion_missingPim())
                {
                        MM_mPim_twoPi_sec[sec - 1]->Fill(_e->MM_mPim(), _e->weight());
                        MM2_mPim_twoPi_sec[sec - 1]->Fill(_e->MM2_mPim(), _e->weight());
                }
                if (_e->TwoPion_missingPip())
                        MM2_twoPi_missingPip_sec[sec - 1]->Fill(_e->MM2_mpip(), _e->weight());
                if (_e->TwoPion_missingProt())
                        MM2_twoPi_missingProt_sec[sec - 1]->Fill(_e->MM2_mprot(), _e->weight());
        }

        short det = _e->det();
        if (det == 1 && _e->W() <= 3.5)
        {
                W_det[0]->Fill(_e->W(), _e->weight());
                WQ2_det[0]->Fill(_e->W(), _e->Q2(), _e->weight());
        }
        else if (det == 2)
        {
                W_det[1]->Fill(_e->W(), _e->weight());
                WQ2_det[1]->Fill(_e->W(), _e->Q2(), _e->weight());
        }
        else
        {
                W_det[2]->Fill(_e->W(), _e->weight());
                WQ2_det[2]->Fill(_e->W(), _e->Q2(), _e->weight());
        }*/
}

// void Histogram::Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<MCReaction> &_e)
// {
//         // short sec = _e->sec();
//         // weight_hist->Fill(_e->weight());

//         W_vs_q2_twoPi_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
//         // W_hist_twoPi_thrown->Fill(_e->W_mc(), _e->weight());
//         // Q2_hist_twoPi_thrown->Fill(_e->Q2_mc(), _e->weight());
//         // W_vs_Q2_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
//         W_thrown->Fill(_e->W_mc(), _e->weight());
//         Q2_thrown->Fill(_e->Q2_mc(), _e->weight());

//         // theta_prot_thrown->Fill(_e->MCprot_theta_thrown(), _e->weight());
//         // theta_pip_thrown->Fill(_e->MCpip_theta_thrown(), _e->weight());
//         // theta_pim_thrown->Fill(_e->MCpim_theta_thrown(), _e->weight());

//         // mc_pid_at_zero->Fill(_d->mc_pid(0));
// }
void Histogram::Write_WvsQ2()
{
        W_thrown->SetXTitle("W_thrown (GeV)");
        if (W_thrown->GetEntries())
                W_thrown->Write();

        Q2_thrown->SetXTitle("Q2_thrown (GeV)");
        Q2_thrown->Write();

        W_vs_q2_twoPi_thrown->SetXTitle("W_thrown (GeV)");
        W_vs_q2_twoPi_thrown->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2_twoPi_thrown->SetOption("COLZ1");
        if (W_vs_q2_twoPi_thrown->GetEntries())
                W_vs_q2_twoPi_thrown->Write();

        W_vs_Q2_thrown->SetXTitle("W (GeV)");
        W_vs_Q2_thrown->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_Q2_thrown->SetOption("COLZ1");
        if (W_vs_Q2_thrown->GetEntries())
                W_vs_Q2_thrown->Write();

        weight_hist->SetXTitle("weight");
        weight_hist->Write();
        // mc_pid_at_zero->SetXTitle("mc pid at zero");
        // mc_pid_at_zero->Write();

        // pid_at_zero->SetXTitle("pid at zero");
        // pid_at_zero->Write();

        dp_prot_cdfd_hist->SetXTitle("dp");
        dp_prot_cdfd_hist->Write();
        dp_pip_cdfd_hist->SetXTitle("dp");
        dp_pip_cdfd_hist->Write();

        dth_prot_cdfd_hist->SetXTitle("dp");
        dth_prot_cdfd_hist->Write();
        dth_pip_cdfd_hist->SetXTitle("dp");
        dth_pip_cdfd_hist->Write();

        dphi_prot_cdfd_hist->SetXTitle("dp");
        dphi_prot_cdfd_hist->Write();
        dphi_pip_cdfd_hist->SetXTitle("dp");
        dphi_pip_cdfd_hist->Write();

        inv_mass_pPip->SetXTitle("Mass (GeV)");
        inv_mass_pPip->Write();
        inv_mass_pPim->SetXTitle("Mass (GeV)");
        inv_mass_pPim->Write();
        inv_mass_pipPim->SetXTitle("Mass (GeV)");
        inv_mass_pipPim->Write();

        theta_Prot_cm->SetXTitle("Theta (deg)");
        theta_Prot_cm->Write();
        theta_Pip_cm->SetXTitle("Theta (deg)");
        theta_Pip_cm->Write();
        theta_Pim_cm->SetXTitle("Theta (deg)");
        theta_Pim_cm->Write();

        phi_Prot_cm->SetXTitle("phi (deg)");
        phi_Prot_cm->Write();
        phi_Pip_cm->SetXTitle("phi (deg)");
        phi_Pip_cm->Write();
        phi_Pim_cm->SetXTitle("phi (deg)");
        phi_Pim_cm->Write();

        alpha_Prot_cm->SetXTitle("alpha (deg)");
        alpha_Prot_cm->Write();
        alpha_Pip_cm->SetXTitle("alpha (deg)");
        alpha_Pip_cm->Write();
        alpha_Pim_cm->SetXTitle("alpha (deg)");
        alpha_Pim_cm->Write();

        inv_mass_pPip_swapped->SetXTitle("Mass (GeV)");
        inv_mass_pPip_swapped->Write();
        inv_mass_pPim_swapped->SetXTitle("Mass (GeV)");
        inv_mass_pPim_swapped->Write();
        inv_mass_pipPim_swapped->SetXTitle("Mass (GeV)");
        inv_mass_pipPim_swapped->Write();

        theta_Prot_cm_swapped->SetXTitle("Theta (deg)");
        theta_Prot_cm_swapped->Write();
        theta_Pip_cm_swapped->SetXTitle("Theta (deg)");
        theta_Pip_cm_swapped->Write();
        theta_Pim_cm_swapped->SetXTitle("Theta (deg)");
        theta_Pim_cm_swapped->Write();

        // phi_Prot_cm_swapped->SetXTitle("phi (deg)");
        // phi_Prot_cm_swapped->Write();
        // phi_Pip_cm_swapped->SetXTitle("phi (deg)");
        // phi_Pip_cm_swapped->Write();
        // phi_Pim_cm_swapped->SetXTitle("phi (deg)");
        // phi_Pim_cm_swapped->Write();

        alpha_Prot_cm_swapped->SetXTitle("alpha (deg)");
        alpha_Prot_cm_swapped->Write();
        alpha_Pip_cm_swapped->SetXTitle("alpha (deg)");
        alpha_Pip_cm_swapped->Write();
        alpha_Pim_cm_swapped->SetXTitle("alpha (deg)");
        alpha_Pim_cm_swapped->Write();

        W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2->SetXTitle("W (GeV)");
        W_vs_q2->SetOption("COLZ");
        W_vs_q2->Write();
        for (short i = 0; i < 3; i++)
        {
                WQ2_det[i]->SetXTitle("W (GeV)");
                WQ2_det[i]->SetYTitle("Q^{2} (GeV^2)");
                WQ2_det[i]->SetOption("COLZ");
                if (WQ2_det[i]->GetEntries())
                        WQ2_det[i]->Write();
                W_det[i]->SetXTitle("W (GeV)");
                if (W_det[i]->GetEntries())
                        W_det[i]->Write();
        }

        // // W_vs_MM->SetXTitle("W (GeV)");
        // // W_vs_MM->SetYTitle("MM (GeV)");
        // // W_vs_MM->SetOption("COLZ1");
        // // if (W_vs_MM->GetEntries())
        // //         W_vs_MM->Write();
        // // W_vs_MM2->SetXTitle("W (GeV)");
        // // W_vs_MM2->SetYTitle("MMSQ (GeV^{2})");
        // // W_vs_MM2->SetOption("COLZ1");
        // // if (W_vs_MM2->GetEntries())
        // //         W_vs_MM2->Write();

        // // W_vs_q2->SetXTitle("W (GeV)");
        // // W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
        // // W_vs_q2->SetOption("COLZ1");
        // // if (W_vs_q2->GetEntries())
        // //         W_vs_q2->Write();

        W_hist->SetXTitle("W (GeV)");
        if (W_hist->GetEntries())
                W_hist->Write();

        // W_P2pi_hist->SetXTitle("W_P2pi (GeV)");
        // if (W_P2pi_hist->GetEntries())
        //         W_P2pi_hist->Write();

        Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
        if (Q2_hist->GetEntries())
                Q2_hist->Write();

        // W_vs_q2_twoPi->SetXTitle("W (GeV)");
        // W_vs_q2_twoPi->SetYTitle("Q^{2} (GeV^{2})");
        // W_vs_q2_twoPi->SetOption("COLZ1");
        // if (W_vs_q2_twoPi->GetEntries())
        //         W_vs_q2_twoPi->Write();

        // // W_hist_twoPi->SetXTitle("W (GeV)");
        // // if (W_hist_twoPi->GetEntries())
        // //         W_hist_twoPi->Write();

        // Q2_hist_twoPi->SetXTitle("Q^{2} (GeV^{2})");
        // if (Q2_hist_twoPi->GetEntries())
        //         Q2_hist_twoPi->Write();

        auto mmsq_4_topology = RootOutputFile->mkdir("mmsq_4_topology");
        mmsq_4_topology->cd();

        // MM2_twoPi_excl->SetXTitle("MMSQ Excl (GeV^{2})");
        // if (MM2_twoPi_excl->GetEntries())
        //         MM2_twoPi_excl->Write();

        // MM_twoPi_excl->SetXTitle("MM Excl (GeV)");
        // if (MM_twoPi_excl->GetEntries())
        //         MM_twoPi_excl->Write();

        MM2_twoPi_mPim->SetXTitle("MMSQ mPim (GeV^{2})");
        if (MM2_twoPi_mPim->GetEntries())
                MM2_twoPi_mPim->Write();

        // MM2_twoPi_missingPip->SetXTitle("MM2 mPip (GeV^{2})");
        // if (MM2_twoPi_missingPip->GetEntries())
        //         MM2_twoPi_missingPip->Write();

        // MM2_twoPi_missingProt->SetXTitle("MM2 mProt (GeV^{2})");
        // if (MM2_twoPi_missingProt->GetEntries())
        //         MM2_twoPi_missingProt->Write();

        // // missing_Energy_hist->SetXTitle("me (GeV)");
        // // if (missing_Energy_hist->GetEntries())
        // //         missing_Energy_hist->Write();

        // // MM_twoPi_mPim->SetXTitle("MM mPim (GeV^{2})");
        // // if (MM_twoPi_mPim->GetEntries())
        // //         MM_twoPi_mPim->Write();

        // // auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
        // // wvsq2_sec->cd();
        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
        // //         W_vs_q2_sec[i]->SetXTitle("W (GeV)");
        // //         W_vs_q2_sec[i]->SetOption("COLZ1");
        // //         W_vs_q2_sec[i]->Write();
        // // }
        // // auto w_sec = RootOutputFile->mkdir("w_sec");
        // // w_sec->cd();
        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         W_sec[i]->SetXTitle("W (GeV)");

        // //         //  W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.85, 1.05);
        // //         // gStyle->SetOptFit(01);
        // //         W_sec[i]->Write();
        // // }

        ///////////////////////// numner of photoelectrons in HTCC //////////////////////////
        auto nphe_htcc_sec = RootOutputFile->mkdir("nphe_htcc_sec");
        nphe_htcc_sec->cd();

        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                // for (short i = 0; i < num_sectors; i++)
                // {
                //         htcc_nphe_sec[c][i]->SetXTitle("Nphe");
                //         htcc_nphe_sec[c][i]->SetYTitle("Count");

                //         // if (htcc_nphe_sec[c][i]->GetEntries())
                //         htcc_nphe_sec[c][i]->Write();
                // }
        }
        ///////////////////////// //////////////////////////
        auto Vz_sec = RootOutputFile->mkdir("Vz_sec");
        Vz_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        vz_sec[c][i]->SetXTitle("Vz (cm)");
                        if (vz_sec[c][i]->GetEntries())
                                vz_sec[c][i]->Write();
                }
        }

        ///////////////////////// //////////////////////////
        auto Pcal_lu_sec = RootOutputFile->mkdir("pcal_lu_sec");
        Pcal_lu_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        pcal_lu_sec[c][i]->SetXTitle("Pcal lu (cm)");
                        if (pcal_lu_sec[c][i]->GetEntries())
                                pcal_lu_sec[c][i]->Write();
                }
        }
        ///////////////////////// //////////////////////////
        auto Pcal_lv_sec = RootOutputFile->mkdir("pcal_lv_sec");
        Pcal_lv_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        pcal_lv_sec[c][i]->SetXTitle("Pcal lv (cm)");
                        if (pcal_lv_sec[c][i]->GetEntries())
                                pcal_lv_sec[c][i]->Write();
                }
        }
        ///////////////////////// //////////////////////////
        auto Pcal_lw_sec = RootOutputFile->mkdir("pcal_lw_sec");
        Pcal_lw_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        pcal_lw_sec[c][i]->SetXTitle("Pcal lw (cm)");
                        if (pcal_lw_sec[c][i]->GetEntries())
                                pcal_lw_sec[c][i]->Write();
                }
        }
        // ///////////////////////// //////////////////////////
        // auto ECAL_VS_PCAL_sec = RootOutputFile->mkdir("ECAL_VS_PCAL_sec");
        // ECAL_VS_PCAL_sec->cd();
        // for (auto &&cut : before_after_cut)
        // {
        //         int c = cut.first;

        //         for (short i = 0; i < num_sectors; i++)
        //         {
        //                 for (short j = 0; j < 9; j++) // mom bins
        //                 {
        //                         ECAL_VS_PCAL[c][i][j]->SetOption("COLZ1");
        //                         ECAL_VS_PCAL[c][i][j]->SetYTitle("SF PCAL ");
        //                         ECAL_VS_PCAL[c][i][j]->SetXTitle("SF ECIN");
        //                         // if (ECAL_VS_PCAL[i][j]->GetEntries())
        //                         ECAL_VS_PCAL[c][i][j]->Write();
        //                 }
        //         }
        // }

        /////////////////////////  //////////////////////////
        auto SF_VS_MOM_sec = RootOutputFile->mkdir("SF_VS_MOM_sec");
        SF_VS_MOM_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        SF_VS_MOM[c][i]->SetOption("COLZ1");
                        SF_VS_MOM[c][i]->SetYTitle("SF ");
                        SF_VS_MOM[c][i]->SetXTitle("MOM (GeV)");
                        // if (SF_VS_MOM[i]->GetEntries())
                        SF_VS_MOM[c][i]->Write();
                }
        }
        /////////////////////////  //////////////////////////
        auto chi2pid_elec_sec = RootOutputFile->mkdir("chi2pid_elec_sec");
        chi2pid_elec_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        elec_Chi2pid_sec[c][i]->SetXTitle("Chi2pid");
                        elec_Chi2pid_sec[c][i]->SetYTitle("Count");

                        if (elec_Chi2pid_sec[c][i]->GetEntries())
                                elec_Chi2pid_sec[c][i]->Write();
                }
        }

        auto elec_theta_vs_mom_fd = RootOutputFile->mkdir("elec_theta_vs_mom_fd");
        elec_theta_vs_mom_fd->cd();
        for (short i = 0; i < num_sectors; i++)
        {
                Theta_fd_elec_lab_vs_mom_elec[i]->SetXTitle("mom_prot (GeV)");
                Theta_fd_elec_lab_vs_mom_elec[i]->SetYTitle("theta_prot (Deg)");
                Theta_fd_elec_lab_vs_mom_elec[i]->SetOption("COLZ1");
                if (Theta_fd_elec_lab_vs_mom_elec[i]->GetEntries())
                        Theta_fd_elec_lab_vs_mom_elec[i]->Write("");
        }

        // auto twoPi_sec = RootOutputFile->mkdir("twoPi_sec");
        // twoPi_sec->cd();
        // for (short i = 0; i < num_sectors; i++)
        // {
        //         W_vs_q2_twoPi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
        //         W_vs_q2_twoPi_sec[i]->SetXTitle("W (GeV)");
        //         W_vs_q2_twoPi_sec[i]->SetOption("COLZ1");
        //         if (W_vs_q2_twoPi_sec[i]->GetEntries())
        //                 W_vs_q2_twoPi_sec[i]->Write();
        // }

        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         //  if (MM_mPim_twoPi_sec[i]->GetEntries()) MM_mPim_twoPi_sec[i]->Fit("gaus", "QMR+",
        // //         //  "QMR+", -0.1, 0.1);
        // //         MM_mPim_twoPi_sec[i]->SetXTitle("MisingMass (GeV)");
        // //         MM_mPim_twoPi_sec[i]->Write();
        // // }

        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         // if (MM2_mPim_twoPi_sec[i]->GetEntries()) MM2_mPim_twoPi_sec[i]->Fit("gaus", "QMR+",
        // //         // "QMR+", -0.1, 0.1);
        // //         MM2_mPim_twoPi_sec[i]->SetXTitle("MM2 (GeV2)");
        // //         MM2_mPim_twoPi_sec[i]->Write();
        // // }

        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         MM2_twoPi_missingPip_sec[i]->SetXTitle("MM2 (GeV2)");
        // //         MM2_twoPi_missingPip_sec[i]->Write();
        // // }

        // // for (short i = 0; i < num_sectors; i++)
        // // {
        // //         MM2_twoPi_missingProt_sec[i]->SetXTitle("MM2 (GeV2)");
        // //         MM2_twoPi_missingProt_sec[i]->Write();
        // // }
}

void Histogram::makeHists_electron_cuts()
{
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                auto type = cut.second.c_str();
                EC_sampling_fraction[c] =
                    std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type),
                                           Form("EC_sampling_fraction%s", type), bins, 1.5,
                                           10.0, bins, zero, 0.5);

                ECin_sf_vs_PCAL_sf[c] = std::make_shared<TH2D>(
                    Form("ECin_sf_vs_PCAL_sf%s", type), Form("ECin_sf_vs_PCAL_sf%s", type), bins,
                    -0., 0.25, bins, 0, 0.25);
                vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type),
                                                        Form("vz_position%s", type), bins, -15, 15);
                momentum[c] = std::make_shared<TH1D>(Form("mom%s", type), Form("mom%s", type), bins, 0, 10);

                pcal_sec[c] =
                    std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s", type), bins, -420, 420, bins, -420, 420);
                pcal_sec_ineff_cuts[c] =
                    std::make_shared<TH2D>(Form("pcal_sec_ineff_cut%s", type), Form("pcal_sec_ineff_cut%s", type), bins, -420, 420, bins, -420, 420);

                pcal_hx_hy_sec[c] =
                    std::make_shared<TH2D>(Form("pcal_hx_hy_sec%s", type), Form("pcal_hx_hy_sec%s", type), bins, -420, 420, bins, -420, 420);
                dcr1_sec[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s", type), bins, -320, 320, bins, -320, 320);

                prot_Delta_vz_cut_fd[c] = std::make_shared<TH1D>(Form("fd_prot_dvz_position%s", type),
                                                                 Form("fd_prot_dvz_position%s", type), bins, -40, 40);

                prot_Chi2pid_cut_fd[c] = std::make_shared<TH1D>(Form("fd_prot_chi2pid%s", type),
                                                                Form("fd_prot_chi2pid%s", type), bins, -20, 20);

                prot_Delta_vz_cut_cd[c] = std::make_shared<TH1D>(Form("cd_prot_dvz_position%s", type),
                                                                 Form("cd_prot_dvz_position%s", type), bins, -40, 40);

                prot_Chi2pid_cut_cd[c] = std::make_shared<TH1D>(Form("cd_prot_chi2pid%s", type),
                                                                Form("cd_prot_chi2pid%s", type), bins, -20, 20);

                pip_Delta_vz_cut_fd[c] = std::make_shared<TH1D>(Form("fd_pip_dvz_position%s", type),
                                                                Form("fd_pip_dvz_position%s", type), bins, -40, 40);

                pip_Chi2pid_cut_fd[c] = std::make_shared<TH1D>(Form("fd_pip_chi2pid%s", type),
                                                               Form("fd_pip_chi2pid%s", type), bins, -20, 20);

                pip_Delta_vz_cut_cd[c] = std::make_shared<TH1D>(Form("cd_pip_dvz_position%s", type),
                                                                Form("cd_pip_dvz_position%s", type), bins, -40, 40);

                pip_Chi2pid_cut_cd[c] = std::make_shared<TH1D>(Form("cd_pip_chi2pid%s", type),
                                                               Form("cd_pip_chi2pid%s", type), bins, -20, 20);

                pim_Delta_vz_cut[c] = std::make_shared<TH1D>(Form("pim_dvz_position%s", type),
                                                             Form("pim_dvz_position%s", type), bins, -40, 40);

                pim_Chi2pid_cut[c] = std::make_shared<TH1D>(Form("pim_chi2pid%s", type),
                                                            Form("pim_chi2pid%s", type), bins, -20, 20);

                phi_vs_mom_prot_fd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_mom_fd_prot%s", type), Form("phi_vs_mom_fd_prot%s", type), bins, -180, 180, bins, 0, 5);
                phi_vs_mom_pip_fd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_mom_fd_pip%s", type), Form("phi_vs_mom_fd_pip%s", type), bins, -180, 180, bins, 0, 5);
                phi_vs_mom_pim_fd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_mom_fd_pim%s", type), Form("phi_vs_mom_fd_pim%s", type), bins, -180, 180, bins, 0, 5);

                theta_prot_fd[c] =
                    std::make_shared<TH1D>(Form("fd_theta_prot%s", type), Form("fd_theta_prot%s", type), bins, 0, 50);
                theta_pip_fd[c] = std::make_shared<TH1D>(Form("fd_theta_pip%s", type), Form("fd_theta_pip%s", type), bins, 0, 50);
                theta_pim_fd[c] = std::make_shared<TH1D>(Form("fd_theta_pim%s", type), Form("fd_theta_pim%s", type), bins, 0, 50);

                Theta_prot_lab_vs_mom_prot_fd[c] = std::make_shared<TH2D>(
                    Form("Theta_fd_prot_lab_vs_mom_prot%s", type), Form("Theta_fd_prot_lab_vs_mom_prot%s", type), bins, zero,
                    7.0, bins, 0, 50);
                Theta_pip_lab_vs_mom_pip_fd[c] = std::make_shared<TH2D>(
                    Form("Theta_fd_pip_lab_vs_mom_pip%s", type), Form("Theta_fd_pip_lab_vs_mom_pip%s", type), bins, zero, 7.0,
                    bins, 0, 50);
                Theta_pim_lab_vs_mom_pim_fd[c] = std::make_shared<TH2D>(
                    Form("Theta_fd_pim_lab_vs_mom_pim%s", type), Form("Theta_fd_pim_lab_vs_mom_pim%s", type), bins, zero, 7.0,
                    bins, 0, 50);

                phi_vs_momT_prot_cd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_momT_cd_prot%s", type), Form("phi_vs_momT_cd_prot%s", type), bins, -180, 180, bins, 0, 2);
                phi_vs_momT_pip_cd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_momT_cd_pip%s", type), Form("phi_vs_momT_cd_pip%s", type), bins, -180, 180, bins, 0, 2);
                phi_vs_momT_pim_cd[c] =
                    std::make_shared<TH2D>(Form("phi_vs_momT_cd_pim%s", type), Form("phi_vs_momT_cd_pim%s", type), bins, -180, 180, bins, 0, 2);

                theta_prot_cd[c] =
                    std::make_shared<TH1D>(Form("cd_theta_prot%s", type), Form("cd_theta_prot%s", type), bins, 30., 150);
                theta_pip_cd[c] = std::make_shared<TH1D>(Form("cd_theta_pip%s", type), Form("cd_theta_pip%s", type), bins, 30, 150);
                theta_pim_cd[c] = std::make_shared<TH1D>(Form("cd_theta_pim%s", type), Form("cd_theta_pim%s", type), bins, 30, 150);

                Theta_prot_lab_vs_mom_prot_cd[c] = std::make_shared<TH2D>(
                    Form("Theta_cd_prot_lab_vs_mom_prot%s", type), Form("Theta_cd_prot_lab_vs_mom_prot%s", type), bins, zero,
                    4.0, bins, 30, 150);
                Theta_pip_lab_vs_mom_pip_cd[c] = std::make_shared<TH2D>(
                    Form("Theta_cd_pip_lab_vs_mom_pip%s", type), Form("Theta_cd_pip_lab_vs_mom_pip%s", type), bins, zero, 4.0,
                    bins, 30, 150);
                Theta_pim_lab_vs_mom_pim_cd[c] = std::make_shared<TH2D>(
                    Form("Theta_cd_pim_lab_vs_mom_pim%s", type), Form("Theta_cd_pim_lab_vs_mom_pim%s", type), bins, zero, 4.0,
                    bins, 30, 150);

                dcr1_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_prot%s", type), Form("dcr1_sec_prot%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_prot%s", type), Form("dcr2_sec_prot%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_prot%s", type), Form("dcr3_sec_prot%s", type), bins, -420, 420, bins, -420, 420);

                dcr1_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_pip%s", type), Form("dcr1_sec_pip%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_pip%s", type), Form("dcr2_sec_pip%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_pip%s", type), Form("dcr3_sec_pip%s", type), bins, -420, 420, bins, -420, 420);

                dcr1_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_pim%s", type), Form("dcr1_sec_pim%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_pim%s", type), Form("dcr2_sec_pim%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_pim%s", type), Form("dcr3_sec_pim%s", type), bins, -360, 360, bins, -360, 360);
        }
}

void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e)
{
        pid_at_zero->Fill(_d->pid(0));
        mc_pid_at_zero->Fill(_d->mc_pid(0));

        auto elec_cuts = std::make_shared<Pass2_Cuts>(_d);

        int sec = _e->sec();
        int ith_part = 0;
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {

                vz_position[before_any_cuts]->Fill(_d->vz(0));
                if (_d->vz(0) > -10 && _d->vz(0) < 5)
                        vz_position[with_one_cut]->Fill(_d->vz(0));
                else
                        vz_position[outside_one_cut]->Fill(_d->vz(0));

                //// pcal
                pcal_sec[before_any_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                if (elec_cuts->EC_hit_position_fiducial_cut_homogeneous(condition_of_cut))
                        pcal_sec[with_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                else
                        pcal_sec[outside_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                //// pcal hx hy
                pcal_hx_hy_sec[before_any_cuts]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0));
                if (elec_cuts->PCAL_fiducial_cut_X_Y(condition_of_cut))
                        pcal_hx_hy_sec[with_one_cut]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0));
                else
                        pcal_hx_hy_sec[outside_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));

                //// pcal ineff
                pcal_sec_ineff_cuts[before_any_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                if (elec_cuts->PCAL_Ineff_cut_X_Y())
                        pcal_sec_ineff_cuts[with_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                else
                        pcal_sec_ineff_cuts[outside_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                // dc
                dcr1_sec[before_any_cuts]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
                dcr2_sec[before_any_cuts]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
                dcr3_sec[before_any_cuts]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));

                if (elec_cuts->DC_fiducial_cut_XY_E(condition_of_cut))
                {
                        dcr1_sec[with_one_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
                        dcr2_sec[with_one_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
                        dcr3_sec[with_one_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
                }
                else
                {
                        dcr1_sec[outside_one_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
                        dcr2_sec[outside_one_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
                        dcr3_sec[outside_one_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
                }

                momentum[before_any_cuts]->Fill(_d->p(0));
                if (_d->p(0) > 1.50)
                        momentum[with_one_cut]->Fill(_d->p(0));
                else
                        momentum[outside_one_cut]->Fill(_d->p(0));

                EC_sampling_fraction[before_any_cuts]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                if (elec_cuts->EC_sampling_fraction_cut(condition_of_cut))

                        EC_sampling_fraction[with_one_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                else
                        EC_sampling_fraction[outside_one_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));

                //
                ECin_sf_vs_PCAL_sf[before_any_cuts]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)));
                if (elec_cuts->EC_sampling_fraction_cut(condition_of_cut))

                        ECin_sf_vs_PCAL_sf[with_one_cut]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)));
                else
                        ECin_sf_vs_PCAL_sf[outside_one_cut]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)));
                /////////
                if (sec > 0 && sec <= 6)
                {
                        /// nphe
                        // htcc_nphe_sec[before_any_cuts][sec - 1]->Fill(_d->cc_htcc_nphe(0));
                        // if (_d->cc_htcc_nphe(0) > 2)
                        //         htcc_nphe_sec[with_one_cut][sec - 1]->Fill(_d->cc_htcc_nphe(0));
                        // else
                        //         htcc_nphe_sec[outside_one_cut][sec - 1]->Fill(_d->cc_htcc_nphe(0));

                        // chi2pid
                        elec_Chi2pid_sec[before_any_cuts][sec - 1]->Fill(_d->chi2pid(0));
                        if (_d->chi2pid(0) < 3)
                                elec_Chi2pid_sec[with_one_cut][sec - 1]->Fill(_d->chi2pid(0));
                        else
                                elec_Chi2pid_sec[outside_one_cut][sec - 1]->Fill(_d->chi2pid(0));
                        // vz cut
                        vz_sec[before_any_cuts][sec - 1]->Fill(_d->vz(0));
                        // if (_d->vz(0) > -(2.78 + 3 * 2.16) && _d->vz(0) < (-2.78 + 3 * 2.16))
                        if (_d->vz(0) > -10 && _d->vz(0) < 5)
                                vz_sec[with_one_cut][sec - 1]->Fill(_d->vz(0));
                        else
                                vz_sec[outside_one_cut][sec - 1]->Fill(_d->vz(0));

                        pcal_lu_sec[before_any_cuts][sec - 1]->Fill(_d->ec_pcal_lu(0));
                        pcal_lv_sec[before_any_cuts][sec - 1]->Fill(_d->ec_pcal_lv(0));
                        pcal_lw_sec[before_any_cuts][sec - 1]->Fill(_d->ec_pcal_lw(0));

                        if (_d->ec_pcal_lu(0) > 14.5)
                                pcal_lu_sec[with_one_cut][sec - 1]->Fill(_d->ec_pcal_lu(0));
                        else
                                pcal_lu_sec[outside_one_cut][sec - 1]->Fill(_d->ec_pcal_lu(0));
                        if (_d->ec_pcal_lv(0) > 14.5)
                                pcal_lv_sec[with_one_cut][sec - 1]->Fill(_d->ec_pcal_lv(0));
                        else
                                pcal_lv_sec[outside_one_cut][sec - 1]->Fill(_d->ec_pcal_lv(0));
                        if (_d->ec_pcal_lw(0) > 14.5)
                                pcal_lw_sec[with_one_cut][sec - 1]->Fill(_d->ec_pcal_lw(0));
                        else
                                pcal_lw_sec[outside_one_cut][sec - 1]->Fill(_d->ec_pcal_lw(0));

                        // ecal vs pcal
                        int momRangeIdx = getMomRange(_d->p(0));
                        // ECAL_VS_PCAL[before_any_cuts][sec - 1][momRangeIdx]->Fill((_d->ec_ecin_energy(0) / _d->p(0)), (_d->ec_pcal_energy(0) / _d->p(0)));
                        // if (elec_cuts->EC_inner_vs_EC_outer())
                        //         ECAL_VS_PCAL[with_one_cut][sec - 1][momRangeIdx]->Fill((_d->ec_ecin_energy(0) / _d->p(0)), (_d->ec_pcal_energy(0) / _d->p(0)));
                        // else
                        //         ECAL_VS_PCAL[outside_one_cut][sec - 1][momRangeIdx]->Fill((_d->ec_ecin_energy(0) / _d->p(0)), (_d->ec_pcal_energy(0) / _d->p(0)));

                        SF_VS_MOM[before_any_cuts][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                        if (elec_cuts->EC_sampling_fraction_cut(condition_of_cut))
                                SF_VS_MOM[with_one_cut][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                        else
                                SF_VS_MOM[outside_one_cut][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                }
        }
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e)
{
        int sec = _e->sec();
        int momRangeIdx1 = getMomRange(_d->p(0));

        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {

                vz_position[after_all_cuts]->Fill(_d->vz(0));
                pcal_sec[after_all_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
                pcal_hx_hy_sec[after_all_cuts]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0));
                pcal_sec_ineff_cuts[after_all_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));

                dcr1_sec[after_all_cuts]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
                dcr2_sec[after_all_cuts]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
                dcr3_sec[after_all_cuts]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
                EC_sampling_fraction[after_all_cuts]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
                ECin_sf_vs_PCAL_sf[after_all_cuts]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)));

                momentum[after_all_cuts]->Fill(_d->p(0));

                if (sec > 0 && sec <= 6)
                {
                        // htcc_nph
                        // e_sec[after_all_cuts][sec - 1]->Fill(_d->cc_htcc_nphe(0));
                        elec_Chi2pid_sec[after_all_cuts][sec - 1]->Fill(_d->chi2pid(0));
                        vz_sec[after_all_cuts][sec - 1]->Fill(_d->vz(0));
                        // ECAL_VS_PCAL[after_all_cuts][sec - 1][momRangeIdx1]->Fill((_d->ec_ecin_energy(0) / _d->p(0)), (_d->ec_pcal_energy(0) / _d->p(0)));
                        SF_VS_MOM[after_all_cuts][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));

                        Theta_fd_elec_lab_vs_mom_elec[sec - 1]->Fill(_e->elec_mom(), _e->elec_th());

                        pcal_lu_sec[after_all_cuts][sec - 1]->Fill(_d->ec_pcal_lu(0));
                        pcal_lv_sec[after_all_cuts][sec - 1]->Fill(_d->ec_pcal_lv(0));
                        pcal_lw_sec[after_all_cuts][sec - 1]->Fill(_d->ec_pcal_lw(0));
                }
        }
}

void Histogram::FillHists_prot_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {
                if (abs(_d->status(i)) < 4000)
                {
                        prot_Delta_vz_cut_fd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        prot_Chi2pid_cut_fd[before_any_cuts]->Fill(_d->chi2pid(i));

                        // theta_prot_fd[before_any_cuts]->Fill(_e->prot_theta_lab());
                        // phi_vs_mom_prot_fd[before_any_cuts]->Fill(_e->prot_Phi_lab(), _e->prot_momentum());
                        // Theta_prot_lab_vs_mom_prot_fd[before_any_cuts]->Fill(_e->prot_momentum(), _e->prot_theta_lab());

                        dcr1_sec_prot[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_prot[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_prot[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                        if (_cuts->Hadron_Delta_vz_cut(i, condition_of_cut))
                                prot_Delta_vz_cut_fd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)));
                        else
                                prot_Delta_vz_cut_fd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)));

                        if (_cuts->Hadron_Chi2pid_cut(i, condition_of_cut))
                                prot_Chi2pid_cut_fd[with_one_cut]->Fill(_d->chi2pid(i));
                        else
                                prot_Chi2pid_cut_fd[outside_one_cut]->Fill(_d->chi2pid(i));

                        if (_cuts->DC_fiducial_cut_XY_PROT(i, 1, condition_of_cut))
                        {
                                dcr1_sec_prot[with_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                                dcr2_sec_prot[with_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                                dcr3_sec_prot[with_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                        }
                        else
                        {

                                dcr1_sec_prot[outside_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                                dcr2_sec_prot[outside_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                                dcr3_sec_prot[outside_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                        }
                }
                else if (abs(_d->status(i)) >= 4000)
                {
                        prot_Delta_vz_cut_cd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        prot_Chi2pid_cut_cd[before_any_cuts]->Fill((_d->chi2pid(i)));
                        // phi_vs_momT_prot_cd[before_any_cuts]->Fill(_e->prot_Phi_lab(), _e->prot_momT());
                        // theta_prot_cd[before_any_cuts]->Fill(_e->prot_theta_lab());
                        // Theta_prot_lab_vs_mom_prot_cd[before_any_cuts]->Fill(_e->prot_momentum(), _e->prot_theta_lab());

                        if (_cuts->Hadron_Delta_vz_cut(i, condition_of_cut))
                                prot_Delta_vz_cut_cd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)));
                        else
                                prot_Delta_vz_cut_cd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)));

                        if (_cuts->Hadron_Chi2pid_cut(i, condition_of_cut))
                                prot_Chi2pid_cut_cd[with_one_cut]->Fill(_d->chi2pid(i));
                        else
                                prot_Chi2pid_cut_cd[outside_one_cut]->Fill(_d->chi2pid(i));

                        // if (_cuts->CD_fiducial_had(i))
                        //         phi_vs_momT_prot_cd[with_one_cut]->Fill(_e->prot_Phi_lab(), _e->prot_momT());
                        // else
                        //         phi_vs_momT_prot_cd[outside_one_cut]->Fill(_e->prot_Phi_lab(), _e->prot_momT());
                }
        }
}

void Histogram::FillHists_prot_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &prot)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0) // && _cuts->HadronsCuts(i))
        {
                int sec = _d->dc_sec(i);

                auto dt = std::make_shared<Delta_T>(_d);
                if (abs(_d->status(i)) < 4000)
                // if (dt->isCtof() == false)
                {
                        dcr1_sec_prot[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_prot[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_prot[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                        prot_Delta_vz_cut_fd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        phi_vs_mom_prot_fd[after_all_cuts]->Fill(_e->prot_Phi_lab(prot), _e->prot_momentum(prot));
                        prot_Chi2pid_cut_fd[after_all_cuts]->Fill((_d->chi2pid(i)));
                        theta_prot_fd[after_all_cuts]->Fill(_e->prot_theta_lab(prot));
                        Theta_prot_lab_vs_mom_prot_fd[after_all_cuts]->Fill(_e->prot_momentum(prot), _e->prot_theta_lab(prot));
                        if (sec > 0 && sec <= 6)
                        {
                                Theta_fd_prot_lab_vs_mom_prot[sec - 1]->Fill(_e->prot_momentum(prot), _e->prot_theta_lab(prot));
                        }
                }
                else if (abs(_d->status(i)) > 4000)
                // else if (dt->isCtof() == true)
                {
                        // std::cout << "  status of ctof particle :  " << _d->status(i) << "  prot lab angle theta " << _e->prot_theta_lab(prot) << std::endl;

                        prot_Delta_vz_cut_cd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        prot_Chi2pid_cut_cd[after_all_cuts]->Fill(_d->chi2pid(i));
                        phi_vs_momT_prot_cd[after_all_cuts]->Fill(_e->prot_Phi_lab(prot), _e->prot_momT(prot));
                        theta_prot_cd[after_all_cuts]->Fill(_e->prot_theta_lab(prot));
                        Theta_prot_lab_vs_mom_prot_cd[after_all_cuts]->Fill(_e->prot_momentum(prot), _e->prot_theta_lab(prot));
                }

                // std::cout << "weighr = " << _e->weight() << '\n';
        }
}
void Histogram::FillHists_pip_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {

                if (abs(_d->status(i)) < 4000)
                {
                        pip_Delta_vz_cut_fd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        pip_Chi2pid_cut_fd[before_any_cuts]->Fill(_d->chi2pid(i));
                        // phi_vs_mom_pip_fd[before_any_cuts]->Fill(_e->pip_Phi_lab(), _e->pip_momentum());
                        // theta_pip_fd[before_any_cuts]->Fill(_e->pip_theta_lab());
                        // Theta_pip_lab_vs_mom_pip_fd[before_any_cuts]->Fill(_e->pip_momentum(), _e->pip_theta_lab());

                        dcr1_sec_pip[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_pip[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_pip[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));

                        if (_cuts->Hadron_Delta_vz_cut(i, condition_of_cut))
                                pip_Delta_vz_cut_fd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)));
                        else
                                pip_Delta_vz_cut_fd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)));

                        if (_cuts->Hadron_Chi2pid_cut(i, condition_of_cut))
                                pip_Chi2pid_cut_fd[with_one_cut]->Fill(_d->chi2pid(i));
                        else
                                pip_Chi2pid_cut_fd[outside_one_cut]->Fill(_d->chi2pid(i));
                }
                else if (abs(_d->status(i)) > 4000)
                {
                        pip_Delta_vz_cut_cd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        pip_Chi2pid_cut_cd[before_any_cuts]->Fill(_d->chi2pid(i));
                        // phi_vs_momT_pip_cd[before_any_cuts]->Fill(_e->pip_Phi_lab(), _e->pip_momT());
                        // theta_pip_cd[before_any_cuts]->Fill(_e->pip_theta_lab());
                        // Theta_pip_lab_vs_mom_pip_cd[before_any_cuts]->Fill(_e->pip_momentum(), _e->pip_theta_lab());

                        // if (_cuts->CD_fiducial_had(i))
                        //         phi_vs_momT_pip_cd[with_one_cut]->Fill(_e->pip_Phi_lab(), _e->pip_momT());
                        // else
                        //         phi_vs_momT_pip_cd[outside_one_cut]->Fill(_e->pip_Phi_lab(), _e->pip_momT());

                        if (_cuts->Hadron_Delta_vz_cut(i, condition_of_cut))
                                pip_Delta_vz_cut_cd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)));
                        else
                                pip_Delta_vz_cut_cd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)));

                        if (_cuts->Hadron_Chi2pid_cut(i, condition_of_cut))
                                pip_Chi2pid_cut_cd[with_one_cut]->Fill(_d->chi2pid(i));
                        else
                                pip_Chi2pid_cut_cd[outside_one_cut]->Fill(_d->chi2pid(i));
                }

                if (_cuts->DC_fiducial_cut_XY_PIP(i, 2, condition_of_cut))
                {
                        dcr1_sec_pip[with_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_pip[with_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_pip[with_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                }
                else
                {
                        dcr1_sec_pip[outside_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_pip[outside_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_pip[outside_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
                }
        }
}
void Histogram::FillHists_pip_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &pip)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0) // && _cuts->HadronsCuts(i))
        {
                short sec = _d->dc_sec(i);
                auto dt = std::make_shared<Delta_T>(_d);

                if (abs(_d->status(i)) < 4000)
                // if (dt->isCtof() == false)
                {
                        pip_Delta_vz_cut_fd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        pip_Chi2pid_cut_fd[after_all_cuts]->Fill(_d->chi2pid(i));
                        phi_vs_mom_pip_fd[after_all_cuts]->Fill(_e->pip_Phi_lab(pip), _e->pip_momentum(pip));
                        theta_pip_fd[after_all_cuts]->Fill(_e->pip_theta_lab(pip));
                        Theta_pip_lab_vs_mom_pip_fd[after_all_cuts]->Fill(_e->pip_momentum(pip), _e->pip_theta_lab(pip));
                        dcr1_sec_pip[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
                        dcr2_sec_pip[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
                        dcr3_sec_pip[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));

                        if (sec > 0 && sec <= 6)
                        {
                                Theta_fd_pip_lab_vs_mom_pip[sec - 1]->Fill(_e->pip_momentum(pip), _e->pip_theta_lab(pip));
                        }
                }
                else if (abs(_d->status(i)) > 4000)
                // else if (dt->isCtof() == true)
                {
                        pip_Delta_vz_cut_cd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)));
                        pip_Chi2pid_cut_cd[after_all_cuts]->Fill(_d->chi2pid(i));
                        phi_vs_momT_pip_cd[after_all_cuts]->Fill(_e->pip_Phi_lab(pip), _e->pip_momT(pip));
                        theta_pip_cd[after_all_cuts]->Fill(_e->pip_theta_lab(pip));
                        Theta_pip_lab_vs_mom_pip_cd[after_all_cuts]->Fill(_e->pip_momentum(pip), _e->pip_theta_lab(pip));
                }
                // std::cout << "weighr = " << _e->weight() << '\n';
        }
}

void Histogram::FillHists_pim_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        pim_Delta_vz_cut[before_any_cuts]->Fill((_d->vz(0) - _d->vz(i)));
        pim_Chi2pid_cut[before_any_cuts]->Fill(_d->chi2pid(i));
        dcr1_sec_pim[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
        dcr2_sec_pim[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
        dcr3_sec_pim[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
        // phi_vs_momT_cd_pim[before_any_cuts]->Fill(_e->pim_Phi_lab_measured(), _e->pim_momT());
        //         theta_pim->Fill(_e->pim_theta_lab_measured());
        //         Theta_pim_lab_vs_mom_pim->Fill(_e->pim_momentum_measured(), _e->pim_theta_lab_measured());
}
void Histogram::FillHists_pim_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        pim_Delta_vz_cut[after_all_cuts]->Fill((_d->vz(0) - _d->vz(i)));
        pim_Chi2pid_cut[after_all_cuts]->Fill(_d->chi2pid(i));
        dcr1_sec_pim[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i));
        dcr2_sec_pim[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i));
        dcr3_sec_pim[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i));
        // phi_vs_momT_cd_pim[after_all_cuts]->Fill(_e->pim_Phi_lab(pim), _e->pim_momT(pim));
}
void Histogram::Write_Electron_cuts()
{
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                // vz_position[c] -> Fit("gaus", "QMR+", "QMR+", -7.089, 2.0);
                // gStyle->SetOptFit(1111);
                vz_position[c]->SetXTitle("vz (cm)");
                if (vz_position[c]->GetEntries())
                        vz_position[c]->Write();
                pcal_sec[c]->SetXTitle("x (cm)");
                pcal_sec[c]->SetYTitle("y (cm)");
                pcal_sec[c]->SetOption("COLZ1");
                if (pcal_sec[c]->GetEntries())
                        pcal_sec[c]->Write();
                pcal_sec_ineff_cuts[c]->SetXTitle("x (cm)");
                pcal_sec_ineff_cuts[c]->SetYTitle("y (cm)");
                pcal_sec_ineff_cuts[c]->SetOption("COLZ1");
                if (pcal_sec_ineff_cuts[c]->GetEntries())
                        pcal_sec_ineff_cuts[c]->Write();

                pcal_hx_hy_sec[c]->SetXTitle("x (cm)");
                pcal_hx_hy_sec[c]->SetYTitle("y (cm)");
                pcal_hx_hy_sec[c]->SetOption("COLZ1");
                if (pcal_hx_hy_sec[c]->GetEntries())
                        pcal_hx_hy_sec[c]->Write();

                dcr1_sec[c]->SetXTitle("x (cm)");
                dcr1_sec[c]->SetYTitle("y (cm)");
                dcr1_sec[c]->SetOption("COLZ1");
                if (dcr1_sec[c]->GetEntries())
                        dcr1_sec[c]->Write();

                dcr2_sec[c]->SetXTitle("x (cm)");
                dcr2_sec[c]->SetYTitle("y (cm)");
                dcr2_sec[c]->SetOption("COLZ1");
                if (dcr2_sec[c]->GetEntries())
                        dcr2_sec[c]->Write();

                dcr3_sec[c]->SetXTitle("x (cm)");
                dcr3_sec[c]->SetYTitle("y (cm)");
                dcr3_sec[c]->SetOption("COLZ1");
                if (dcr3_sec[c]->GetEntries())
                        dcr3_sec[c]->Write();

                EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
                EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
                EC_sampling_fraction[c]->SetOption("COLZ1");
                if (EC_sampling_fraction[c]->GetEntries())
                        EC_sampling_fraction[c]->Write();

                ECin_sf_vs_PCAL_sf[c]->SetYTitle("ECin SF (GeV)");
                ECin_sf_vs_PCAL_sf[c]->SetXTitle("(0.2 - PCAL) SF");
                ECin_sf_vs_PCAL_sf[c]->SetOption("COLZ1");
                if (ECin_sf_vs_PCAL_sf[c]->GetEntries())
                        ECin_sf_vs_PCAL_sf[c]->Write();

                momentum[c]->SetXTitle("Momentum (GeV)");
                if (momentum[c]->GetEntries())
                        momentum[c]->Write();
        }
}

void Histogram::Write_Hadrons_cuts()
{

        auto proton_cuts_fd = RootOutputFile->mkdir("proton_cuts_fd");
        proton_cuts_fd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                prot_Delta_vz_cut_fd[c]->SetXTitle("#Delta vz (cm)");
                if (prot_Delta_vz_cut_fd[c]->GetEntries())
                        prot_Delta_vz_cut_fd[c]->Write();

                prot_Chi2pid_cut_fd[c]->SetXTitle("Chi2pid");
                if (prot_Chi2pid_cut_fd[c]->GetEntries())
                        prot_Chi2pid_cut_fd[c]->Write();
                theta_prot_fd[c]->SetXTitle("theta_prot (Deg)");
                if (theta_prot_fd[c]->GetEntries())
                        theta_prot_fd[c]->Write();

                Theta_prot_lab_vs_mom_prot_fd[c]->SetXTitle("mom_prot (GeV)");
                Theta_prot_lab_vs_mom_prot_fd[c]->SetYTitle("theta_prot (Deg)");
                Theta_prot_lab_vs_mom_prot_fd[c]->SetOption("COLZ1");
                if (Theta_prot_lab_vs_mom_prot_fd[c]->GetEntries())
                        Theta_prot_lab_vs_mom_prot_fd[c]->Write("");

                phi_vs_mom_prot_fd[c]->SetXTitle("phi (deg)");
                phi_vs_mom_prot_fd[c]->SetYTitle("Mom (GeV)");
                phi_vs_mom_prot_fd[c]->SetOption("COLZ1");
                if (phi_vs_mom_prot_fd[c]->GetEntries())
                        phi_vs_mom_prot_fd[c]->Write();

                dcr1_sec_prot[c]->SetXTitle("x (cm)");
                dcr1_sec_prot[c]->SetYTitle("y (cm)");
                dcr1_sec_prot[c]->SetOption("COLZ1");
                if (dcr1_sec_prot[c]->GetEntries())
                        dcr1_sec_prot[c]->Write();

                dcr2_sec_prot[c]->SetXTitle("x (cm)");
                dcr2_sec_prot[c]->SetYTitle("y (cm)");
                dcr2_sec_prot[c]->SetOption("COLZ1");
                if (dcr2_sec_prot[c]->GetEntries())
                        dcr2_sec_prot[c]->Write();

                dcr3_sec_prot[c]->SetXTitle("x (cm)");
                dcr3_sec_prot[c]->SetYTitle("y (cm)");
                dcr3_sec_prot[c]->SetOption("COLZ1");
                if (dcr3_sec_prot[c]->GetEntries())
                        dcr3_sec_prot[c]->Write();
        }
        auto proton_theta_vs_mom_fd_sec = RootOutputFile->mkdir("proton_theta_vs_mom_fd_sec");
        proton_theta_vs_mom_fd_sec->cd();
        for (short i = 0; i < num_sectors; i++)
        {
                Theta_fd_prot_lab_vs_mom_prot[i]->SetXTitle("mom_prot (GeV)");
                Theta_fd_prot_lab_vs_mom_prot[i]->SetYTitle("theta_prot (Deg)");
                Theta_fd_prot_lab_vs_mom_prot[i]->SetOption("COLZ1");
                if (Theta_fd_prot_lab_vs_mom_prot[i]->GetEntries())
                        Theta_fd_prot_lab_vs_mom_prot[i]->Write("");
        }
        auto proton_cuts_cd = RootOutputFile->mkdir("proton_cuts_cd");
        proton_cuts_cd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                prot_Delta_vz_cut_cd[c]->SetXTitle("#Delta vz (cm)");
                if (prot_Delta_vz_cut_cd[c]->GetEntries())
                        prot_Delta_vz_cut_cd[c]->Write();

                prot_Chi2pid_cut_cd[c]->SetXTitle("Chi2pid");
                if (prot_Chi2pid_cut_cd[c]->GetEntries())
                        prot_Chi2pid_cut_cd[c]->Write();
                theta_prot_cd[c]->SetXTitle("theta_prot (Deg)");
                if (theta_prot_cd[c]->GetEntries())
                        theta_prot_cd[c]->Write();

                Theta_prot_lab_vs_mom_prot_cd[c]->SetXTitle("mom_prot (GeV)");
                Theta_prot_lab_vs_mom_prot_cd[c]->SetYTitle("theta_prot (Deg)");
                Theta_prot_lab_vs_mom_prot_cd[c]->SetOption("COLZ1");
                if (Theta_prot_lab_vs_mom_prot_cd[c]->GetEntries())
                        Theta_prot_lab_vs_mom_prot_cd[c]->Write("");

                phi_vs_momT_prot_cd[c]->SetXTitle("phi (deg)");
                phi_vs_momT_prot_cd[c]->SetYTitle("Mom (GeV)");
                phi_vs_momT_prot_cd[c]->SetOption("COLZ1");
                if (phi_vs_momT_prot_cd[c]->GetEntries())
                        phi_vs_momT_prot_cd[c]->Write();
        }

        auto pip_cuts_fd = RootOutputFile->mkdir("pip_cuts_fd");
        pip_cuts_fd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pip_Delta_vz_cut_fd[c]->SetXTitle("#Delta vz (cm)");
                if (pip_Delta_vz_cut_fd[c]->GetEntries())
                        pip_Delta_vz_cut_fd[c]->Write();

                pip_Chi2pid_cut_fd[c]->SetXTitle("Chi2pid");
                if (pip_Chi2pid_cut_fd[c]->GetEntries())
                        pip_Chi2pid_cut_fd[c]->Write();

                phi_vs_mom_pip_fd[c]->SetXTitle("phi (deg)");
                phi_vs_mom_pip_fd[c]->SetYTitle("Mom (GeV)");
                phi_vs_mom_pip_fd[c]->SetOption("COLZ1");
                if (phi_vs_mom_pip_fd[c]->GetEntries())
                        phi_vs_mom_pip_fd[c]->Write();
                theta_pip_fd[c]->SetXTitle("theta_prot (Deg)");
                if (theta_pip_fd[c]->GetEntries())
                        theta_pip_fd[c]->Write();

                Theta_pip_lab_vs_mom_pip_fd[c]->SetXTitle("mom_prot (GeV)");
                Theta_pip_lab_vs_mom_pip_fd[c]->SetYTitle("theta_prot (Deg)");
                Theta_pip_lab_vs_mom_pip_fd[c]->SetOption("COLZ1");
                if (Theta_pip_lab_vs_mom_pip_fd[c]->GetEntries())
                        Theta_pip_lab_vs_mom_pip_fd[c]->Write("");

                dcr1_sec_pip[c]->SetXTitle("x (cm)");
                dcr1_sec_pip[c]->SetYTitle("y (cm)");
                dcr1_sec_pip[c]->SetOption("COLZ1");
                if (dcr1_sec_pip[c]->GetEntries())
                        dcr1_sec_pip[c]->Write();

                dcr2_sec_pip[c]->SetXTitle("x (cm)");
                dcr2_sec_pip[c]->SetYTitle("y (cm)");
                dcr2_sec_pip[c]->SetOption("COLZ1");
                if (dcr2_sec_pip[c]->GetEntries())
                        dcr2_sec_pip[c]->Write();

                dcr3_sec_pip[c]->SetXTitle("x (cm)");
                dcr3_sec_pip[c]->SetYTitle("y (cm)");
                dcr3_sec_pip[c]->SetOption("COLZ1");
                if (dcr3_sec_pip[c]->GetEntries())
                        dcr3_sec_pip[c]->Write();
        }
        auto pip_th_vs_mom_fd_sec = RootOutputFile->mkdir("pip_th_vs_mom_fd_sec");
        pip_th_vs_mom_fd_sec->cd();
        for (short i = 0; i < num_sectors; i++)
        {
                Theta_fd_pip_lab_vs_mom_pip[i]->SetXTitle("mom_prot (GeV)");
                Theta_fd_pip_lab_vs_mom_pip[i]->SetYTitle("theta_prot (Deg)");
                Theta_fd_pip_lab_vs_mom_pip[i]->SetOption("COLZ1");
                if (Theta_fd_pip_lab_vs_mom_pip[i]->GetEntries())
                        Theta_fd_pip_lab_vs_mom_pip[i]->Write("");
        }
        auto pip_cuts_cd = RootOutputFile->mkdir("pip_cuts_cd");
        pip_cuts_cd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pip_Delta_vz_cut_cd[c]->SetXTitle("#Delta vz (cm)");
                if (pip_Delta_vz_cut_cd[c]->GetEntries())
                        pip_Delta_vz_cut_cd[c]->Write();

                pip_Chi2pid_cut_cd[c]->SetXTitle("Chi2pid");
                if (pip_Chi2pid_cut_cd[c]->GetEntries())
                        pip_Chi2pid_cut_cd[c]->Write();

                phi_vs_momT_pip_cd[c]->SetXTitle("phi (deg)");
                phi_vs_momT_pip_cd[c]->SetYTitle("Mom (GeV)");
                phi_vs_momT_pip_cd[c]->SetOption("COLZ1");
                if (phi_vs_momT_pip_cd[c]->GetEntries())
                        phi_vs_momT_pip_cd[c]->Write();
                theta_pip_cd[c]->SetXTitle("theta_prot (Deg)");
                if (theta_pip_cd[c]->GetEntries())
                        theta_pip_cd[c]->Write();

                Theta_pip_lab_vs_mom_pip_cd[c]->SetXTitle("mom_prot (GeV)");
                Theta_pip_lab_vs_mom_pip_cd[c]->SetYTitle("theta_prot (Deg)");
                Theta_pip_lab_vs_mom_pip_cd[c]->SetOption("COLZ1");
                if (Theta_pip_lab_vs_mom_pip_cd[c]->GetEntries())
                        Theta_pip_lab_vs_mom_pip_cd[c]->Write("");
        }

        auto pim_cuts = RootOutputFile->mkdir("pim_cuts");
        pim_cuts->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pim_Delta_vz_cut[c]->SetXTitle("#Delta vz (cm)");
                if (pim_Delta_vz_cut[c]->GetEntries())
                        pim_Delta_vz_cut[c]->Write();

                pim_Chi2pid_cut[c]->SetXTitle("Chi2pid");
                if (pim_Chi2pid_cut[c]->GetEntries())
                        pim_Chi2pid_cut[c]->Write();

                // phi_vs_momT_cd_pim[c]->SetXTitle("phi (deg)");
                // phi_vs_momT_cd_pim[c]->SetYTitle("Mom (GeV)");
                // phi_vs_momT_cd_pim[c]->SetOption("COLZ1");
                // if (phi_vs_momT_cd_pim[c]->GetEntries())
                //         phi_vs_momT_cd_pim[c]->Write();

                dcr1_sec_pim[c]->SetXTitle("x (cm)");
                dcr1_sec_pim[c]->SetYTitle("y (cm)");
                dcr1_sec_pim[c]->SetOption("COLZ1");
                if (dcr1_sec_pim[c]->GetEntries())
                        dcr1_sec_pim[c]->Write();

                dcr2_sec_pim[c]->SetXTitle("x (cm)");
                dcr2_sec_pim[c]->SetYTitle("y (cm)");
                dcr2_sec_pim[c]->SetOption("COLZ1");
                if (dcr2_sec_pim[c]->GetEntries())
                        dcr2_sec_pim[c]->Write();

                dcr3_sec_pim[c]->SetXTitle("x (cm)");
                dcr3_sec_pim[c]->SetYTitle("y (cm)");
                dcr3_sec_pim[c]->SetOption("COLZ1");
                if (dcr3_sec_pim[c]->GetEntries())
                        dcr3_sec_pim[c]->Write();
        }
}

void Histogram::makeHists_sector()
{

        for (short i = 0; i < 3; i++)
        {
                W_det[i] = std::make_shared<TH1D>(Form("W_det_%d", i + 1),
                                                  Form("W detector: %d", i + 1), bins, w_min,
                                                  w_max);
                if (i == 0)
                        WQ2_det[i] = std::make_shared<TH2D>(
                            Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1),
                            bins, w_min, w_max, bins, zero, 0.5);
                else
                        WQ2_det[i] = std::make_shared<TH2D>(
                            Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1),
                            bins, w_min, w_max, bins, zero, q2_max);
        }

        for (short i = 0; i < num_sectors; i++)
        {

                W_vs_q2_sec[i] = std::make_shared<TH2D>(
                    Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                    w_min, w_max, bins, zero, q2_max);

                W_sec[i] =
                    std::make_shared<TH1D>(Form("w_sec_%d", i + 1),
                                           Form("W Sector: %d", i + 1), bins, w_min, w_max);

                W_vs_q2_twoPi_sec[i] =
                    std::make_shared<TH2D>(Form("wvsq2_sec_twoPi_%d", i + 1),
                                           Form("W vs Q^{2} W_twoPi Sector: %d", i + 1),
                                           bins, w_min, w_max, bins, zero, q2_max);

                MM_mPim_twoPi_sec[i] = std::make_shared<TH1D>(
                    Form("MM_missingPim_sec_%d", i + 1), Form("MM missingPim Sector: %d", i + 1), bins,
                    -0.4, 0.4);
                MM2_mPim_twoPi_sec[i] = std::make_shared<TH1D>(
                    Form("MM_SQ_missingPim_sec_%d", i + 1), Form("MM_SQ missingPim Sector: %d", i + 1),
                    bins, -0.4, 0.4);
                MM2_twoPi_missingPip_sec[i] = std::make_shared<TH1D>(
                    Form("MM_SQ_missingPip_sec_%d", i + 1),
                    Form("MM_SQ missingPiP Sector: %d", i + 1), bins, -0.4, 0.4);
                MM2_twoPi_missingProt_sec[i] = std::make_shared<TH1D>(
                    Form("MM_SQt_missingProt_sec_%d", i + 1),
                    Form("MM_SQ missingProt Sector: %d", i + 1), bins, 0.6, 1.2);

                Theta_fd_prot_lab_vs_mom_prot[i] = std::make_shared<TH2D>(
                    Form("Theta_fd_prot_lab_vs_mom_prot_sec_%d", i + 1), Form("Theta_fd_prot_lab_vs_mom_prot_%d", i + 1), bins,
                    0, 10.0, bins, 0, 50);

                Theta_fd_pip_lab_vs_mom_pip[i] = std::make_shared<TH2D>(
                    Form("Theta_fd_pip_lab_vs_mom_pip_sec_%d", i + 1), Form("Theta_fd_pip_lab_vs_mom_pip_%d", i + 1), bins,
                    0, 5.0, bins, 0, 50);
                Theta_fd_elec_lab_vs_mom_elec[i] = std::make_shared<TH2D>(
                    Form("Theta_fd_elec_lab_vs_mom_elec_sec_%d", i + 1), Form("Theta_fd_elec_lab_vs_mom_elec_%d", i + 1), bins,
                    3, 10.0, bins, 0, 30);

                for (auto &&cut : before_after_cut)
                {
                        int c = cut.first;
                        auto type = cut.second.c_str();
                        // std::cout << " cut number = " << c << "  cut type is : " << type << std::endl;

                        htcc_nphe_sec[c][i] = std::make_shared<TH1D>(Form("htcc_nphe_sec%d%s", i + 1, type),
                                                                     Form("htcc nphe sec: %d%s", i + 1, type),
                                                                     200, 0, 70);
                        elec_Chi2pid_sec[c][i] = std::make_shared<TH1D>(Form("elec_Chi2pid_sec%d%s", i + 1, type),
                                                                        Form("elec_Chi2pid sec: %d%s", i + 1, type),
                                                                        200, -10, 5);
                        vz_sec[c][i] = std::make_shared<TH1D>(Form("vz_sec%d%s", i + 1, type),
                                                              Form("vz sec: %d%s", i + 1, type),
                                                              200, -20, 20);

                        pcal_lu_sec[c][i] = std::make_shared<TH1D>(Form("pcal_lu_sec%d%s", i + 1, type),
                                                                   Form("pcal_lu_sec sec: %d%s", i + 1, type),
                                                                   200, 0, 250);

                        pcal_lv_sec[c][i] = std::make_shared<TH1D>(Form("pcal_lv_sec%d%s", i + 1, type),
                                                                   Form("pcal_lv_sec sec: %d%s", i + 1, type),
                                                                   200, 0, 250);

                        pcal_lw_sec[c][i] = std::make_shared<TH1D>(Form("pcal_lw_sec%d%s", i + 1, type),
                                                                   Form("pcal_lw_sec sec: %d%s", i + 1, type),
                                                                   200, 0, 250);

                        // for (short j = 0; j < 9; j++) // mom bins
                        // {
                        //         auto mom_bin = ECIN_ECOUT_MOM_NAME[j].c_str();

                        //         ECAL_VS_PCAL[c][i][j] = std::make_shared<TH2D>(
                        //             Form("ECAL_VS_PCAL_%d%s%s", i + 1, type, mom_bin), Form("ECAL_VS_PCAL_%d%s%s", i + 1, type, mom_bin), bins,
                        //             0, 0.2, bins, 0, 0.25);
                        // }
                        SF_VS_MOM[c][i] = std::make_shared<TH2D>(
                            Form("SF_VS_MOM_%d%s", i + 1, type), Form("SF_VS_MOM_%d%s", i + 1, type), bins,
                            0, 10.0, bins, 0, 0.5);
                }
        }
}

void Histogram::makeHists_deltat()
{
        for (int i = 0; i < FDmomArraySize; ++i)
        {
                float dt_lim = 0.0;
                if (i < 2)
                {
                        dt_lim = 2.0;
                }
                else if (i < 10)
                {
                        dt_lim = 1.5;
                }
                else if (i < 20)
                {
                        dt_lim = 1.0;
                }
                else
                {
                        dt_lim = 0.5;
                }
                auto name_dt_prot_fd = Form("dt_prot_fd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));
                auto name_dt_pip_fd = Form("dt_pip_fd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));
                auto name_dt_pim_fd = Form("dt_pim_fd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));

                dt_prot_fd_hist[i] = std::make_shared<TH1D>(name_dt_prot_fd, name_dt_prot_fd, 100, -dt_lim, dt_lim);
                dt_pip_fd_hist[i] = std::make_shared<TH1D>(name_dt_pip_fd, name_dt_pip_fd, 100, -dt_lim, dt_lim);
                dt_pim_fd_hist[i] = std::make_shared<TH1D>(name_dt_pim_fd, name_dt_pim_fd, 100, -dt_lim, dt_lim);

                // std::cout << " name is " << name_dt_prot_fd << " " << i << ": " << FDmomArray[i] << std::endl;
        }
        for (int i = 0; i < CDmomArraySize; ++i)
        {
                float dt_lim = 1.0;
                float offset = 0.0;

                // if (i < 3)
                // {
                //         dt_lim = 1.0;
                //         offset = 0.0;
                // }
                // else if (i < 5)
                // {
                //         dt_lim = 0.75;
                //         offset = 0.25;
                // }
                // else if (i < 20)
                // {
                //         dt_lim = 1.0;
                // }
                // else
                // {
                //         dt_lim = 0.5;
                // }
                auto name_dt_prot_cd = Form("dt_prot_cd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));
                auto name_dt_pip_cd = Form("dt_pip_cd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));
                auto name_dt_pim_cd = Form("dt_pim_cd_%.2f<=Mom<%.2f GeV", (0.2 + 0.25 * i), (0.2 + 0.25 * i + 0.25));

                dt_prot_cd_hist[i] = std::make_shared<TH1D>(name_dt_prot_cd, name_dt_prot_cd, 100, -dt_lim, dt_lim);
                dt_pip_cd_hist[i] = std::make_shared<TH1D>(name_dt_pip_cd, name_dt_pip_cd, 100, -dt_lim, dt_lim);
                dt_pim_cd_hist[i] = std::make_shared<TH1D>(name_dt_pim_cd, name_dt_pim_cd, 100, -dt_lim, dt_lim);

                // std::cout << " name is " << name_dt_prot_cd << " " << i << ": " << CDmomArray[i] << std::endl;
        }

        {
                std::string tof = "";
                static const short p_num = 3; // 0-P 1-Pip 2-Pim
                std::string p_name[p_num] = {"P", "pip", "Pim"};
                static const short cut_num = 2; // 0-no cuts, 1- with cuts
                std::string cut_name[cut_num] = {"before cut", "after cut"};

                for (short p = 0; p < p_num; p++)
                {
                        for (short c = 0; c < cut_num; c++)
                        {

                                tof = "both";
                                delta_t_hist[p][c][0] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, Dt_min, Dt_max);

                                tof = "ftof";
                                delta_t_hist[p][c][1] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, Dt_min, Dt_max);
                                tof = "ctof";
                                delta_t_hist[p][c][2] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, -6, 6);
                        }
                }
        }
        // }
        // delta_t_hist_pip[0][0] = std::make_shared<TH2D>("pip_both_bc ", "#Dt pip both bc ", bins, p_min, p_max, bins, Dt_min, Dt_max);
        // delta_t_hist_pip[1][0] = std::make_shared<TH2D>("pip_both_ac ", "#Dt pip both ac", bins, p_min, p_max, bins, Dt_min, Dt_max);
        // delta_t_hist_pip[0][1] = std::make_shared<TH2D>("pip_FD_bc", "#Dt pip FD bc ", bins, p_min, p_max, bins, Dt_min, Dt_max);
        // delta_t_hist_pip[1][1] = std::make_shared<TH2D>("pip_FD_ac", "#Dt pip FD ac ", bins, p_min, p_max, bins, Dt_min, Dt_max);
        // delta_t_hist_pip[0][2] = std::make_shared<TH2D>("pip_CD_bc", "#Dt pip CD bc ", bins, p_min, p_max, bins, Dt_min, Dt_max);
        // delta_t_hist_pip[1][2] = std::make_shared<TH2D>("pip_CD_ac", "#Dt pip CD ac ", bins, p_min, p_max, bins, Dt_min, Dt_max);
}

void Histogram::Fill_deltat_before_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                       const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {
                // auto _cuts = std::make_unique<Cuts>(data, dt);
                int charge = data->charge(part);
                // bool fc = dt->ctof();
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time[] = {NAN, NAN, NAN};
                float time_cd[] = {NAN, NAN, NAN};
                float time_fd[] = {NAN, NAN, NAN};
                // std::cout<<" charge before =  "<< charge<<std::endl;

                // if (_e->TwoPion_exclusive())// || _e->TwoPion_missingPim() || _e->TwoPion_missingPip() || _e->TwoPion_missingProt())

                // if (_e->TwoPion_exclusive())
                //         {
                if (fd_part)
                {

                        if (charge == 1)
                        {
                                if (pid == PIP)
                                {
                                        time[1] = dt->dt_Pi();
                                        time_fd[1] = dt->dt_Pi();
                                        delta_t_hist[1][0][1]->Fill(mom, time_fd[1]);
                                }
                                // std::cout << " dt_Pi() .... =  " << dt->dt_Pi() << std::endl;
                                // std::cout << " dt_Pi(int i) .... =  " << dt->dt_Pi(pid) << std::endl;

                                else if (pid == PROTON)
                                {

                                        time[0] = dt->dt_P();
                                        time_fd[0] = dt->dt_P();
                                        delta_t_hist[0][0][1]->Fill(mom, time_fd[0]);
                                }
                        }
                        else if (charge == -1)
                        {
                                if (pid == PIM)
                                // if (pid != ELECTRON)
                                {

                                        time[2] = dt->dt_Pi();
                                        time_fd[2] = dt->dt_Pi();
                                        delta_t_hist[2][0][1]->Fill(mom, time_fd[2]);
                                }
                        }

                        if (mom < 7.7 && mom >= 0.2)
                        {
                                // std::cout << mom << " int mom =  " << int((mom - 0.2) / 0.25) << std::endl;

                                dt_prot_fd_hist[int((mom - 0.2) / 0.25)]->Fill(time_fd[0]);
                                dt_pip_fd_hist[int((mom - 0.2) / 0.25)]->Fill(time_fd[1]);
                                dt_pim_fd_hist[int((mom - 0.2) / 0.25)]->Fill(time_fd[2]);
                        }
                }
                else if (cd_part)
                {
                        if (charge == 1)
                        {
                                // std::cout << " charge after cd =  " << charge << std::endl;
                                if (pid == PROTON)
                                {

                                        time[0] = dt->dt_P();
                                        time_cd[0] = dt->dt_P();
                                        delta_t_hist[0][0][2]->Fill(mom, time_cd[0]);
                                }

                                else if (pid == PIP)
                                {
                                        time[1] = dt->dt_Pi();
                                        time_cd[1] = dt->dt_Pi();
                                        delta_t_hist[1][0][2]->Fill(mom, time_cd[1]);
                                }
                        }
                        else if (charge == -1)
                        {
                                if (pid == PIM)
                                // if (pid != ELECTRON)
                                {
                                        time[2] = dt->dt_Pi();
                                        time_cd[2] = dt->dt_Pi();
                                        delta_t_hist[2][0][2]->Fill(mom, time_cd[2]);
                                }
                        }

                        if (mom < 3.95 && mom >= 0.2)
                        {
                                // std::cout << mom << " int mom =  " << int((mom - 0.2) / 0.25) << std::endl;

                                dt_prot_cd_hist[int((mom - 0.2) / 0.25)]->Fill(time_cd[0]);
                                dt_pip_cd_hist[int((mom - 0.2) / 0.25)]->Fill(time_cd[1]);
                                dt_pim_cd_hist[int((mom - 0.2) / 0.25)]->Fill(time_cd[2]);
                        }
                }

                delta_t_hist[0][0][0]->Fill(mom, time[0]);
                delta_t_hist[1][0][0]->Fill(mom, time[1]);
                delta_t_hist[2][0][0]->Fill(mom, time[2]);
        }
}
void Histogram::Fill_deltat_prot_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                           const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {
                // auto _cuts = std::make_unique<Cuts>(data, dt);
                int charge = data->charge(part);
                // bool fc = dt->ctof();
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (_e->TwoPion_exclusive())
                // {

                // if (pid == PROTON)
                {

                        if (charge == 1)

                        {
                                if (fd_part)
                                {
                                        time = dt->dt_P();
                                        time_fd = dt->dt_P();
                                        delta_t_hist[0][1][1]->Fill(mom, time_fd);
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_P();
                                        time_cd = dt->dt_P();
                                        delta_t_hist[0][1][2]->Fill(mom, time_cd);
                                }

                                delta_t_hist[0][1][0]->Fill(mom, time);
                        }
                }
        }
}
void Histogram::Fill_deltat_pip_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                          const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {
                // auto _cuts = std::make_unique<Cuts>(data, dt);
                int charge = data->charge(part);
                // bool fc = dt->ctof();
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (_e->TwoPion_exclusive() || _e->TwoPion_missingPim() || _e->TwoPion_missingPip() || _e->TwoPion_missingProt())

                //         if (_e->TwoPion_exclusive())
                // {
                // if (pid == PIP)
                {

                        if (charge == 1)
                        {
                                if (fd_part)
                                {
                                        time = dt->dt_Pi();
                                        time_fd = dt->dt_Pi();
                                        delta_t_hist[1][1][1]->Fill(mom, time_fd);
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_Pi();
                                        time_cd = dt->dt_Pi();
                                        delta_t_hist[1][1][2]->Fill(mom, time_cd);
                                }

                                delta_t_hist[1][1][0]->Fill(mom, time);
                        }
                }
        }
}
void Histogram::Fill_deltat_pim_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                          const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.35 && _e->W() <= 2.15 && _e->Q2() > 1.95 && _e->Q2() <= 9.0)
        {
                // auto _cuts = std::make_unique<Cuts>(data, dt);
                int charge = data->charge(part);
                // bool fc = dt->ctof();
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (pid == PIM)
                // if (_e->TwoPion_exclusive() || _e->TwoPion_missingPim() || _e->TwoPion_missingPip() || _e->TwoPion_missingProt())
                {

                        if (charge == -1)
                        {
                                if (fd_part)
                                {
                                        time = dt->dt_Pi();
                                        time_fd = dt->dt_Pi();
                                        delta_t_hist[2][1][1]->Fill(mom, time_fd);
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_Pi();
                                        time_cd = dt->dt_Pi();
                                        delta_t_hist[2][1][2]->Fill(mom, time_cd);
                                }

                                delta_t_hist[2][1][0]->Fill(mom, time);
                        }
                }
        }
}
void Histogram::Write_deltat()
{

        for (short i = 0; i < 3; i++)
        {
                for (short j = 0; j < 2; j++)
                {
                        for (short k = 0; k < 3; k++)
                        {
                                delta_t_hist[i][j][k]->SetXTitle("Momentum (GeV)");
                                delta_t_hist[i][j][k]->SetYTitle("#Deltat");
                                delta_t_hist[i][j][k]->SetOption("COLZ1");
                                // if (delta_t_hist[i][j][k]->GetEntries() > 1)
                                delta_t_hist[i][j][k]->Write();
                        }
                }
        }

        auto dt_hist_prot_fd = RootOutputFile->mkdir("dt_hist_prot_fd");
        dt_hist_prot_fd->cd();
        for (int i = 0; i < FDmomArraySize; ++i)
        {
                dt_prot_fd_hist[i]->Write();
        }
        auto dt_hist_pip_fd = RootOutputFile->mkdir("dt_hist_pip_fd");
        dt_hist_pip_fd->cd();
        for (int i = 0; i < FDmomArraySize; ++i)
        {
                dt_pip_fd_hist[i]->Write();
        }

        auto dt_hist_pim_fd = RootOutputFile->mkdir("dt_hist_pim_fd");
        dt_hist_pim_fd->cd();
        for (int i = 0; i < FDmomArraySize; ++i)
        {
                dt_pim_fd_hist[i]->Write();
        }

        //////////////////////////////////////////////////////

        auto dt_hist_prot_cd = RootOutputFile->mkdir("dt_hist_prot_cd");
        dt_hist_prot_cd->cd();
        for (int i = 0; i < CDmomArraySize; ++i)
        {
                dt_prot_cd_hist[i]->Write();
        }
        auto dt_hist_pip_cd = RootOutputFile->mkdir("dt_hist_pip_cd");
        dt_hist_pip_cd->cd();
        for (int i = 0; i < CDmomArraySize; ++i)
        {
                dt_pip_cd_hist[i]->Write();
        }

        auto dt_hist_pim_cd = RootOutputFile->mkdir("dt_hist_pim_cd");
        dt_hist_pim_cd->cd();
        for (int i = 0; i < CDmomArraySize; ++i)
        {
                dt_pim_cd_hist[i]->Write();
        }
}
void Histogram::makeHists_MomVsBeta()
{
        for (short p = 0; p < particle_num; p++)
        {
                for (short c = 0; c < charge_num; c++)
                {
                        for (short i = 0; i < with_id_num; i++)
                        {
                                momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
                                    Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
                                    Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
                                    bins, p_min, p_max, bins, zero, 1.2);
                        }
                }
        }
}
void Histogram::Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part, const std::shared_ptr<Reaction> &_e)
{
        int good_ID = 0;
        float beta = data->beta(part);
        float mom = data->p(part);
        int charge = data->charge(part);
        int pid = data->pid(part);
        if (beta != 0)
        {
                for (short p = 0; p < particle_num; p++)
                {
                        switch (p)
                        {
                        case 0:
                                good_ID = ELECTRON;
                                break;
                        case 1:
                                good_ID = PIP;
                                break;
                        case 2:
                                good_ID = PROTON;
                                break;
                        case 3:
                                good_ID = KP;
                                break;
                        }
                        if (charge == -1)
                        {
                                momvsbeta_hist[p][1][0]->Fill(mom, beta);
                                if (good_ID == 11)
                                        momvsbeta_hist[0][1][1]->Fill(mom, beta);
                                else if (-good_ID == pid)
                                        momvsbeta_hist[p][1][1]->Fill(mom, beta);
                        }
                        else if (charge == 1)
                        {
                                momvsbeta_hist[p][0][0]->Fill(mom, beta);
                                if (good_ID == pid)
                                {
                                        momvsbeta_hist[p][0][1]->Fill(mom, beta);
                                }
                        }
                }
        }
}

void Histogram::Write_MomVsBeta()
{

        for (short p = 0; p < particle_num; p++)
        {
                for (short c = 0; c < charge_num; c++)
                {
                        for (short i = 0; i < with_id_num; i++)
                        {
                                momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
                                momvsbeta_hist[p][c][i]->SetYTitle("#beta");
                                momvsbeta_hist[p][c][i]->SetOption("COLZ1");
                                momvsbeta_hist[p][c][i]->Write();
                        }
                }
        }
}

void Histogram::Fill_deltaP_prot(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_prot_hist->Fill(dp);
}
void Histogram::Fill_deltaP_pip(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_pip_hist->Fill(dp);
}

void Histogram::Fill_deltaP_pip_for_prot(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_pip_for_prot_hist->Fill(dp);
}
void Histogram::Fill_deltaP_prot_for_pip(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_prot_for_pip_hist->Fill(dp);
}

void Histogram::Fill_deltaP_sum(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_sum_hist->Fill(dp);
}
void Histogram::Fill_deltaP_sum_twoPi(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_sum_hist_twoPi->Fill(dp);
}

void Histogram::Fill_deltaP_ambi_prot(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_ambi_prot_all_hist->Fill(dp);
}
void Histogram::Fill_deltaP_ambi_pip(const std::shared_ptr<Reaction> &_e, double dp)
{
        dp_ambi_pip_all_hist->Fill(dp);
}

void Histogram::Fill_Entries(int num_entries)
{
        entries_in_each_event->Fill(num_entries);
}
void Histogram::Fill_Entries_prot(int num_entries)
{
        entries_prot->Fill(num_entries);
}
void Histogram::Fill_Entries_pip(int num_entries)
{
        entries_pip->Fill(num_entries);
}

void Histogram::Fill_all_Combi(const std::shared_ptr<Reaction> &_e)
{
        // MM2_mPim_all_comb->Fill(_e->MM2_mPim(), _e->weight());
        MM2_mPim_all_comb->Fill(_e->MM2_mPim());
}
void Histogram::Fill_1_Combi(const std::shared_ptr<Reaction> &_e)
{
        MM2_mPim_1_comb->Fill(_e->MM2_mPim(), _e->weight());
}
void Histogram::Fill_2_Combi(const std::shared_ptr<Reaction> &_e)
{
        // MM2_mPim_2_comb->Fill(_e->MM2_mPim(), _e->weight());
        MM2_mPim_2_comb->Fill(_e->MM2_mPim_swapped(), _e->weight());
}
// void Histogram::Fill_3_Combi(const std::shared_ptr<Reaction> &_e)
// {
//         MM2_mPim_3_comb->Fill(_e->MM2_mPim() - _e->MM2_mPim_swapped(), _e->weight());
// }

void Histogram::Fill_3_Combi(float dv2, const std::shared_ptr<Reaction> &_e)
{
        MM2_mPim_3_comb->Fill(dv2, _e->weight());
}
void Histogram::Fill_4_or_more_Combi(float dv2, const std::shared_ptr<Reaction> &_e)
{
        // MM2_mPim_4_or_more_comb->Fill(_e->MM2_mPim(), _e->weight());
        MM2_mPim_4_or_more_comb->Fill(dv2, _e->weight());
}

void Histogram::Write_deltaP()
{
        dp_prot_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_prot_hist->Write();

        dp_pip_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_pip_hist->Write();

        dp_prot_for_pip_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_prot_for_pip_hist->Write();

        dp_pip_for_prot_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_pip_for_prot_hist->Write();

        dp_sum_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_sum_hist->Write();

        dp_sum_hist_twoPi->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_sum_hist_twoPi->Write();

        dp_ambi_prot_all_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_ambi_prot_all_hist->Write();

        dp_ambi_pip_all_hist->SetXTitle(" (Gen - Rec ) Mom (GeV)");
        dp_ambi_pip_all_hist->Write();

        entries_in_each_event->Write();
        entries_prot->Write();
        entries_pip->Write();

        MM2_mPim_all_comb->SetXTitle("MMSQ (GeV^2)");
        MM2_mPim_all_comb->Write();
        MM2_mPim_1_comb->SetXTitle("MMSQ (GeV^2)");
        MM2_mPim_1_comb->Write();
        MM2_mPim_2_comb->SetXTitle("MMSQ (GeV^2)");
        MM2_mPim_2_comb->Write();
        MM2_mPim_3_comb->SetXTitle("MMSQ (GeV^2)");
        MM2_mPim_3_comb->Write();
        MM2_mPim_4_or_more_comb->SetXTitle("MMSQ (GeV^2)");
        MM2_mPim_4_or_more_comb->Write();
};
