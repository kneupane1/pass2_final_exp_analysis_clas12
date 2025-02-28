
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "QADB.h"

template <class CutType>
//size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id)
 size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, const std::shared_ptr<QA::QADB> &_qa = nullptr, int thread_id=1)
{        if (_mc)
        {
                // Monte Carlo data processing (no QA)
                std::cout << "Processing without QA (MC)" << std::endl;
        }
        else
        {
                // Real data processing with QA checks
                std::cout << "Processing with QA" << std::endl;
        }

        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();

        float beam_energy = 10.6041;

        if (std::is_same<CutType, Pass2_Cuts>::value)
        {
                if (thread_id == 0)
                        std::cout << BLUE << "Using Pass2 RGA Cuts" << DEF << std::endl;

                beam_energy = 10.6041;
        }

        if (getenv("CLAS12_E") != NULL)
                beam_energy = atof(getenv("CLAS12_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // // // Make a data object which all the branches can be accessed from
        // auto data = std::make_shared<Branches12>(_chain, true);
        auto data = std::make_shared<Branches12>(_chain, _mc);

        // if (!_qa->Golden(data->getRun(), data->getEvent()))
        //         continue;

        // Total number of events "Processed"
        size_t total = 0;
        float no_of_events = 0, miss_prot = 0, miss_pim = 0, miss_pip = 0, other = 0, excl_events = 0,
              twopi = 0;
        int two_pion_mPim_events = 0, two_pion_Excl_events = 0;
        int prot = 0, pip = 0, pim = 0, elec = 0, failed_prot = 0, pid_zero_elec = 0;
        int no_prot_pip = 0, no_prot_pip_mc = 0;
        int Pip_pid_mc = -9999;
        int Prot_pid_mc = -9999;
        int Pim_pid_mc = -9999;
        int Pip_pid_rec = -9999;
        int Prot_pid_rec = -9999;
        int Pim_pid_rec = -9999;
        double dp_prot1 = NAN;
        double dp_pip1 = NAN;
        double dp_prot2 = NAN;
        double dp_pip2 = NAN;
        double dp_prot3 = NAN;
        double dp_pip3 = NAN;
        double dp_pim1 = NAN;
        double dp_pim2 = NAN;
        double dp_pim3 = NAN;
        // For each event

        int first_entries = 0;
        int second_entries = 0;
        int third_entries = 0;
        int four_or_more_entries = 0;
        int events_with_non_zero_wt = 0, events_with_zero_wt = 0, events_passes_w_q2_cuts = 0;

        for (size_t current_event = 0; current_event < num_of_events; current_event++)
        {
                prot = 0;
                pip = 0;
                int entries_in_this_event = 0;

                // Get current event
                _chain->GetEntry(current_event);
                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

                ///////// This part is only for mc events /////////////
                ///////// This part is only for mc events /////////////
                ///////// This part is only for mc events /////////////
                ///////// This part is only for mc events /////////////
                ///////// This part is only for mc events /////////////
/*                if (_mc)
                {
                        if (data->mc_weight() <= 0)
                                events_with_zero_wt++;
                        if (data->mc_npart() < 1 || data->mc_weight() <= 0)
                                continue;
                        events_with_non_zero_wt++;

                        // Make a reaction class from the data given
                        auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

                        for (int part = 0; part < data->mc_npart(); part++)
                        {
                                // Check particle ID's and fill the reaction class
                                if (data->mc_pid(part) == PIP)
                                {
                                        mc_event->SetMCPip(part);
                                }
                                if (data->mc_pid(part) == PROTON)
                                {
                                        mc_event->SetMCProton(part);
                                }
                                if (data->mc_pid(part) == PIM)
                                {
                                        mc_event->SetMCPim(part);
                                        // } else {
                                        //   mc_event->SetMCOther(part);
                                }
                        }

                        // if (mc_event->W_mc() < 3.0 && mc_event->W_mc() > 1.0 && mc_event->Q2_mc() < 12.0 && mc_event->Q2_mc() > 1.0)
                        if (mc_event->W_mc() <= 2.15 && mc_event->W_mc() > 1.35 && mc_event->Q2_mc() < 9.0 && mc_event->Q2_mc() > 1.95)
                        {

                                // // Retrieve the number of protons and pions in the event
                                size_t num_protons_mc = mc_event->GetMcProtons().size();
                                size_t num_pips_mc = mc_event->GetMcPips().size();
                                size_t num_pims_mc = mc_event->GetMcPims().size();

                                for (size_t i = 0; i < num_protons_mc; ++i)
                                {
                                        for (size_t j = 0; j < num_pips_mc; ++j)
                                        {

                                                if (mc_event->GetProtonMcIndices()[i] == mc_event->GetPipMcIndices()[j])
                                                        no_prot_pip_mc++;

                                                // Exclude the case where the same particle is assigned as both proton and pip
                                                if (mc_event->GetProtonMcIndices()[i] != mc_event->GetPipMcIndices()[j])
                                                {

                                                        for (size_t k = 0; k < num_pims_mc; ++k)
                                                        {

                                                                mc_event->boost_mc(*mc_event->GetMcProtons()[i], *mc_event->GetMcPips()[j], *mc_event->GetMcPims()[k]);
                                                                // if (mc_event->MCinv_Ppip() > 1.22 && mc_event->MCinv_Ppip() < 1.32)
                                                                {
                                                                        _hists->Fill_WvsQ2_twoPi_thrown(data, mc_event);
                                                                        _hists->Fill_histSevenD_thrown_pim(mc_event);
                                                                        _hists->Fill_histSevenD_thrown_prot(mc_event);
                                                                        _hists->Fill_histSevenD_thrown_pip(mc_event);
                                                                        _hists->Fill_histSevenD_thrown_pim_evt(mc_event);
                                                                        _hists->Fill_histSevenD_thrown_prot_evt(mc_event);
                                                                        _hists->Fill_histSevenD_thrown_pip_evt(mc_event);
                                                                        //         // // // _hists->Fill_W_bin_check_th(mc_event);

                                                                        // // // ///// bin centering corr
                                                                        _hists->Fill_hist1D_thrown_w_q2(mc_event);
                                                                        _hists->Fill_hist1D_thrown_inv_mass(mc_event);
                                                                        _hists->Fill_hist1D_thrown_theta(mc_event);
                                                                        _hists->Fill_hist1D_thrown_alpha(mc_event);
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }

*/
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////

                auto event = std::make_shared<Reaction>(data, beam_energy);
                no_of_events++;
                // if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() > 1.95)
                if (data->charge(0) == -1)
                        _hists->FillHists_electron_cuts(data, event);

                ///// auto cuts = std::make_unique<rga_Cuts>(data);
                auto cuts = std::make_unique<Pass2_Cuts>(data);

                if (!cuts->ElectronCuts())
                        continue;
                // event->SetMomCorrElec();

                 if (!_qa->Golden(data->getRun(), data->getEvent())) continue;

                // if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() > 1.95)
                _hists->FillHists_electron_with_cuts(data, event);

                // If we pass electron cuts the event is processed
                elec++;
                total++;
                int statusPim = -99;
                int statusPip = -99;
                int statusProt = -99;
                int sectorPim = -99;
                int sectorPip = -99;
                int sectorProt = -99;
                int prot_idx = -1;
                int pip_idx = -1;

                // // Make a reaction class from the data
                auto dt = std::make_shared<Delta_T>(data);
                auto dt_proton = std::make_shared<Delta_T>(data);
                auto dt_pip = std::make_shared<Delta_T>(data);

                // ///////////////////////////////////   Original method ///////////////////////////
                // ///////////////////////////////////   Original method ///////////////////////////

                for (int part = 1; part < data->gpart(); part++)
                {
                        if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() > 1.95)
                        {

                                dt->dt_calc(part);

                                if (data->charge(part) != 0)
                                {

                                        dt->dt_calc(part);
                                        // // dt->dt_calc_ctof(part);
                                        _hists->Fill_MomVsBeta(data, part, event);
                                        _hists->Fill_deltat_before_cut(data, dt, part, event);

                                        if (data->charge(part) > 0)
                                        {
                                                _hists->FillHists_prot_pid_cuts(data, event, part);
                                                _hists->FillHists_pip_pid_cuts(data, event, part);

                                                /*   if ((cuts->IsProton(part)) && (cuts->IsPip(part)))
                                                   {
                                                           {
                                                                   {
                                                                           event->SetProton(part);
                                                                           prot_idx++;
                                                                           statusProt = abs(data->status(part));
                                                                           sectorProt = data->dc_sec(part);
                                                                           // _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                                           // _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                                                           prot++;
                                                                   }
                                                                   {
                                                                           {
                                                                                   pip++;

                                                                                   // _hists->(event, dp_pip1);
                                                                                   event->SetPip(part);
                                                                                   pip_idx++;
                                                                                   statusPip = abs(data->status(part));
                                                                                   sectorPip = data->dc_sec(part);
                                                                                   // _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                                                   // _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                                                           }
                                                                   }
                                                           }
                                                   }*/

                                                if (cuts->IsProton(part))
                                                {

                                                        prot++;
                                                        event->SetProton(part);
                                                        // event->SetPip(part);
                                                        prot_idx++;

                                                        statusProt = abs(data->status(part));
                                                        sectorProt = data->dc_sec(part);

                                                        _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                        _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                                }

                                                if (cuts->IsPip(part))
                                                {
                                                        {
                                                                pip++;
                                                                event->SetPip(part);
                                                                // event->SetProton(part);
                                                                pip_idx++;

                                                                statusPip = abs(data->status(part));
                                                                sectorPip = data->dc_sec(part);

                                                                _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                                _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                                        }
                                                }
                                        }
                                        else
                                        {
                                                _hists->FillHists_pim_pid_cuts(data, event, part);

                                                if (cuts->IsPim(part))

                                                {
                                                        event->SetPim(part);

                                                        pim++;

                                                        statusPim = abs(data->status(part));
                                                        sectorPim = data->dc_sec(part);
                                                        _hists->FillHists_pim_pid_cuts(data, event, part);

                                                        {
                                                                _hists->Fill_deltat_pim_after_cut(data, dt, part, event);
                                                                _hists->FillHists_pim_pid_with_cuts(data, event, part);
                                                        }
                                                }
                                        }
                                }
                        }
                }
                /////////////////////////////////////  For  dp cut method ////////////////////
                ///////////////////////////////////////  For  dp cut method ////////////////////
                /*
                                // // Define vectors to store dp values and indices for protons and pions
                                std::vector<std::pair<int, double>> proton_dps; // Pair of (index, dp_prot)
                                std::vector<std::pair<int, double>> pip_dps;    // Pair of (index, dp_pip)

                                for (int part = 1; part < data->gpart(); part++)
                                {

                                        if (data->charge(part) != 0)
                                        {
                                                dt->dt_calc(part);
                                                // _hists->Fill_MomVsBeta(data, part, event);
                                                _hists->Fill_deltat_before_cut(data, dt, part, event);

                                                if (data->charge(part) > 0)
                                                {

                                                        _hists->FillHists_prot_pid_cuts(data, event, part);
                                                        _hists->FillHists_pip_pid_cuts(data, event, part);

                                                        // Check if the particle satisfies proton and/or pion conditions
                                                        if (cuts->IsProton(part))
                                                        {
                                                                prot++;
                                                                double dp_prot = pow(mc_event->prot_momX_mc_gen() - data->px(part), 2) +
                                                                                 pow(mc_event->prot_momY_mc_gen() - data->py(part), 2) +
                                                                                 pow(mc_event->prot_momZ_mc_gen() - data->pz(part), 2);
                                                                proton_dps.push_back(std::make_pair(part, dp_prot)); // Store index and dp value for proton

                                                                _hists->Fill_deltaP_ambi_prot(event, dp_prot);
                                                        }

                                                        if (cuts->IsPip(part))
                                                        {
                                                                pip++;
                                                                double dp_pip = pow(mc_event->pip_momX_mc_gen() - data->px(part), 2) +
                                                                                pow(mc_event->pip_momY_mc_gen() - data->py(part), 2) +
                                                                                pow(mc_event->pip_momZ_mc_gen() - data->pz(part), 2);
                                                                pip_dps.push_back(std::make_pair(part, dp_pip)); // Store index and dp value for proton

                                                                _hists->Fill_deltaP_ambi_pip(event, dp_pip);
                                                        }
                                                }

                                                else
                                                {
                                                        _hists->FillHists_pim_pid_cuts(data, event, part);

                                                        if (cuts->IsPim(part))

                                                        {
                                                                event->SetPim(part);

                                                                pim++;

                                                                statusPim = abs(data->status(part));
                                                                sectorPim = data->dc_sec(part);

                                                                _hists->Fill_deltat_pim_after_cut(data, dt, part, event);
                                                                _hists->FillHists_pim_pid_with_cuts(data, event, part);
                                                        }
                                                }
                                        }
                                }

                                // // // Now, find the pair of proton and pip with the minimum dp_prot + dp_pip
                                double min_dp_sum = std::numeric_limits<double>::max();
                                int best_proton_index = -1;
                                int best_pip_index = -1;
                                std::vector<std::pair<int, int>> non_minimum_pairs; // Stores all non-minimum proton-pip pairs

                                // Loop over all combinations of protons and pions to find the minimum dp_prot + dp_pip
                                for (size_t i = 0; i < proton_dps.size(); i++)
                                {
                                        int prot_index = proton_dps[i].first;
                                        double dp_prot = proton_dps[i].second;

                                        for (size_t j = 0; j < pip_dps.size(); j++)
                                        {

                                                int pip_index = pip_dps[j].first;
                                                double dp_pip = pip_dps[j].second;

                                                double dp_sum = dp_prot + dp_pip;
                                                if (dp_sum < min_dp_sum)
                                                {
                                                        min_dp_sum = dp_sum;
                                                        best_proton_index = prot_index;
                                                        best_pip_index = pip_index;
                                                }
                                                // if (event->TwoPion_missingPim())
                                                // {
                                                _hists->Fill_deltaP_prot_for_pip(event, dp_prot);
                                                _hists->Fill_deltaP_pip_for_prot(event, dp_pip);
                                                _hists->Fill_deltaP_sum(event, dp_prot + dp_pip);

                                                // std::cout << "  .......................   Outsside .................... " << std::endl;
                                                // std::cout << "  event number  " << current_event << "  Number of particles: " << data->gpart()
                                                //           << "  prot size  " << proton_dps.size() << "   pip size  " << pip_dps.size()
                                                //           << std::endl;
                                                // std::cout << "proton_dps " << proton_dps[i].second << std::endl;
                                                // std::cout << "pip_dps " << pip_dps[j].second << "   sum dp  " << dp_sum << std::endl;
                                                // }
                                        }
                                }

                                // // Set the proton and pip with the minimum dp_sum for further processing
                                // if (best_proton_index != -1 && best_pip_index != -1)
                                // {
                                //         event->SetProton(best_proton_index);
                                //         event->SetPip(best_pip_index);
                                // }

                                // Overlapped loop over all combinations of protons and pions
                                for (size_t i = 0; i < proton_dps.size(); i++)
                                {
                                        int prot_index = proton_dps[i].first;
                                        event->SetProton(prot_index);
                                        // event->SetPip(prot_index); // for swap
                                }

                                for (size_t j = 0; j < pip_dps.size(); j++)
                                {
                                        int pip_index = pip_dps[j].first;
                                        event->SetPip(pip_index);
                                        // event->SetProton(pip_index); // for swap
                                }
                */
                ///////////////////////////////////  Now start cutting on our kinamatics and selecting twoPion events for analysis //////////////
                // if (event->W() > 1. && event->W() <= 3.0 && event->Q2() <= 12.0 && event->Q2() >= 1.)
                if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() > 1.95)
                {
                        events_passes_w_q2_cuts++;
                        {

                                // // if (event->TwoPion_exclusive())
                                // // {
                                // //         for (size_t i = 0; i < event->GetProtons().size(); ++i)
                                // //         {
                                // //                 for (size_t j = 0; j < event->GetPips().size(); ++j)
                                // //                 {
                                // for (size_t k = 0; k < event->GetPims().size(); ++k)
                                // {
                                //         //                                 event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

                                //         //                                 two_pion_Excl_events++;
                                //         //                                 _hists->Fill_WvsQ2(event);

                                //         // // You should have a similar method for π⁻ if applicable
                                //         dt->dt_calc(event->GetPimIndices()[k]);
                                //         _hists->Fill_deltat_pim_after_cut(data, dt, event->GetPimIndices()[k], event);
                                //         _hists->FillHists_pim_pid_with_cuts(data, event, event->GetPimIndices()[k]);
                                // }
                                // //                 }
                                // //         }
                                // // }

                                if (event->TwoPion_missingPim())
                                { // // Retrieve the number of protons and pions in the event

                                        int num_protons = event->GetProtons().size();
                                        int num_pips = event->GetPips().size();
                                        int num_combinations = 0;
                                        // std::cout << "  .......................   Insside .................... " << std::endl;
                                        // std::cout << "  event number  " << current_event << "   Number of particles: " << data->gpart() << std::endl;
                                        // std::cout << std::endl;
                                        // // Calculate the total number of combinations
                                        // size_t num_combinations = num_protons * num_pips;

                                        // // Get the vector of protons
                                        // const auto &protons = event->GetProtons();
                                        // const auto &pip = event->GetPips();

                                        ////////////  CONTROL OVER HAOW MANY Prot/Pip PER EVENT /////////
                                        ////////////  CONTROL OVER HAOW MANY Prot/Pip PER EVENT /////////
                                        ////////////  CONTROL OVER HAOW MANY Prot/Pip PER EVENT /////////
                                        _hists->Fill_WvsQ2(event);

                                        // First loop: count valid combinations
                                        for (size_t i = 0; i < num_protons; ++i)
                                        {
                                                // std::cout << " prot dp inside loop : " << proton_dps[i].second << std::endl;

                                                for (size_t j = 0; j < num_pips; ++j)
                                                {
                                                        // std::cout << " pip dp inside loop   : " << pip_dps[j].second << std::endl;

                                                        // if (event->GetProtonIndices()[i] != event->GetPipIndices()[j])
                                                        {
                                                                num_combinations++;
                                                        }
                                                }
                                        }

                                        ///////////////////////////////////////////////////

                                        bool proton_cdfd_cut = false;
                                        bool pip_cdfd_cut = false;
                                        double dp_Prot = NAN;
                                        double dtheta_Prot = NAN;
                                        double dphi_Prot = NAN;
                                        double dp_Pip = NAN;
                                        double dtheta_Pip = NAN;
                                        double dphi_Pip = NAN;

                                        // Loop over particles to check proton and pip
                                        for (int part1 = 0; part1 < data->gpart(); part1++)
                                        {
                                                // std::cout << "  current event  " << current_event << "  no of parts = " << data->gpart()
                                                //           << " part at " << part1 << " is  " << data->pid(part1) << std::endl;
                                                dt->dt_calc(part1);

                                                // Proton Block
                                                if (cuts->IsProton(part1))
                                                {
                                                        bool isFD1 = ((data->status(part1) > 2000) && (data->status(part1) < 4000));
                                                        bool isCD1 = (data->status(part1) > 4000);
                                                        // if (isFD1)
                                                        //   std::cout << "   part1 " << part1 << "  is in FD with status " << data->status(part1) << "   with id "
                                                        //             << data->pid(part1) << std::endl;
                                                        // if (isCD1)
                                                        //   std::cout << "   part1 " << part1 << "  is in CD with status " << data->status(part1) << "   with id "
                                                        //             << data->pid(part1) << std::endl;

                                                        if (isFD1 || isCD1)
                                                        { // Ensure part1 is in FD or CD
                                                                TLorentzVector track1;

                                                                track1.SetXYZM(data->px(part1), data->py(part1), data->pz(part1), MASS_P);

                                                                for (int part2 = part1 + 1; part2 < data->gpart(); part2++)
                                                                {
                                                                        if (cuts->IsProton(part2))
                                                                        {
                                                                                bool isFD2 = ((data->status(part2) > 2000) && (data->status(part2) < 4000));
                                                                                bool isCD2 = (data->status(part2) > 4000);

                                                                                // if (isFD2)
                                                                                //   std::cout << "   part2  " << part2 << "  is in FD with status " << data->status(part2)
                                                                                //             << "   with id " << data->pid(part2) << std::endl;
                                                                                // if (isCD2)
                                                                                //   std::cout << "   part2 " << part2 << "  is in CD with status " << data->status(part2)
                                                                                //             << "   with id " << data->pid(part2) << std::endl;

                                                                                if ((isFD1 && isCD2) || (isCD1 && isFD2))
                                                                                {
                                                                                        TLorentzVector track2;
                                                                                        track2.SetXYZM(data->px(part2), data->py(part2), data->pz(part2), MASS_P);

                                                                                        TLorentzVector trackFD = isFD1 ? track1 : track2;
                                                                                        TLorentzVector trackCD = isCD1 ? track1 : track2;

                                                                                        dp_Prot = (trackFD.P() - trackCD.P());
                                                                                        dtheta_Prot = (trackFD.Theta() - trackCD.Theta()) * 180 / PI;
                                                                                        dphi_Prot = (trackFD.Phi() - trackCD.Phi()) * 180 / PI;
                                                                                        // h_dp_prot->Fill(dp_Prot, event->weight());
                                                                                        // h_dtheta_prot->Fill(dtheta_Prot, event->weight());
                                                                                        // h_dphi_prot->Fill(dphi_Prot, event->weight());
                                                                                        // Apply proton cuts
                                                                                        // if (dp_Prot > -0.6 && dp_Prot < 0.2 && dtheta_Prot > -7 && dphi_Prot > -20 && dphi_Prot < 5)
                                                                                        if (dp_Prot > -0.3 && dp_Prot < 0.0 && dtheta_Prot > -3 && dtheta_Prot < 3 && dphi_Prot > -5 && dphi_Prot < 2)
                                                                                        {
                                                                                                proton_cdfd_cut = true;
                                                                                                //  // } else {  // Fill histograms
                                                                                                // h_dp_prot->Fill(dp_Prot, event->weight());
                                                                                                // h_dtheta_prot->Fill(dtheta_Prot, event->weight());
                                                                                                // h_dphi_prot->Fill(dphi_Prot, event->weight());
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }

                                                ///////////////////////////
                                                // Pip Block
                                                if (cuts->IsPip(part1))
                                                {
                                                        bool isFD1 = ((data->status(part1) > 2000) && (data->status(part1) < 4000));
                                                        bool isCD1 = (data->status(part1) > 4000);

                                                        if (isFD1 || isCD1)
                                                        { // Ensure part1 is in FD or CD
                                                                TLorentzVector track1;
                                                                track1.SetXYZM(data->px(part1), data->py(part1), data->pz(part1), MASS_PIP);

                                                                for (int part2 = part1 + 1; part2 < data->gpart(); part2++)
                                                                {
                                                                        if (cuts->IsPip(part2))
                                                                        {
                                                                                bool isFD2 = ((data->status(part2) > 2000) && (data->status(part2) < 4000));
                                                                                bool isCD2 = (data->status(part2) > 4000);

                                                                                if ((isFD1 && isCD2) || (isCD1 && isFD2))
                                                                                {
                                                                                        TLorentzVector track2;
                                                                                        track2.SetXYZM(data->px(part2), data->py(part2), data->pz(part2), MASS_PIP);

                                                                                        TLorentzVector trackFD = isFD1 ? track1 : track2;
                                                                                        TLorentzVector trackCD = isCD1 ? track1 : track2;

                                                                                        dp_Pip = (trackFD.P() - trackCD.P());
                                                                                        dtheta_Pip = (trackFD.Theta() - trackCD.Theta()) * 180 / PI;
                                                                                        dphi_Pip = (trackFD.Phi() - trackCD.Phi()) * 180 / PI;
                                                                                        // h_dp_pip->Fill(dp_Pip, event->weight());
                                                                                        // h_dtheta_pip->Fill(dtheta_Pip, event->weight());
                                                                                        // h_dphi_pip->Fill(dphi_Pip, event->weight());
                                                                                        // Apply pip cuts
                                                                                        // if (dp_Pip > -0.4 && dp_Pip < 0.2 && dtheta_Pip > -10 && dtheta_Pip < 10 && dphi_Pip > -20 &&
                                                                                        //     dphi_Pip < 5)
                                                                                        if (dp_Prot > -0.25 && dp_Prot < 0.0 && dtheta_Prot > -3 && dtheta_Prot < 3 && dphi_Prot > -5 && dphi_Prot < 2)

                                                                                        {
                                                                                                pip_cdfd_cut = true;
                                                                                                // // } else {  // Fill histograms
                                                                                                // h_dp_pip->Fill(dp_Pip, event->weight());
                                                                                                // h_dtheta_pip->Fill(dtheta_Pip, event->weight());
                                                                                                // h_dphi_pip->Fill(dphi_Pip, event->weight());
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }

                                        for (size_t i = 0; i < num_protons; ++i)
                                        {
                                                for (size_t j = 0; j < num_pips; ++j)
                                                {
/////                                                        if (!(proton_cdfd_cut == true || pip_cdfd_cut == true))
                                                        {
                                                                // Exclude the case where the same particle is assigned as both proton and pip
                                                                if (event->GetProtonIndices()[i] != event->GetPipIndices()[j])
                                                                {
                                                                        int proton_part_idx = event->GetProtonIndices()[i];
                                                                        int pip_part_idx = event->GetPipIndices()[j];

                                                                        // if ((best_proton_index != event->GetProtonIndices()[i]) || (best_pip_index != event->GetPipIndices()[j]))
                                                                        {

                                                                                two_pion_mPim_events++;
                                                                                entries_in_this_event++;
                                                                                event->SetSwappedProton(pip_part_idx);
                                                                                event->SetSwappedPip(proton_part_idx);
                                                                                // Extract velocity components with energy normalization
                                                                                double v_original_x_Prot = event->GetProtons()[i]->Px() / event->GetProtons()[i]->E();
                                                                                double v_original_y_Prot = event->GetProtons()[i]->Py() / event->GetProtons()[i]->E();
                                                                                double v_original_z_Prot = event->GetProtons()[i]->Pz() / event->GetProtons()[i]->E();
                                                                                double v_swapped_x_Prot = event->GetProtonsSwapped()->Px() / event->GetProtonsSwapped()->E();
                                                                                double v_swapped_y_Prot = event->GetProtonsSwapped()->Py() / event->GetProtonsSwapped()->E();
                                                                                double v_swapped_z_Prot = event->GetProtonsSwapped()->Pz() / event->GetProtonsSwapped()->E();

                                                                                // Calculate delta_V^2
                                                                                double dv2_Prot = pow(v_swapped_x_Prot - v_original_x_Prot, 2) +
                                                                                                  pow(v_swapped_y_Prot - v_original_y_Prot, 2) +
                                                                                                  pow(v_swapped_z_Prot - v_original_z_Prot, 2);
                                                                                // std::cout << "Delta V^2: " << dv2 << std::endl;

                                                                                // // Extract velocity components with energy normalization
                                                                                // double v_original_x_Pip = event->GetPips()[j]->Px() / event->GetPips()[j]->E();
                                                                                // double v_original_y_Pip = event->GetPips()[j]->Py() / event->GetPips()[j]->E();
                                                                                // double v_original_z_Pip = event->GetPips()[j]->Pz() / event->GetPips()[j]->E();
                                                                                // double v_swapped_x_Pip = event->GetPipsSwapped()->Px() / event->GetPipsSwapped()->E();
                                                                                // double v_swapped_y_Pip = event->GetPipsSwapped()->Py() / event->GetPipsSwapped()->E();
                                                                                // double v_swapped_z_Pip = event->GetPipsSwapped()->Pz() / event->GetPipsSwapped()->E();

                                                                                // // Calculate delta_V^2
                                                                                // double dv2_Pip = pow(v_swapped_x_Pip - v_original_x_Pip, 2) +
                                                                                //                  pow(v_swapped_y_Pip - v_original_y_Pip, 2) +
                                                                                //                  pow(v_swapped_z_Pip - v_original_z_Pip, 2);

                                                                                // std::cout << "  event->GetProtonIndices()[i]  : " << event->GetProtonIndices()[i] << "  Pip index  : " << event->GetPipIndices()[j] << std::endl;

                                                                                event->CalcMissMassPim(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                                event->boost(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                                event->CalcMissMassPimSwapped();
                                                                                ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                                                                                ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                                                                                ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                                                                                // if (num_combinations == 2)
                                                                                // if (event->MM2_mPim() < -0.1)
                                                                                // if (dv2_Prot >= 0.0005 && dv2_Prot < 0.001)
                                                                                // if (event->MM2_mPim() > -0.1 && event->MM2_mPim() < 0.1)
                                                                                if (_hists->MM_cut(event->W(), event->Q2(), event->MM2_mPim()))
                                                                                {
                                                                                        // if (dv2_Prot > 0.01)
                                                                                        {
                                                                                                // if (num_combinations == 1)
                                                                                                // {
                                                                                                _hists->Fill_all_Combi(event);
                                                                                                _hists->Fill_MMSQ_mPim(event);

                                                                                                // Fill the appropriate histogram based on num_combinations
                                                                                                if (num_combinations == 1)
                                                                                                {
                                                                                                        _hists->Fill_1_Combi(event);
                                                                                                        _hists->Fill_MMSQ_mPim_1_comb(event);
                                                                                                }
                                                                                                // if (num_combinations == 2)
                                                                                                // {
                                                                                                //         _hists->Fill_2_Combi(event);
                                                                                                //         _hists->Fill_MMSQ_mPim_2_comb(event);
                                                                                                // }
                                                                                                // if (num_combinations == 3)
                                                                                                // {
                                                                                                //         _hists->Fill_3_Combi(event);
                                                                                                //         _hists->Fill_MMSQ_mPim_3_comb(event);
                                                                                                // }
                                                                                                // if (num_combinations >= 4)
                                                                                                // {
                                                                                                //         _hists->Fill_4_or_more_Combi(event);
                                                                                                //         _hists->Fill_MMSQ_mPim_4_or_more_comb(event);
                                                                                                // }

                                                                                                //         // if (event->Fixed_MM_cut())
                                                                                                //         // if (MM_cut(event->W(), event->Q2(), event->MM2_mPim()))
                                                                                                //         // if ((data->p(event->GetProtonIndices()[i]) > 3.0) || (data->p(event->GetPipIndices()[j]) > 3.0))

                                                                                                //         {
                                                                                                // _hists->Fill_WvsQ2(event);

                                                                                                // _hists->FillHists_electron_with_cuts(data, event);

                                                                                                // two_pion_mPim_events++;
                                                                                                // {
                                                                                                // _hists->Fill_MMSQ_mPim(event);
                                                                                                // if (entries_in_this_event == 1)
                                                                                                // _hists->Fill_WvsQ2(event);

                                                                                                _hists->Fill_histSevenD_prot(event);
                                                                                                _hists->Fill_histSevenD_pip(event);
                                                                                                _hists->Fill_histSevenD_pim(event);

                                                                                                _hists->Fill_histSevenD_prot_evt(event);
                                                                                                _hists->Fill_histSevenD_pip_evt(event);
                                                                                                _hists->Fill_histSevenD_pim_evt(event);

                                                                                                // //         }
                                                                                                // // }
                                                                                                // //         else if (event->TwoPion_missingPip())
                                                                                                // //                 miss_pip++;
                                                                                                // //         else if (event->TwoPion_missingProt())
                                                                                                // //                 miss_prot++;

                                                                                                // //         //{
                                                                                                // // //         twopi++;
                                                                                                // // First, check if the index is for proton or pip, then use it as needed

                                                                                                // dt_proton->dt_calc(proton_part_idx);

                                                                                                // _hists->Fill_MomVsBeta(data, proton_part_idx, event);
                                                                                                // _hists->Fill_deltat_prot_after_cut(data, dt_proton, proton_part_idx, event);
                                                                                                // _hists->FillHists_prot_pid_with_cuts(data, event, proton_part_idx, *event->GetProtons()[i]);
                                                                                                // // }

                                                                                                // // std::cout << "   pip_part_idx  " << pip_part_idx << std::endl;
                                                                                                // dt_pip->dt_calc(pip_part_idx);
                                                                                                // _hists->Fill_MomVsBeta(data, pip_part_idx, event);
                                                                                                // _hists->Fill_deltat_pip_after_cut(data, dt_pip, pip_part_idx, event);
                                                                                                // _hists->FillHists_pip_pid_with_cuts(data, event, pip_part_idx, *event->GetPips()[j]);

                                                                                                // _hists->Fill_2_Combi(event);
                                                                                                // _hists->Fill_MMSQ_mPim_2_comb(event);

                                                                                                // _hists->Fill_3_Combi(dv2_Pip, event);
                                                                                                // _hists->Fill_MMSQ_mPim_3_comb(dv2_Pip, event);

                                                                                                // ///////////
                                                                                                // _hists->Fill_WvsQ2(event);
                                                                                                // // // std::cout << " proton velocity : " << event->GetProtons()[i]->Px() << std::endl;

                                                                                                // _hists->Fill_4_or_more_Combi(dv2_Prot, event);
                                                                                                // _hists->Fill_MMSQ_mPim_4_or_more_comb(dv2_Prot, event);

                                                                                                // _hists->Fill_deltaP_prot(event, proton_dps[i].second);
                                                                                                // _hists->Fill_deltaP_pip(event, pip_dps[j].second);
                                                                                                // _hists->Fill_deltaP_sum_twoPi(event, proton_dps[i].second + pip_dps[j].second);

                                                                                                // }

                                                                                                // }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }

                                        // if (num_combinations == 1)
                                        // {
                                        //         first_entries++;
                                        // }
                                        // if (num_combinations == 2)
                                        // {
                                        //         second_entries++;
                                        // }
                                        // if (num_combinations == 3)
                                        // {
                                        //         third_entries++;
                                        // }
                                        // if (num_combinations > 3)
                                        // {
                                        //         four_or_more_entries++;
                                        // }
                                        _hists->Fill_Entries_prot(num_protons);
                                        _hists->Fill_Entries_pip(num_pips);
                                        _hists->Fill_Entries(num_combinations);
                                }
                        }
                }
        }
        // std::cout.precision(3);

        std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
        std::cout << "   no of total events  " << num_of_events << std::endl;
        // std::cout << " elec " << elec << "  electron as pid(0)  " << pid_zero_elec << " prot " << prot << " pip " << pip << " pim " << pim << '\n';
        // std::cout << "   nonzero wt events   " << events_with_non_zero_wt << " , " << events_with_non_zero_wt / (float)num_of_events * 100 << std::endl;
        // std::cout << "   zero wt events " << events_with_zero_wt << "  ,  " << events_with_zero_wt / (float)num_of_events * 100 << std::endl;
        // std::cout << "   events passing electron cuts  " << elec << ",  " << elec / (float)(no_of_events) * 100 << std::endl;
        // std::cout << "   events passing w-q2 cuts " << events_passes_w_q2_cuts << "  ,  " << events_passes_w_q2_cuts / (float)(no_of_events) * 100 << std::endl;
        // std::cout << "   no of twoPion events (mPim topo ) = " << two_pion_mPim_events << "  ; " << float(two_pion_mPim_events) / float(num_of_events) * 100 << std::endl;
        // // std::cout << "  first entry only " << first_entries << "  % is : " << (first_entries) / (float)(two_pion_mPim_events) * 100 << std::endl;
        // // std::cout << "  second entry only " << second_entries << "  % is : " << (second_entries) / (float)(two_pion_mPim_events) * 100 << std::endl;
        // // std::cout << "  third entry only " << third_entries << "  % is : " << (third_entries) / (float)(two_pion_mPim_events) * 100 << std::endl;
        // // std::cout << "  four or more entries only " << four_or_more_entries << "  % is : " << (four_or_more_entries) / (float)(two_pion_mPim_events) * 100 << std::endl;

        // // Return the total number of events
        return num_of_events;
}
#endif
