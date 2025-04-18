
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "cuts.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id)
{
        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();

        float beam_energy = 10.6041;
        // float beam_energy = 5.754;

        if (std::is_same<CutType, Pass2_Cuts>::value)
        {
                beam_energy = 10.6041;
                // beam_energy = 5.754;
        }
        // Get the number of events in this thread
        // size_t num_of_events = (int)_chain->GetEntries();
        // float beam_energy = NAN;
        if (getenv("CLAS12_E") != NULL)
                beam_energy = atof(getenv("CLAS12_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // Make a data object which all the branches can be accessed from
        auto data = std::make_shared<Branches12>(_chain, true);

        // Total number of events "Processed"
        size_t total = 0;
        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++)
        {
                // Get current event
                _chain->GetEntry(current_event);
                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
                // std::cout << "  event number " << current_event << '\n';

                // if (data->mc_npart() < 1)
                //         continue;

                // If we pass electron cuts the event is processed
                total++;
                int status_pim = -9999;
                int status_pip = -9999;
                int status_prot = -9999;

                bool proton_1 = false;
                bool pip_1 = false;
                bool pim_1 = false;

                // Make a reaction class from the data given
                auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

                for (int part = 1; part < data->mc_npart(); part++)
                {

                        // Check particle ID's and fill the reaction class
                        if (data->mc_pid(part) == PIP)
                        {
                                mc_event->SetMCPip(part);
                        }
                        else if (data->mc_pid(part) == PROTON)
                        {
                                mc_event->SetMCProton(part);
                        }
                        else if (data->mc_pid(part) == PIM)
                        {
                                mc_event->SetMCPim(part);
                                // } else {
                                //   mc_event->SetMCOther(part);
                        }
                }
                // if (mc_event->weight() > 0.0)
                // {
                //         // if (mc_event->W_mc() < 2.15 && mc_event->W_mc() > 1.3 && mc_event->Q2_mc() < 9.0 && mc_event->Q2_mc() > 1.5 && mc_event->weight() > 0.0)
                //         // {

                //         // // _hists->Fill_WvsQ2_twoPi_thrown(mc_event);
                //         // _hists->Fill_histSevenD_thrown_pim(mc_event);
                //         // _hists->Fill_histSevenD_thrown_prot(mc_event);
                //         // _hists->Fill_histSevenD_thrown_pip(mc_event);
                //         // _hists->Fill_histSevenD_thrown_pim_evt(mc_event);
                //         // _hists->Fill_histSevenD_thrown_prot_evt(mc_event);
                //         // _hists->Fill_histSevenD_thrown_pip_evt(mc_event);
                //         // //         // // // _hists->Fill_W_bin_check_th(mc_event);

                //         // // // ///// bin centering corr
                //         // _hists->Fill_hist1D_thrown_w_q2(mc_event);
                //         // _hists->Fill_hist1D_thrown_inv_mass(mc_event);
                //         // _hists->Fill_hist1D_thrown_theta(mc_event);
                //         // _hists->Fill_hist1D_thrown_alpha(mc_event);
                // }

                auto event = std::make_shared<Reaction>(data, beam_energy);

                if (data->gpart() == 0)
                        continue;
                bool elec = true;
                elec &= (data->charge(0) == NEGATIVE);
                elec &= (data->pid(0) == 11);
                if (!elec)
                        continue;
                // if (data->charge(0) == -1)
                // _hists->FillHists_electron_cuts(data, event);
                auto cuts = std::make_shared<Pass2_Cuts>(data);
                if (!cuts->ElectronCuts("mid"))
                        continue;
                // _hists->FillHists_electron_with_cuts(data, event);

                auto dt = std::make_shared<Delta_T>(data);
                // For each particle in the event

                for (int part = 1; part < data->gpart(); part++)
                {
                        dt->dt_calc(part);
                        // //        dt->dt_calc_ctof(part);
                        // _hists->Fill_MomVsBeta(data, part, event);

                        // _hists->Fill_deltat_before_cut(data, dt, part, event);

                        if (cuts->IsProton(part, "mid"))
                        {
                                // _hists->FillHists_prot_pid_cuts(data, event, part);

                                // if (cuts->HadronsCuts(part))
                                {
                                        event->SetProton(part);
                                        status_prot = abs(data->status(part));
                                        // _hists->FillHists_prot_pid_with_cuts(data, event, part);
                                        // _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                }
                        }
                        else if (cuts->IsPip(part, "mid"))
                        {
                                // _hists->FillHists_pip_pid_cuts(data, event, part);
                                // if (cuts->HadronsCuts(part))
                                {
                                        pip_1 = true;

                                        // _hists->FillHists_pip_pid_with_cuts(data, event, part);
                                        // _hists->Fill_deltat_pip_after_cut(data, dt, part, event);

                                        event->SetPip(part);
                                        status_pip = abs(data->status(part));
                                }
                        }
                        else if (cuts->IsPim(part, "mid"))
                        {
                                // _hists->FillHists_pim_pid_cuts(data, event, part);

                                // if (cuts->HadronsCuts(part))
                                {
                                        pim_1 = true;

                                        // if (event->MM_cut())
                                        event->SetPim(part);
                                        status_pim = abs(data->status(part));
                                        // _hists->FillHists_pim_pid_with_cuts(data, event, part);
                                        // _hists->Fill_deltat_pim_after_cut(data, dt, part, event);

                                        // } else if (cuts->IsmissingPim(part)) {
                                        //   event->SetmissingPim(part);
                                }
                        }

                        else
                        {
                                event->SetOther(part);
                        }

                        // // if (proton_1 && pip_1 )
                        // // {
                        // _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                        // _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                        // _hists->Fill_deltat_pim_after_cut(data, dt, part, event);
                        // // }
                }

                if (event->W() < 2.5 && event->W() > 1.3 && event->Q2() < 10.5 && event->Q2() > 1.5 && event->weight() > 0.0)
                {
                        // _hists->Fill_MMSQ_mPim(event);
                        _hists->Fill_WvsQ2(event);

                        // // // // if (event->TwoPion_missingPim())
                        // // // // {
                        // // // // //         // if (event->TwoPion_exclusive()) {
                        // // // //         // if (event->Fixed_MM_cut())
                        // // // //         // {
                        // // // // //                 // _hists->Fill_WvsQ2(event);
                        // // // // //                 // _hists->Fill_x_mu(event);
                        // // // // //                 // _hists->Fill_WvsQ2_twoPi(event);
                        // _hists->Fill_histSevenD_pim(event);
                        // _hists->Fill_histSevenD_pip(event);
                        // _hists->Fill_histSevenD_prot(event);
                        // _hists->Fill_histSevenD_pim_evt(event);
                        // _hists->Fill_histSevenD_pip_evt(event);
                        // _hists->Fill_histSevenD_prot_evt(event);
                        // // // // //         }
                        // // // // }
                }
        }

        // Return the total number of events
        return num_of_events;
}

#endif
