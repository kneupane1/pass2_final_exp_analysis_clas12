#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12> &data, float beam_energy)
{
        _data = data;
        _beam = std::make_unique<TLorentzVector>();
        _beam_energy = 10.6041; //
        // _beam_energy = atof(getenv("CLAS12_E"));

        _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

        _gamma = std::make_unique<TLorentzVector>();
        _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
        _elecUnSmear = std::make_unique<TLorentzVector>();

        _elec = std::make_unique<TLorentzVector>();
        // // mom corr
        // _mom_corr_elec = std::make_unique<TLorentzVector>();
        this->SetElec();

        _prot = std::vector<std::unique_ptr<TLorentzVector>>(); // Initialize _protons as an empty vector
        _pip = std::vector<std::unique_ptr<TLorentzVector>>();
        _pim = std::vector<std::unique_ptr<TLorentzVector>>();

        _mom_corr_prot = std::vector<std::unique_ptr<TLorentzVector>>();
        _mom_corr_pip = std::vector<std::unique_ptr<TLorentzVector>>();
        _mom_corr_pim = std::vector<std::unique_ptr<TLorentzVector>>();

        _prot_indices = std::vector<int>();
        _pip_indices = std::vector<int>();
        _pim_indices = std::vector<int>();

        _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
        _Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
        _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();

        _other = std::make_unique<TLorentzVector>();
        _neutron = std::make_unique<TLorentzVector>();

        _protUnSmear = std::make_unique<TLorentzVector>();
        _pipUnSmear = std::make_unique<TLorentzVector>();
        _pimUnSmear = std::make_unique<TLorentzVector>();

        _swapped_prot = std::make_unique<TLorentzVector>();
        _swapped_pip = std::make_unique<TLorentzVector>();

        // _weight = _data->mc_weight();   //
        //          1.0;
}

Reaction::~Reaction() {}
auto objMomCorr = std::make_shared<mom_corr>();
auto objEffCorr = std::make_shared<EffCorr>();

void Reaction::SetElec()
{
        _numPart++;
        _hasE = true;
        // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
        // ////////_elec->SetPxPyPzE(_data->px(0), _data->py(0), _data->pz(0), sqrt(_data->px(0) * _data->px(0) + _data->py(0) * _data->py(0) + _data->pz(0) * _data->pz(0) + MASS_E * MASS_E));

        // *_gamma += *_beam - *_elec;
        // // // // // Can calculate W and Q2 here
        // _W = physics::W_calc(*_beam, *_elec);
        // _Q2 = physics::Q2_calc(*_beam, *_elec);
        // // _mom_corr_elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E); // unsmeared
        // // // // //
        // // // // // // //////////////////////////////////////////////////////////////
        // // // // //////////////////////////////////////////////////////////////
        _sectorElec = _data->dc_sec(0);
        _elec_status = abs(_data->status(0));

        if (_mc)
        {

                _elecUnSmear->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

                double W_unsmear = physics::W_calc(*_beam, *_elecUnSmear);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear;

                pUnSmear = _elecUnSmear->P();

                thetaUnSmear = _elecUnSmear->Theta() * 180 / PI;

                if (_elecUnSmear->Phi() > 0)
                        phiUnSmear = _elecUnSmear->Phi() * 180 / PI;
                else if (_elecUnSmear->Phi() < 0)
                        phiUnSmear = (_elecUnSmear->Phi() + 2 * PI) * 180 / PI;

                ////////////////////////////////////////////////////////////////

                // Generate new values
                Reaction::SmearingFunc(ELECTRON, _elec_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _elecUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _elecUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _elecUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _elecSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
                _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E); // smeared
                // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);  // unsmeared

                *_gamma += *_beam - *_elec; // be careful you are commenting this only to include the momentum correction

                // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
                _W = physics::W_calc(*_beam, *_elec);
                _Q2 = physics::Q2_calc(*_beam, *_elec);

                _P_elec = _elec->P();
                _E_elec = _elec->E();
                _theta_e = _elec->Theta() * 180 / PI;

                // if (_elec->Phi() > 0)
                //         _phi_elec = _elec->Phi() * 180 / PI;
                // else if (_elec->Phi() < 0)
                //         _phi_elec = (_elec->Phi() + 2 * PI) * 180 / PI;
        }
        else
        {
                fe = objMomCorr->dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) * 0.9 + 1;
                _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                               MASS_E); // this is new electron mom corrections aug 2022
                // _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                //                MASS_E); // elec and mom corr elec are SAME !!!!!
                // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E); // unsmeared

                *_gamma += *_beam - *_elec;
                _W = physics::W_calc(*_beam, *_elec);
                _Q2 = physics::Q2_calc(*_beam, *_elec);

                _P_elec = _elec->P();
                _E_elec = _elec->E();
                _theta_e = _elec->Theta() * 180 / PI;
        }
}
// // // ///////////////////////////// MOM CORR /////////////////////////////////
// // //////////////////////////////// MOM CORR //////////////////////////////

// void Reaction::SetMomCorrElec()
// { // New electron momentum corrections

//         if (!_mc)
//         {
//                 fe = objMomCorr->dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1;
//                 _mom_corr_elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
//                                         MASS_E); // this is new electron mom corrections aug 2022
//                 // _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
//                 //                MASS_E); // elec and mom corr elec are SAME !!!!!

//                 *_gamma += *_beam - *_mom_corr_elec;
//                 // // // _W_after = physics::W_calc(*_beam, *_mom_corr_elec);
//                 // _W = physics::W_calc(*_beam, *_mom_corr_elec);
//                 // _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);

//                 _P_elec = _mom_corr_elec->P();
//                 _E_elec = _mom_corr_elec->E();
//                 _theta_e = _mom_corr_elec->Theta() * 180 / PI;
//         }
// }

void Reaction::SetProton(int i)
{
        _numProt++;
        _numPos++;
        _hasP = true;
        auto proton = std::make_unique<TLorentzVector>();
        auto mom_corr_proton = std::make_unique<TLorentzVector>();

        _prot_status = abs(_data->status(i));
        _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

        //  proton->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

        _sectorProt = _data->dc_sec(i);
        _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
        _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
        if (_Energy_loss_uncorr_prot->Phi() > 0)
                _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_prot->Phi() < 0)
                _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

        _is_FD_Prot = objMomCorr->is_FD(_prot_status);
        _is_CD_Prot = objMomCorr->is_CD(_prot_status);
        if (_is_CD_Prot)
        {
                _prot_mom_tmt = _prot_mom_uncorr;
                // _prot_theta_tmt = _prot_theta_uncorr;
                // _prot_phi_tmt = _prot_phi_uncorr;

                // _prot_mom_tmt = objMomCorr->CD_prot_Emom_corr(_prot_mom_uncorr, _prot_theta_uncorr);
                // _prot_theta_tmt = objMomCorr->CD_prot_Eth_corr(_prot_mom_uncorr, _prot_theta_uncorr);
                // _prot_phi_tmt = objMomCorr->CD_prot_Eph_corr(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
        }
        if (_is_FD_Prot)
        {

                // // these are Andrey's corrections
                if (_prot_theta_uncorr < 27)
                {
                        // _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
                        if (_prot_mom_tmt < 2.4)
                                _prot_mom_tmt = _prot_mom_uncorr + (0.001046) * pow(_prot_mom_uncorr, 4) +
                                                (-0.010446) * pow(_prot_mom_uncorr, 3) + (0.036945) * pow(_prot_mom_uncorr, 2) +
                                                (-0.055368) * _prot_mom_uncorr + 0.034539;
                        else
                                _prot_mom_tmt = _prot_mom_uncorr + 0.004741;
                }
                else
                {
                        // _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
                        if (_prot_mom_tmt < 2.4)
                                _prot_mom_tmt = _prot_mom_uncorr + (0.005519) * pow(_prot_mom_uncorr, 4) +
                                                (-0.046289) * pow(_prot_mom_uncorr, 3) + (0.137504) * pow(_prot_mom_uncorr, 2) +
                                                (-0.177027) * _prot_mom_uncorr + 0.094555;
                        else
                                _prot_mom_tmt = _prot_mom_uncorr + 0.004899;
                }
        }

        _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

        // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);
        /////// _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected
        /////// _mom_corr_prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss
        /// corrected

        if (_mc)
        {
                // /////////////////// SMEARING PART ////////////////////////////////////////////////////////////////////////////

                _protUnSmear->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected

                //////////////////////////////////////////////////////////////
                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _protUnSmear->P();

                thetaUnSmear = _protUnSmear->Theta() * 180 / PI;

                if (_protUnSmear->Phi() > 0)
                        phiUnSmear = _protUnSmear->Phi() * 180 / PI;
                else if (_protUnSmear->Phi() < 0)
                        phiUnSmear = (_protUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values

                Reaction::SmearingFunc(PROTON, _prot_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _protUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _protUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _protUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _protSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared
                // _prot->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P); // smeared

                proton->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P); // smeared

                _prot.push_back(std::move(proton)); // Add proton to the vector
                _prot_indices.push_back(i);         // Store the index
        }
        // // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
        else
        {
                // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
                if (_is_FD_Prot)
                {
                        fpro = 1.0;
                        //// // fpro = objMomCorr->dppC(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, _data->dc_sec(i), 3) + 1;// pass1 corr
                }
                else
                {
                        fpro = 1.0;
                }

                mom_corr_proton->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);
                /////////     proton->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);

                _prot_mom = mom_corr_proton->P();

                if (mom_corr_proton->Phi() > 0)
                        _prot_phi = mom_corr_proton->Phi() * 180 / PI;
                else if (mom_corr_proton->Phi() < 0)
                        _prot_phi = (mom_corr_proton->Phi() + 2 * PI) * 180 / PI;

                if (_is_CD_Prot)
                {
                        _prot_mom_prime = objMomCorr->CD_prot_Hmom_corr(_prot_mom, _prot_phi);
                }
                if (_is_FD_Prot)
                {
                        _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr(_prot_mom, _sectorProt);
                }

                _px_prime_prot_mom = mom_corr_proton->Px() * (_prot_mom_prime) / (_prot_mom);
                _py_prime_prot_mom = mom_corr_proton->Py() * (_prot_mom_prime) / (_prot_mom);
                _pz_prime_prot_mom = mom_corr_proton->Pz() * (_prot_mom_prime) / (_prot_mom);

                // proton->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);
                proton->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);

                _prot.push_back(std::move(proton)); // Add proton to the vector
                // _mom_corr_prot.push_back(std::move(mom_corr_proton)); // Add proton to the vector
                _prot_indices.push_back(i);

        } // Store the index
}

void Reaction::SetPip(int i)
{
        _numPip++;
        _numPos++;
        _hasPip = true;
        auto pip = std::make_unique<TLorentzVector>();
        auto mom_corr_pip = std::make_unique<TLorentzVector>();
        /////        pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

        _pip_status = abs(_data->status(i));
        _sectorPip = _data->dc_sec(i);
        _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
        // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
        // _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

        _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
        _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
        if (_Energy_loss_uncorr_pip->Phi() > 0)
                _pip_phi_uncorr = _Energy_loss_uncorr_pip->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_pip->Phi() < 0)
                _pip_phi_uncorr = (_Energy_loss_uncorr_pip->Phi() + 2 * PI) * 180 / PI;

        _is_FD_Pip = objMomCorr->is_FD(_pip_status);
        _is_CD_Pip = objMomCorr->is_CD(_pip_status);

        // if (_is_CD_Pip)
        // {
        //         _pip_mom_tmt = _pip_mom_uncorr;
        // }
        // if (_is_FD_Pip)
        // {
        //         // _pim_mom_tmt = _pim_mom_uncorr;
        //         if (_pip_theta_uncorr < 27)
        //         {
        //                 _pip_mom_tmt = _pip_mom_uncorr + 0.0002468543 * _pip_mom_uncorr + 0.00324120;
        //         }
        //         else
        //         {
        //                 _pip_mom_tmt = _pip_mom_uncorr + -0.0004140691 * _pip_mom_uncorr + 0.007524105;
        //         }
        // }

        // eloss used by stefan
        if (_is_CD_Pip)
        {
                // _pip_mom_tmt = _pip_mom_uncorr;
                _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipCD(_pip_mom_uncorr, _pip_theta_uncorr);
        }
        if (_is_FD_Pip)
        {
                // _pip_mom_tmt = _pip_mom_uncorr;
                _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipFD(_pip_mom_uncorr, _pip_theta_uncorr);
        }
        _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
        _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
        _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

        // /////////////////////////////////     SMEARING PART  /////////////////////////////
        if (_mc)
        {
                _pipUnSmear->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _pipUnSmear->P();

                thetaUnSmear = _pipUnSmear->Theta() * 180 / PI;

                if (_pipUnSmear->Phi() > 0)
                        phiUnSmear = _pipUnSmear->Phi() * 180 / PI;
                else if (_pipUnSmear->Phi() < 0)
                        phiUnSmear = (_pipUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values
                Reaction::SmearingFunc(PIP, _pip_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _pipUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _pipUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear = _pipUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _pipSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
                // _pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP); // smeared

                pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP); // smeared
                _pip.push_back(std::move(pip));
                _pip_indices.push_back(i); // Store the index
        }
        else
        {
                if (_is_FD_Pip)
                {
                        // fpip = 1.0;
                        fpip = objMomCorr->dppC(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, _data->dc_sec(i), 1) * 0.9 + 1;
                }
                else
                {
                        fpip = 1.0;
                }

                //// unique ptr method for one pip
                // _pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
                mom_corr_pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
                /////// vector method for many pip
                /////////            pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

                _pip_mom = mom_corr_pip->P();
                if (mom_corr_pip->Phi() > 0)
                        _pip_phi = mom_corr_pip->Phi() * 180 / PI;
                else if (mom_corr_pip->Phi() < 0)
                        _pip_phi = (mom_corr_pip->Phi() + 2 * PI) * 180 / PI;

                if (_is_CD_Pip)
                {
                        _pip_mom_prime = objMomCorr->CD_pip_Hmom_corr(_pip_mom, _pip_phi);
                }
                if (_is_FD_Pip)
                {
                        _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr(_pip_mom, _sectorPip);
                }
                _px_prime_pip_mom = mom_corr_pip->Px() * ((_pip_mom_prime) / (_pip_mom));
                _py_prime_pip_mom = mom_corr_pip->Py() * ((_pip_mom_prime) / (_pip_mom));
                _pz_prime_pip_mom = mom_corr_pip->Pz() * ((_pip_mom_prime) / (_pip_mom));

                pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
                /////// vector method for many pip
                // pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

                _pip.push_back(std::move(pip)); // Add pip to the vector
                // _mom_corr_pip.push_back(std::move(mom_corr_pip)); // Add pip to the vector
                _pip_indices.push_back(i); // Store the index}
        }
}

void Reaction::SetSwappedProton(int i)
{
        ////_swapped_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);  // smeared

        _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
        // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

        _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();

        _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

        if (_Energy_loss_uncorr_prot->Phi() > 0)
                _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_prot->Phi() < 0)
                _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

        _is_FD_Prot = objMomCorr->is_FD(_prot_status);
        _is_CD_Prot = objMomCorr->is_CD(_prot_status);

        // _is_lower_band = mom_corr::is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

        if (_is_CD_Prot)
        {
                _prot_mom_tmt = _prot_mom_uncorr;
        }

        if (_is_FD_Prot)
        {
                _prot_mom_tmt = _prot_mom_uncorr;
        }

        _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

        /// // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected

        ////////////////////////// Mom CORR /////////////////////////////////
        if (!_mc)
        {
                // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
                if (_is_FD_Prot)
                {
                        fpro = 1.0;
                        ////   // fpro = objMomCorr->dppC(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, _data->dc_sec(i), 3) + 1;///pass1
                }
                else
                {
                        fpro = 1.0;
                }
                _swapped_prot->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);
                // Store the index
        }
        /////////////////// SMEARING PART ////////////////////////////////////////////////////////////////////////////
        if (_mc)
        {
                _protUnSmear->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected

                //////////////////////////////////////////////////////////////
                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _protUnSmear->P();

                thetaUnSmear = _protUnSmear->Theta() * 180 / PI;

                if (_protUnSmear->Phi() > 0)
                        phiUnSmear = _protUnSmear->Phi() * 180 / PI;
                else if (_protUnSmear->Phi() < 0)
                        phiUnSmear = (_protUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values

                Reaction::SmearingFunc(PROTON, _prot_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _protUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _protUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _protUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _protSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared

                _swapped_prot->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P); // smeared
                                                                                             // _swapped_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
                                                                                             // std::cout << "   swapped prot E " << _swapped_prot->E() << std::endl;
        }
}
void Reaction::SetSwappedPip(int i)
{
        ///// _swapped_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

        _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

        _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
        _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;

        _is_FD_Pip = objMomCorr->is_FD(_pip_status);
        _is_CD_Pip = objMomCorr->is_CD(_pip_status);

        // eloss used by stefan
        if (_is_CD_Pip)
        {
                _pip_mom_tmt = _pip_mom_uncorr;
                // _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipCD(_pip_mom_uncorr, _pip_theta_uncorr);
        }
        if (_is_FD_Pip)
        {
                _pip_mom_tmt = _pip_mom_uncorr;
                // _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipFD(_pip_mom_uncorr, _pip_theta_uncorr);
        }
        _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
        _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
        _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
        //// _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);  // energy loss corrected

        ////////////////////////// Mom CORR /////////////////////////////////
        if (!_mc)
        {
                if (_is_FD_Pip)
                {
                        // fpip = 1.0;
                        fpip = objMomCorr->dppC(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, _data->dc_sec(i), 1) + 1;
                }
                else
                {
                        fpip = 1.0;
                }

                /////// vector method for many pip
                _swapped_pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
        } // Store the index

        // /////////////////////////////////     SMEARING PART  /////////////////////////////
        if (_mc)
        {
                _pipUnSmear->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP); // energy loss corrected

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _pipUnSmear->P();

                thetaUnSmear = _pipUnSmear->Theta() * 180 / PI;

                if (_pipUnSmear->Phi() > 0)
                        phiUnSmear = _pipUnSmear->Phi() * 180 / PI;
                else if (_pipUnSmear->Phi() < 0)
                        phiUnSmear = (_pipUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values
                Reaction::SmearingFunc(PIP, _pip_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _pipUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _pipUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _pipUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                _swapped_pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP); // smeared

                // _swapped_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
                // std::cout << "   swapped pip E " << _swapped_pip->E() << std::endl;
        }
}

void Reaction::SetPim(int i)
{
        _numPim++;
        _numNeg++;
        _hasPim = true;
        auto pim = std::make_unique<TLorentzVector>();
        auto mom_corr_pim = std::make_unique<TLorentzVector>();

        _pim_status = abs(_data->status(i));
        _sectorPim = _data->dc_sec(i);

        pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
        //// // _pim->SetPxPyPzE(_data->px(i), _data->py(i), _data->pz(i), sqrt(_data->px(i) * _data->px(i) + _data->py(i) * _data->py(i) + _data->pz(i) * _data->pz(i) + MASS_PIM * MASS_PIM));
        // _thetaDC_r1_Pim = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)),
        // _data->dc_r1_z(i)));

        _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
        // _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
        // _mom_corr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

        _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
        _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
        if (_Energy_loss_uncorr_pim->Phi() > 0)
                _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_pim->Phi() < 0)
                _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

        _is_FD_Pim = objMomCorr->is_FD(_pim_status);
        _is_CD_Pim = objMomCorr->is_CD(_pim_status);
        // _is_lower_band = objMomCorr->is_lower_band(_pim_mom_uncorr, _thetaDC_r1_Pim, _pim_status);

        if (_is_CD_Pim)
        {
                _pim_mom_tmt = _pim_mom_uncorr;
        }
        if (_is_FD_Pim)
        {
                // _pim_mom_tmt = _pim_mom_uncorr;
                if (_pim_theta_uncorr < 27)
                {
                        _pim_mom_tmt = _pim_mom_uncorr + 0.00046571 * _pim_mom_uncorr + 0.00322164;
                }
                else
                {
                        if (_pim_mom_uncorr < 1.7)
                                _pim_mom_tmt = _pim_mom_uncorr + (-0.0024313) * pow(_pim_mom_uncorr, 3) +
                                               (0.0094416) * pow(_pim_mom_uncorr, 2) + (-0.01257967) * pow(_pim_mom_uncorr, 1) + 0.0122432;
                        else
                                _pim_mom_tmt = _pim_mom_uncorr + 0.006199071;
                }
        }
        _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
        _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
        _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

        // /////////////////////////////////     SMEARING PART  /////////////////////////////
        if (_mc)
        {
                _pimUnSmear->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _pimUnSmear->P();

                thetaUnSmear = _pimUnSmear->Theta() * 180 / PI;

                if (_pimUnSmear->Phi() > 0)
                        phiUnSmear = _pimUnSmear->Phi() * 180 / PI;
                else if (_pimUnSmear->Phi() < 0)
                        phiUnSmear = (_pimUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values
                Reaction::SmearingFunc(PIM, _pim_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _pimUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _pimUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear = _pimUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _pimSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM);  // smeared
                // _pim->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM); // smeared

                pim->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM); // smeared

                _pim.push_back(std::move(pim)); // Add pim to the vector
                _pim_indices.push_back(i);      // Store the index
        }
        else
        {
                if (_is_FD_Pim)
                {
                        fpim = 1.0;
                        /// // fpim = objMomCorr->dppC(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, _data->dc_sec(i), 2) + 1; ///pass1
                }
                else
                {
                        fpim = 1.0;
                }

                /////// vector method for many pim

                ///////////                 pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E , MASS_PIM);

                mom_corr_pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

                _pim_mom = mom_corr_pim->P();

                if (mom_corr_pim->Phi() > 0)
                        _pim_phi = mom_corr_pim->Phi() * 180 / PI;
                else if (mom_corr_pim->Phi() < 0)
                        _pim_phi = (mom_corr_pim->Phi() + 2 * PI) * 180 / PI;

                if (_is_CD_Pim)
                {
                        _pim_mom_prime = objMomCorr->CD_pim_Hmom_corr(_pim_mom, _pim_phi);
                        // if (_pim_mom > 2.0) std::cout << "  pim mom cd  " << _pim_mom << "   dp val  " << _pim_mom_prime << std::endl;
                }
                if (_is_FD_Pim)
                {
                        _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr(_pim_mom, _sectorPim);
                }

                _px_prime_pim_mom = mom_corr_pim->Px() * ((_pim_mom_prime) / (_pim_mom));
                _py_prime_pim_mom = mom_corr_pim->Py() * ((_pim_mom_prime) / (_pim_mom));
                _pz_prime_pim_mom = mom_corr_pim->Pz() * ((_pim_mom_prime) / (_pim_mom));

                pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
                // /////pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

                _pim.push_back(std::move(pim)); // Add pim to the vector
                // _mom_corr_pim.push_back(std::move(mom_corr_pim)); // Add pim to the vector
                _pim_indices.push_back(i); // Store the index}
        }
}
void Reaction::SetNeutron(int i)
{
        _numNeutral++;
        _hasNeutron = true;
        _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i)
{
        if (_data->pid(i) == NEUTRON)
        {
                SetNeutron(i);
        }
        else
        {
                _numOther++;
                _hasOther = true;
                _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
        }
}

/////////////////// new added ////////////////
void Reaction::CalcMissMassPim(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto mm_mpim = std::make_unique<TLorentzVector>();

        *mm_mpim += (*_gamma + *_target);
        *mm_mpim -= prot;
        *mm_mpim -= pip;

        _MM_mPim = mm_mpim->M();
        _MM2_mPim = mm_mpim->M2();
}
/////////////////// new added ////////////////
void Reaction::CalcMissMassPimSwapped()
{
        auto mm_mpim_swapped = std::make_unique<TLorentzVector>();

        *mm_mpim_swapped += (*_gamma + *_target);
        *mm_mpim_swapped -= *_swapped_prot;
        *mm_mpim_swapped -= *_swapped_pip;
        // std::cout << "   swapped prot E " << _swapped_prot->E() << std::endl;
        // std::cout << "   swapped pip E " << _swapped_pip->E() << std::endl;

        _MM2_mPim_swapped = mm_mpim_swapped->M2();
}

void Reaction::CalcMissMassExcl(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim)
//     void Reaction::CalcMissMass()
{
        auto mm_mpim = std::make_unique<TLorentzVector>();
        auto mm_mpip = std::make_unique<TLorentzVector>();
        auto mm_mprot = std::make_unique<TLorentzVector>();
        auto mm_excl = std::make_unique<TLorentzVector>();

        if (TwoPion_exclusive())
        {

                *mm_excl += (*_gamma + *_target);
                *mm_excl -= prot;
                *mm_excl -= pip;
                *mm_excl -= pim;
                _MM_exclusive = mm_excl->M();
                _MM2_exclusive = mm_excl->M2();
                _excl_Energy = mm_excl->E();
                // }

                // //         if (TwoPion_missingPim())
                // // {
                *mm_mpim += (*_gamma + *_target);
                *mm_mpim -= prot;
                *mm_mpim -= pip;

                _MM_mPim = mm_mpim->M();
                _MM2_mPim = mm_mpim->M2();
                // // }
                // // if (TwoPion_missingPip())
                // // {
                *mm_mpip += (*_gamma + *_target);
                *mm_mpip -= prot;
                *mm_mpip -= pim;
                _MM2_mpip = mm_mpip->M2();
                // }
                // if (TwoPion_missingProt())
                // {
                *mm_mprot += (*_gamma + *_target);
                *mm_mprot -= pip;
                *mm_mprot -= pim;
                _MM2_mprot = mm_mprot->M2();
                // }
        }
}

float Reaction::MM_mPim()
{
        // if (_MM_mPim != _MM_mPim)
        //         CalcMissMass(*_prot[0], *_pip[0]); // This is just a default case for first proton/pion pair
        return _MM_mPim;
}

float Reaction::MM2_mPim()
{
        // if (_MM2_mPim != _MM2_mPim)
        //         CalcMissMass(*_prot[0], *_pip[0]); // This is just a default case for first proton/pion pair
        return _MM2_mPim;
}

float Reaction::MM2_mPim_swapped()
{
        if (_MM2_mPim_swapped != _MM2_mPim_swapped)
                CalcMissMassPimSwapped();
        return _MM2_mPim_swapped;
}
float Reaction::MM2_exclusive()
{
        // if (_MM2_exclusive != _MM2_exclusive)
        // CalcMissMass();
        return _MM2_exclusive;
}
float Reaction::MM_exclusive()
{
        // if (_MM_exclusive != _MM_exclusive)
        //         CalcMissMass();
        return _MM_exclusive;
}
float Reaction::MM2_mpip()
{
        // if (_MM2_mpip != _MM2_mpip)
        //         CalcMissMass();
        return _MM2_mpip;
}
float Reaction::MM2_mprot()
{
        // if (_MM2_mprot != _MM2_mprot)
        //         CalcMissMass();
        return _MM2_mprot;
}
float Reaction::Energy_excl()
{
        // if (_excl_Energy != _excl_Energy)
        //         CalcMissMass();
        //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
        //  if (_x_mu_E > 0)
        return _excl_Energy;
        // else
        // return NAN;
}

float Reaction::elec_momentum()
{
        if (TwoPion_missingPim())
                return _elec->P();
        else
                return NAN;
}

float Reaction::prot_momentum(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.P();
        // return _prot->P();
        else
                return NAN;
}

float Reaction::pip_momentum(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.P();
        // return _pip->P();
        else
                return NAN;
}
float Reaction::pim_momentum(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->P();
        }
        else
                return NAN;
}

float Reaction::pim_E(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->E();
        }
        else
                return NAN;
}
float Reaction::pim_momentum_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
                return pim.P();
        else
                return NAN;
}
float Reaction::pim_E_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
                return pim.E();
        else
                return NAN;
}
float Reaction::theta_elec()
{ /// lab theta mattrai hunchha electron ko case ma
        // if (TwoPion_missingPim())
        if (_elec->Theta() > -500)
                return _elec->Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::Phi_elec()
{ /// lab theta mattrai hunchha electron ko case ma
        // if (TwoPion_missingPim()) {
        //         if ((_elec->Phi() * 180.0 / PI) > -150)
        //                 return _elec->Phi() * 180.0 / PI;
        //         else
        //                 return ((_elec->Phi() * 180.0 / PI) + 360);
        // }
        // else
        //         return NAN;

        if (_elec->Phi() > -500)
        {
                if (_elec->Phi() > 0)
                        return _elec->Phi() * 180 / PI;
                else if (_elec->Phi() < 0)
                        return (_elec->Phi() + 2 * PI) * 180 / PI;
                else
                        return NAN;
        }
        else
                return NAN;
}
float Reaction::prot_theta_lab(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::pip_theta_lab(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::pim_theta_lab(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->Theta() * 180.0 / PI;
        }
        else
                return NAN;
}
float Reaction::pim_theta_lab_measured(const TLorentzVector &pim)
{ /////////////////////////////////////work here

        if (TwoPion_exclusive())
                return pim.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::prot_Phi_lab(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        // {
        if (_hasP)
        {
                // if (prot.Phi() > 0)
                return prot.Phi() * 180 / PI;
                // else if (prot.Phi() < 0)
                // return (prot.Phi() + 2 * PI) * 180 / PI;
                // else return NAN;
        }
        else
                return NAN;
}
float Reaction::prot_momT(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.Perp();
        else
                return NAN;
}
float Reaction::pip_momT(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.Perp();
        else
                return NAN;
}
float Reaction::pim_momT(const TLorentzVector &pim)
{
        // if (TwoPion_missingPim())
        return pim.Perp();
        // else
        //         return NAN;
}
float Reaction::pip_Phi_lab(const TLorentzVector &pip)
{

        if (_hasPip)
        {
                // if (pip.Phi() > 0)
                return pip.Phi() * 180 / PI;
                // else if (pip.Phi() <= 0)
                //         return (pip.Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}
float Reaction::pim_Phi_lab(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto missingpim_ = std::make_unique<TLorentzVector>();
        *missingpim_ += *_gamma + *_target - prot - pip;
        if (TwoPion_missingPim())
        {
                // if (missingpim_->Phi() > 0)
                return missingpim_->Phi() * 180 / PI;
                // else if (missingpim_->Phi() <= 0)
                //         return (missingpim_->Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}

float Reaction::pim_Phi_lab_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
        {
                // if (pim.Phi() > 0)
                return pim.Phi() * 180 / PI;
                // else if (pim.Phi() <= 0)
                //         return (pim.Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}

void Reaction::EffCorrFactor(const TLorentzVector &prot, const TLorentzVector &pip)
{
        _is_eff_corrected = true;

        _pr_p = prot.P();

        _pr_th = prot.Theta() * 180 / PI;

        // if (prot.Phi() >= 0)
        _pr_ph_eff = prot.Phi() * 180 / PI;
        // else if (prot.Phi() < 0)
        //         _pr_ph_eff = (prot.Phi() + 2 * PI) * 180 / PI;

        _pip_p = pip.P();

        _pip_th = pip.Theta() * 180 / PI;

        if (pip.Phi() >= 0)
                _pip_ph_eff = pip.Phi() * 180 / PI;
        else if (pip.Phi() < 0)
                _pip_ph_eff = (pip.Phi() + 2 * PI) * 180 / PI;
        // std::cout << "  pip_phi 1.  :  " << _pip_ph_eff << std::endl;

        // _eff_corr_fact_Excl = objEffCorr->EFF_CORR_FACT(_pr_p, _pr_th, _pr_ph_eff, _pip_p, _pip_th, _pip_ph_eff, _pim_p, _pim_th, _pim_ph_eff);
        _eff_corr_fact_mPim = objEffCorr->EFF_CORR_FACT1(_pr_p, _pr_th, _pr_ph_eff, _pip_p, _pip_th, _pip_ph_eff);

        // if (_eff_corr_fact_mPim == 1.0 )
        // {
        //         std::cout << "  pr_p :  " << _pr_p << "  pr_th :  " << _pr_th << "  pr_phi :  " << _pr_ph_eff << std::endl;
        //         std::cout << "  pip_p :  " << _pip_p << "  pip_th :  " << _pip_th << "  pip_phi 2 :  " << _pip_ph_eff << std::endl;
        //         std::cout << "  pim_p :  " << _pim_p << "  pim_th :  " << _pim_th << "  pim_phi 2 :  " << _pim_ph_eff << std::endl;

        // std::cout << "  final eff corr factor is  :  " << _eff_corr_fact_mPim << std::endl;
        //         // std::cout << "  pip_fact :  " << objEffCorr->PIP_EFF_CORR_FACT(_pip_p, _pip_th, _pip_ph_eff) << std::endl;
        //         // std::cout << "  pr_fact * pip_fact :  " << objEffCorr->PROT_EFF_CORR_FACT(_pr_p, _pr_th, _pr_ph_eff) * objEffCorr->PIP_EFF_CORR_FACT(_pip_p, _pip_th, _pip_ph_eff) << std::endl;
        // }

        // std::cout << " the factor is:  " << 1 / _eff_corr_fact_mPim << std::endl;

        // return (1 / _eff_corr_fact_mPim);
}

float Reaction::weight()
{
        // std::cout << " the factor is:  " << 1 / _eff_corr_fact_mPim << std::endl;

        // if (_mc)
        //         return _data->mc_weight();
        // else
        return (1 / _eff_corr_fact_mPim);
        // return 1.0;
};

// // Boost the particles to the center-of-mass system
// // void Reaction::boost(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim)
void Reaction::boost(const TLorentzVector &prot, const TLorentzVector &pip)
{
        _is_boosted = true;
        // std::unique_ptr<TLorentzVector> pim = std::make_unique<TLorentzVector>(*_gamma + *_target - prot - pip);
        const TLorentzVector pim = *_gamma + *_target - prot - pip;

        // Boost all particles to the center of mass frame
        _boosted_gamma = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(*_gamma, *_gamma, *_elec, _Q2));
        _boosted_prot = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(prot, *_gamma, *_elec, _Q2));
        _boosted_pip = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pip, *_gamma, *_elec, _Q2));
        _boosted_pim = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pim, *_gamma, *_elec, _Q2));
}
// // // Calculate invariant masse
float Reaction::inv_Ppip()
{
        return boost_cms::calculateInvariantMass(*_boosted_prot, *_boosted_pip);
}
float Reaction::inv_Ppim()
{
        return boost_cms::calculateInvariantMass(*_boosted_prot, *_boosted_pim);
}
float Reaction::inv_pip_pim()
{
        return boost_cms::calculateInvariantMass(*_boosted_pip, *_boosted_pim);
}

// float Reaction::w_P2pi_rec()
// {
//         if (_W_P2pi != _W_P2pi)
//                 W_2pi_P();
//         return _W_P2pi;
// }

// //////////////
float Reaction::prot_theta()
{
        return boost_cms::calculateTheta(*_boosted_prot);
}
float Reaction::pip_theta()
{
        return boost_cms::calculateTheta(*_boosted_pip);
}
float Reaction::pim_theta()
{
        return boost_cms::calculateTheta(*_boosted_pim);
}

// //////////////
float Reaction::gamma_Phi()
{
        return boost_cms::calculatePhi(*_boosted_gamma);
}

float Reaction::prot_Phi()
{
        return boost_cms::calculatePhi(*_boosted_prot);
}

float Reaction::pip_Phi()
{
        return boost_cms::calculatePhi(*_boosted_pip);
}
float Reaction::pim_Phi()
{
        return boost_cms::calculatePhi(*_boosted_pim);
}

// //////////////
float Reaction::alpha_ppip_pipim()
{
        return boost_cms::calculateAlpha(
            _boosted_pim->Vect().Unit(),
            _boosted_pip->Vect().Unit(),
            _boosted_pim->Vect());
}

float Reaction::alpha_pippim_pipf()
{
        return boost_cms::calculateAlpha(
            _boosted_prot->Vect().Unit(),
            _boosted_pip->Vect().Unit(),
            _boosted_prot->Vect());
}

float Reaction::alpha_ppim_pipip()
{
        return boost_cms::calculateAlpha(
            _boosted_pip->Vect().Unit(),
            _boosted_pim->Vect().Unit(),
            _boosted_pip->Vect());
}

///////////////////////// Trying boost for swapped particle /////////////////////////
void Reaction::boost_swapped()
{
        _is_boosted_swapped = true;

        _boosted_prot_swapped = std::make_unique<TLorentzVector>(*_swapped_prot);
        _boosted_pip_swapped = std::make_unique<TLorentzVector>(*_swapped_pip);
        _boosted_pim_swapped = std::make_unique<TLorentzVector>(*_gamma + *_target - *_swapped_prot - *_swapped_pip);
        _boosted_gamma_swapped = std::make_unique<TLorentzVector>(*_gamma);

        TVector3 uz = _boosted_gamma_swapped->Vect().Unit();         // uit vector along virtual photon
        TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit(); // unit vector along e cross e'

        TRotation rot;
        ux.Rotate(3. * PI / 2, uz); // rotating ux by 3pi/2 with uz as axis of roration

        rot.SetZAxis(uz, ux).Invert(); // setting TRotation rot

        _boosted_gamma_swapped->Transform(rot);
        _boosted_prot_swapped->Transform(rot);
        _boosted_pip_swapped->Transform(rot);
        _boosted_pim_swapped->Transform(rot);

        // note beta is calculated only after transforming gamma

        float_t beta_1 =
            ((sqrt(_boosted_gamma_swapped->E() * _boosted_gamma_swapped->E() + _Q2)) / (_boosted_gamma_swapped->E() + MASS_P));

        _boosted_prot_swapped->Boost(0, 0, -beta_1);
        _boosted_pip_swapped->Boost(0, 0, -beta_1);
        _boosted_pim_swapped->Boost(0, 0, -beta_1);
        _boosted_gamma_swapped->Boost(0, 0, -beta_1);
        // std::cout << "   boosted swapped prot E  " << _boosted_prot_swapped->E() << std::endl;
        // std::cout << "   boosted pip swapped  E  " << _boosted_pip_swapped->E() << std::endl;
}
// // // // Calculate invariant masse
void Reaction::invMassPpim_swapped()
{

        if (!_is_boosted_swapped)
                boost_swapped();
        auto inv_Ppim = std::make_unique<TLorentzVector>();
        *inv_Ppim += *_boosted_prot_swapped;
        *inv_Ppim += *_boosted_pim_swapped;

        // *inv_Ppim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
        // if (TwoPion_missingPim())
        _inv_Ppim_swapped = inv_Ppim->M();
}
void Reaction::invMasspippim_swapped()
{
        if (!_is_boosted_swapped)
                boost_swapped();
        auto inv_pip_pim = std::make_unique<TLorentzVector>();
        *inv_pip_pim += *_boosted_pip_swapped;
        *inv_pip_pim += *_boosted_pim_swapped;

        // *inv_pip_pim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
        // if (TwoPion_missingPim())
        _inv_pip_pim_swapped = inv_pip_pim->M();
}

void Reaction::invMassPpip_swapped()
{
        if (!_is_boosted_swapped)
                boost_swapped();
        auto inv_Ppip = std::make_unique<TLorentzVector>();
        *inv_Ppip += *_boosted_prot_swapped;
        *inv_Ppip += *_boosted_pip_swapped;

        // if (TwoPion_missingPim())
        _inv_Ppip_swapped = inv_Ppip->M();
}

float Reaction::inv_Ppip_swapped()
{
        // return boost_cms::calculateInvariantMass(*_boosted_prot, *_boosted_pip);

        if (_inv_Ppip_swapped != _inv_Ppip_swapped)
                invMassPpip_swapped();
        // std::cout << "   _inv_Ppip_swapped  " << _inv_Ppip_swapped << std::endl;
        return _inv_Ppip_swapped;
}
float Reaction::inv_Ppim_swapped()
{
        if (_inv_Ppim_swapped != _inv_Ppim_swapped)
                invMassPpim_swapped();
        return _inv_Ppim_swapped;
}
float Reaction::inv_pip_pim_swapped()
{
        if (_inv_pip_pim_swapped != _inv_pip_pim_swapped)
                invMasspippim_swapped();
        return _inv_pip_pim_swapped;
}

// //////////////
float Reaction::alpha_ppip_pipim_swapped()
{
        return boost_cms::calculateAlpha(
            _boosted_pim_swapped->Vect().Unit(),
            _boosted_pip_swapped->Vect().Unit(),
            _boosted_pim_swapped->Vect());
}

float Reaction::alpha_pippim_pipf_swapped()
{
        return boost_cms::calculateAlpha(
            _boosted_prot_swapped->Vect().Unit(),
            _boosted_pip_swapped->Vect().Unit(),
            _boosted_prot_swapped->Vect());
}

float Reaction::alpha_ppim_pipip_swapped()
{
        return boost_cms::calculateAlpha(
            _boosted_pip_swapped->Vect().Unit(),
            _boosted_pim_swapped->Vect().Unit(),
            _boosted_pip_swapped->Vect());
}

// //////////////
float Reaction::prot_theta_swapped()
{
        return _boosted_prot_swapped->Theta() * (180 / PI);
}
float Reaction::pip_theta_swapped()
{
        return _boosted_pip_swapped->Theta() * (180 / PI);
}
float Reaction::pim_theta_swapped()
{
        return _boosted_pim_swapped->Theta() * (180 / PI);
}

// // //////////////
// float Reaction::gamma_Phi()
// {
//         return boost_cms::calculatePhi(*_boosted_gamma);
// }

// float Reaction::prot_Phi()
// {
//         return boost_cms::calculatePhi(*_boosted_prot);
// }

// float Reaction::pip_Phi()
// {
//         return boost_cms::calculatePhi(*_boosted_pip);
// }
// float Reaction::pim_Phi()
// {
//         return boost_cms::calculatePhi(*_boosted_pim);
// }

// //////////////
// float Reaction::alpha_ppip_pipim()
// {
//         return boost_cms::calculateAlpha(
//             _boosted_pim->Vect().Unit(),
//             _boosted_pip->Vect().Unit(),
//             _boosted_pim->Vect());
// }

// float Reaction::alpha_pippim_pipf()
// {
//         return boost_cms::calculateAlpha(
//             _boosted_prot->Vect().Unit(),
//             _boosted_pip->Vect().Unit(),
//             _boosted_prot->Vect());
// }

// float Reaction::alpha_ppim_pipip()
// {
//         return boost_cms::calculateAlpha(
//             _boosted_pip->Vect().Unit(),
//             _boosted_pim->Vect().Unit(),
//             _boosted_pip->Vect());
// }

///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////
///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////
///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////
///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////
///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////
///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////

MCReaction::MCReaction(const std::shared_ptr<Branches12> &data,
                       float beam_enrgy)
{
        _data = data;
        if (!_data->mc())
                _data->mc_branches();
        _beam_mc = std::make_unique<TLorentzVector>();
        _beam_energy = 10.6041;
        // _beam_energy = atof(getenv("CLAS12_E"));
        _weight_mc = _data->mc_weight();
        _beam_mc->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);
        _gamma_mc = std::make_unique<TLorentzVector>();
        _target = std::make_unique<TLorentzVector>();
        // _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
        _target->SetXYZM(0.0, 0.0, 0.0, MASS_P);
        _elec_mc = std::make_unique<TLorentzVector>();
        this->SetMCElec();
        // _prot_mc = std::make_unique<TLorentzVector>();
        // _pip_mc = std::make_unique<TLorentzVector>();
        // _pim_mc = std::make_unique<TLorentzVector>();
        _other_mc = std::make_unique<TLorentzVector>();
        //_neutron = std::make_unique<TLorentzVector>();

        _prot_mc = std::vector<std::unique_ptr<TLorentzVector>>(); // Initialize _protons as an empty vector
        _pip_mc = std::vector<std::unique_ptr<TLorentzVector>>();
        _pim_mc = std::vector<std::unique_ptr<TLorentzVector>>();

        _prot_mc_indices = std::vector<int>();
        _pip_mc_indices = std::vector<int>();
        _pim_mc_indices = std::vector<int>();
}
// Reaction::~Reaction() {} // why this is not here
void MCReaction::SetMCElec()
{
        //  _hasE = true;  //??
        _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);
        // _elec_mc->SetPxPyPzE(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), sqrt(_data->mc_px(0) * _data->mc_px(0)
        // + _data->mc_py(0) * _data->mc_py(0) + _data->mc_pz(0) * _data->mc_pz(0)  + MASS_E * MASS_E));

        *_gamma_mc += *_beam_mc - *_elec_mc;

        // Can calculate W and Q2 here
        _W_mc = physics::W_calc(*_beam_mc, *_elec_mc);
        _Q2_mc = physics::Q2_calc(*_beam_mc, *_elec_mc);
}

void MCReaction::SetMCProton(int i)
{
        auto proton_mc = std::make_unique<TLorentzVector>();
        proton_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P);
        _prot_mc.push_back(std::move(proton_mc));
        _prot_mc_indices.push_back(i);
        // std::cout << "mc_prot  index  " << i << std::endl;
}

void MCReaction::SetMCPip(int i)
{

        auto pionp_mc = std::make_unique<TLorentzVector>();
        pionp_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP);
        _pip_mc.push_back(std::move(pionp_mc));
        _pip_mc_indices.push_back(i);
        // std::cout << "mc_pip  index  " << i << std::endl;
}

void MCReaction::SetMCPim(int i)
{
        auto pionm_mc = std::make_unique<TLorentzVector>();
        pionm_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM);
        _pim_mc.push_back(std::move(pionm_mc));
        _pim_mc_indices.push_back(i);
        // std::cout << "mc_pim  index  " << i << std::endl;
}
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
// }

float MCReaction::pim_mom_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->P();
}
float MCReaction::pip_mom_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->P();
}
float MCReaction::prot_mom_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->P();
}

float MCReaction::pim_momX_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->Px();
}
float MCReaction::pip_momX_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->Px();
}
float MCReaction::prot_momX_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->Px();
}

float MCReaction::pim_momY_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->Py();
}
float MCReaction::pip_momY_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->Py();
}
float MCReaction::prot_momY_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->Py();
}

float MCReaction::pim_momZ_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->Pz();
}
float MCReaction::pip_momZ_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->Pz();
}
float MCReaction::prot_momZ_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->Pz();
}

float MCReaction::pim_theta_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->Theta() * 180 / PI;
}
float MCReaction::pip_theta_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->Theta() * 180 / PI;
}
float MCReaction::prot_theta_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->Theta() * 180 / PI;
}

float MCReaction::pim_phi_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        if (pim_mc->Phi() >= 0)
                return (pim_mc->Phi() * 180 / PI);
        else if (pim_mc->Phi() < 0)
                return ((pim_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}
float MCReaction::pip_phi_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        if (pip_mc->Phi() >= 0)
                return (pip_mc->Phi() * 180 / PI);
        else if (pip_mc->Phi() < 0)
                return ((pip_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}
float MCReaction::prot_phi_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        if (prot_mc->Phi() >= 0)
                return (prot_mc->Phi() * 180 / PI);
        else if (prot_mc->Phi() < 0)
                return ((prot_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}
////////////////////////////////  BOOST TO CM SYSTEM ///////////////////////////////
////////////////////////////////  BOOST TO CM SYSTEM ///////////////////////////////
////////////////////////////////  BOOST TO CM SYSTEM ///////////////////////////////
////////////////////////////////  BOOST TO CM SYSTEM ///////////////////////////////

void MCReaction::boost_mc(const TLorentzVector &prot_mc, const TLorentzVector &pip_mc, const TLorentzVector &pim_mc)
{
        _is_boosted_mc = true;

        // Boost all particles to the center of mass frame
        _boosted_gamma_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(*_gamma_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
        _boosted_prot_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(prot_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
        _boosted_pip_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pip_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
        _boosted_pim_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pim_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
}
// // // Calculate invariant masse
float MCReaction::MCinv_Ppip()
{
        return boost_cms::calculateInvariantMass(*_boosted_prot_mc, *_boosted_pip_mc);
}
float MCReaction::MCinv_Ppim()
{
        return boost_cms::calculateInvariantMass(*_boosted_prot_mc, *_boosted_pim_mc);
}
float MCReaction::MCinv_pip_pim()
{
        return boost_cms::calculateInvariantMass(*_boosted_pip_mc, *_boosted_pim_mc);
}

// //////////////
float MCReaction::MCprot_theta_thrown()
{
        return boost_cms::calculateTheta(*_boosted_prot_mc);
}
float MCReaction::MCpip_theta_thrown()
{
        return boost_cms::calculateTheta(*_boosted_pip_mc);
}
float MCReaction::MCpim_theta_thrown()
{
        return boost_cms::calculateTheta(*_boosted_pim_mc);
}

// //////////////
float MCReaction::MCgamma_Phi_thrown()
{
        return boost_cms::calculatePhi(*_boosted_gamma_mc);
}

float MCReaction::MCprot_Phi_thrown()
{
        return boost_cms::calculatePhi(*_boosted_prot_mc);
}

float MCReaction::MCpip_Phi_thrown()
{
        return boost_cms::calculatePhi(*_boosted_pip_mc);
}
float MCReaction::MCpim_Phi_thrown()
{
        return boost_cms::calculatePhi(*_boosted_pim_mc);
}

// //////////////
float MCReaction::MCalpha_ppip_pipim_thrown()
{
        return boost_cms::calculateAlpha(
            _boosted_pim_mc->Vect().Unit(),
            _boosted_pip_mc->Vect().Unit(),
            _boosted_pim_mc->Vect());
}

float MCReaction::MCalpha_pippim_pipf_thrown()
{
        return boost_cms::calculateAlpha(
            _boosted_prot_mc->Vect().Unit(),
            _boosted_pip_mc->Vect().Unit(),
            _boosted_prot_mc->Vect());
}

float MCReaction::MCalpha_ppim_pipip_thrown()
{
        return boost_cms::calculateAlpha(
            _boosted_pip_mc->Vect().Unit(),
            _boosted_pim_mc->Vect().Unit(),
            _boosted_pip_mc->Vect());
}
