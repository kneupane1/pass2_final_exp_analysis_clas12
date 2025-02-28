
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"

class Cuts
{
protected:
    std::shared_ptr<Branches12> _data = nullptr;
    std::shared_ptr<Delta_T> _dt = nullptr;

public:
    Cuts(const std::shared_ptr<Branches12> &data);
    Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt);
    ~Cuts();

    bool ElectronCuts();
    bool FiducialCuts();
    bool IsPip(int i);
    bool IsProton(int i);
    bool IsPim(int i);
};

// class rga_Cuts : public Cuts
// {
// public:
//     rga_Cuts(const std::shared_ptr<Branches12> &data) : Cuts(data)
//     {
//     }
//     rga_Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : Cuts(data, dt){};
// };

class Pass2_Cuts : public Cuts
{
public:
    Pass2_Cuts(const std::shared_ptr<Branches12> &data) : Cuts(data)
    {
    }
    Pass2_Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : Cuts(data, dt) {};
    bool ElectronCuts();
    bool IsPip(int i);
    bool IsProton(int i);
    bool IsPim(int i);
    // bool CC_nphe_cut(double nphe);
    bool CC_nphe_cut();
    bool PCAL_minimum_energy();
    bool EC_sampling_fraction_cut();
    bool EC_hit_position_fiducial_cut_homogeneous();
    bool DC_fiducial_cut_XY(int i, int pid);
    bool DC_z_vertex_cut();
    bool EC_inner_vs_EC_outer();
    bool PCAL_fiducial_cut_X_Y();
    bool PCAL_Ineff_cut_X_Y();
    bool DC_Ineff_cut_X_Y(int i, int pid);

    bool HadronsCuts(int i);
    bool DC_fiducial_cut_theta_phi(int i);
    bool Hadron_Delta_vz_cut(int i);
    bool Hadron_Chi2pid_cut(int i);
    bool CD_fiducial_had(int i);
    // Function to get the momentum range index based on the value of p

    int getMomRangeIndex(double p)
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
