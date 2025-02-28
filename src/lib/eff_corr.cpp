#include "eff_corr.hpp"
#include <iostream>
EffCorr::~EffCorr() {} /// Since your class does not involve dynamic memory allocation (e.g., using new),
// there's no need to define a custom destructor (e.g., ~EffCorr() {}). The compiler-generated destructor
//  will work fine in this case, and you can omit the empty destructor.

////////////////////////////////// PROTON EFF FACTORS

float EffCorr::PROT_EFF_CORR_FACT(float mom_mes, float theta_mes, float phi_mes)
{
    float theta = NAN;
    if (theta_mes < 42)
        theta = theta_mes / ((ProtThMap[0]) * pow(theta_mes, 2) + (ProtThMap[1]) * theta_mes + (ProtThMap[2]));
    if (theta_mes >= 42)
        theta = theta_mes / 1.01;
    // std::cout << " \n theta prot miss " << theta << " theta p mes  " << theta_mes << "  dividing by " << ((ProtThMap[0]) * pow(theta_mes, 2) + (ProtThMap[1]) * theta_mes + (ProtThMap[2])) << std::endl;

    if (theta >= 0 && theta < 37)
    {
        float phi = phi_mes / ((ProtPhiMap[0][0]) * pow(phi_mes, 2) + (ProtPhiMap[0][1]) * phi_mes + (ProtPhiMap[0][2]));

        for (int p = 0; p < FD_SEC; p++)
        {

            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                /////////// Map the measured momentum and missing momentum /////////////
                float mom = mom_mes / ((FDProtMap[p][0]) * pow(mom_mes, 2) + (FDProtMap[p][1]) * mom_mes + (FDProtMap[p][2]));
                // std::cout << " theta p mes  " << theta_mes << "   phi_mes  " << phi_mes << "   phi_miss  " << phi << "  mom_mes is : " << mom_mes << "  mom_miss is : " << mom << "  division fact  " << (FDProtMap[p][0]) * pow(mom_mes, 2) + (FDProtMap[p][1]) * mom_mes + (FDProtMap[p][2]) << std::endl;

                // std::cout << " theta p mes  " << theta_mes << "   phi_mes  " << phi_mes << "  mom_mes is : " << mom_mes << "  mom_miss is : " << mom << "  division fact  " << (FDProtMap[p][0]) * pow(mom_mes, 2) + (FDProtMap[p][1]) * mom_mes + (FDProtMap[p][2]) << std::endl;

                return (FDProtCoef[p][0]) * pow(mom, 2) + (FDProtCoef[p][1]) * mom + (FDProtCoef[p][2]);
            }
        }
    }
    else if (theta >= 37 && theta <= 180)
    {
        float phi = phi_mes; // / ((ProtPhiMap[1][0]) * pow(phi_mes, 2) + (ProtPhiMap[1][1]) * phi_mes + (ProtPhiMap[1][2]));

        for (int p = 0; p < CD_SEC; p++)
        {

            // std::cout << " phi : " << phi << "  ProtPhiMap[1][0] "
            //           << ProtPhiMap[1][0] << "  1st coeff =  " << (CDProtCoef[p][0]) << std::endl;
            if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((CDProtMap[p][0]) * pow(mom_mes, 2) + (CDProtMap[p][1]) * mom_mes + (CDProtMap[p][2]));
                // std::cout << " theta p mes  " << theta_mes << "   phi_mes  " << phi_mes << "   phi_miss  " << phi << "  mom_mes is : " << mom_mes << "  mom_miss is : " << mom << "  division fact  " << (FDProtMap[p][0]) * pow(mom_mes, 2) + (FDProtMap[p][1]) * mom_mes + (FDProtMap[p][2]) << std::endl;

                // // std::cout << " Prot mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " (ProtThMap[0]  " << ProtThMap[0] << " phi : " << phi << "  ProtPhiMap[1][0] "
                // //           << ProtPhiMap[1][0] << std::endl;
                // std::cout << " cd sec is : " << p + 1 << " CD Prot eff corr val =  "
                //           << ((CDProtCoef[p][0]) * pow(mom, 2) + (CDProtCoef[p][1]) * mom + (CDProtCoef[p][2])) << std::endl;

                return (CDProtCoef[p][0]) * pow(mom, 2) + (CDProtCoef[p][1]) * mom + (CDProtCoef[p][2]);
            }
        }
    }

    return 1.0; // to test if you are getting good result or not
}

////////////////////////////////// PIP EFF FACTORS

float EffCorr::PIP_EFF_CORR_FACT(float mom_mes, float theta_mes, float phi_mes)
{
    float theta = NAN;
    if (theta_mes < 40)
        theta = theta_mes / ((PipThMap[0]) * pow(theta_mes, 3) + (PipThMap[1]) * pow(theta_mes, 2) + (PipThMap[2] * pow(theta_mes, 1)) + PipThMap[3]);
    if (theta_mes >= 40)
        theta = theta_mes / 1.02347;
    // std::cout << " theta pip miss " << theta << " theta pip mes  " << theta_mes << "  dividing by " << (PipThMap[0]) * pow(theta_mes, 3) + (PipThMap[1]) * pow(theta_mes, 2) + (PipThMap[2] * pow(theta_mes, 1)) + PipThMap[3] << std::endl;

    if (theta >= 0 && theta < 37)
    {
        float phi = phi_mes / ((PipPhiMap[0][0]) * pow(phi_mes, 2) + (PipPhiMap[0][1]) * phi_mes + (PipPhiMap[0][2]));

        for (int p = 0; p < FD_SEC; p++)
        {

            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((FDPipMap[p][0]) * pow(mom_mes, 2) + (FDPipMap[p][1]) * mom_mes + (FDPipMap[p][2]));

                if (mom < 1.6)
                {
                    return (FDPipCoef[p][0]) * pow(mom, 2) + (FDPipCoef[p][1]) * pow(mom, 1) + (FDPipCoef[p][2]);
                    // std::cout << " fd sec is : " << p + 1 << " FD pip eff corr val =  "
                    //           << (FDPipCoef[p][0]) * pow(mom, 2) + (FDPipCoef[p][1]) * pow(mom, 1) + (FDPipCoef[p][2]) << std::endl;
                }

                else if (mom >= 1.6)
                {
                    return (FDPipCoef1[p][0]) * pow(mom, 1) + (FDPipCoef1[p][1]);
                    // std::cout << " fd sec is : " << p + 1 << " FD pip  eff corr val =  "
                    //           << (FDPipCoef1[p][0]) * pow(mom, 1) + (FDPipCoef1[p][1]) << std::endl;
                }
            }
        }
    }
    else if (theta >= 37 && theta <= 180)
    {
        float phi = phi_mes / ((PipPhiMap[1][0]) * pow(phi_mes, 2) + (PipPhiMap[1][1]) * phi_mes + (PipPhiMap[1][2]));

        for (int p = 0; p < CD_SEC; p++)
        {

            if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((CDPipMap[p][0]) * pow(mom_mes, 2) + (CDPipMap[p][1]) * mom_mes + (CDPipMap[p][2]));

                // std::cout << " Pip mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (CDPipCoef[p][0]) << std::endl;
                // std::cout << " cd is : " << p + 1 << "  eff corr val =  "
                //           << ((CDPipCoef[p][0]) * pow(mom, 2) + (CDPipCoef[p][1]) * mom + (CDPipCoef[p][2])) << std::endl;

                return (CDPipCoef[p][0]) * pow(mom, 2) + (CDPipCoef[p][1]) * mom + (CDPipCoef[p][2]);
            }
        }
    }
    return 1.0; // to test if you are getting good result or not
}

// ////////////////////////////////// PIM EFF FACTORS

// float EffCorr::PIM_EFF_CORR_FACT(float mom_mes, float theta_mes, float phi_mes)
// {
//     if (theta >= 0 && theta < 20)
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 float mom = mom_mes / ((FDPimMap[0][p][0]) * pow(mom_mes, 2) + (FDPimMap[0][p][1]) * mom_mes + (FDPimMap[0][p][2]));

//                 return (FDPimCoef[0][p][0]) * pow(mom, 2) + (FDPimCoef[0][p][1]) * mom + (FDPimCoef[0][p][2]);
//             }
//         }

//     else if (theta >= 20 && theta < 40)
//     {

//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 float mom = mom_mes / ((FDPimMap[1][p][0]) * pow(mom_mes, 2) + (FDPimMap[1][p][1]) * mom_mes + (FDPimMap[1][p][2]));

//                 // std::cout << "mom is : " << mom << " theta : " << theta << " phi : " << phi
//                 //           << "  1st coeff =  " << (FDPimCoef[1][p][0]) << std::endl;
//                 // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
//                 //           << ((FDPimCoef[1][p][0]) * pow(mom, 2) + (FDPimCoef[1][p][1]) * mom + (FDPimCoef[1][p][2])) << std::endl;

//                 return (FDPimCoef[1][p][0]) * pow(mom, 2) + (FDPimCoef[1][p][1]) * mom + (FDPimCoef[1][p][2]);
//             }
//         }
//     }
//     else if (theta >= 40 && theta <= 180)
//     {
//         for (int p = 0; p < CD_SEC; p++)
//         {
//             if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
//             {
//                 float mom = mom_mes / ((CDPimMap[p][0]) * pow(mom_mes, 2) + (CDPimMap[p][1]) * mom_mes + (CDPimMap[p][2]));

//                 return (CDPimCoef[p][0]) * pow(mom, 2) + (CDPimCoef[p][1]) * mom + (FDPimCoef[1][p][2]);
//             }
//         }
//     }
//     return 1.0; // to test if you are getting good result or not
// }
// float EffCorr::EFF_CORR_FACT(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip, float mom_pim, float theta_pim, float phi_pim)
// {
//     return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
//             (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)) *
//             (EffCorr::PIM_EFF_CORR_FACT(mom_pim, theta_pim, phi_pim)));
// };
float EffCorr::EFF_CORR_FACT1(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip)
{

    // std::cout
    //     << "  pip theta  = " << theta_pip << "   phi " << phi_pip << "  eff corr fact = " << PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip) << std::endl;
    // std::cout << "  prot theta  = " << theta_p << "   phi " << phi_p << "  eff corr fact = " << PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p) << std::endl;

    return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
            (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)));
};
