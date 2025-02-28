
#include "mom_corr.hpp"
#include <cstdlib>
#include <ctime>
#include "iostream"

mom_corr::~mom_corr() {}

bool mom_corr::is_FD(int part_status) {
  // if (dc_sec >= 1 && dc_sec <= 6)
  if (part_status > 2000 && part_status <= 4000)
    return true;
  else
    return false;
}

// bool mom_corr::is_AllFD(int part1_status, int part2_status, int part3_status) {
//   // if (dc_sec >= 1 && dc_sec <= 6)
//   if ((part1_status > 2000 && part1_status <= 4000) && (part2_status > 2000 && part2_status <= 4000) &&
//       (part3_status > 2000 && part3_status <= 4000))
//     return true;
//   else
//     return false;
// }

bool mom_corr::is_CD(int part_status) {
  // if (dc_sec < 1 || dc_sec > 6)
  if (part_status > 4000 && part_status <= 6000)
    return true;
  else
    return false;
}

//////////////////// new mom correction start pass2 2024
double mom_corr::dppC(float Px, float Py, float Pz, int sec, int ivec) {
  // auto dppC = [&](float Px, float Py, float Pz, int sec, int ivec) {
  // ivec = 0 --> Electron Corrections
  // ivec = 1 --> Pi+ Corrections
  // ivec = 2 --> Pi- Corrections
  // ivec = 3 --> Proton Corrections

  // Momentum Magnitude
  double pp = sqrt(Px * Px + Py * Py + Pz * Pz);

  // Initializing the correction factor
  double dp = 0;

  // Defining Phi Angle
  double Phi = (180 / 3.1415926) * atan2(Py, Px);

  // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
  if (((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)) {
    Phi += 360;
  }

  // Getting Local Phi Angle
  double PhiLocal = Phi - (sec - 1) * 60;

  // Applying Shift Functions to Phi Angles (local shifted phi = phi)
  double phi = PhiLocal;

  // For Electron Shift
  if (ivec == 0) {
    phi = PhiLocal - 30 / pp;
  }

  // For Pi+ Pion/Proton Shift
  if (ivec == 1 || ivec == 3) {
    phi = PhiLocal + (32 / (pp - 0.05));
  }

  // For Pi- Pion Shift
  if (ivec == 2) {
    phi = PhiLocal - (32 / (pp - 0.05));
  }

  //==========//  PARTICLE = ELECTRON  //==========//

  if (ivec == 0) {
    if (sec == 1) {
      dp = ((1.57e-06) * phi * phi + (5.021e-05) * phi + (-1.74089e-03)) * pp * pp +
           ((-2.192e-05) * phi * phi + (-1.12528e-03) * phi + (0.0146476)) * pp +
           ((8.504e-05) * phi * phi + (2.08012e-03) * phi + (-0.0122501));
    }
    if (sec == 2) {
      dp = ((-3.98e-06) * phi * phi + (1.66e-05) * phi + (-1.55918e-03)) * pp * pp +
           ((2.136e-05) * phi * phi + (-5.7373e-04) * phi + (0.0143591)) * pp +
           ((2.4e-06) * phi * phi + (1.6656e-03) * phi + (-0.0218711));
    }

    if (sec == 3) {
      dp = ((5.57e-06) * phi * phi + (2.3e-07) * phi + (-2.26999e-03)) * pp * pp +
           ((-7.761e-05) * phi * phi + (4.1437e-04) * phi + (0.0152985)) * pp +
           ((2.2542e-04) * phi * phi + (-9.442e-04) * phi + (-0.0231432));
    }

    if (sec == 4) {
      dp = ((3.48e-06) * phi * phi + (2.166e-05) * phi + (-2.29e-04)) * pp * pp +
           ((-2.758e-05) * phi * phi + (7.226e-05) * phi + (-3.38e-03)) * pp +
           ((3.166e-05) * phi * phi + (6.93e-05) * phi + (0.04767));
    }

    if (sec == 5) {
      dp = ((1.19e-06) * phi * phi + (-2.286e-05) * phi + (-1.6332e-04)) * pp * pp +
           ((-1.05e-06) * phi * phi + (7.04e-05) * phi + (-5.0754e-03)) * pp +
           ((-7.22e-06) * phi * phi + (4.1748e-04) * phi + (0.04441));
    }

    if (sec == 6) {
      dp = ((-5.97e-06) * phi * phi + (-3.689e-05) * phi + (5.782e-05)) * pp * pp +
           ((6.573e-05) * phi * phi + (2.1376e-04) * phi + (-9.54576e-03)) * pp +
           ((-1.7732e-04) * phi * phi + (-8.62e-04) * phi + (0.0618975));
    }
  }

  //==========//  PARTICLE = ELECTRON (END)  //==========//

  //==========//  PARTICLE = PI+ PION  //==========//

  if (ivec == 1) {
    if (sec == 1) {
      dp = ((-5.2e-07) * phi * phi + (-1.383e-05) * phi + (4.7179e-04)) * pp * pp +
           ((8.33e-06) * phi * phi + (3.8849e-04) * phi + (-6.81319e-03)) * pp +
           ((-1.645e-05) * phi * phi + (-5.0057e-04) * phi + (1.9902e-02));
    }

    if (sec == 2) {
      dp = ((-1.88e-06) * phi * phi + (3.303e-05) * phi + (1.1331e-03)) * pp * pp +
           ((1.569e-05) * phi * phi + (-3.974e-05) * phi + (-1.25869e-02)) * pp +
           ((-2.903e-05) * phi * phi + (-1.0638e-04) * phi + (2.61529e-02));
    }
    if (sec == 3) {
      dp = ((2.4e-07) * phi * phi + (-1.04e-05) * phi + (7.0864e-04)) * pp * pp +
           ((8.0e-06) * phi * phi + (-5.156e-05) * phi + (-8.12169e-03)) * pp +
           ((-2.42e-05) * phi * phi + (8.928e-05) * phi + (2.13223e-02));
    }
    if (sec == 4) {
      dp = ((-4.0e-08) * phi * phi + (-3.59e-05) * phi + (1.32146e-03)) * pp * pp +
           ((1.023e-05) * phi * phi + (2.2199e-04) * phi + (-1.33043e-02)) * pp +
           ((-2.801e-05) * phi * phi + (-1.576e-04) * phi + (3.27995e-02));
    }
    if (sec == 5) {
      dp = ((2.7e-06) * phi * phi + (5.03e-06) * phi + (1.59668e-03)) * pp * pp +
           ((-1.28e-05) * phi * phi + (-1.99e-06) * phi + (-1.71578e-02)) * pp +
           ((2.091e-05) * phi * phi + (-4.14e-05) * phi + (3.25434e-02));
    }
    if (sec == 6) {
      dp = ((2.13e-06) * phi * phi + (-7.49e-05) * phi + (1.75565e-03)) * pp * pp +
           ((-7.37e-06) * phi * phi + (5.8222e-04) * phi + (-1.27969e-02)) * pp +
           ((4.9e-07) * phi * phi + (-7.2253e-04) * phi + (3.11499e-02));
    }
  }

  //==========//  PARTICLE = PI+ PION (END)  //==========//

  //==========//  PARTICLE = PI- PION  //==========//

  if (ivec == 2) {
    if (sec == 1) {
      dp = ((-4.0192658422317425e-06) * phi * phi - (2.660222128967742e-05) * phi + 0.004774434682983547) * pp * pp;
      dp = dp + ((1.9549520962477972e-05) * phi * phi - 0.0002456062756770577 * phi - 0.03787692408323466) * pp;
      dp = dp + (-2.128953094937459e-05) * phi * phi + 0.0002461708852239913 * phi + 0.08060704449822174 - 0.01;
    }

    if (sec == 2) {
      dp = ((1.193010521758372e-05) * phi * phi - (5.996221756031922e-05) * phi + 0.0009093437955814359) * pp * pp;
      dp = dp + ((-4.89113824430594e-05) * phi * phi + 0.00021676479488147118 * phi - 0.01861892053916726) * pp;
      dp = dp + (4.446394152208071e-05) * phi * phi - (3.6592784167335244e-05) * phi + 0.05498710249944096 - 0.01;
    }

    if (sec == 3) {
      dp = ((-1.6596664895992133e-07) * phi * phi + (6.317189710683516e-05) * phi + 0.0016364212312654086) * pp * pp;
      dp = dp + ((-2.898409777520318e-07) * phi * phi - 0.00014531513577533802 * phi - 0.025456145839203827) * pp;
      dp = dp + (2.6432552410603506e-06) * phi * phi + 0.00018447151306275443 * phi + 0.06442602664627255 - 0.01;
    }

    if (sec == 4) {
      dp = ((2.4035259647558634e-07) * phi * phi - (8.649647351491232e-06) * phi + 0.004558993439848128) * pp * pp;
      dp = dp + ((-5.981498144060984e-06) * phi * phi + 0.00010582131454222416 * phi - 0.033572004651981686) * pp;
      dp = dp + (8.70140266889548e-06) * phi * phi - 0.00020137414379966883 * phi + 0.07258774523336173 - 0.01;
    }

    if (sec == 5) {
      dp = ((2.5817024702834863e-06) * phi * phi + 0.00010132810066914441 * phi + 0.003397314538804711) * pp * pp;
      dp = dp + ((-1.5116941263931812e-05) * phi * phi - 0.00040679799541839254 * phi - 0.028144285760769876) * pp;
      dp = dp + (1.4701931057951464e-05) * phi * phi + 0.0002426350390593454 * phi + 0.06781682510174941 - 0.01;
    }

    if (sec == 6) {
      dp = ((-8.196823669099362e-07) * phi * phi - (5.280412421933636e-05) * phi + 0.0018457238328451137) * pp * pp;
      dp = dp + ((5.2675062282094536e-06) * phi * phi + 0.0001515803461044587 * phi - 0.02294371578470564) * pp;
      dp = dp + (-9.459454671739747e-06) * phi * phi - 0.0002389523716779765 * phi + 0.06428970810739926 - 0.01;
    }
  }

  //==========//  PARTICLE = PI- PION (END)  //==========//

  //==========//  PARTICLE = PROTON  //==========//

  if (ivec == 3) {
    // The following lines should be added up in the order given for the full correction
    // Applying this code as given will give the exact corrections of this analysis
    // These parameters will be combined into a single line at a later point

    if (sec == 1) {
      dp = (5.415e-04) * pp * pp + (-1.0262e-02) * pp + (7.78075e-03);
      dp = dp + ((1.2129e-04) * pp * pp + (1.5373e-04) * pp + (-2.7084e-04));
    }
    if (sec == 2) {
      dp = (-9.5439e-04) * pp * pp + (-2.86273e-03) * pp + (3.38149e-03);
      dp = dp + ((-1.6890e-03) * pp * pp + (4.3744e-03) * pp + (-2.1218e-03));
    }
    if (sec == 3) {
      dp = (-5.5541e-04) * pp * pp + (-7.69739e-03) * pp + (5.7692e-03);
      dp = dp + ((7.6422e-04) * pp * pp + (-1.5425e-03) * pp + (5.4255e-04));
    }
    if (sec == 4) {
      dp = (-1.944e-04) * pp * pp + (-5.77104e-03) * pp + (3.42399e-03);
      dp = dp + ((1.1174e-03) * pp * pp + (-3.2747e-03) * pp + (2.3687e-03));
    }
    if (sec == 5) {
      dp = (1.54009e-03) * pp * pp + (-1.69437e-02) * pp + (1.04656e-02);
      dp = dp + ((-2.1067e-04) * pp * pp + (1.2266e-03) * pp + (-1.0553e-03));
    }
    if (sec == 6) {
      dp = (2.38182e-03) * pp * pp + (-2.07301e-02) * pp + (1.72325e-02);
      dp = dp + ((-3.6002e-04) * pp * pp + (8.9582e-04) * pp + (-1.0093e-03));
    }
  }

  //==========//  PARTICLE = PROTON (END)  //==========//

  return dp / pp;
}

// Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
// auto fe = dppC(ex, ey, ez, esec, 0) + 1;
// auto fpip = dppC(pipx, pipy, pipz, pipsec, 1) + 1;
// auto fpim = dppC(pimx, pimy, pimz, pimsec, 2) + 1;
// auto fpro = dppC(prox, proy, proz, prosec, 3) + 1;

// auto eleC = ROOT::Math::PxPyPzMVector(ex * fe, ey* fe, ez* fe, 0);
// auto pipC = ROOT::Math::PxPyPzMVector(pipx * fpip, pipy* fpip, pipz* fpip, 0.13957);
// auto pimC = ROOT::Math::PxPyPzMVector(pimx * fpim, pimy* fpim, pimz* fpim, 0.13957);
// auto proC = ROOT::Math::PxPyPzMVector(prox * fpro, proy* fpro, proz* fpro, 0.938);
////////////////// Eloss corr pip ////////////////////////////////////

// double eloss_pip(double pim_p, double pim_theta, double status_pim) {
//   double dp_pim = 0.0;

//   // INBENDING
//   if (is_FD(status_pion)) {  // Forward Detector
//     if (pim_theta < 27) {
//       dp_pim = 0.00044836 * pim_p + 0.00325965;
//     } else if (pim_theta >= 27) {
//       dp_pim = -0.00208368 * pim_p + 0.00908514;
//     }
//   }
// }
////////////////// Eloss corr pip ////////////////////////////////////

double mom_corr::elossPipFD(double pion_p, double pip_theta) {
  // momentum loss correction for low momentum pions:
  // input: p = pion momentum in GeV, pip_theta = pion theta in degree,
  //        pion_det = pion detector (2 = FD, 3 = CD),  outbending = torus polarity
  // output: dp_pion_fd = generated momentum - reconstructed momentum = momentum loss (+) / gain (-)

  double dp_pion_fd = 0.0;

  // INBENDING
  if (pip_theta < 27) {
    dp_pion_fd = 0.00342646 + (-0.00282934) * pion_p + (0.00205983) * pow(pion_p, 2) + (-0.00043158) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta < 27 && pion_p >= 2.5) {
    dp_pion_fd =
        0.00342646 + (-0.00282934) * 2.5 + (0.00205983) * pow(2.5, 2) + (-0.00043158) * pow(2.5, 3) + (0) * pow(2.5, 4);
  }
  if (pip_theta > 27 && pip_theta < 28) {
    dp_pion_fd = 0.00328565 + (-0.00376042) * pion_p + (0.00433886) * pow(pion_p, 2) + (-0.00141614) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 27 && pip_theta < 28 && pion_p >= 1.83) {
    dp_pion_fd = 0.00328565 + (-0.00376042) * 1.83 + (0.00433886) * pow(1.83, 2) + (-0.00141614) * pow(1.83, 3) +
                 (0) * pow(1.83, 4);
  }
  if (pip_theta > 28 && pip_theta < 29) {
    dp_pion_fd = 0.00328579 + (-0.00281121) * pion_p + (0.00342749) * pow(pion_p, 2) + (-0.000932614) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 28 && pip_theta < 29 && pion_p >= 2) {
    dp_pion_fd =
        0.00328579 + (-0.00281121) * 2 + (0.00342749) * pow(2, 2) + (-0.000932614) * pow(2, 3) + (0) * pow(2, 4);
  }
  if (pip_theta > 29 && pip_theta < 30) {
    dp_pion_fd = 0.00167358 + (0.00441871) * pion_p + (-0.000834667) * pow(pion_p, 2) +
                 (-0.000137968) * pow(pion_p, 3) + (0) * pow(pion_p, 4);
  }
  if (pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9) {
    dp_pion_fd = 0.00167358 + (0.00441871) * 1.9 + (-0.000834667) * pow(1.9, 2) + (-0.000137968) * pow(1.9, 3) +
                 (0) * pow(1.9, 4);
  }
  if (pip_theta > 30 && pip_theta < 31) {
    dp_pion_fd = 0.00274159 + (0.00635686) * pion_p + (-0.00380977) * pow(pion_p, 2) + (0.00071627) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9) {
    dp_pion_fd =
        0.00274159 + (0.00635686) * 1.9 + (-0.00380977) * pow(1.9, 2) + (0.00071627) * pow(1.9, 3) + (0) * pow(1.9, 4);
  }
  if (pip_theta > 31 && pip_theta < 32) {
    dp_pion_fd = 0.00450241 + (0.00248969) * pion_p + (-0.00336795) * pow(pion_p, 2) + (0.00111193) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00450241 + (0.00248969) * 1.8 + (-0.00336795) * pow(1.8, 2) + (0.00111193) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 32 && pip_theta < 33) {
    dp_pion_fd = 0.00505593 + (-0.00246203) * pion_p + (0.00172984) * pow(pion_p, 2) + (-0.000406701) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 32 && pip_theta < 33 && pion_p >= 1.8) {
    dp_pion_fd = 0.00505593 + (-0.00246203) * 1.8 + (0.00172984) * pow(1.8, 2) + (-0.000406701) * pow(1.8, 3) +
                 (0) * pow(1.8, 4);
  }
  if (pip_theta > 33 && pip_theta < 34) {
    dp_pion_fd = 0.00273402 + (0.00440449) * pion_p + (-0.00373488) * pow(pion_p, 2) + (0.000996612) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 33 && pip_theta < 34 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00273402 + (0.00440449) * 1.8 + (-0.00373488) * pow(1.8, 2) + (0.000996612) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 34 && pip_theta < 35) {
    dp_pion_fd = 0.00333542 + (0.00439874) * pion_p + (-0.00397776) * pow(pion_p, 2) + (0.00105586) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 34 && pip_theta < 35 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00333542 + (0.00439874) * 1.8 + (-0.00397776) * pow(1.8, 2) + (0.00105586) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 35 && pip_theta < 36) {
    dp_pion_fd = 0.00354663 + (0.00565397) * pion_p + (-0.00513503) * pow(pion_p, 2) + (0.00153346) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 35 && pip_theta < 36 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00354663 + (0.00565397) * 1.8 + (-0.00513503) * pow(1.8, 2) + (0.00153346) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 36 && pip_theta < 37) {
    dp_pion_fd = 0.00333909 + (0.00842367) * pion_p + (-0.0077321) * pow(pion_p, 2) + (0.0022489) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00333909 + (0.00842367) * 1.8 + (-0.0077321) * pow(1.8, 2) + (0.0022489) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 37 && pip_theta < 38) {
    dp_pion_fd = 0.00358828 + (0.0112108) * pion_p + (-0.0133854) * pow(pion_p, 2) + (0.00486924) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 37 && pip_theta < 38 && pion_p >= 1.4) {
    dp_pion_fd =
        0.00358828 + (0.0112108) * 1.4 + (-0.0133854) * pow(1.4, 2) + (0.00486924) * pow(1.4, 3) + (0) * pow(1.4, 4);
  }
  if (pip_theta > 38 && pip_theta < 39) {
    dp_pion_fd = 0.00354343 + (0.0117121) * pion_p + (-0.0129649) * pow(pion_p, 2) + (0.00455602) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 38 && pip_theta < 39 && pion_p >= 1.3) {
    dp_pion_fd =
        0.00354343 + (0.0117121) * 1.3 + (-0.0129649) * pow(1.3, 2) + (0.00455602) * pow(1.3, 3) + (0) * pow(1.3, 4);
  }
  if (pip_theta > 39 && pip_theta < 40) {
    dp_pion_fd = -0.00194951 + (0.0409713) * pion_p + (-0.0595861) * pow(pion_p, 2) + (0.0281588) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 39 && pip_theta < 40 && pion_p >= 0.9) {
    dp_pion_fd =
        -0.00194951 + (0.0409713) * 0.9 + (-0.0595861) * pow(0.9, 2) + (0.0281588) * pow(0.9, 3) + (0) * pow(0.9, 4);
  }
  if (pip_theta > 40 && pip_theta < 41) {
    dp_pion_fd = -0.0099217 + (0.0808096) * pion_p + (-0.119836) * pow(pion_p, 2) + (0.0559553) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 40 && pip_theta < 41 && pion_p >= 0.75) {
    dp_pion_fd =
        -0.0099217 + (0.0808096) * 0.75 + (-0.119836) * pow(0.75, 2) + (0.0559553) * pow(0.75, 3) + (0) * pow(0.75, 4);
  }
  if (pip_theta > 41 && pip_theta < 42) {
    dp_pion_fd = 0.00854898 + (0.00025037) * pion_p + (-0.0113992) * pow(pion_p, 2) + (0.0145178) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 41 && pip_theta < 42 && pion_p >= 0.65) {
    dp_pion_fd = 0.00854898 + (0.00025037) * 0.65 + (-0.0113992) * pow(0.65, 2) + (0.0145178) * pow(0.65, 3) +
                 (0) * pow(0.65, 4);
  }
  if (pip_theta > 42) {
    dp_pion_fd = 0.00564818 + (0.00706606) * pion_p + (0.0042602) * pow(pion_p, 2) + (-0.01141) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 42 && pion_p >= 0.65) {
    dp_pion_fd =
        0.00564818 + (0.00706606) * 0.65 + (0.0042602) * pow(0.65, 2) + (-0.01141) * pow(0.65, 3) + (0) * pow(0.65, 4);
  }
  return dp_pion_fd;
}
double mom_corr::elossPipCD(double pion_p, double pip_theta) {
  double dp_pion_cd = 0.0;

  if (pip_theta < 39) {
    dp_pion_cd =
        -0.045 + (-0.102652) + (0.455589) * pion_p + (-0.671635) * pow(pion_p, 2) + (0.303814) * pow(pion_p, 3);
  }
  if (pip_theta < 39 && pion_p >= 0.7) {
    dp_pion_cd = -0.045 + (-0.102652) + (0.455589) * 0.7 + (-0.671635) * pow(0.7, 2) + (0.303814) * pow(0.7, 3);
  }
  if (pip_theta > 39 && pip_theta < 40) {
    dp_pion_cd = 0.0684552 + (-0.766492) * pion_p + (1.73092) * pow(pion_p, 2) + (-1.46215) * pow(pion_p, 3) +
                 (0.420127) * pow(pion_p, 4);
  }
  if (pip_theta > 39 && pip_theta < 40 && pion_p >= 1.4) {
    dp_pion_cd =
        0.0684552 + (-0.766492) * 1.4 + (1.73092) * pow(1.4, 2) + (-1.46215) * pow(1.4, 3) + (0.420127) * pow(1.4, 4);
  }
  if (pip_theta > 40 && pip_theta < 41) {
    dp_pion_cd = 0.751549 + (-7.4593) * pion_p + (26.8037) * pow(pion_p, 2) + (-47.1576) * pow(pion_p, 3) +
                 (43.8527) * pow(pion_p, 4) + (-20.7039) * pow(pion_p, 5) + (3.90931) * pow(pion_p, 6);
  }
  if (pip_theta > 40 && pip_theta < 41 && pion_p >= 1.45) {
    dp_pion_cd = 0.751549 + (-7.4593) * 1.45 + (26.8037) * pow(1.45, 2) + (-47.1576) * pow(1.45, 3) +
                 (43.8527) * pow(1.45, 4) + (-20.7039) * pow(1.45, 5) + (3.90931) * pow(1.45, 6);
  }
  if (pip_theta > 41 && pip_theta < 42) {
    dp_pion_cd = -1.35043 + (10.0788) * pion_p + (-30.4829) * pow(pion_p, 2) + (47.7792) * pow(pion_p, 3) +
                 (-40.996) * pow(pion_p, 4) + (18.2662) * pow(pion_p, 5) + (-3.30449) * pow(pion_p, 6);
  }
  if (pip_theta > 41 && pip_theta < 42 && pion_p >= 1.2) {
    dp_pion_cd = -1.35043 + (10.0788) * 1.2 + (-30.4829) * pow(1.2, 2) + (47.7792) * pow(1.2, 3) +
                 (-40.996) * pow(1.2, 4) + (18.2662) * pow(1.2, 5) + (-3.30449) * pow(1.2, 6);
  }
  if (pip_theta > 42 && pip_theta < 43) {
    dp_pion_cd = -0.0231195 + (0.0744589) * pion_p + (-0.0807029) * pow(pion_p, 2) + (0.0264266) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 42 && pip_theta < 43 && pion_p >= 1.3) {
    dp_pion_cd =
        -0.0231195 + (0.0744589) * 1.3 + (-0.0807029) * pow(1.3, 2) + (0.0264266) * pow(1.3, 3) + (0) * pow(1.3, 4);
  }
  if (pip_theta > 43 && pip_theta < 44) {
    dp_pion_cd = -0.00979928 + (0.0351043) * pion_p + (-0.0365865) * pow(pion_p, 2) + (0.00977218) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 43 && pip_theta < 44 && pion_p >= 1.1) {
    dp_pion_cd =
        -0.00979928 + (0.0351043) * 1.1 + (-0.0365865) * pow(1.1, 2) + (0.00977218) * pow(1.1, 3) + (0) * pow(1.1, 4);
  }
  if (pip_theta > 44 && pip_theta < 45) {
    dp_pion_cd = 0.00108491 + (-0.00924885) * pion_p + (0.0216431) * pow(pion_p, 2) + (-0.0137762) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 44 && pip_theta < 45 && pion_p >= 1.1) {
    dp_pion_cd =
        0.00108491 + (-0.00924885) * 1.1 + (0.0216431) * pow(1.1, 2) + (-0.0137762) * pow(1.1, 3) + (0) * pow(1.1, 4);
  }
  if (pip_theta > 45 && pip_theta < 55) {
    dp_pion_cd = 0.0092263 + (-0.0676178) * pion_p + (0.168778) * pow(pion_p, 2) + (-0.167463) * pow(pion_p, 3) +
                 (0.05661) * pow(pion_p, 4);
  }
  if (pip_theta > 45 && pip_theta < 55 && pion_p >= 1.3) {
    dp_pion_cd =
        0.0092263 + (-0.0676178) * 1.3 + (0.168778) * pow(1.3, 2) + (-0.167463) * pow(1.3, 3) + (0.05661) * pow(1.3, 4);
  }
  if (pip_theta > 55 && pip_theta < 65) {
    dp_pion_cd = 0.00805642 + (-0.0670962) * pion_p + (0.188536) * pow(pion_p, 2) + (-0.20571) * pow(pion_p, 3) +
                 (0.0765) * pow(pion_p, 4);
  }
  if (pip_theta > 55 && pip_theta < 65 && pion_p >= 1.05) {
    dp_pion_cd = 0.00805642 + (-0.0670962) * 1.05 + (0.188536) * pow(1.05, 2) + (-0.20571) * pow(1.05, 3) +
                 (0.0765) * pow(1.05, 4);
  }
  if (pip_theta > 65 && pip_theta < 75) {
    dp_pion_cd = 0.00312202 + (-0.0269717) * pion_p + (0.0715236) * pow(pion_p, 2) + (-0.0545622) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 65 && pip_theta < 75 && pion_p >= 0.75) {
    dp_pion_cd = 0.00312202 + (-0.0269717) * 0.75 + (0.0715236) * pow(0.75, 2) + (-0.0545622) * pow(0.75, 3) +
                 (0) * pow(0.75, 4);
  }
  if (pip_theta > 75 && pip_theta < 85) {
    dp_pion_cd = 0.00424971 + (-0.0367683) * pion_p + (0.10417) * pow(pion_p, 2) + (-0.0899651) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 75 && pip_theta < 85 && pion_p >= 0.65) {
    dp_pion_cd =
        0.00424971 + (-0.0367683) * 0.65 + (0.10417) * pow(0.65, 2) + (-0.0899651) * pow(0.65, 3) + (0) * pow(0.65, 4);
  }
  if (pip_theta > 85 && pip_theta < 95) {
    dp_pion_cd = 0.00654123 + (-0.0517915) * pion_p + (0.147888) * pow(pion_p, 2) + (-0.14253) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 85 && pip_theta < 95 && pion_p >= 0.5) {
    dp_pion_cd =
        0.00654123 + (-0.0517915) * 0.5 + (0.147888) * pow(0.5, 2) + (-0.14253) * pow(0.5, 3) + (0) * pow(0.5, 4);
  }
  if (pip_theta > 95 && pip_theta < 105) {
    dp_pion_cd = -0.00111721 + (0.00478119) * pion_p + (0.0158753) * pow(pion_p, 2) + (-0.052902) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 95 && pip_theta < 105 && pion_p >= 0.45) {
    dp_pion_cd = -0.00111721 + (0.00478119) * 0.45 + (0.0158753) * pow(0.45, 2) + (-0.052902) * pow(0.45, 3) +
                 (0) * pow(0.45, 4);
  }
  if (pip_theta > 105 && pip_theta < 115) {
    dp_pion_cd = -0.00239839 + (0.00790738) * pion_p + (0.0311713) * pow(pion_p, 2) + (-0.104157) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 105 && pip_theta < 115 && pion_p >= 0.35) {
    dp_pion_cd = -0.00239839 + (0.00790738) * 0.35 + (0.0311713) * pow(0.35, 2) + (-0.104157) * pow(0.35, 3) +
                 (0) * pow(0.35, 4);
  }
  if (pip_theta > 115 && pip_theta < 125) {
    dp_pion_cd = -0.00778793 + (0.0256774) * pion_p + (0.0932503) * pow(pion_p, 2) + (-0.32771) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 115 && pip_theta < 125 && pion_p >= 0.35) {
    dp_pion_cd =
        -0.00778793 + (0.0256774) * 0.35 + (0.0932503) * pow(0.35, 2) + (-0.32771) * pow(0.35, 3) + (0) * pow(0.35, 4);
  }
  if (pip_theta > 125 && pip_theta < 135) {
    dp_pion_cd = -0.00292778 + (-0.00536697) * pion_p + (-0.00414351) * pow(pion_p, 2) + (0.0196431) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35) {
    dp_pion_cd = -0.00292778 + (-0.00536697) * 0.35 + (-0.00414351) * pow(0.35, 2) + (0.0196431) * pow(0.35, 3) +
                 (0) * pow(0.35, 4);
  }

  return dp_pion_cd;
}

//// Our dp mom corrections USING pass2 data
double CDProt[3][3] = {{0.006386, 0.00122, -0.04572}, {0.011795, -0.02394, -0.01157}, {0.002148, -0.002533, -0.014534}};

float mom_corr::CD_prot_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - (alpha_CD[0][0] * (CDProt[0][0] * pow(mom_, 2) + CDProt[0][1] * pow(mom_, 1) + CDProt[0][2]));
  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - (alpha_CD[0][1] * (CDProt[1][0] * pow(mom_, 2) + CDProt[1][1] * pow(mom_, 1) + CDProt[1][2]));
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - (alpha_CD[0][2] * (CDProt[2][0] * pow(mom_, 2) + CDProt[2][1] * pow(mom_, 1) + CDProt[2][2]));
  } else
    return NAN;
}

float FDProt[6][2] = {{0.02086, -0.05795},  {0.02266, -0.04367}, {0.01221, -0.02475},
                      {0.012794, -0.03397}, {0.01819, -0.04025}, {0.00882, -0.03973}};

float mom_corr::FD_prot_Hmom_corr(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_FD[0] * (FDProt[0][0] * pow(mom_, 1) + FDProt[0][1]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_FD[0] * (FDProt[1][0] * pow(mom_, 1) + FDProt[1][1]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_FD[0] * (FDProt[2][0] * pow(mom_, 1) + FDProt[2][1]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_FD[0] * (FDProt[3][0] * pow(mom_, 1) + FDProt[3][1]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_FD[0] * (FDProt[4][0] * pow(mom_, 1) + FDProt[4][1]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_FD[0] * (FDProt[5][0] * pow(mom_, 1) + FDProt[5][1]);
  } else
    return NAN;
}

///// pip
double CDPip[3][4] = {{0.04364, -0.0612, -0.001775}, {0.02333, -0.03705, 0.001743}, {0.0001144, 0.009926, -0.01627}};

float mom_corr::CD_pip_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_CD[1][0] * (CDPip[0][0] * pow(mom_, 2) + CDPip[0][1] * pow(mom_, 1) + CDPip[0][2]);
  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_CD[1][1] * (CDPip[1][0] * pow(mom_, 2) + CDPip[1][1] * pow(mom_, 1) + CDPip[1][2]);
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_CD[1][2] * (CDPip[2][0] * pow(mom_, 2) + CDPip[2][1] * pow(mom_, 1) + CDPip[2][2]);
  } else
    return NAN;
}

float FDPip[6][3] = {{0.001547, 0.01054, -0.03177},    {-0.006725, 0.03934, -0.0363}, {-0.00738, 0.03766, -0.03226},
                     {-0.0009484, 0.00803, -0.007244}, {-0.00301, 0.0178, -0.02011},  {0.002329, -0.003761, -0.01123}};

float mom_corr::FD_pip_Hmom_corr(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_FD[1] * (FDPip[0][0] * pow(mom_, 2) + FDPip[0][1] * mom_ + FDPip[0][2]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_FD[1] * (FDPip[1][0] * pow(mom_, 2) + FDPip[1][1] * mom_ + FDPip[1][2]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_FD[1] * (FDPip[2][0] * pow(mom_, 2) + FDPip[2][1] * mom_ + FDPip[2][2]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_FD[1] * (FDPip[3][0] * pow(mom_, 2) + FDPip[3][1] * mom_ + FDPip[3][2]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_FD[1] * (FDPip[4][0] * pow(mom_, 2) + FDPip[4][1] * mom_ + FDPip[4][2]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_FD[1] * (FDPip[5][0] * pow(mom_, 2) + FDPip[5][1] * mom_ + FDPip[5][2]);
  } else
    return NAN;
}
///// pim
double CDPim[3][3] = {{0.02408, -0.02676, 0.003197}, {0.006775, 0.0002553, 0.003998}, {0.01688, -0.04636, 0.012726}};

float mom_corr::CD_pim_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    // if (mom_ > 1.5 && mom_ < 1.6)
    // std::cout << "  pim mom cd  " << mom_ << "   dp val  "
    //           << alpha_CD[2][0] *
    //                  (CDPim[0][0] * pow(mom_, 2) + CDPim[0][1] * pow(mom_, 1) + CDPim[0][2] * pow(mom_, 1))
    //           << std::endl;
    return mom_ -
           alpha_CD[2][0] * (CDPim[0][0] * pow(mom_, 2) + CDPim[0][1] * pow(mom_, 1) + CDPim[0][2] * pow(mom_, 1));

  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ -
           alpha_CD[2][1] * (CDPim[1][0] * pow(mom_, 2) + CDPim[1][1] * pow(mom_, 1) + CDPim[1][2] * pow(mom_, 1));
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ -
           alpha_CD[2][2] * (CDPim[2][0] * pow(mom_, 2) + CDPim[2][1] * pow(mom_, 1) + CDPim[2][2] * pow(mom_, 1));
  } else
    return NAN;
}

float FDPim[6][3] = {{0.00857, -0.02242, -0.007652}, {0.002808, 0.00238, -0.02382},  {-0.000836, 0.0139, -0.03366},
                     {-0.001736, 0.01741, -0.04303}, {0.007187, -0.004913, -0.0368}, {0.0067, -0.0194, -0.02003}};

float mom_corr::FD_pim_Hmom_corr(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    // std::cout << "  pim mom cd  " << mom_ << "   dp val  "
    //           << alpha_FD[2] * (FDPim[0][0] * pow(mom_, 2) + FDPim[0][1] * mom_ + FDPim[0][2]) << std::endl;
    return mom_ - alpha_FD[2] * (FDPim[0][0] * pow(mom_, 2) + FDPim[0][1] * mom_ + FDPim[0][2]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_FD[2] * (FDPim[1][0] * pow(mom_, 2) + FDPim[1][1] * mom_ + FDPim[1][2]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_FD[2] * (FDPim[2][0] * pow(mom_, 2) + FDPim[2][1] * mom_ + FDPim[2][2]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_FD[2] * (FDPim[3][0] * pow(mom_, 2) + FDPim[3][1] * mom_ + FDPim[3][2]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_FD[2] * (FDPim[4][0] * pow(mom_, 2) + FDPim[4][1] * mom_ + FDPim[4][2]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_FD[2] * (FDPim[5][0] * pow(mom_, 2) + FDPim[5][1] * mom_ + FDPim[5][2]);
  } else
    return NAN;
}