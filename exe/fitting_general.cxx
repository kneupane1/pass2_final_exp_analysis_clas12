
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <string>
using TH1D_ptr = std::shared_ptr<TH1D>;
Double_t Estimates[18];
// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}
// double gaussian(double *x, double *par) {
//   const double amp = par[0];
//   const double mu = par[1];
//   const double sigma = par[2];
//   return amp / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma) *
//          TMath::Exp(-0.5 * TMath::Power((x[0] - mu) / sigma, 2.0));
// }
// // Sum of background and peak function
// Double_t fitFunction(Double_t *x, Double_t *par) {
//   return background(x, par) + gaussian(x, &par[3]);
// }
void fitting_general(std::string inFileName1 =
                         "/home/krishna/Desktop/May-2020/may_19"
                         "/rga_missingPim_3d_05_24_2020.root ") {
  //   // gStyle->SetHistMinimumZero();
  gStyle->SetOptStat("e");
  //  gStyle->SetOptFit(1011);
  TFile *exp_data1 = new TFile(inFileName1.c_str());
  //   //
  THnSparseD *h_3d = (THnSparseD *)exp_data1->Get("THnSparse_3D/threeD_hist");
  h_3d->GetAxis(0)->SetRange(24, 25);
  h_3d->GetAxis(1)->SetRange(2, 3);

  TCanvas *MMSQ_FITS = new TCanvas("MMSQ", "MMSQ", 980, 720);
  MMSQ_FITS->Divide(6, 3); // column , row
  TH1D *mmsq_missingPim[18];
  for (int q2 = 0; q2 < 18; q2 = q2 + 1) {
    MMSQ_FITS->cd(q2 + 1);
    if (q2 != 17)
      h_3d->GetAxis(1)->SetRange(q2, q2 + 1);
    else if (q2 == 17)
      h_3d->GetAxis(1)->SetRange(q2, q2 + 5);
    mmsq_missingPim[q2] = h_3d->Projection(2);
    mmsq_missingPim[q2]->SetTitle("MMSQ ");
    mmsq_missingPim[q2]->GetXaxis()->SetRangeUser(-0.2, 0.2);
    mmsq_missingPim[q2]->SetXTitle("MMSQ (GeV2)");
    // create a TF1 with the range from 0 to 3 and 6 parameters
    TF1 *fitFcn = new TF1("fitFcn", background, -0.15, 0.15, 3);
    fitFcn->SetLineWidth(4);
    fitFcn->SetLineColor(kMagenta);
    fitFcn->SetParameters(-206.0, 4.12, 49.5);
    mmsq_missingPim[q2]->Fit("fitFcn", "0");
    fitFcn->SetParLimits(0, 2.0, 100.0);
    mmsq_missingPim[q2]->Fit("fitFcn", "Q", "ep");

    Double_t par[3];

    // writes the fit results into the par array
    fitFcn->GetParameters(par);

    TF1 *extrapolation =
        new TF1("extrapolation", background, -0.2, 0.2, 3); // extrapolation
    extrapolation->SetParameters(fitFcn->GetParameters());
    extrapolation->SetLineColor(kRed);
    extrapolation->Draw("same");
    TH1F *hnew = (TH1F *)mmsq_missingPim[q2]->Clone("hnew");
    hnew->Add(extrapolation, -1);
    hnew->Draw("same");
    // hnew->SetLineWidth();
    hnew->SetLineColor(kGreen);
    // The final result looks like that
    Double_t estimate = 0;
    estimate = hnew->Integral(25, 175);

    Estimates[q2] = estimate; //* (200 / 0.4);
  }
  for (size_t i = 0; i < 18; i++) {
    std::cout << Estimates[i] << '\n';

    // << "   ratio = " << Ratio[i]
    //           << " Entries " << Entries[i]
    //           << "  Estimates_by_ratio =  " << Estimates_by_ratio[i]
    //           << "  Diff =  " << Estimates[i] - Estimates_by_ratio[i] <<
    //           '\n';
  }
}
