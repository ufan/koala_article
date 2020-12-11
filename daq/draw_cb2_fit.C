#include "get_hist.hpp"
#include "CB2Shape.cxx"
#include "KoaColors.h"

using namespace KoaUtility;
using namespace RooFit;

// MIP Mean: Si1/Si2 = 0.35, Ge1 = 2. , Ge2 = 7
// Exponential: Si1/Si2 = -4.5/-0.7, Ge1 = -2.8/-0.2, Ge2 = -0.6/-0.07
void draw_cb2_fit(const char* sensor = "Ge1",
                  const int channel = 5,
                  double range_low = 0.8,
                  double range_high = 20,
                  double param_cb_mean = 10.7,
                  double param_cb_sigma = 0.15,
                  double param_cb_alpha_ld = 1.6,
                  double param_cb_alpha_tl = 2.6,
                  double param_cb_n_ld = 4,
                  double param_cb_n_tl = 9,
                  double param_expo1 = -2.7,
                  double param_expo2 = -0.2,
                  double param_mip_mpv = 2.2,
                  double param_mip_sigma = 0.3,
                  double frac_expo1 = 0.9,
                  double frac_mip = 0.2,
                  double frac_elatic = 0.4)
{
  init_KoaColors();
  int color = kTVibBlue;

  /*********************************************************************************************************/
  // Define the x-axis: recoil energy in MeV [0, 20]
  /*********************************************************************************************************/
  RooRealVar energy("energy", "Energy (MeV)", 0, 20);
  energy.setBins(1000);

  /*********************************************************************************************************/
  // Read binned data from file
  /*********************************************************************************************************/
  #include "get_hist.hpp"
  auto h1 = get_hist("./trigger_seed/P_2.6_calib_cluster_purification_multiple_smear_result.root", sensor, channel);
  h1->Rebin(4);
  RooDataHist dh("dh", "dh", energy, Import(*h1));

  /*********************************************************************************************************/
  // Double Exponential: as continuous background
  /*********************************************************************************************************/
  RooRealVar et1("bkg_expo1_param", "Background Exponential_1 Parameter", param_expo1, -10., 0.);
  RooExponential bkg_expo1("bkg_expo1", "Exponential Background 1", energy, et1);

  RooRealVar et2("bkg_expo2_param", "Background Exponential_2 Parameter", param_expo2, -2, 0.);
  RooExponential bkg_expo2("bkg_expo2", "Exponential Background 2", energy, et2);

  // Sum up the two components 
  RooRealVar f_bkg_expo1("f_bkg_expo1", "Fraction of exponential_1 in continuous background", frac_expo1, 0., 1.);
  RooAddPdf dbl_expo_bkg("dbl_expo_bkg", "Continuous Background (double exponential)", RooArgList(bkg_expo1, bkg_expo2), RooArgList(f_bkg_expo1), kTRUE);

  /*********************************************************************************************************/
  // MIP Background: Landau * Gaussian convolution
  /*********************************************************************************************************/
  RooRealVar ml("landau_mpv", "Landau MPV", param_mip_mpv, param_mip_mpv-0.5, param_mip_mpv+0.5);
  RooRealVar sl("landau_sigma", "Landau Sigma", param_mip_sigma, 0.0, 1.);
  RooLandau  landau("lx", "Landau", energy, ml, sl);

  RooRealVar landau_mg("landau_mg", "Landau gaussian mean", 0);
  RooRealVar landau_sg("landau_sg", "Landau gaussian sigma", 0.2, 0.0, 1);
  RooGaussian landau_gauss("landau_gauss", "Landau gauss", energy, landau_mg, landau_sg);

  energy.setBins(10000, "cache");
  RooFFTConvPdf mip_bkg("mip_bkg", "landau (X) gauss", energy, landau, landau_gauss);

  /*********************************************************************************************************/
  // Composite background model: (1 - f_mip) * dbl_expo_bkg(energy, et1, et2) + f_mip * mip_bkg(energy, ml, sl, landau_sg)
  /*********************************************************************************************************/
  RooRealVar f_mip("f_mip", "fraction of background", frac_mip, 0., 1.);
  RooAddPdf  bkg_model("bkg_model", "mip_bkg + dbl_expo_bkg", RooArgList(mip_bkg, dbl_expo_bkg), f_mip);

  /*********************************************************************************************************/
  // Elastic peak: (flat gaussian convolution)
  /*********************************************************************************************************/

  // Convolution with Crystal Ball
  RooRealVar cb_m0("cb_m0", "CrystalBall mean", param_cb_mean, param_cb_mean-0.5, param_cb_mean+0.5);
  RooRealVar cb_sigma("cb_sigma", "CrystallBall sigma", param_cb_sigma, 0., 1);

  RooRealVar cb_alpha_ld("cb_alpha_ld", "CrystallBall head turning point", param_cb_alpha_ld, 1, 3);
  RooRealVar cb_n_ld("cb_n_ld", "CrystallBall head power index", param_cb_n_ld, 1, 5);
  RooRealVar cb_alpha_tl("cb_alpha_tl", "CrystallBall tail turning point", param_cb_alpha_tl, 2, 5);
  RooRealVar cb_n_tl("cb_n_tl", "CrystallBall tail power index", param_cb_n_tl, 4, 10);

  CB2Shape elastic_model("elastic_model", "CrystallBall Pdf (Double Side)", energy, cb_m0, cb_sigma, cb_alpha_ld, cb_n_ld, cb_alpha_tl, cb_n_tl);

  /*********************************************************************************************************/
  // Final model for fitting: (1 - f_elastic) * bkg_model + f_elastic * elastic_model
  /*********************************************************************************************************/
  RooRealVar f_elastic("f_elastic", "Fraction of elastic components", frac_elatic, 0., 1.);
  RooAddPdf  model("model", "elastic_model + bkg_model", RooArgList(elastic_model, bkg_model), f_elastic);


  /*********************************************************************************************************/
  // Fitting
  /*********************************************************************************************************/
  model.fitTo(dh, Range(range_low, range_high));

  /*********************************************************************************************************/
  // Drawing
  /*********************************************************************************************************/
  RooPlot *frame = energy.frame(Title("Fitting"));
  dh.plotOn(frame, MarkerSize(0.5));
  model.plotOn(frame, Range(range_low, 20.), LineColor(kTVibRed));
  // model.plotOn(frame, Components(bkg_model), LineStyle(kSolid), LineColor(kTVibRed), Range(range_low, 20.));

  /* model.paramOn(frame, Layout(0.55)); */
  /* dh.statOn(frame, Layout(0.55, 0.99, 0.8)); */

  /* model.Print("t"); */

  // chi2 and residual
  auto chi2 = frame->chiSquare();
  auto chi2_over_ndf = frame->chiSquare(10);

  // add texts
  // TText *txt = new TText(2, 100, Form("chi^2/ndf = %.2f", chi2));
  // txt->SetTextSize(0.04);
  // txt->SetTextColor(kBlack);
  // frame->addObject(txt);

  TText *txt2 = new TText(2, 200, Form("%s_%d", sensor, channel));
  txt2->SetTextSize(0.04);
  txt2->SetTextColor(kBlack);
  frame->addObject(txt2);

  // get residual
  RooHist *hpull = frame->pullHist();
  hpull->SetMarkerSize(0.5);

  // model.plotOn(frame, Components(bkg_expo1), LineStyle(kDashed), LineColor(kTVibTeal), Range(range_low, 20.));
  // model.plotOn(frame, Components(bkg_expo2), LineStyle(kDashed), LineColor(kTVibTeal), Range(range_low, 20.));
  model.plotOn(frame, Components(mip_bkg), LineStyle(kDashed), LineColor(kTVibCyan), Range(range_low, 20.));
  model.plotOn(frame, Components(dbl_expo_bkg), LineStyle(kDashed), LineColor(kTVibTeal), Range(range_low, 20.));
  model.plotOn(frame, Components(elastic_model), LineStyle(kSolid), LineColor(kTVibBlue), Range(range_low, 20.));

  RooPlot *frame2 = energy.frame(Title("Pull Distribution"));
  frame2->addPlotable(hpull, "P");

  // Draw all frames on a canvas
  auto c = new TCanvas("c","c",900,900);
  auto pad1 = new TPad("p1","p1",0.0, 0.30, 0.999, 0.999);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.005);
  pad1->SetTopMargin(0.005);
  pad1->SetBottomMargin(0.00);
  pad1->Draw();

  auto pad2 = new TPad("p2","p2",0.0, 0.0, 0.999, 0.3);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.005);
  pad2->SetTopMargin(0.000);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();

  pad1->cd();
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->Draw();

  pad2->cd();
  frame2->GetYaxis()->SetTitleOffset(1.4);
  frame2->Draw();

  cout << "chi2 = " << chi2 << endl;
  cout << "chi2/ndf = " << chi2_over_ndf << endl;
}
