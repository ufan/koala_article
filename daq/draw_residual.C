void draw_residual(const char* filename_bg="calib_pulser_beamGain_2019_Sep_23_21_38_16_calib.root",
                   const char* filename_cg="calib_pulser_commGain_2019_Sep_23_22_39_27_calib.roo",
                   const char* detector="Ge1",
                   int channel=21)
{
  auto fbg = TFile::Open(filename_bg);
  auto* gdir = fbg->GetDirectory("graphs");
  TGraphErrors* gr_bg = nullptr;
  /* gdir->GetObject(Form("g_%s_%d",detector,channel), gr_bg); */
  gdir->GetObject(Form("g_inverse_%s_%d",detector,channel), gr_bg);

  // TF1* f1 = gr->GetFunction("fpol1");
  auto n = gr_bg->GetN();
  auto x = gr_bg->GetX();
  auto y = gr_bg->GetY();
  auto erry = gr_bg->GetEY();

  auto gr_bg_new = new TGraph(n-1, x+1, y+1);
  gr_bg_new->SetName(Form("g_bg_new_%s_%d", detector, channel));
  gr_bg_new->SetMarkerStyle(20);
  gr_bg_new->SetMarkerSize(1.1);
  gr_bg_new->Fit("pol1");
  TF1* f1 = gr_bg_new->GetFunction("pol1");

  auto gr_bg_residual = new TGraph();
  gr_bg_residual->SetName(Form("g_bg_residual_%s_%d", detector, channel));
  gr_bg_residual->SetMarkerStyle(20);
  gr_bg_residual->SetMarkerSize(1.1);
  for(int i=1;i<n;i++){
    gr_bg_residual->SetPoint(i-1, x[i], f1->Eval(x[i])-y[i]);
    /* gr_bg_residual->SetPoint(i-1, x[i], (f1->Eval(x[i])-y[i])/y[i]); */
    // gr_bg_residual->SetPoint(i-1, x[i], (f1->Eval(x[i])-y[i])/erry[i]);
    // std::cout << x[i] << "\t" << y[i] << "\t" << f1->Eval(x[i])-y[i] << std::endl;
  }
  delete fbg;

  auto fcg = TFile::Open(filename_cg);
  gdir = fcg->GetDirectory("graphs");
  TGraphErrors* gr_cg = nullptr;
  /* gdir->GetObject(Form("g_%s_%d",detector,channel), gr_cg); */
  gdir->GetObject(Form("g_inverse_%s_%d",detector,channel), gr_cg);

  // TF1* f1 = gr->GetFunction("fpol1");
  auto ncg = gr_cg->GetN();
  auto xcg = gr_cg->GetX();
  auto ycg = gr_cg->GetY();
  auto erry_cg = gr_cg->GetEY();

  auto gr_cg_new = new TGraph(ncg-1, xcg+1, ycg+1);
  gr_cg_new->SetName(Form("g_cg_new_%s_%d", detector, channel));
  gr_cg_new->SetMarkerStyle(22);
  gr_cg_new->SetMarkerSize(1.3);
  gr_cg_new->SetMarkerColor(kBlue);
  gr_cg_new->Fit("pol1");
  TF1* f1_cg = gr_cg_new->GetFunction("pol1");

  auto gr_cg_residual = new TGraph();
  gr_cg_residual->SetName(Form("g_cg_residual_%s_%d", detector, channel));
  gr_cg_residual->SetMarkerStyle(22);
  gr_cg_residual->SetMarkerSize(1.3);
  gr_cg_residual->SetMarkerColor(kBlue);
  for(int i=1;i<ncg;i++){
    gr_cg_residual->SetPoint(i-1, xcg[i], f1_cg->Eval(xcg[i])-ycg[i]);
    /* gr_cg_residual->SetPoint(i-1, xcg[i], (f1_cg->Eval(xcg[i])-ycg[i])/ycg[i]); */
    // gr_cg_residual->SetPoint(i-1, xcg[i], (f1_cg->Eval(xcg[i])-ycg[i])/erry_cg[i]);
  }
  delete fcg;

  // auto haxis = new TH2("haxis", "haxis", 10, 0, x[n-1]+100, 10, );
  auto mg = new TMultiGraph();
  mg->SetName("mg");
  mg->Add(gr_bg_new,"P");
  mg->Add(gr_cg_new,"P");

  auto mg_residual = new TMultiGraph();
  mg_residual->SetName("mg_residual");
  mg_residual->Add(gr_bg_residual,"P");
  mg_residual->Add(gr_cg_residual,"P");

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

  // auto pad3 = new TPad("p3","p3",0.1, 0.1,  0.998, 0.25);
  // pad3->SetTopMargin(0.001);
  // pad3->Draw();
  //


  ///////////////////////////////////////////
  pad1->cd();
  mg->Draw("AP");
  mg->GetXaxis()->SetRangeUser(0, 1100);
  mg->GetYaxis()->SetRangeUser(-300, 7700);

  mg->GetYaxis()->SetTitle("ADC counts");
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleOffset(1.0);
  mg->GetYaxis()->SetLabelSize(0.05);

  pad2->cd();
  mg_residual->Draw("AP");
  // gr_bg_residual->GetXaxis()->SetRangeUser(0, 1100);
  mg_residual->GetXaxis()->SetRangeUser(0, 1100);

  //Int_t font=72;
  mg_residual->GetXaxis()->SetTitle("Pulser Amplitude (mV)");
  mg_residual->GetXaxis()->SetTitleSize(0.12);
  mg_residual->GetXaxis()->SetTitleOffset(1.2);
  mg_residual->GetXaxis()->SetLabelSize(0.12);

  mg_residual->GetYaxis()->SetTitle("Residual");
  mg_residual->GetYaxis()->SetLabelSize(0.12);
  mg_residual->GetYaxis()->SetTitleSize(0.14);
  mg_residual->GetYaxis()->SetTitleOffset(0.42);
  /* mg_residual->GetXaxis()->SetLabelOffset(0); */
  /* mg_residual->GetXaxis()->SetLabelFont(font,"y"); */
  /* mg_residual->GetXaxis()->SetTitleFont(font,"y"); */
  /* mg_residual->GetXaxis()->SetLabelFont(font,"z"); */
  /* mg_residual->GetXaxis()->SetTitleFont(font,"z"); */
  /* mg_residual->GetXaxis()->SetLabelSize(tsize,"x"); */
  /* mg_residual->GetXaxis()->SetTitleSize(tsize,"x"); */
  /* mg_residual->GetXaxis()->SetLabelSize(tsize,"y"); */
  /* mg_residual->GetXaxis()->SetTitleSize(tsize,"y"); */
  /* mg_residual->GetXaxis()->SetLabelSize(tsize,"z"); */
  /* mg_residual->GetXaxis()->SetTitleSize(tsize,"z"); */

}
