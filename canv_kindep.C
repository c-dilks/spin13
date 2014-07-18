{
//=========Macro generated from canvas: canv_kindep/canv_kindep
//=========  (Fri Jul 18 17:14:38 2014) by ROOT version5.34/14
   TCanvas *canv_kindep = new TCanvas("canv_kindep", "canv_kindep",974,23,800,1000);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   canv_kindep->Range(0,0,1,1);
   canv_kindep->SetFillColor(0);
   canv_kindep->SetBorderMode(0);
   canv_kindep->SetBorderSize(2);
   canv_kindep->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: kpad1
   TPad *kpad1 = new TPad("kpad1", "kpad1",0.05,0.6633333,0.95,0.95);
   kpad1->Draw();
   kpad1->cd();
   kpad1->Range(-34.375,-0.2362198,309.375,0.05983118);
   kpad1->SetFillColor(0);
   kpad1->SetBorderMode(0);
   kpad1->SetBorderSize(0);
   kpad1->SetGridx();
   kpad1->SetGridy();
   kpad1->SetTopMargin(0);
   kpad1->SetBottomMargin(0);
   kpad1->SetFrameBorderMode(0);
   kpad1->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(6);
   gre->SetName("Graph0");
   gre->SetTitle("A_{L}^{Y} #pm #sigma A_{L}^{Y} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff9900");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.3);
   gre->SetPoint(0,12.5,-0.0002202982);
   gre->SetPointError(0,12.5,0.0005727711);
   gre->SetPoint(1,37.5,0.0001839363);
   gre->SetPointError(1,12.5,0.000543869);
   gre->SetPoint(2,62.5,0.0002735636);
   gre->SetPointError(2,12.5,0.00113479);
   gre->SetPoint(3,87.5,-0.003034532);
   gre->SetPointError(3,12.5,0.003622449);
   gre->SetPoint(4,125,0.02304185);
   gre->SetPointError(4,25,0.01211841);
   gre->SetPoint(5,200,-0.1033239);
   gre->SetPointError(5,50,0.108225);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","A_{L}^{Y} #pm #sigma A_{L}^{Y} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)",100,0,275);
   Graph_Graph1->SetMinimum(-0.2362198);
   Graph_Graph1->SetMaximum(0.05983118);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("A_{L}^{Y}");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   gre->Draw("ape");
   TLine *line = new TLine(0,0,275,0);

   ci = TColor::GetColor("#006666");
   line->SetLineColor(ci);
   line->SetLineStyle(2);
   line->SetLineWidth(3);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.2415355,0.8961509,0.7584645,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("A_{L}^{Y} #pm #sigma A_{L}^{Y} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   pt->Draw();
   kpad1->Modified();
   canv_kindep->cd();
  
// ------------>Primitives in pad: kpad2
   kpad2 = new TPad("kpad2", "kpad2",0.05,0.3766667,0.95,0.6633333);
   kpad2->Draw();
   kpad2->cd();
   kpad2->Range(-34.375,-0.2357746,309.375,0.03365139);
   kpad2->SetFillColor(0);
   kpad2->SetBorderMode(0);
   kpad2->SetBorderSize(0);
   kpad2->SetGridx();
   kpad2->SetGridy();
   kpad2->SetTopMargin(0);
   kpad2->SetBottomMargin(0);
   kpad2->SetFrameBorderMode(0);
   kpad2->SetFrameBorderMode(0);
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph0");
   gre->SetTitle("A_{L}^{B} #pm #sigma A_{L}^{B} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.3);
   gre->SetPoint(0,12.5,-0.0008192614);
   gre->SetPointError(0,12.5,0.0005955166);
   gre->SetPoint(1,37.5,-0.0003339422);
   gre->SetPointError(1,12.5,0.0005662829);
   gre->SetPoint(2,62.5,-0.0001936806);
   gre->SetPointError(2,12.5,0.00118145);
   gre->SetPoint(3,87.5,-0.0009344188);
   gre->SetPointError(3,12.5,0.003770878);
   gre->SetPoint(4,125,-0.01629984);
   gre->SetPointError(4,25,0.01261188);
   gre->SetPoint(5,200,-0.1010616);
   gre->SetPointError(5,50,0.1122608);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","A_{L}^{B} #pm #sigma A_{L}^{B} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)",100,0,275);
   Graph_Graph2->SetMinimum(-0.2357746);
   Graph_Graph2->SetMaximum(0.03365139);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetTitle("A_{L}^{B}");
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph2);
   
   gre->Draw("ape");
   line = new TLine(0,0,275,0);

   ci = TColor::GetColor("#006666");
   line->SetLineColor(ci);
   line->SetLineStyle(2);
   line->SetLineWidth(3);
   line->Draw();
   
   pt = new TPaveText(0.2415355,0.8961509,0.7584645,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("A_{L}^{B} #pm #sigma A_{L}^{B} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   pt->Draw();
   kpad2->Modified();
   canv_kindep->cd();
  
// ------------>Primitives in pad: kpad3
   kpad3 = new TPad("kpad3", "kpad3",0.05,0.05,0.95,0.3766667);
   kpad3->Draw();
   kpad3->cd();
   kpad3->Range(-34.375,-0.3273398,309.375,0.2123854);
   kpad3->SetFillColor(0);
   kpad3->SetBorderMode(0);
   kpad3->SetBorderSize(0);
   kpad3->SetGridx();
   kpad3->SetGridy();
   kpad3->SetTopMargin(0);
   kpad3->SetBottomMargin(0.122449);
   kpad3->SetFrameBorderMode(0);
   kpad3->SetFrameBorderMode(0);
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph0");
   gre->SetTitle("A_{LL} #pm #sigma A_{LL} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.3);
   gre->SetPoint(0,12.5,-2.486712e-05);
   gre->SetPointError(0,12.5,0.001046049);
   gre->SetPoint(1,37.5,0.0002990677);
   gre->SetPointError(1,12.5,0.0009945243);
   gre->SetPoint(2,62.5,-0.003084279);
   gre->SetPointError(2,12.5,0.00207396);
   gre->SetPoint(3,87.5,0.008418852);
   gre->SetPointError(3,12.5,0.006619331);
   gre->SetPoint(4,125,0.00917514);
   gre->SetPointError(4,25,0.02212718);
   gre->SetPoint(5,200,-0.02443281);
   gre->SetPointError(5,50,0.1973485);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","A_{LL} #pm #sigma A_{LL} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)",100,0,275);
   Graph_Graph3->SetMinimum(-0.261251);
   Graph_Graph3->SetMaximum(0.2123854);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetTitle("A_{LL}");
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph3->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph3);
   
   gre->Draw("ape");
   line = new TLine(0,0,275,0);

   ci = TColor::GetColor("#006666");
   line->SetLineColor(ci);
   line->SetLineStyle(2);
   line->SetLineWidth(3);
   line->Draw();
   
   pt = new TPaveText(0.1968677,0.9157979,0.8031323,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("A_{LL} #pm #sigma A_{LL} vs. E_{#gamma#gamma} for p_{T}#in[0.00,15.00) and #eta#in[2.50,4.00)");
   pt->Draw();
   kpad3->Modified();
   canv_kindep->cd();
   
   pt = new TPaveText(0.3178392,0.9589322,0.6545226,0.9917864,"br");
   text = pt->AddText("hello world");
   pt->Draw();
   canv_kindep->Modified();
   canv_kindep->cd();
   canv_kindep->SetSelected(canv_kindep);
   canv_kindep->ToggleToolBar();
}
