
#include <exception>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCutG.h"
#include "../../AnalysisCode/tdrstyle.C"
#include "../../AnalysisCode/Analysis_PlotFunction.h"

using namespace std;

bool Plot_WFake = true;

void Plot (void){

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.125);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1); 
   gStyle->SetNdivisions(510,"X");

   FILE* fin = fopen (Plot_WFake?"Superpose_WFake.txt":"Superpose.txt", "r");
   size_t N = 13;
   double* Mass     = new double [N];
   double* harm2Eff = new double [N];
   double* hybr2Eff = new double [N];
   double* ratioEff = new double [N];

   double MinRange = 1,
          MaxRange = 0;
   for (size_t i = 0; i < N; i++){
      fscanf (fin, "%lf %lf %lf", &Mass[i], &harm2Eff[i], &hybr2Eff[i]);

      ratioEff[i] = hybr2Eff[i]/harm2Eff[i];

      MinRange = MinRange>harm2Eff[i]?harm2Eff[i]:MinRange;
      MinRange = MinRange>hybr2Eff[i]?hybr2Eff[i]:MinRange;

      if (harm2Eff[i] < 1) MaxRange = MaxRange<harm2Eff[i]?harm2Eff[i]:MaxRange;
      if (hybr2Eff[i] < 1) MaxRange = MaxRange<hybr2Eff[i]?hybr2Eff[i]:MaxRange;
   }

   TGraph* harm2 = new TGraph (N, Mass, harm2Eff);
   TGraph* hybr2 = new TGraph (N, Mass, hybr2Eff);
   TGraph* ratio = new TGraph (N, Mass, ratioEff);

   TCanvas*  c1 = new TCanvas ("c1", "c1", 600, 600);
   TLegend* leg = new TLegend (0.50, 0.70, 0.8, 0.85);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);

   TPad* t1 = new TPad ("t1", "t1", 0.0, 0.20, 1.0, 1.0);
   t1->Draw();
   t1->cd();
   t1->SetTopMargin (0.06);

   TH1D h ("temp1", "temp1", 1, 0, Mass[N-1]);
   h.SetTitle("");
   h.SetStats(0);
   h.GetXaxis()->SetRangeUser (Mass[0], Mass[N-1]);
   h.GetYaxis()->SetRangeUser (MinRange, MaxRange*1.2);
   h.GetXaxis()->SetTitle ("Mass (GeV)");
   h.GetYaxis()->SetTitle ("Signal efficiency");
   h.Draw();

   harm2->SetLineColor(kBlue);
   harm2->SetMarkerColor(kBlue);
   harm2->SetMarkerStyle(20);
   hybr2->SetLineColor(kRed);
   hybr2->SetMarkerColor(kRed);
   hybr2->SetMarkerStyle(22);

   harm2->Draw("same P");
   hybr2->Draw("same P");

   leg->AddEntry (harm2, "Harmonic-2 an.", "LP");
   leg->AddEntry (hybr2, "Hybrid-2-15 an.", "LP");
   leg->Draw();

   c1->cd();
   TPad* t2 = new TPad ("t2", "t2", 0.0, 0.0, 1.0, 0.2);
   t2->Draw();
   t2->cd();
   t2->SetGridy(true);
   t2->SetTopMargin(0);
   t2->SetBottomMargin(0.5);

   TH1D h2 ("temp2", "temp2", 1, 0, Mass[N-1]);
   h2.SetTitle("");
   h2.SetStats(0);
   h2.GetXaxis()->SetRangeUser (Mass[0], Mass[N-1]);
   h2.GetXaxis()->SetTitle ("");
   h2.GetXaxis()->SetLabelFont (43);
   h2.GetXaxis()->SetLabelSize (15);
   h2.GetXaxis()->SetTitleFont (43);
   h2.GetXaxis()->SetTitleSize (15);
   if (!Plot_WFake) h2.GetYaxis()->SetRangeUser (0.98, 1.02);
   else             h2.GetYaxis()->SetRangeUser (0.80, 1.20);
   h2.GetYaxis()->SetTitle ("Hybr / Harm");
   h2.GetYaxis()->SetLabelFont (43);
   h2.GetYaxis()->SetLabelSize (15);
   h2.GetYaxis()->SetTitleFont (43);
   h2.GetYaxis()->SetTitleSize (15);
   h2.GetYaxis()->SetNdivisions (205);
   h2.Draw("same");

   ratio->SetMarkerColor (kRed);
   ratio->SetMarkerStyle (22);
   ratio->Draw("same P");

   c1->cd();
   DrawPreliminary ("Tracker - Only", 13.0, "Gluino, f = 10");
   SaveCanvas (c1, "Efficiency", Plot_WFake?"Plot_WRatio_WFake":"Plot_WRatio");
   delete t1;
   delete t2;
   delete c1;
   delete leg;
   delete harm2;
   delete hybr2;
   delete ratio;
   delete [] Mass;
   delete [] ratioEff;
   delete [] harm2Eff;
   delete [] hybr2Eff;
   fclose (fin);
}
