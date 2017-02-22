
#include <exception>
#include <vector>
#include <algorithm>
#include <unordered_map>

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
#include "TCutG.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TPaveText.h"


namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


//#include "../../AnalysisCode_NewSyst_Hybr_WOverflow/Analysis_Step1_EventLoop.C"
//#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"
#include "../../AnalysisCode/tdrstyle.C"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "../../AnalysisCode/Analysis_CommonFunction.h"


#endif


double DistToHSCP (const reco::TrackRef& track, const std::vector<reco::GenParticle>& genColl);
bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl);
double GetMass (double P, double I, double K, double C);


const double P_Min               = 1   ;
const double P_Max               = 16  ; // 1 + 14 + 1; final one is for pixel!
const int    P_NBins             = 15  ; // 15th bin = pixel; 0 is underflow
const double Path_Min            = 0.2 ;
const double Path_Max            = 1.6 ;
const int    Path_NBins          = 42  ;
const double Charge_Min          = 0   ;
const double Charge_Max          = 5000;
const int    Charge_NBins        = 500 ;

struct perTrackHistos 
{
    TH1F * chi2                 = new TH1F ("chi2","chi2;#Chi^{2}/ndof;;",50,0,5);

    TH1I * dEdXHits             = new TH1I ("dEdXHits","dEdXHits;dE/dX clusters;;",40,0,40);
    TH1I * pixelHits            = new TH1I ("pixelHits","pixelHits;pixel clusters;;",40,0,40);
    TH1I * phase2sHits          = new TH1I ("phase2sHits","phase2sHits;phase2 strip clusters;;",40,0,40);
    TH1I * phase2sHitsHoT       = new TH1I ("phase2sHitsHoT","phase2sHitsHoT;phase2 strip clusters w/ HoT;;",40,0,40);

    TH2F * dEdXHitsVsEta        = new TH2F ("dEdXHitsVsEta","dEdXHitsVsEta;dE/dX clusters;#eta;",40,-0.5,39.5,40,-4,4);
    TH2F * pixelHitsVsEta       = new TH2F ("pixelHitsVsEta","pixelHitsVsEta;pixel clusters;#eta;",40,-0.5,39.5,40,-4,4);
    TH2F * phase2sHitsVsEta     = new TH2F ("phase2sHitsVsEta","phase2sHitsVsEta;phase2 strip clusters;#eta;",40,-0.5,39.5,40,-4,4);
    TH2F * phase2sHitsHoTVsEta  = new TH2F ("phase2sHitsHoTVsEta","phase2sHitsHoTVsEta;phase2 strip clusters w/ HoT;#eta;",40,-0.5,39.5,40,-4,4);

    TH2F * dEdXHitsVsPt         = new TH2F ("dEdXHitsVsPt","dEdXHitsVsPt;dE/dX clusters;p_{T} (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * pixelHitsVsPt        = new TH2F ("pixelHitsVsPt","pixelHitsVsPt;pixel clusters;p_{T} (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * phase2sHitsVsPt      = new TH2F ("phase2sHitsVsPt","phase2sHitsVsPt;phase2 strip clusters;p_{T} (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * phase2sHitsHoTVsPt   = new TH2F ("phase2sHitsHoTVsPt","phase2sHitsHoTVsPt;phase2 strip clusters w/ HoT;p_{T} (GeV);",40,-0.5,39.5,100,0,100);

    TH2F * dEdXHitsVsP          = new TH2F ("dEdXHitsVsP","dEdXHitsVsP;dE/dX clusters;p (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * pixelHitsVsP         = new TH2F ("pixelHitsVsP","pixelHitsVsP;pixel clusters;p (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * phase2sHitsVsP       = new TH2F ("phase2sHitsVsP","phase2sHitsVsP;phase2 strip clusters;p (GeV);",40,-0.5,39.5,100,0,100);
    TH2F * phase2sHitsHoTVsP    = new TH2F ("phase2sHitsHoTVsP","phase2sHitsHoTVsP;phase2 strip clusters w/ HoT;p (GeV);",40,-0.5,39.5,100,0,100);
    
    TH2I * dEdXHitsVsTrkHits        = new TH2I ("dEdXHitsVsTrkHits","dEdXHitsVsTrkHits;dE/dX clusters;track valid clusters;",40,0,40,40,0,40);
    TH2I * pixelHitsVsTrkHits       = new TH2I ("pixelHitsVsTrkHits","pixelHitsVsTrkHits;pixel clusters;track valid clusters;",40,0,40,40,0,40);
    TH2I * phase2sHitsVsTrkHits     = new TH2I ("phase2sHitsVsTrkHits","phase2sHitsVsTrkHits;phase2 strip clusters;track valid clusters;",40,0,40,40,0,40);
    TH2I * phase2sHitsHoTVsTrkHits  = new TH2I ("phase2sHitsHoTVsTrkHits","phase2sHitsHoTVsTrkHits;phase2 strip clusters w/ HoT;track valid clusters;",40,0,40,40,0,40);

    TH2F * pixelLayVsEta            = new TH2F ("pixelLayVsEta","pixelLayVsEta;pixel layer;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * pixelBarLayVsEta         = new TH2F ("pixelBarLayVsEta","pixelBarLayVsEta;pixel layer;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * pixelEndLayVsEta         = new TH2F ("pixelEndLayVsEta","pixelEndLayVsEta;pixel layer;#eta;",16,-0.5,15.5,40,-4,4);

    TH2F * phase2sLayVsEta          = new TH2F ("phase2sLayVsEta","phase2sLayVsEta;phase2 strip layer;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * phase2sBarLayVsEta       = new TH2F ("phase2sBarLayVsEta","phase2sBarLayVsEta;phase2 strip layer;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * phase2sEndLayVsEta       = new TH2F ("phase2sEndLayVsEta","phase2sEndLayVsEta;phase2 strip layer;#eta;",16,-0.5,15.5,40,-4,4);

    TH2F * phase2sHoTLayVsEta       = new TH2F ("phase2sHoTLayVsEta","phase2sHoTLayVsEta;phase2 strip layer w/ HoT;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * phase2sHoTBarLayVsEta    = new TH2F ("phase2sHoTBarLayVsEta","phase2sHoTBarLayVsEta;phase2 strip layer w/ HoT;#eta;",16,-0.5,15.5,40,-4,4);
    TH2F * phase2sHoTEndLayVsEta    = new TH2F ("phase2sHoTEndLayVsEta","phase2sHoTEndLayVsEta;phase2 strip layer w/ HoT;#eta;",16,-0.5,15.5,40,-4,4);
    
};

struct dEdxStudyObj
{
   string Name;
   bool isDiscrim;
   bool isEstim;
   bool useTrunc;
   bool isHit;

   bool usePixel;
   bool useStrip;

   bool mustBeInside;
   bool removeCosmics;
   bool correctFEDSat;
   bool useClusterCleaning;
   bool fakeHIP;

   double dropLowerDeDxValue;
   bool trimPixel;

   double Kconst;
   double Cconst;

   int crossTalkInvAlgo; // 0  -- do not use crossTalkInversion
                         // 1  -- use existing algorithm developed by Claude

   TH3D* Charge_Vs_Path;
   TH1D* HdedxMIP;
   TH2D* HdedxVsP;
   TH2D* HdedxVsPSyst;
//   TH2D* HdedxVsQP;
//   TProfile2D* HdedxVsP_NS;
   TProfile* HdedxVsPProfile;
   TProfile* HdedxVsEtaProfile;
   TProfile* HdedxVsNOH;
   TProfile* HNOMVsdEdxProfile;
   TH2D* HdedxVsEta;
   TH2D* HNOMVsdEdx;
   TProfile* HNOSVsEtaProfile;
   TProfile* HNOMVsEtaProfile;
   TProfile* HNOMSVsEtaProfile;
   TH1D* HMass;
   TH1D* HMassHSCP;
   TH1D* HP;
   TH1D* HHit;
   TProfile* HHitProfile; 
   TProfile* HHitProfile_U; 
   TH1D* HProtonHitSO; 
   TH1D* HProtonHitPO; 
   TProfile* Charge_Vs_FS[16];
   TH2D* Charge_Vs_XYH[16];
   TH2D* Charge_Vs_XYL[16];
   TH2D* Charge_Vs_XYHN[16];
   TH2D* Charge_Vs_XYLN[16];


   TH3F* dEdxTemplates = NULL;
   std::unordered_map<unsigned int,double>* TrackerGains = NULL;

   dEdxStudyObj(string Name_, int type_, int subdet_, TH3F* dEdxTemplates_=NULL, double K_=2.7, double C_=3.2, std::unordered_map<unsigned int,double>* TrackerGains_=NULL, bool mustBeInside_=false, bool removeCosmics_=false, bool correctFEDSat_=false, int useClusterCleaning_=0, int crossTalkInvAlgo_=0, double dropLowerDeDxValue_ = 0, bool trimPixel_ = false, bool fakeHIP_=true){
      Name = Name_;

      if     (type_==0){ isHit=true;  isEstim= false; isDiscrim = false; useTrunc = false;} // hit level only
      else if(type_==1){ isHit=false; isEstim= true;  isDiscrim = false; useTrunc = false;} // harm2
      else if(type_==2){ isHit=false; isEstim= false; isDiscrim = true;  useTrunc = false;} // Ias via harm2
      else if(type_==3){ isHit=false; isEstim= true;  isDiscrim = false; useTrunc = true; } // trunc40
      else             { isHit=false; isEstim= false; isDiscrim = false;}

           if(subdet_==1){ usePixel = true;  useStrip = false;}
      else if(subdet_==2){ usePixel = false; useStrip = true; }
      else               { usePixel = true;  useStrip = true; }

      dEdxTemplates      = dEdxTemplates_;
      TrackerGains       = TrackerGains_;
      mustBeInside       = mustBeInside_;
      removeCosmics      = removeCosmics_; 
      correctFEDSat      = correctFEDSat_;
      useClusterCleaning = static_cast<bool> (useClusterCleaning_);
      crossTalkInvAlgo   = crossTalkInvAlgo_;
      dropLowerDeDxValue = dropLowerDeDxValue_;
      trimPixel          = trimPixel_;
      fakeHIP            = fakeHIP_;

      Kconst = K_;
      Cconst = C_;

      string HistoName;
      //HitLevel plot      
      if(isHit){ 
         HistoName = Name + "_Hit";               HHit                  = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20); 
         HistoName = Name + "_HitProfile";        HHitProfile           = new TProfile(  HistoName.c_str(), HistoName.c_str(),  50, 0, 100); 
         HistoName = Name + "_HitProfile_U";      HHitProfile_U         = new TProfile(  HistoName.c_str(), HistoName.c_str(),  50, 0, 100);
         if(usePixel && useStrip){ 
            HistoName = Name + "_ChargeVsPath";      Charge_Vs_Path        = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
            for(unsigned int g=0;g<16;g++){
               char Id[255]; sprintf(Id, "%02i", g);
               HistoName = Name + "_ChargeVsFS"+Id;       Charge_Vs_FS[g]       = new TProfile  ( HistoName.c_str(), HistoName.c_str(),  769, 0, 769);
               HistoName = Name + "_ChargeVsXYH"+Id;      Charge_Vs_XYH[g]      = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -15, 15, 250, -15, 15);
               HistoName = Name + "_ChargeVsXYL"+Id;      Charge_Vs_XYL[g]      = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -15, 15, 250, -15, 15);
               HistoName = Name + "_ChargeVsXYHN"+Id;     Charge_Vs_XYHN[g]     = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -1.5, 1.5, 250, -1.5, 1.5);
               HistoName = Name + "_ChargeVsXYLN"+Id;     Charge_Vs_XYLN[g]     = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -1.5, 1.5, 250, -1.5, 1.5);
            }
         }
      }

      //Track Level plots
      if(isEstim || isDiscrim){
         HistoName = Name + "_MIP";               HdedxMIP              = new TH1D(      HistoName.c_str(), HistoName.c_str(), 1000, 0, isDiscrim?1.0:25);
         HistoName = Name + "_dedxVsP";           HdedxVsP              = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsPSyst";       HdedxVsPSyst          = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_Profile";           HdedxVsPProfile       = new TProfile(  HistoName.c_str(), HistoName.c_str(),   50, 0,100);
         HistoName = Name + "_Eta";               HdedxVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_dedxVsNOH";         HdedxVsNOH            = new TProfile(  HistoName.c_str(), HistoName.c_str(),   80, 0, 80);
         HistoName = Name + "_NOMVsdEdxProfile";  HNOMVsdEdxProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   200, 0, isDiscrim?1.0:25);
         HistoName = Name + "_NOMVsdEdx";         HNOMVsdEdx            = new TH2D(      HistoName.c_str(), HistoName.c_str(), 200, 0, isDiscrim?1.0:25, 30, 0, 30);
         HistoName = Name + "_Eta2D";             HdedxVsEta            = new TH2D(      HistoName.c_str(), HistoName.c_str(),   60,-3,  3, 100,0, isDiscrim?1.0:5);
         HistoName = Name + "_NOS";               HNOSVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_NOM";               HNOMVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_NOMS";              HNOMSVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_P";                 HP                    = new TH1D(      HistoName.c_str(), HistoName.c_str(),   50, 0, 100);  
      }

      //estimator plot only
      if(isEstim){
         HistoName = Name + "_Mass";              HMass                 = new TH1D(      HistoName.c_str(), HistoName.c_str(),  250, 0, 10);
         HistoName = Name + "_MassHSCP";          HMassHSCP             = new TH1D(      HistoName.c_str(), HistoName.c_str(),  300, 0, 3000);
         HistoName = Name + "_ProtonHitPO";       HProtonHitPO          = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20);
         HistoName = Name + "_ProtonHitSO";       HProtonHitSO          = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20);
      }

   }
};

void DeDxStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
{
  if(DIRNAME=="COMPILE") return;

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
   TH1::AddDirectory(kTRUE);

   TH3F* dEdxTemplates      = NULL;
   bool isSignal            = true;

   std::vector<string> FileName;
   if(INPUT.find(".root")<std::string::npos){
      char* pch=strtok(&INPUT[0],",");
      while (pch!=NULL){
         FileName.push_back(pch);    
         pch=strtok(NULL,",");
      }
   }

   double dEdxSF [2];
   bool SuppressFakeHIP = false;
       dEdxSF [0]      = 1.09711;
       dEdxSF [1]      = 1.09256;

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");  //File must be opened before the histogram are created

   std::vector<dEdxStudyObj*> results;
   results.push_back(new dEdxStudyObj("hit_PO"           , 0, 1, NULL, 2.7, 3.2, NULL, false, false, false, 1, 1) );
   results.push_back(new dEdxStudyObj("hit_SO_in_noC_CCC", 0, 2, NULL, 2.7, 3.2, NULL,  true,  true, false, 1, 1) );
   results.push_back(new dEdxStudyObj("hit_SP_in_noC_CCC", 0, 3, NULL, 2.7, 3.2, NULL,  true,  true, false, 1, 1) );

    perTrackHistos hists;
   
   // FILE LOOP
   for(unsigned int f=0;f<FileName.size();f++){
     TFile* file = TFile::Open(FileName[f].c_str() );
     fwlite::Event ev(file);
     // EVENT LOOP
     for(ev.toBegin(); !ev.atEnd(); ++ev){

         fwlite::Handle<DeDxHitInfoAss> dedxCollH;
         dedxCollH.getByLabel(ev, "dedxHitInfo");
         if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

         fwlite::Handle< std::vector<reco::Track> > trackCollHandle;
         trackCollHandle.getByLabel(ev,"RefitterForDeDx");
         if(!trackCollHandle.isValid()){
            trackCollHandle.getByLabel(ev,"generalTracks");
               if (!trackCollHandle.isValid()){
                  printf("Invalid trackCollHandle\n");
                  continue;
               }
         }

         fwlite::Handle < std::vector<reco::Vertex> > vertexCollHandle;
         vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
         if(!vertexCollHandle.isValid()){printf("Vertex Collection not found!\n"); continue;}
         const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
         if(vertexColl.size()<1){printf("NO VERTICES\n"); continue;}


         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         if(isSignal){
            //get the collection of generated Particles
            genCollHandle.getByLabel(ev, "genParticlesSkimmed");
            if(!genCollHandle.isValid()){
               genCollHandle.getByLabel(ev, "genParticles");
               if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
            }
         }

         std::cerr << "Event contains " << trackCollHandle->size() << " tracks." << std::endl;
         // TEST TRACK LOOP
         for(unsigned int c=0;c<trackCollHandle->size();c++){
            reco::TrackRef track = reco::TrackRef( trackCollHandle.product(), c );
            if(track.isNull())continue;
            if(fabs(track->eta())>4)continue; // should be useless
            if(track->pt()<55)continue; // might be a quite tight cut
            if( (fabs(track->dz(vertexColl[0].position())) > 0.5) || (fabs(track->dxy(vertexColl[0].position())) > 0.5) ) continue; // possibly related to cosmics more than tracks
            if(track->ptError()>0.25*track->pt()) continue;        
            if(track->chi2()/track->ndof()>5 )continue;

            const std::vector<reco::GenParticle>& genColl = *genCollHandle;
            if (DistToHSCP (track, genColl)>0.03) continue;

            
            //load the track dE/dx hit info
            const DeDxHitInfo* dedxHits = NULL;
            DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
            if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
            if(!dedxHits)continue;

            unsigned int StripHits = 0,
                         PixelHits = 0;
            for(unsigned int h=0;h<dedxHits->size();h++){
                DetId detId (dedxHits->detId(h));
                if (detId.subdetId() < SiStripDetId::TIB)
                    PixelHits++;
                else StripHits++;
                
                switch (detId.subdetId()){
                    case 1: hists.pixelBarLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                            hists.pixelLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                            break;//pixel barrel
                    case 2: hists.pixelEndLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                            hists.pixelLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                            break;//pixel endcap
                    case 5: hists.phase2sBarLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                            hists.phase2sLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                            if (dedxHits->phase2cluster(h)->threshold()>0){
                                hists.phase2sHoTBarLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                                hists.phase2sHoTLayVsEta->Fill(int(detId>>20 & 0xF),track->eta());
                            }
                            break;//tracker outer barrel
                    case 4: hists.phase2sEndLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                            hists.phase2sLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                            if (dedxHits->phase2cluster(h)->threshold()>0){
                                hists.phase2sHoTEndLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                                hists.phase2sHoTLayVsEta->Fill(int(detId>>18 & 0xF),track->eta());
                            }
                            break;//tracker outer endcap
                    default : std::cerr << "This should not happen: not in the tracker...\n"; 
                              break;
                }
                
            }
            const std::vector <SiPixelCluster> pixels = dedxHits->pixelClusters();
            const std::vector <SiStripCluster> strips = dedxHits->stripClusters();
            const std::vector <Phase2TrackerCluster1D> phase2s = dedxHits->phase2TrackerCluster1D();

            unsigned int phase2sHoT = 0;
            for (auto & it : phase2s ){
                if ( it.threshold() > 0 ) ++phase2sHoT;
            }

            fprintf (stderr, "Track %d: %lu/%lu hits: %lu P, %lu S, %lu P2\n", c, dedxHits->size(), track->recHitsSize(),
                    pixels.size(), strips.size(), phase2s.size());

            hists.chi2->Fill(track->chi2()/track->ndof());
            
            hists.dEdXHits->Fill(dedxHits->size());
            hists.pixelHits->Fill(pixels.size());
            hists.phase2sHits->Fill(phase2s.size());
            hists.phase2sHitsHoT->Fill(phase2sHoT);
            
            hists.dEdXHitsVsEta->Fill(dedxHits->size(),track->eta());
            hists.pixelHitsVsEta->Fill(pixels.size(),track->eta());
            hists.phase2sHitsVsEta->Fill(phase2s.size(),track->eta());
            hists.phase2sHitsHoTVsEta->Fill(phase2sHoT,track->eta());
            
            hists.dEdXHitsVsPt->Fill(dedxHits->size(),track->pt());
            hists.pixelHitsVsPt->Fill(pixels.size(),track->pt());
            hists.phase2sHitsVsPt->Fill(phase2s.size(),track->pt());
            hists.phase2sHitsHoTVsPt->Fill(phase2sHoT,track->pt());
            
            hists.dEdXHitsVsP->Fill(dedxHits->size(),track->p());
            hists.pixelHitsVsP->Fill(pixels.size(),track->p());
            hists.phase2sHitsVsP->Fill(phase2s.size(),track->p());
            hists.phase2sHitsHoTVsP->Fill(phase2sHoT,track->p());
            
            hists.dEdXHitsVsTrkHits->Fill(dedxHits->size(),track->numberOfValidHits());
            hists.pixelHitsVsTrkHits->Fill(pixels.size(),track->numberOfValidHits());
            hists.phase2sHitsVsTrkHits->Fill(phase2s.size(),track->numberOfValidHits());            
            hists.phase2sHitsHoTVsTrkHits->Fill(phase2sHoT,track->numberOfValidHits());
         } // END TEST TRACK LOOP

         // LOOP OVER TRACKS
         
         for(unsigned int c=0;c<trackCollHandle->size();c++){
            //basic track quality cuts
            reco::TrackRef track = reco::TrackRef( trackCollHandle.product(), c );
            if(track.isNull())continue;
            if(track->chi2()/track->ndof()>5 )continue;
            if(track->found()<8)continue;

//             if(isSignal){
//                if(track->pt()<45)continue;
//                const std::vector<reco::GenParticle>& genColl = *genCollHandle;
//                if (DistToHSCP (track, genColl)>0.03) continue;
//             }else{
//                if(track->pt()>=45)continue;
//             }

            //load dEdx information
            const DeDxHitInfo* dedxHits = NULL;
            DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
            if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
            if(!dedxHits)continue;

            //hit level dEdx information (only done for MIPs)
            if(track->pt() > 5){
               for(unsigned int h=0;h<dedxHits->size();h++){
                   DetId detid(dedxHits->detId(h));
                   int moduleGeometry = 0; // underflow bin -- debug purposes
                   if(detid.subdetId()>=3){ SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry(); if (moduleGeometry==15) {cerr << "ERROR! There is no SiStrip geometry 15!" << endl; exit (EXIT_FAILURE);}}
                   else if(detid.subdetId()<3){moduleGeometry = 15;} // 15 is for pixel

                   for(unsigned int R=0;R<results.size();R++){
//                      if (results[R]->Name.find("newCCC")!=string::npos){dEdxSF[0] = dEdxSF_NewCC[0]; dEdxSF[1] = dEdxSF_NewCC[1];}
//                      else {dEdxSF[0] = dEdxSF_OldCC[0]; dEdxSF[1] = dEdxSF_OldCC[1];}
                      double scaleFactor = dEdxSF[0];
                      if (detid.subdetId()<3) scaleFactor *= dEdxSF[1];
                      double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;


                      if(!results[R]->isHit) continue; //only consider results related to hit info here
                      if(!results[R]->usePixel && detid.subdetId() <3)continue; // skip pixels
                      if(!results[R]->useStrip && detid.subdetId()>=3)continue; // skip strips
// // //                       if(results[R]->mustBeInside && !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->phase2cluster(h):NULL) )continue; /// NOTE -> JACOPO
                      if(results[R]->removeCosmics){ if (isCompatibleWithCosmic(track, vertexColl))continue;} //don't consider hits, which belong to cosmic tracks

// // //                       if(results[R]->useClusterCleaning && detid.subdetId()>=3 && !clusterCleaning(dedxHits->phase2cluster(h), results[R]->crossTalkInvAlgo)) continue; //if it fails clusterCleaning, skip it! /// NOTE -> JACOPO
//                      if(results[R]->CCFunction && detid.subdetId()>=3 && !results[R]->CCFunction(dedxHits->phase2cluster(h), results[R]->crossTalkInvAlgo, NULL)) continue; //if it fails clusterCleaning, skip it!

                      int charge = dedxHits->charge(h);
                      double ChargeOverPathlength   = 0.;
                      double ChargeOverPathlength_U = 0.;
                      if (detid.subdetId()>=3 ){
                        ChargeOverPathlength   = dedxHits->phase2cluster(h)->threshold(); // FIXME -> not dividing for pathlength
                        ChargeOverPathlength_U = dedxHits->phase2cluster(h)->threshold(); // FIXME -> not dividing for pathlength
                      } else {                          
                        ChargeOverPathlength   = scaleFactor*Norm*charge/dedxHits->pathlength(h);
                        ChargeOverPathlength_U = 1.0*Norm*charge/dedxHits->pathlength(h);
                      }
                      results[R]->HHit->Fill(ChargeOverPathlength);

//                       if (fabs(track->eta())<0.4){
//                          results[R]->HHitProfile->Fill(track->p(), ChargeOverPathlength);
//                          results[R]->HHitProfile_U->Fill(track->p(), ChargeOverPathlength_U);
//                       }

                      if(results[R]->usePixel && results[R]->useStrip){
                         
                         results[R]->Charge_Vs_Path->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1))); 
                         if(detid.subdetId()>=3)results[R]->Charge_Vs_FS[moduleGeometry]->Fill(dedxHits->phase2cluster(h)->firstStrip(), charge); 
                         results[R]->Charge_Vs_XYH[moduleGeometry]->Fill(dedxHits->pos(h).x(), dedxHits->pos(h).y()); 
                         if(ChargeOverPathlength<1.6)results[R]->Charge_Vs_XYL[moduleGeometry]->Fill(dedxHits->pos(h).x(), dedxHits->pos(h).y()); 
    
//                          if(moduleGeometry>=1 && moduleGeometry<=14){ // FIXME we don't have the geometry information for Pixels yet (TkModGeom* arrays) !!!
//                             double nx, ny;
//                             if(moduleGeometry<=4){
//                                ny = dedxHits->pos(h).y() /  TkModGeomLength[moduleGeometry];
//                                nx = dedxHits->pos(h).x() /  TkModGeomWidthT[moduleGeometry];
//                             }else{
//                                double  offset = TkModGeomLength[moduleGeometry] * (TkModGeomWidthT[moduleGeometry]+TkModGeomWidthB[moduleGeometry]) / (TkModGeomWidthT[moduleGeometry]-TkModGeomWidthB[moduleGeometry]);  // check sign if GeomWidthT[moduleGeometry] < TkModGeomWidthB[moduleGeometry] !!! 
//                                double  tan_a = TkModGeomWidthT[moduleGeometry] / std::abs(offset + TkModGeomLength[moduleGeometry]);
//                                ny = dedxHits->pos(h).y() /  TkModGeomLength[moduleGeometry];
//                                nx = dedxHits->pos(h).x() / (tan_a*std::abs(dedxHits->pos(h).y()+offset));
//                             }
//                             //printf("%i - %f - %f --> %f - %f\n", moduleGeometry, dedxHits->pos(h).x(), dedxHits->pos(h).y(), nx, ny);
//                             results[R]->Charge_Vs_XYHN[moduleGeometry]->Fill(nx, ny); 
//                             if(ChargeOverPathlength<1.6)results[R]->Charge_Vs_XYLN[moduleGeometry]->Fill(nx, ny);
//                          }
                      }
                   }
                }
             }

//              bool isCosmic = isCompatibleWithCosmic(track, vertexColl);
//              bool lockOnTrack=false;
//              double dEdxDebug = 0;
//              for(unsigned int R=0;R<results.size();R++){
//                 if(!results[R]->isEstim and !results[R]->isDiscrim) continue; //only consider results related to estimator/discriminator variables here
//                 if(results[R]->removeCosmics && isCosmic)continue; //don't consider cosmic tracks
// 
//                 DeDxData dedxObj   = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, NULL);
//  
// 
// 		if (isSignal) results[R]->HdedxVsP->SetBins(1000, 0, 2400, results[R]->isDiscrim?1000:2000, 0, results[R]->isDiscrim?1.0:30); // if it's signal sample increase axis range
//                 results[R]->HdedxVsP    ->Fill(track->p(), dedxObj.dEdx() );
// 
//                 if(track->pt()>10 && track->pt()<45 && dedxObj.numberOfMeasurements()>=(results[R]->useStrip?7:3) ){
//                   results[R]->HdedxVsEtaProfile->Fill(track->eta(), dedxObj.dEdx() );
//                   results[R]->HdedxVsEta->Fill(track->eta(), dedxObj.dEdx() );
//                   results[R]->HNOMVsEtaProfile->Fill(track->eta(),dedxObj.numberOfMeasurements() );
//                   results[R]->HNOSVsEtaProfile->Fill(track->eta(),dedxObj.numberOfSaturatedMeasurements() );
//                   results[R]->HNOMSVsEtaProfile->Fill(track->eta(),dedxObj.numberOfMeasurements() - dedxObj.numberOfSaturatedMeasurements() );
//                 }
// 
//                 if(fabs(track->eta())>2.1) continue;
//                 if((int)dedxObj.numberOfMeasurements()<(results[R]->useStrip?10:3))continue;
// //                if(track->found()<10) continue; // we cut on total number of hits instead of valid measurements
// 
//                 if(track->pt()>5){
//                    results[R]->HdedxVsNOH->Fill(track->found(), dedxObj.dEdx());
//                    results[R]->HNOMVsdEdxProfile->Fill(dedxObj.dEdx(), (int)dedxObj.numberOfMeasurements());
//                    results[R]->HNOMVsdEdx->Fill(dedxObj.dEdx(), (int)dedxObj.numberOfMeasurements());
//                    results[R]->HdedxMIP  ->Fill(dedxObj.dEdx());
//                    results[R]->HP->Fill(track->p());
//                 }
//                 if(fabs(track->eta())<0.4){
//                    results[R]->HdedxVsPProfile  ->Fill(track->p(), dedxObj  .dEdx() );
//                 }
// 
//                 DeDxData dedxObjEstim   = results[R]->isEstim?dedxObj:computedEdx(dedxHits, dEdxSF, NULL                     , results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, NULL);
//                 double Mass = GetMass(track->p(),dedxObjEstim.dEdx(), results[R]->Kconst, results[R]->Cconst);
//                 if (Mass > 0.938-0.20 && Mass < 0.938+0.20 && dedxObjEstim.dEdx() > 5){// proton candidates
//                    results[R]->HdedxVsPSyst->Fill(track->p(), dedxObj.dEdx() );
// 		   
//                    if (results[R]->isEstim && dedxObj.dEdx() > 7){
//                       for(unsigned int h=0;h<dedxHits->size();h++){
//                          DetId detid(dedxHits->detId(h));
//                          double scaleFactor = dEdxSF[0];
//                          if (detid.subdetId()<3) scaleFactor *= dEdxSF[1];
//                          double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
//                        
//                          int moduleGeometry = 0; // underflow bin -- debug purposes
//                          if(detid.subdetId()>=3){ SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry(); if (moduleGeometry==15) {cerr << "ERROR! There is no SiStrip geometry 15!" << endl; exit (EXIT_FAILURE);}}
//                          else if(detid.subdetId()<3){moduleGeometry = 15;} // 15 is for pixel
//                          if(!results[R]->usePixel && detid.subdetId() <3)continue; // skip pixels
//                          if(!results[R]->useStrip && detid.subdetId()>=3)continue; // skip strips
//                          if(results[R]->mustBeInside && !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->phase2cluster(h):NULL) )continue;
//                          if(results[R]->removeCosmics){ if (isCompatibleWithCosmic(track, vertexColl))continue;} //don't consider hits, which belong to cosmic tracks
//                          if(results[R]->useClusterCleaning && detid.subdetId()>=3 && !clusterCleaning(dedxHits->phase2cluster(h), results[R]->crossTalkInvAlgo)) continue; //if it fails clusterCleaning, skip it!
//                        
//                          int charge = dedxHits->charge(h);
//                          if (detid.subdetId()>=3 && results[R]->crossTalkInvAlgo!=0){
//                             vector <int> amps = CrossTalkInv(convert(dedxHits->phase2cluster(h)->amplitudes()), 0.10, 0.04, true);
//                             charge = std::accumulate(amps.begin(), amps.end(), 0);
//                          }
//                          double ChargeOverPathlength = scaleFactor*Norm*charge/dedxHits->pathlength(h);
// 
// 			 if (detid.subdetId()<3) results[R]->HProtonHitPO->Fill(ChargeOverPathlength);
// 			 else                    results[R]->HProtonHitSO->Fill(ChargeOverPathlength);
//                       }
//                    }
//                 }
// 
//                 if(results[R]->isEstim && dedxObj.dEdx()>results[R]->Cconst + 3.0){  //mass can only be computed for dEdx estimators
//                    if(track->p()<3.0){      results[R]->HMass->Fill(Mass);
//                    }else{                   results[R]->HMassHSCP->Fill(Mass);
//                    }
//                 }
//              }
         } // END TRACK LOOP
      } // END EVENT LOOP
      printf("\n");
      delete file;
   } // END FILE LOOP

   OutputHisto->Write();
   OutputHisto->Close();  
}

double DistToHSCP (const reco::TrackRef& track, const std::vector<reco::GenParticle>& genColl){
   if(track.isNull())return false; // FIXME does this make sense? returning false to a double function?

   double RMin = 9999;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000 && AbsPdg!=17)continue;
      double dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
      if(dR<RMin)RMin=dR;
   }
   return RMin;
}

bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl){
   for (unsigned int vertex_i=0;vertex_i<vertexColl.size();vertex_i++){
      if(fabs(track->dz (vertexColl[vertex_i].position())) < 0.5 && fabs(track->dxy(vertexColl[vertex_i].position())) < 0.2)return false;
   }
   return true;
}

double GetMass (double P, double I, double K, double C){
   if (I-C<0) return -1;
   return sqrt((I-C)/K)*P;
}

