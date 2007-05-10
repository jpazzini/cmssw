#ifndef SiStripDigitizerAlgorithm_h
#define SiStripDigitizerAlgorithm_h

/** \class SiStripDigitizerAlgorithm
 *
 * SiStripDigitizerAlgorithm converts hits to digis
 *
 * \author Andrea Giammanco
 *
 * \version   1st Version Sep. 29, 2005
 *
 ************************************************************/

#include <string>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/SiStripDigitizer/interface/SiTrivialZeroSuppress.h"
#include "SimTracker/SiStripDigitizer/interface/SiTrivialDigitalConverter.h"
#include "SimTracker/SiStripDigitizer/interface/SiGaussianTailNoiseAdder.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "CommonTools/SiStripZeroSuppression/interface/SiStripNoiseService.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

namespace CLHEP {
  class HepRandomEngine;
}

class SiStripDigitizerAlgorithm 
{
 public:

  typedef  SiDigitalConverter::DigitalMapType DigitalMapType;
  typedef  SiPileUpSignals::HitToDigisMapType HitToDigisMapType;
  typedef std::map< int, float, std::less<int> > hit_map_type;
  typedef float Amplitude;

  //digisimlink
  std::vector<StripDigiSimLink> link_coll;
  std::vector<StripDigiSimLink> make_link(){ return link_coll;}

  
  SiStripDigitizerAlgorithm(const edm::ParameterSet& conf, StripGeomDetUnit *det, uint32_t& idForNoise, SiStripNoiseService*,CLHEP::HepRandomEngine&);

  ~SiStripDigitizerAlgorithm();

  // Runs the algorithm
  edm::DetSet<SiStripDigi>::collection_type  run(const std::vector<PSimHit> &input, StripGeomDetUnit *det,GlobalVector,float langle);

  void setParticleDataTable(const ParticleDataTable * pdt);
  
 private:
  int ndigis; 
  std::vector<short int> adcVec;

  edm::ParameterSet conf_;
  //-- make_digis 
  float theElectronPerADC;     // Gain, number of electrons per adc count.
  int theAdcFullScale;         // Saturation count, 255=8bit.
  float theNoiseInElectrons;   // Noise (RMS) in units of electrons.
  float theStripThreshold;     // Strip threshold in units of noise.
  float theStripThresholdInE;  // Strip noise in electorns.
  bool peakMode;
  bool noNoise;
  float tofCut;             // Cut on the particle TOF
  float theThreshold;          // ADC threshold

  double depletionVoltage;
  double appliedVoltage;
  double chargeMobility;
  double temperature;
  bool noDiffusion;
  double chargeDistributionRMS;

  SiChargeDivider* theSiChargeDivider;
  SiGaussianTailNoiseAdder* theSiNoiseAdder;
  SiPileUpSignals* theSiPileUpSignals;
  SiHitDigitizer* theSiHitDigitizer;
  SiTrivialZeroSuppress* theSiZeroSuppress;
  SiTrivialDigitalConverter* theSiDigitalConverter;
  SiStripNoiseService* SiStripNoiseService_;
  int theStripsInChip;  // num of columns per APV (for strip ineff.)

  int numStrips;  // number of strips in the module
  int strip;  // number used for noise calculation
  float moduleThickness; // sensor thickness 


  void push_digis(const DigitalMapType&,
		  const HitToDigisMapType&,
		  const SiPileUpSignals::signal_map_type&,
		  unsigned int);
 
  //-- calibration smearing
  bool doMissCalibrate;         // Switch on the calibration smearing
  float theGainSmearing;        // The sigma of the gain fluctuation (around 1)
  float theOffsetSmearing;      // The sigma of the offset fluct. (around 0)

  //-- charge fluctuation
  double tMax;  // The delta production cut, should be as in OSCAR = 30keV
                //                                           cmsim = 100keV
  std::vector<const PSimHit*> ss;
  void fluctuateEloss(int particleId, float momentum, float eloss, 
		      float length, int NumberOfSegments,
		      float elossVector[]);
  
  GeomDetType::SubDetector stripPart;            // is it barrel on forward
  const StripGeomDetUnit* _detp;
  const StripTopology* topol;
  std::vector<SiStripDigi> digis;
  CLHEP::HepRandomEngine& rndEngine;
};

#endif
