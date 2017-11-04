// -*- C++ -*-
//
//
// Description: psi(2S)-> gamme p p
//
// Original Author:  SHI Xiaodong <wherenpc@mail.ustc.edu.cn>
//

//#define TEST

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/LoadFactoryEntries.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"

#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"

#include "ParticleID/ParticleID.h"
#include "McTruth/McParticle.h"


#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <vector>

class GammaPP: public Algorithm {
  
public:
  GammaPP(const std::string&, ISvcLocator*);
  ~GammaPP(); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
	NTuple::Tuple* m_tuple1;  // 
  // declare some cuts for charged tracks
	double m_vz0cut;
	double m_vr0cut;
	double m_cha_costheta_cut;
	double m_prob_pion_min;
	double m_ep_electron_min;
	double m_cha_momentum_max_cut;
  // declare some cuts for neutral tracks
	double m_min_emctime;
	double m_max_emctime;
	double m_costheta_barrel_max;
	double m_costheta_endcap_min;
	double m_costheta_endcap_max;
	double m_energy_barrel_min;
	double m_energy_endcap_min;
	double m_gammaCosCut;
	double m_photon_iso_angle_min;





  NTuple::Item<int>  m_run;
  NTuple::Item<int>  m_rec;
  NTuple::Item<int>  m_ngch;
  NTuple::Item<int>  m_nGam;
  NTuple::Item<int>  m_nCharge;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_trkidx;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
  NTuple::Item<int>  m_np;
  NTuple::Item<int>  m_nap;
  NTuple::Item<int>  m_ext_pi0;
  NTuple::Item<int>  m_ext_eta;
  NTuple::Item<double> m_mc_costheta_p;
  NTuple::Item<double> m_mc_phi_p;
  NTuple::Item<double> m_mc_px_p;
  NTuple::Item<double> m_mc_py_p;
  NTuple::Item<double> m_mc_pz_p;
  NTuple::Item<double> m_mc_rho_p;
  NTuple::Item<double> m_mc_pxy_p;

  NTuple::Item<double> m_mc_costheta_ap;
  NTuple::Item<double> m_mc_phi_ap;
  NTuple::Item<double> m_mc_px_ap;
  NTuple::Item<double> m_mc_py_ap;
  NTuple::Item<double> m_mc_pz_ap;
  NTuple::Item<double> m_mc_rho_ap;
  NTuple::Item<double> m_mc_pxy_ap;

  NTuple::Item<double> m_mc_costheta_gamma;
  NTuple::Item<double> m_mc_phi_gamma;
  NTuple::Item<double> m_mc_px_gamma;
  NTuple::Item<double> m_mc_py_gamma;
  NTuple::Item<double> m_mc_pz_gamma;
  NTuple::Item<double> m_mc_rho_gamma;
  NTuple::Item<double> m_mc_pxy_gamma;


  NTuple::Item<double> m_mc_u_p;
  NTuple::Item<double> m_mc_u_ap;


  NTuple::Item<double> m_4c_oksq;
  NTuple::Item<double> m_4c_chisq;

  NTuple::Item<double> m_4c_p_mass;
  NTuple::Item<double> m_4c_ap_mass;
  NTuple::Item<double> m_4c_pap_mass;
  NTuple::Item<double> m_4c_pgamma_mass;
  NTuple::Item<double> m_4c_apgamma_mass;
  NTuple::Item<double> m_4c_p_pxy;
  NTuple::Item<double> m_4c_ap_pxy;
  NTuple::Item<double> m_4c_pap_pxy;
  NTuple::Item<double> m_4c_pgamma_pxy;
  NTuple::Item<double> m_4c_apgamma_pxy;
  NTuple::Item<double> m_4c_gamma_pxy;
  NTuple::Item<double> m_4c_p_p;
  NTuple::Item<double> m_4c_ap_p;
  NTuple::Item<double> m_4c_pap_p;
  NTuple::Item<double> m_4c_pgamma_p;
  NTuple::Item<double> m_4c_apgamma_p;
  NTuple::Item<double> m_4c_gamma_p;
  NTuple::Item<double> m_4c_p_costheta;
  NTuple::Item<double> m_4c_ap_costheta;
  NTuple::Item<double> m_4c_pap_costheta;
  NTuple::Item<double> m_4c_pgamma_costheta;
  NTuple::Item<double> m_4c_apgamma_costheta;
  NTuple::Item<double> m_4c_gamma_costheta;
  NTuple::Item<double> m_4c_p_phi;
  NTuple::Item<double> m_4c_ap_phi;
  NTuple::Item<double> m_4c_pap_phi;
  NTuple::Item<double> m_4c_pgamma_phi;
  NTuple::Item<double> m_4c_apgamma_phi;
  NTuple::Item<double> m_4c_gamma_phi;


  NTuple::Item<double> 	m_kal_p_p;
  NTuple::Item<double> 	m_kal_p_costheta;
  NTuple::Item<double> 	m_kal_p_e;
  NTuple::Item<double> 	m_kal_p_mass;
  NTuple::Item<double>   m_kal_p_pxy;
  NTuple::Item<double> 	m_kal_ap_p;
  NTuple::Item<double> 	m_kal_ap_costheta;
  NTuple::Item<double> 	m_kal_ap_e;
  NTuple::Item<double> 	m_kal_ap_mass;
  NTuple::Item<double>   m_kal_ap_pxy;
  NTuple::Item<double> 	m_kal_pap_p;
  NTuple::Item<double> 	m_kal_pap_costheta;
  NTuple::Item<double> 	m_kal_pap_e;
  NTuple::Item<double> 	m_kal_pap_mass;
  NTuple::Item<double>   m_kal_pap_pxy;
  NTuple::Item<double> 	m_kal_gamma_p;
  NTuple::Item<double> 	m_kal_gamma_pxy;
  NTuple::Item<double> 	m_kal_gamma_costheta;
  NTuple::Item<double> 	m_kal_pgamma_p;
  NTuple::Item<double> 	m_kal_pgamma_costheta;
  NTuple::Item<double> 	m_kal_pgamma_e;
  NTuple::Item<double> 	m_kal_pgamma_mass;
  NTuple::Item<double>   m_kal_pgamma_pxy;
  NTuple::Item<double> 	m_kal_apgamma_p;
  NTuple::Item<double> 	m_kal_apgamma_costheta;
  NTuple::Item<double> 	m_kal_apgamma_e;
  NTuple::Item<double> 	m_kal_apgamma_mass;
  NTuple::Item<double>   m_kal_apgamma_pxy;
  NTuple::Item<double> 	m_kal_jpsi_p;
  NTuple::Item<double> 	m_kal_jpsi_costheta;
  NTuple::Item<double> 	m_kal_jpsi_e;
  NTuple::Item<double> 	m_kal_jpsi_mass;
  NTuple::Item<double> 	m_kal_jpsi_pxy;



   NTuple::Item<double> 	m_bak_p_p;
   NTuple::Item<double> 	m_bak_p_costheta;
   NTuple::Item<double> 	m_bak_p_e;
   NTuple::Item<double> 	m_bak_p_mass;
   NTuple::Item<double>   m_bak_p_pxy;
   NTuple::Item<double> 	m_bak_ap_p;
   NTuple::Item<double> 	m_bak_ap_costheta;
   NTuple::Item<double> 	m_bak_ap_e;
   NTuple::Item<double> 	m_bak_ap_mass;
   NTuple::Item<double>   m_bak_ap_pxy;
   NTuple::Item<double> 	m_bak_pap_p;
   NTuple::Item<double> 	m_bak_pap_costheta;
   NTuple::Item<double> 	m_bak_pap_e;
   NTuple::Item<double> 	m_bak_pap_mass;
   NTuple::Item<double>   m_bak_pap_pxy;
   NTuple::Item<double> 	m_bak_gamma_p;
   NTuple::Item<double> 	m_bak_gamma_costheta;
   NTuple::Item<double>   m_bak_gamma_pxy;
   NTuple::Item<double>   m_bak_gamma_u;
   NTuple::Item<double>   m_bak_gamma_pt2;
   NTuple::Item<double>   m_bak_gamma_e;
   NTuple::Item<double> 	m_bak_pgamma_p;
   NTuple::Item<double> 	m_bak_pgamma_costheta;
   NTuple::Item<double> 	m_bak_pgamma_e;
   NTuple::Item<double> 	m_bak_pgamma_mass;
   NTuple::Item<double>   m_bak_pgamma_pxy;
   NTuple::Item<double> 	m_bak_apgamma_p;
   NTuple::Item<double> 	m_bak_apgamma_costheta;
   NTuple::Item<double> 	m_bak_apgamma_e;
   NTuple::Item<double> 	m_bak_apgamma_mass;
   NTuple::Item<double>   m_bak_apgamma_pxy;

  NTuple::Item<double> m_diff_costheta_ap;
  NTuple::Item<double> m_diff_costheta_p;
  NTuple::Item<double> m_diff_costheta_gamma;
  NTuple::Item<double> m_diff_p_ap;
  NTuple::Item<double> m_diff_p_p;
  NTuple::Item<double> m_diff_p_gamma;
  NTuple::Item<double> m_diff_pxy_ap;
  NTuple::Item<double> m_diff_pxy_p;
  NTuple::Item<double> m_diff_pxy_gamma;
  NTuple::Item<double> m_diff_phi_p;
  NTuple::Item<double> m_diff_phi_ap;
  NTuple::Item<double> m_diff_phi_gamma;


     NTuple::Item<double>  m_chie_p; 
     NTuple::Item<double>  m_chimu_p; 
     NTuple::Item<double>  m_chipi_p ;
     NTuple::Item<double>  m_chik_p ;
     NTuple::Item<double>  m_chip_p ;
     NTuple::Item<double>  m_ghit_p;
     NTuple::Item<double>  m_thit_p;
     NTuple::Item<double>  m_probPH_p;
     NTuple::Item<double>  m_normPH_p;

     NTuple::Item<double>  m_chie_ap; 
     NTuple::Item<double>  m_chimu_ap; 
     NTuple::Item<double>  m_chipi_ap ;
     NTuple::Item<double>  m_chik_ap ;
     NTuple::Item<double>  m_chip_ap ;
     NTuple::Item<double>  m_ghit_ap;
     NTuple::Item<double>  m_thit_ap;
     NTuple::Item<double>  m_probPH_ap;
     NTuple::Item<double>  m_normPH_ap;

     NTuple::Item<double>  m_prob_e_p;
     NTuple::Item<double>  m_prob_pi_p;
     NTuple::Item<double>  m_prob_k_p;
     NTuple::Item<double>  m_prob_p_p;
     NTuple::Item<double>  m_prob_e_m;
     NTuple::Item<double>  m_prob_pi_m;
     NTuple::Item<double>  m_prob_k_m;
     NTuple::Item<double>  m_prob_p_m;


  //  MC truth info
	const int kMaxId; 

	//const ecms_4p
	HepLorentzVector jpsi_cms;

  //bool m_isMonteCarlo; 

  // define Histograms
  TH1F* h_evtflw; 
	int N1,N2,N3,N4,N5,N6;
  
  // define Trees
  TTree* m_tree;



  // functions		
  void book_tree(); 
  void before_execute(); 
  bool buildGammaPP();
  bool saveGenInfo(SmartDataPtr<Event::McParticleCol>); 
  bool passVertexSelection(CLHEP::Hep3Vector,
			   RecMdcKalTrack* ); 
  CLHEP::Hep3Vector getOrigin();
  int selectChargedTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
			  std::vector<int> &,
			  std::vector<int> &); 
	bool ExtraPi0(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int>);
	bool ExtraEta(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int>);
  void calcTrackPID(EvtRecTrackIterator,
		    double& ,
		    double&	,
				double& ,
				double&);
  bool PionPid(SmartDataPtr<EvtRecTrackCol>,
			      int,
			      int);
  bool PPid(SmartDataPtr<EvtRecTrackCol>,
			      int,
			      int);
	void selectNeutralTracks(SmartDataPtr<EvtRecEvent>,
						SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int> &);
	void calInvi(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int> ,
						HepLorentzVector,
						HepLorentzVector &,
						HepLorentzVector &);
	bool GoodNeutral(SmartDataPtr<EvtRecTrackCol>,
			      int);
  //void storeShowersInfo(SmartDataPtr<EvtRecTrackCol>,
	//		      int );
  void storeTracksInfo(HepLorentzVector ,
						HepLorentzVector,
						HepLorentzVector);
  bool KinematicFit(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int> ,
			      int,
			      int,
						HepLorentzVector &,
						HepLorentzVector &,
						HepLorentzVector &,
						int &);



 void GetKalTrackInfo( SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
                      int iPCharge, 
                      int iMCharge,
                      int ig,
                      HepLorentzVector &pP4, 
                      HepLorentzVector &mP4,
                      HepLorentzVector &gP4); 

void GetBakTrackInfo( HepLorentzVector &pP4,
                      HepLorentzVector &mP4,
                      HepLorentzVector &gP4);

bool saveAllGenInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol);

bool saveGenParInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol);


void storeDedxInfo(
           SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
           const int iPCharge,
           const int iMCharge);


void store4CTracksInfo(HepLorentzVector pP4,
                    HepLorentzVector mP4,
                    HepLorentzVector gP4);

void storeKalTracksInfo(HepLorentzVector pP4,
                    HepLorentzVector mP4,
                    HepLorentzVector gP4);

void storeBakTracksInfo(HepLorentzVector pP4, 
                    HepLorentzVector mP4, 
                    HepLorentzVector gP4);


}; 



//
// module declare
//


DECLARE_ALGORITHM_FACTORY( GammaPP )

DECLARE_FACTORY_ENTRIES( GammaPPAlg ) { 
  DECLARE_ALGORITHM(GammaPP);
}

LOAD_FACTORY_ENTRIES( GammaPPAlg )

//
// constants
//

const double PION_MASS = 0.139570;
const double E_MASS = 0.000511;
const double P_MASS = 0.9382723;
const double PION0_MASS = 0.134977;
const double JPSI_MASS = 3.096916;

typedef std::vector<int> Vint;
//
// member functions
//

GammaPP::GammaPP(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
	kMaxId(500)
	//m_pdgid(0),
	//m_motherindex(0)
	{			
  //declareProperty("IsMonteCarlo", m_isMonteCarlo);
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=10.0);
  declareProperty("ChaCosthetaCut", m_cha_costheta_cut=0.93);
  declareProperty("MinEstCut", m_min_emctime=0.0);
  declareProperty("MaxEstCut", m_max_emctime=14.0);
  declareProperty("GammaCosCut",  m_gammaCosCut= 0.8); 
  declareProperty("CosthetaBarrelMax", m_costheta_barrel_max=0.8);
  declareProperty("CosthetaEndcapMin", m_costheta_endcap_min=0.86);
  declareProperty("CosthetaEndcapMax", m_costheta_endcap_max=0.92);
  declareProperty("EnergyBarrelMin", m_energy_barrel_min=0.025); 
  declareProperty("EnergyEndcapMin", m_energy_endcap_min=0.050); 
  declareProperty("PhotonIsoAngleMin", m_photon_iso_angle_min=20);
  declareProperty("ProbPionMin", m_prob_pion_min=0.001);
  declareProperty("EPRation", m_ep_electron_min=0.8);
  declareProperty("RhoPionMax", m_cha_momentum_max_cut=0.450);  //GeV
}

StatusCode GammaPP::initialize(){
	//jpsi_cms.set(3.097*sin(0.011),0.,0.,3.097);     
	jpsi_cms.set(3.097*sin(0.011),0.,0.,3.097);     
	N1=0;
	N2=0;
	N3=0;
	N4=0;
	N5=0;
	N6=0;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;


  book_tree(); 

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;
}

StatusCode GammaPP::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

	before_execute();
  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader) return StatusCode::FAILURE;

  m_run = eventHeader->runNumber();
  m_rec = eventHeader->eventNumber();
  if(m_rec%1000==0) cout<<"event = "<<m_rec<<endl;

  
  if(buildGammaPP() == true) m_tuple1->write();

  return StatusCode::SUCCESS; 
}

StatusCode GammaPP::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
	h_evtflw->SetBinContent(1,N1);
	h_evtflw->SetBinContent(2,N2);
	h_evtflw->SetBinContent(3,N3);
	h_evtflw->SetBinContent(4,N4);
	h_evtflw->SetBinContent(5,N5);
	h_evtflw->SetBinContent(6,N6);

  h_evtflw->GetXaxis()->SetBinLabel(1, "raw");
  h_evtflw->GetXaxis()->SetBinLabel(2, "nCha");
  h_evtflw->GetXaxis()->SetBinLabel(3, "Pid");
  h_evtflw->GetXaxis()->SetBinLabel(4, "nGam");
  h_evtflw->GetXaxis()->SetBinLabel(5, "extpi0");
  h_evtflw->GetXaxis()->SetBinLabel(6, "kitfit");

  return StatusCode::SUCCESS;
}

GammaPP::~GammaPP() {
}

void GammaPP::before_execute(){
}

void GammaPP::book_tree() {
  StatusCode status;
  NTuplePtr nt1(ntupleSvc(), "FILE/GammaPP");
  if ( nt1 ) m_tuple1 = nt1;
  else {
   m_tuple1 = ntupleSvc()->book ("FILE1/GammaPP", CLID_ColumnWiseTuple, " N-Tuple example");
    if ( m_tuple1 )    {
			//MC info 
      status = m_tuple1->addItem ("run",   m_run);
      status = m_tuple1->addItem ("rec", m_rec);
      status = m_tuple1->addItem ("NGch",   m_ngch, 0, 40); 
      status = m_tuple1->addItem ("NCharge",   m_nCharge);
      status = m_tuple1->addItem ("NGam",   m_nGam);
      status = m_tuple1->addItem ("indexmc",   m_idxmc,0,100);
      status = m_tuple1->addIndexedItem ("pdgid",   m_idxmc,m_pdgid);
      status = m_tuple1->addIndexedItem ("motheridx",   m_idxmc,m_motheridx);
      status = m_tuple1->addIndexedItem("trkidx",    m_idxmc, m_trkidx);
      status = m_tuple1->addItem ("mc_costheta_p",   m_mc_costheta_p);
      status = m_tuple1->addItem ("mc_phi_p",   m_mc_phi_p);
      status = m_tuple1->addItem ("mc_px_p",   m_mc_px_p);
      status = m_tuple1->addItem ("mc_py_p",   m_mc_py_p);
      status = m_tuple1->addItem ("mc_pz_p",   m_mc_pz_p);
      status = m_tuple1->addItem ("mc_rho_p",   m_mc_rho_p);
      status = m_tuple1->addItem ("mc_pxy_p",   m_mc_pxy_p);

      status = m_tuple1->addItem ("mc_costheta_ap",   m_mc_costheta_ap);
      status = m_tuple1->addItem ("mc_phi_ap",   m_mc_phi_ap);
      status = m_tuple1->addItem ("mc_px_ap",   m_mc_px_ap);
      status = m_tuple1->addItem ("mc_py_ap",   m_mc_py_ap);
      status = m_tuple1->addItem ("mc_pz_ap",   m_mc_pz_ap);
      status = m_tuple1->addItem ("mc_rho_ap",   m_mc_rho_ap);
      status = m_tuple1->addItem ("mc_pxy_ap",   m_mc_pxy_ap);

      status = m_tuple1->addItem ("mc_costheta_gamma",   m_mc_costheta_gamma);
      status = m_tuple1->addItem ("mc_phi_gamma",   m_mc_phi_gamma);
      status = m_tuple1->addItem ("mc_px_gamma",   m_mc_px_gamma);
      status = m_tuple1->addItem ("mc_py_gamma",   m_mc_py_gamma);
      status = m_tuple1->addItem ("mc_pz_gamma",   m_mc_pz_gamma);
      status = m_tuple1->addItem ("mc_rho_gamma",   m_mc_rho_gamma);
      status = m_tuple1->addItem ("mc_pxy_gamma",   m_mc_pxy_gamma);

      status = m_tuple1->addItem ("mc_u_p",   m_mc_u_p);
      status = m_tuple1->addItem ("mc_u_ap",   m_mc_u_ap);
      status = m_tuple1->addItem ("4c_chisq",   m_4c_chisq);

      status = m_tuple1->addItem ("4c_p_mass",   m_4c_p_mass);
      status = m_tuple1->addItem ("4c_ap_mass",   m_4c_ap_mass);
      status = m_tuple1->addItem ("4c_p_pxy",   m_4c_p_pxy);
      status = m_tuple1->addItem ("4c_ap_pxy",   m_4c_ap_pxy);
      status = m_tuple1->addItem ("4c_gamma_pxy",   m_4c_gamma_pxy);
      status = m_tuple1->addItem ("4c_p_p",   m_4c_p_p);
      status = m_tuple1->addItem ("4c_ap_p",   m_4c_ap_p);
      status = m_tuple1->addItem ("4c_gamma_p",   m_4c_gamma_p);
      status = m_tuple1->addItem ("4c_p_costheta",   m_4c_p_costheta);
      status = m_tuple1->addItem ("4c_ap_costheta",   m_4c_ap_costheta);
      status = m_tuple1->addItem ("4c_gamma_costheta",   m_4c_gamma_costheta);
      status = m_tuple1->addItem ("4c_pap_mass",   m_4c_pap_mass);
      status = m_tuple1->addItem ("4c_pap_pxy",   m_4c_pap_pxy);
      status = m_tuple1->addItem ("4c_pap_costheta",   m_4c_pap_costheta);
      status = m_tuple1->addItem ("4c_pap_p",   m_4c_pap_p);
      status = m_tuple1->addItem ("4c_pgamma_mass",   m_4c_pgamma_mass);
      status = m_tuple1->addItem ("4c_pgamma_pxy",   m_4c_pgamma_pxy);
      status = m_tuple1->addItem ("4c_pgamma_costheta",   m_4c_pgamma_costheta);
      status = m_tuple1->addItem ("4c_pgamma_p",   m_4c_pgamma_p);
      status = m_tuple1->addItem ("4c_apgamma_mass",   m_4c_apgamma_mass);
      status = m_tuple1->addItem ("4c_apgamma_pxy",   m_4c_apgamma_pxy);
      status = m_tuple1->addItem ("4c_apgamma_costheta",   m_4c_apgamma_costheta);
      status = m_tuple1->addItem ("4c_apgamma_p",   m_4c_apgamma_p);
      status = m_tuple1->addItem ("kal_gamma_p",   m_kal_gamma_p);
      status = m_tuple1->addItem ("kal_gamma_pxy",   m_kal_gamma_pxy);
      status = m_tuple1->addItem ("kal_gamma_costheta",   m_kal_gamma_costheta);
      status = m_tuple1->addItem ("kal_p_mass",   m_kal_p_mass);
      status = m_tuple1->addItem ("kal_p_p",   m_kal_p_p);
      status = m_tuple1->addItem ("kal_p_pxy",   m_kal_p_pxy);
      status = m_tuple1->addItem ("kal_p_p_costheta",   m_kal_p_costheta);
      status = m_tuple1->addItem ("kal_p_e",   m_kal_p_e);
      status = m_tuple1->addItem ("kal_ap_p",   m_kal_ap_p);
      status = m_tuple1->addItem ("kal_ap_mass",   m_kal_ap_mass);
      status = m_tuple1->addItem ("kal_ap_pxy",   m_kal_ap_pxy);
      status = m_tuple1->addItem ("kal_ap_costheta",   m_kal_ap_costheta);
      status = m_tuple1->addItem ("kal_ap_e",   m_kal_ap_e);
      status = m_tuple1->addItem ("kal_pap_p",   m_kal_pap_p);
      status = m_tuple1->addItem ("kal_pap_mass",   m_kal_pap_mass);
      status = m_tuple1->addItem ("kal_pap_pxy",   m_kal_pap_pxy);
      status = m_tuple1->addItem ("kal_pap_costheta",   m_kal_pap_costheta);
      status = m_tuple1->addItem ("kal_pap_e",   m_kal_pap_e);
      status = m_tuple1->addItem ("kal_pgamma_p",   m_kal_pgamma_p);
      status = m_tuple1->addItem ("kal_pgamma_mass",   m_kal_pgamma_mass);
      status = m_tuple1->addItem ("kal_pgamma_pxy",   m_kal_pgamma_pxy);
      status = m_tuple1->addItem ("kal_pgamma_costheta",   m_kal_pgamma_costheta);
      status = m_tuple1->addItem ("kal_pgamma_e",   m_kal_pgamma_e);
      status = m_tuple1->addItem ("kal_apgamma_p",   m_kal_apgamma_p);
      status = m_tuple1->addItem ("kal_apgamma_mass",   m_kal_apgamma_mass);
      status = m_tuple1->addItem ("kal_apgamma_pxy",   m_kal_apgamma_pxy);
      status = m_tuple1->addItem ("kal_apgamma_costheta",   m_kal_apgamma_costheta);
      status = m_tuple1->addItem ("kal_apgamma_e",   m_kal_apgamma_e);
      status = m_tuple1->addItem ("kal_jpsi_mass",   m_kal_jpsi_mass);
      status = m_tuple1->addItem ("kal_jpsi_costheta",   m_kal_jpsi_costheta);
      status = m_tuple1->addItem ("kal_jpsi_pxy",   m_kal_jpsi_pxy);
      status = m_tuple1->addItem ("kal_jpsi_p",   m_kal_jpsi_p);
      status = m_tuple1->addItem ("kal_jpsi_e",   m_kal_jpsi_e);
      status = m_tuple1->addItem ("bak_p_mass",   m_bak_p_mass);
      status = m_tuple1->addItem ("bak_p_p",   m_bak_p_p);
      status = m_tuple1->addItem ("bak_p_pxy",   m_bak_p_pxy);
      status = m_tuple1->addItem ("bak_p_p_costheta",   m_bak_p_costheta);
      status = m_tuple1->addItem ("bak_p_e",   m_bak_p_e);
      status = m_tuple1->addItem ("bak_ap_p",   m_bak_ap_p);
      status = m_tuple1->addItem ("bak_ap_mass",   m_bak_ap_mass);
      status = m_tuple1->addItem ("bak_ap_pxy",   m_bak_ap_pxy);
      status = m_tuple1->addItem ("bak_ap_costheta",   m_bak_ap_costheta);
      status = m_tuple1->addItem ("bak_ap_e",   m_bak_ap_e);
      status = m_tuple1->addItem ("bak_pap_p",   m_bak_pap_p);
      status = m_tuple1->addItem ("bak_pap_mass",   m_bak_pap_mass);
      status = m_tuple1->addItem ("bak_pap_pxy",   m_bak_pap_pxy);
      status = m_tuple1->addItem ("bak_pap_costheta",   m_bak_pap_costheta);
      status = m_tuple1->addItem ("bak_pap_e",   m_bak_pap_e);
      status = m_tuple1->addItem ("bak_pgamma_p",   m_bak_pgamma_p);
      status = m_tuple1->addItem ("bak_pgamma_mass",   m_bak_pgamma_mass);
      status = m_tuple1->addItem ("bak_pgamma_pxy",   m_bak_pgamma_pxy);
      status = m_tuple1->addItem ("bak_pgamma_costheta",   m_bak_pgamma_costheta);
      status = m_tuple1->addItem ("bak_pgamma_e",   m_bak_pgamma_e);
      status = m_tuple1->addItem ("bak_apgamma_p",   m_bak_apgamma_p);
      status = m_tuple1->addItem ("bak_apgamma_mass",   m_bak_apgamma_mass);
      status = m_tuple1->addItem ("bak_apgamma_pxy",   m_bak_apgamma_pxy);
      status = m_tuple1->addItem ("bak_apgamma_costheta",   m_bak_apgamma_costheta);
      status = m_tuple1->addItem ("bak_apgamma_e",   m_bak_apgamma_e);
      status = m_tuple1->addItem ("bak_gamma_p",   m_bak_gamma_p);
      status = m_tuple1->addItem ("bak_gamma_pxy",   m_bak_gamma_pxy);
      status = m_tuple1->addItem ("bak_gamma_costheta",   m_bak_gamma_costheta);
      status = m_tuple1->addItem ("bak_gamma_u",   m_bak_gamma_u);
      status = m_tuple1->addItem ("bak_gamma_pt2",   m_bak_gamma_pt2);
      status = m_tuple1->addItem ("bak_gamma_e",   m_bak_gamma_e);
      status = m_tuple1->addItem ("diff_costheta_p",   m_diff_costheta_p);
      status = m_tuple1->addItem ("diff_pxy_p",   m_diff_pxy_p);
      status = m_tuple1->addItem ("diff_costheta_p",   m_diff_costheta_p);
      status = m_tuple1->addItem ("diff_pxy_p",   m_diff_pxy_p);
      status = m_tuple1->addItem ("diff_p_p",   m_diff_p_p);
      status = m_tuple1->addItem ("diff_costheta_ap",   m_diff_costheta_ap);
      status = m_tuple1->addItem ("diff_pxy_ap",   m_diff_pxy_ap);
      status = m_tuple1->addItem ("diff_p_ap",   m_diff_p_ap);
      status = m_tuple1->addItem ("diff_costheta_gamma",   m_diff_costheta_gamma);
      status = m_tuple1->addItem ("diff_pxy_gamma",   m_diff_pxy_gamma);
      status = m_tuple1->addItem ("diff_p_gamma",   m_diff_p_gamma);
      status = m_tuple1->addItem ("chie_p",   m_chie_p);
      status = m_tuple1->addItem ("chimu_p",   m_chimu_p);
      status = m_tuple1->addItem ("chipi_p",   m_chipi_p);
      status = m_tuple1->addItem ("chik_p",   m_chik_p);
      status = m_tuple1->addItem ("chip_p",   m_chip_p);
      status = m_tuple1->addItem ("probPH_p",   m_probPH_p);
      status = m_tuple1->addItem ("normPH_p",   m_normPH_p);
      status = m_tuple1->addItem ("thit_p",   m_thit_p);
      status = m_tuple1->addItem ("ghit_p",   m_ghit_p);
      status = m_tuple1->addItem ("chie_ap",   m_chie_ap);
      status = m_tuple1->addItem ("chimu_ap",   m_chimu_ap);
      status = m_tuple1->addItem ("chipi_ap",   m_chipi_ap);
      status = m_tuple1->addItem ("chik_ap",   m_chik_ap);
      status = m_tuple1->addItem ("chip_ap",   m_chip_ap);
      status = m_tuple1->addItem ("probPH_ap",   m_probPH_ap);
      status = m_tuple1->addItem ("normPH_ap",   m_normPH_ap);
      status = m_tuple1->addItem ("ghit_ap",   m_ghit_ap);
      status = m_tuple1->addItem ("thit_ap",   m_thit_ap);
      status = m_tuple1->addItem ("ext_pi0",   m_ext_pi0);
      status = m_tuple1->addItem ("ext_eta",   m_ext_eta);
      status = m_tuple1->addItem ("prob_e_p",   m_prob_e_p);
      status = m_tuple1->addItem ("prob_pi_p",   m_prob_pi_p);
      status = m_tuple1->addItem ("prob_k_p",   m_prob_k_p);
      status = m_tuple1->addItem ("prob_p_p",   m_prob_p_p);
      status = m_tuple1->addItem ("prob_e_m",   m_prob_e_m);
      status = m_tuple1->addItem ("prob_pi_m",   m_prob_pi_m);
      status = m_tuple1->addItem ("prob_k_m",   m_prob_k_m);
      status = m_tuple1->addItem ("prob_p_m",   m_prob_p_m);

		}
    else    {
      return ;
    }
	}
	h_evtflw= new TH1F("hevtflw", "hevtflw", 6, 0, 6);
  return ;
}


bool GammaPP::buildGammaPP() {
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
 	if (m_run < 0) {
		if(!saveAllGenInfo(mcParticleCol)) return false;
		if(!saveGenParInfo(mcParticleCol)) return false;
		}
	N1++; //raw

  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return false;
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return false;

	std::vector<int> iPGood,iMGood;
  if(selectChargedTracks(evtRecEvent, evtRecTrkCol, iPGood, iMGood) != 2) return false;

	if (iPGood.size() != 1) return false;
	if (iMGood.size() != 1) return false;
	int iPCharge,iMCharge;
	iPCharge=iPGood[0];
	iMCharge=iMGood[0];
	N2++; // N_{Good} == 2 && N_{P} == 1 && N_{M} == 1

	if(PPid(evtRecTrkCol, iPCharge, iMCharge) == false) return false; // N_{p} == 2
	N3++; 
	std::vector<int> iGam;
	selectNeutralTracks(evtRecEvent, evtRecTrkCol, iGam);
	if(iGam.size()==0) return false;
	N4++; 
	m_ext_pi0=0;  m_ext_eta=0;
	if(ExtraPi0(evtRecTrkCol, iGam ) == true) m_ext_pi0 = 1;
	if(ExtraEta(evtRecTrkCol, iGam ) == true) m_ext_eta = 1;
	N5++; 
	int ig;
	HepLorentzVector pP4, mP4, gP4;
	if(KinematicFit(evtRecTrkCol, iGam, iPCharge, iMCharge, pP4, mP4, gP4, ig) == false) return false;
	N6++; 
	store4CTracksInfo(pP4, mP4, gP4);
	m_nGam = iGam.size();

  GetKalTrackInfo(evtRecTrkCol,iPCharge, iMCharge,ig,pP4, mP4, gP4);
  storeKalTracksInfo(pP4, mP4, gP4);

  GetBakTrackInfo(pP4, mP4, gP4);
  storeBakTracksInfo(pP4, mP4, gP4);
  
  storeDedxInfo(evtRecTrkCol,iPCharge, iMCharge);

	return true;
}

void GammaPP::GetKalTrackInfo( SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
                      int iPCharge, 
                      int iMCharge,
                      int ig,
                      HepLorentzVector &pP4, 
                      HepLorentzVector &mP4,
                      HepLorentzVector &gP4){
  MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "begin GetKalTrackInfo" << endmsg;
  EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPCharge;
  EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMCharge;
  EvtRecTrackIterator itTrk_g = evtRecTrkCol->begin() + ig;

if(!(*itTrk_p)->isMdcKalTrackValid()) return;
if(!(*itTrk_m)->isMdcKalTrackValid()) return;


    RecMdcKalTrack* mdcKalTrk_p = (*itTrk_p)->mdcKalTrack();
    pP4.setPx(mdcKalTrk_p->px());
    pP4.setPy(mdcKalTrk_p->py());
    pP4.setPz(mdcKalTrk_p->pz());
    double p3 = pP4.mag();
    pP4.setE(sqrt(p3*p3+P_MASS*P_MASS));

    RecMdcKalTrack* mdcKalTrk_m = (*itTrk_m)->mdcKalTrack();
    mP4.setPx(mdcKalTrk_m->px());
    mP4.setPy(mdcKalTrk_m->py());
    mP4.setPz(mdcKalTrk_m->pz());
    p3 = mP4.mag();
    mP4.setE(sqrt(p3*p3+P_MASS*P_MASS));

    RecEmcShower* emcTrk = (*itTrk_g)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    gP4.setPx(eraw*sin(the)*cos(phi));
    gP4.setPy(eraw*sin(the)*sin(phi));
    gP4.setPz(eraw*cos(the));
    gP4.setE(eraw);

    log << MSG::DEBUG << "end GetKalTrackInfo" << endmsg;
}

void GammaPP::GetBakTrackInfo( HepLorentzVector &pP4, 
                      HepLorentzVector &mP4,
                      HepLorentzVector &gP4){
      MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin GetBakTrackInfo" << endmsg;
      HepLorentzVector temp_pP4,temp_mP4,temp_gP4;

      temp_pP4 = jpsi_cms - mP4 - gP4;
      temp_mP4 = jpsi_cms - pP4 - gP4;
      temp_gP4 = jpsi_cms - mP4 - pP4;
      
      pP4 = temp_pP4;
      mP4 = temp_mP4;
      gP4 = temp_gP4;
          log << MSG::DEBUG << "end GetBakTrackInfo" << endmsg;
}

bool GammaPP::saveAllGenInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol) {
      MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin saveAllGenInfo" << endmsg;
	if (!mcParticleCol){
		std::cout << "Could not retrieve McParticelCol" << std::endl;
		return false;
	}

      bool jpsiDecay = false;
      int rootIndex = -1;
  Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
int m_numParticle = 0;
  for (; iter_mc != mcParticleCol->end(); iter_mc++){
    if ((*iter_mc)->primaryParticle()) continue;
    if (!(*iter_mc)->decayFromGenerator()) continue;

    if ((*iter_mc)->particleProperty()==443)
    {   
        jpsiDecay = true;
        rootIndex = (*iter_mc)->trackIndex();
    }   
    if (!jpsiDecay) continue;
    int trkidx = (*iter_mc)->trackIndex() - rootIndex;
    int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
    int pdgid = (*iter_mc)->particleProperty();
    m_trkidx[m_numParticle] = trkidx;
    m_pdgid[m_numParticle] = pdgid;
    m_motheridx[m_numParticle] = mcidx;
    m_numParticle += 1;
		}
          log << MSG::DEBUG << "end saveAllGenInfo" << endmsg;
  return true;
}

bool GammaPP::saveGenParInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol){
  MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin saveGenParInfo" << endmsg;
	if (!mcParticleCol){
		std::cout << "Could not retrieve McParticelCol" << std::endl;
		return false;
	}
Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
  for(;iter_mc!=mcParticleCol->end();iter_mc++){
      HepLorentzVector initialFourMomentum=(*iter_mc)->initialFourMomentum();
      HepLorentzVector finalPosition=(*iter_mc)->finalPosition();

      if ((*iter_mc)->particleProperty()==2212){
    m_mc_costheta_p=finalPosition.cosTheta();
    m_mc_phi_p=finalPosition.phi();
    m_mc_px_p=initialFourMomentum.px();
    m_mc_py_p=initialFourMomentum.py();
    m_mc_pz_p=initialFourMomentum.pz();
    m_mc_rho_p=initialFourMomentum.rho();
    m_mc_pxy_p=sqrt(m_mc_px_p*m_mc_px_p+m_mc_py_p*m_mc_py_p);
      }

      if ((*iter_mc)->particleProperty()==-2212){
    m_mc_costheta_ap=finalPosition.cosTheta();
    m_mc_phi_ap=finalPosition.phi();
    m_mc_px_ap=initialFourMomentum.px();
    m_mc_py_ap=initialFourMomentum.py();
    m_mc_pz_ap=initialFourMomentum.pz();
    m_mc_rho_ap=initialFourMomentum.rho();
    m_mc_pxy_p=sqrt(m_mc_px_p*m_mc_px_p+m_mc_py_p*m_mc_py_p);
    m_mc_pxy_ap=sqrt(m_mc_px_ap*m_mc_px_ap+m_mc_py_ap*m_mc_py_ap);
      }

      if ((*iter_mc)->particleProperty()==22){
    m_mc_costheta_gamma=finalPosition.cosTheta();
    m_mc_phi_gamma=finalPosition.phi();
    m_mc_px_gamma=initialFourMomentum.px();
    m_mc_py_gamma=initialFourMomentum.py();
    m_mc_pz_gamma=initialFourMomentum.pz();
   m_mc_rho_gamma=initialFourMomentum.rho();
    m_mc_pxy_p=sqrt(m_mc_px_p*m_mc_px_p+m_mc_py_p*m_mc_py_p);
    m_mc_pxy_gamma=sqrt(m_mc_px_gamma*m_mc_px_gamma+m_mc_py_gamma*m_mc_py_gamma);
      }
}


          log << MSG::DEBUG << "end saveGenParInfo" << endmsg;
  return true;
}


CLHEP::Hep3Vector GammaPP::getOrigin() {
  CLHEP::Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double *dbv = vtxsvc->PrimaryVertex(); 
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  return xorigin; 
}

bool GammaPP::passVertexSelection(CLHEP::Hep3Vector xorigin,
				    RecMdcKalTrack* mdcTrk ) {
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);
  VFHelix helixip(point0,a,Ea);
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  
  double vz0 = vecipa[3];
  double vr0 = vecipa[0];
  
  if(fabs(vz0) >= m_vz0cut) return false;
  if(fabs(vr0) >= m_vr0cut) return false;
  
  return true;
}


void GammaPP::storeDedxInfo(
           SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
           const int iPCharge,
           const int iMCharge) {
  MsgStream log(msgSvc(), name());
      log << MSG::DEBUG << "begin sotreDedxInfo" << endmsg;
      EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPCharge;
      if(!(*itTrk_p)->isMdcDedxValid()) return ;
      RecMdcDedx* dedxTrk_p = (*itTrk_p)->mdcDedx();
    m_chie_p = dedxTrk_p->chiE();
    m_chimu_p = dedxTrk_p->chiMu();
    m_chipi_p = dedxTrk_p->chiPi();
    m_chik_p = dedxTrk_p->chiK();
    m_chip_p = dedxTrk_p->chiP();
    m_ghit_p = dedxTrk_p->numGoodHits();
    m_thit_p = dedxTrk_p->numTotalHits();
    m_probPH_p = dedxTrk_p->probPH();
    m_normPH_p = dedxTrk_p->normPH();


      EvtRecTrackIterator itTrk_ap = evtRecTrkCol->begin() + iMCharge;
      if(!(*itTrk_ap)->isMdcDedxValid()) return ;
      RecMdcDedx* dedxTrk_ap = (*itTrk_ap)->mdcDedx();
    m_chie_ap = dedxTrk_ap->chiE();
    m_chimu_ap = dedxTrk_ap->chiMu();
    m_chipi_ap = dedxTrk_ap->chiPi();
    m_chik_ap = dedxTrk_ap->chiK();
    m_chip_ap = dedxTrk_ap->chiP();
    m_ghit_ap = dedxTrk_ap->numGoodHits();
    m_thit_ap = dedxTrk_ap->numTotalHits();
    m_probPH_ap = dedxTrk_ap->probPH();
    m_normPH_ap = dedxTrk_ap->normPH();

      log << MSG::DEBUG << "end sotreDedxInfo" << endmsg;
}


int GammaPP::selectChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				   std::vector<int> & iPGood,
				   std::vector<int> & iMGood) {

  CLHEP::Hep3Vector xorigin = getOrigin(); 
  std::vector<int> iGood;
  iGood.clear();
  iPGood.clear();
  iMGood.clear();
  
  // loop through charged tracks 
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    // get mdcTrk 
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    // Good Kalman Track 
    if(!(*itTrk)->isMdcKalTrackValid()) continue;
    if(!(*itTrk)->isMdcTrackValid()) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();
    // Good Vertex 
    if (!passVertexSelection(xorigin, mdcTrk) ) continue; 
    // Polar angle cut
    if(fabs(cos(mdcTrk->theta())) > m_cha_costheta_cut) continue;
    iGood.push_back((*itTrk)->trackId());
    if(mdcTrk->charge()>0) iPGood.push_back(i);
    if(mdcTrk->charge()<0) iMGood.push_back(i);
  } // end charged tracks

  return iGood.size(); 
}

bool GammaPP::PPid(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
										int iPCharge,
										int iMCharge){
  EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPCharge;
  RecMdcTrack* mdcTrk_p = (*itTrk_p)->mdcTrack();
  double prob_pip, prob_kp, prob_pp, prob_ep;
  calcTrackPID(itTrk_p, prob_pip, prob_kp, prob_pp, prob_ep);  
  if(! (prob_pp > prob_kp &&
		prob_pp > m_prob_pion_min &&
		prob_pp > prob_pip &&
		prob_pp > prob_ep) ) return false;

  EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMCharge;
  RecMdcTrack* mdcTrk_m = (*itTrk_m)->mdcTrack();
	double prob_pim, prob_km, prob_pm, prob_em; 
  calcTrackPID(itTrk_m, prob_pim, prob_km, prob_pm, prob_em);
  if(! (prob_pm > prob_km &&
		prob_pm > m_prob_pion_min &&
		prob_pm > prob_pim &&
		prob_pm > prob_em) ) return false;

  m_prob_e_p = prob_ep;
  m_prob_pi_p = prob_pip;
  m_prob_k_p = prob_kp;
  m_prob_p_p = prob_pp;
  m_prob_e_m = prob_em;
  m_prob_pi_m = prob_pim;
  m_prob_k_m = prob_km;
  m_prob_p_m = prob_pm;

	return true;
}

bool GammaPP::PionPid(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
										int iPCharge,
										int iMCharge){
  EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPCharge;
  RecMdcTrack* mdcTrk_p = (*itTrk_p)->mdcTrack();
  double prob_pip, prob_kp, prob_pp, prob_ep;
  calcTrackPID(itTrk_p, prob_pip, prob_kp, prob_pp, prob_ep);  
  if(! (prob_pip > prob_kp &&
		prob_pip > m_prob_pion_min) ) return false;

  EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMCharge;
  RecMdcTrack* mdcTrk_m = (*itTrk_m)->mdcTrack();
	double prob_pim, prob_km, prob_pm, prob_em; 
  calcTrackPID(itTrk_m, prob_pim, prob_km, prob_pm, prob_em);
  if(! (prob_pim > prob_km &&
		prob_pim > m_prob_pion_min) ) return false;

	return true;
}

void GammaPP::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> & iGam) {
	iGam.clear();
  for(int i=evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i ;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    
    // TDC window
    if ( !(emcTrk->time() >= m_min_emctime && emcTrk->time() <= m_max_emctime) ) continue; 

    double abs_costheta(fabs(cos(emcTrk->theta())));
    bool barrel = (abs_costheta < m_costheta_barrel_max); 
    bool endcap = (abs_costheta > m_costheta_endcap_min
		   && abs_costheta < m_costheta_endcap_max);
    double eraw = emcTrk->energy();
    
    if ( !((barrel && eraw > m_energy_barrel_min)) )  continue; 

    CLHEP::Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    // EMC costheta cut 
    double costhe = cos(emcpos.theta());
    if ( fabs(costhe) >= m_gammaCosCut) continue;
    
    // find the nearest charged track
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      CLHEP::Hep3Vector extpos = extTrk->emcPosition();
      double angd = extpos.angle(emcpos);
      if(angd < dang) dang = angd;	    
    }

    //if(dang>=200) continue;
    dang = dang * 180 / (CLHEP::pi);
		if(dang < 20) continue;

    iGam.push_back(i);
  } 
  //return iGam.size(); 
}
/*
void GammaPP::storeShowersInfo(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
										int ig){
	EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ig;
	RecEmcShower* emcTrk = (*itTrk)->emcShower();
	m_sh_energy = emcTrk->energy();
	m_sh_phi = emcTrk->phi();
	m_sh_theta = emcTrk->theta();
	m_sh_e3 = emcTrk->e3x3();
	m_sh_e5 = emcTrk->e5x5();
	m_sh_status	= emcTrk->status();
	m_sh_numHits = emcTrk->numHits();
	m_sh_secondMoment = emcTrk->secondMoment();
	m_sh_dE = emcTrk->dE();
	m_sh_latMoment = emcTrk->latMoment();
	m_sh_a42Moment = emcTrk->a42Moment();
	m_sh_a20Moment = emcTrk->a20Moment();
	m_sh_seed =emcTrk->eSeed();
	m_sh_px = m_sh_energy*sin(m_sh_theta)*cos(m_sh_phi);
	m_sh_py = m_sh_energy*sin(m_sh_theta)*sin(m_sh_phi);
	m_sh_pz = m_sh_energy*cos(m_sh_theta);
	return;
}
*/
void GammaPP::store4CTracksInfo(HepLorentzVector pP4,
										HepLorentzVector mP4,
										HepLorentzVector gP4){
      MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin store4CTracksInfo" << endmsg;
	m_4c_p_p = pP4.rho();
	m_4c_p_costheta = cos(pP4.theta());
	m_4c_p_mass = pP4.m();
  m_4c_p_pxy = sqrt(pP4.px()*pP4.px()+pP4.py()*pP4.py());

	m_4c_ap_p = mP4.rho();
	m_4c_ap_costheta = cos(mP4.theta());
	m_4c_ap_mass = mP4.m();
  m_4c_ap_pxy = sqrt(mP4.px()*mP4.px()+mP4.py()*mP4.py());

	m_4c_pap_p = (pP4+mP4).rho();
	m_4c_pap_costheta = cos((pP4+mP4).theta());
	m_4c_pap_mass = (pP4+mP4).m();
  m_4c_pap_pxy = sqrt((pP4+mP4).px()*(pP4+mP4).px()+(pP4+mP4).py()*(pP4+mP4).py());

	m_4c_gamma_p = gP4.rho();
	m_4c_gamma_costheta = cos(gP4.theta());
  m_4c_gamma_pxy = sqrt(gP4.px()*gP4.px()+gP4.py()*gP4.py());

	m_4c_pgamma_p = (pP4+gP4).rho();
	m_4c_pgamma_costheta = cos((pP4+gP4).theta());
	m_4c_pgamma_mass = (pP4+gP4).m();
  m_4c_pgamma_pxy = sqrt((pP4+gP4).px()*(pP4+gP4).px()+(pP4+gP4).py()*(pP4+gP4).py());

	m_4c_apgamma_p = (gP4+mP4).rho();
	m_4c_apgamma_costheta = cos((gP4+mP4).theta());
	m_4c_apgamma_mass = (gP4+mP4).m();
  m_4c_apgamma_pxy = sqrt((gP4+mP4).px()*(gP4+mP4).px()+(gP4+mP4).py()*(gP4+mP4).py());
          log << MSG::DEBUG << "end store4CTracksInfo" << endmsg;
}



void GammaPP::storeKalTracksInfo(HepLorentzVector pP4,
										HepLorentzVector mP4,
										HepLorentzVector gP4){
      MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin storeKalTrackInfo" << endmsg;
	m_kal_p_p = pP4.rho();
	m_kal_p_costheta = cos(pP4.theta());
	m_kal_p_e = pP4.e();
	m_kal_p_mass = pP4.m();
  m_kal_p_pxy = sqrt(pP4.px()*pP4.px()+pP4.py()*pP4.py());

	m_kal_ap_p = mP4.rho();
	m_kal_ap_costheta = cos(mP4.theta());
	m_kal_ap_e = mP4.e();
	m_kal_ap_mass = mP4.m();
  m_kal_ap_pxy = sqrt(mP4.px()*mP4.px()+mP4.py()*mP4.py());

	m_kal_pap_p = (pP4+mP4).rho();
	m_kal_pap_costheta = cos((pP4+mP4).theta());
	m_kal_pap_e = (pP4+mP4).e();
	m_kal_pap_mass = (pP4+mP4).m();
  m_kal_pap_pxy = sqrt((pP4+mP4).px()*(pP4+mP4).px()+(pP4+mP4).py()*(pP4+mP4).py());

	m_kal_gamma_p = gP4.rho();
	m_kal_gamma_costheta = cos(gP4.theta());
  m_kal_gamma_pxy = sqrt(gP4.px()*gP4.px()+gP4.py()*gP4.py());

	m_kal_pgamma_p = (pP4+gP4).rho();
	m_kal_pgamma_costheta = cos((pP4+gP4).theta());
	m_kal_pgamma_e = (pP4+gP4).e();
	m_kal_pgamma_mass = (pP4+gP4).m();
  m_kal_pgamma_pxy = sqrt((pP4+gP4).px()*(pP4+gP4).px()+(pP4+gP4).py()*(pP4+gP4).py());

	m_kal_apgamma_p = (gP4+mP4).rho();
	m_kal_apgamma_costheta = cos((gP4+mP4).theta());
	m_kal_apgamma_e = (gP4+mP4).e();
	m_kal_apgamma_mass = (gP4+mP4).m();
  m_kal_apgamma_pxy = sqrt((gP4+mP4).px()*(gP4+mP4).px()+(gP4+mP4).py()*(gP4+mP4).py());

	m_kal_jpsi_p = (pP4+mP4+gP4).rho();
	m_kal_jpsi_costheta = cos((pP4+mP4+gP4).theta());
	m_kal_jpsi_e = (pP4+mP4+gP4).e();
	m_kal_jpsi_mass = (pP4+mP4+gP4).m();
  m_kal_jpsi_pxy = sqrt((pP4+mP4+gP4).px()*(pP4+mP4+gP4).px()+(pP4+mP4+gP4).py()*(pP4+mP4+gP4).py());
          log << MSG::DEBUG << "end storeKalTrackInfo" << endmsg;
}

void GammaPP::storeBakTracksInfo(HepLorentzVector pP4,
										HepLorentzVector mP4,
										HepLorentzVector gP4){
      MsgStream log(msgSvc(), name());
          log << MSG::DEBUG << "begin storeBakTrackInfo" << endmsg;
	m_bak_p_p = pP4.rho();
	m_bak_p_costheta = cos(pP4.theta());
	m_bak_p_e = pP4.e();
	m_bak_p_mass = pP4.m();
  m_bak_p_pxy = sqrt(pP4.px()*pP4.px()+pP4.py()*pP4.py());

	m_bak_ap_p = mP4.rho();
	m_bak_ap_costheta = cos(mP4.theta());
	m_bak_ap_e = mP4.e();
	m_bak_ap_mass = mP4.m();
  m_bak_ap_pxy = sqrt(mP4.px()*mP4.px()+mP4.py()*mP4.py());

	m_bak_pap_p = (pP4+mP4).rho();
	m_bak_pap_costheta = cos((pP4+mP4).theta());
	m_bak_pap_e = (pP4+mP4).e();
	m_bak_pap_mass = (pP4+mP4).m();
  m_bak_pap_pxy = sqrt((pP4+mP4).px()*(pP4+mP4).px()+(pP4+mP4).py()*(pP4+mP4).py());

	m_bak_gamma_p = gP4.rho();
	m_bak_gamma_costheta = cos(gP4.theta());
  m_bak_gamma_pxy = sqrt(gP4.px()*gP4.px()+gP4.py()*gP4.py());
  m_bak_gamma_u = gP4.e() - gP4.rho();
  m_bak_gamma_e = gP4.e();

	m_bak_pgamma_p = (pP4+gP4).rho();
	m_bak_pgamma_costheta = cos((pP4+gP4).theta());
	m_bak_pgamma_e = (pP4+gP4).e();
	m_bak_pgamma_mass = (pP4+gP4).m();
  m_bak_pgamma_pxy = sqrt((pP4+gP4).px()*(pP4+gP4).px()+(pP4+gP4).py()*(pP4+gP4).py());

	m_bak_apgamma_p = (gP4+mP4).rho();
	m_bak_apgamma_costheta = cos((gP4+mP4).theta());
	m_bak_apgamma_e = (gP4+mP4).e();
	m_bak_apgamma_mass = (gP4+mP4).m();
  m_bak_apgamma_pxy = sqrt((gP4+mP4).px()*(gP4+mP4).px()+(gP4+mP4).py()*(gP4+mP4).py());
          log << MSG::DEBUG << "begin storeBakTrackInfo" << endmsg;
}

void GammaPP::calcTrackPID(EvtRecTrackIterator itTrk_p,
			     double& prob_pip,
			     double& prob_kp,
					 double& prob_pp,
					 double& prob_ep) {
  prob_pip = 999.; 
  prob_kp = 999.; 
  ParticleID * pidp = ParticleID::instance();
  pidp->init();
  pidp->setMethod(pidp->methodProbability());
  pidp->setChiMinCut(4);
  pidp->setRecTrack(*itTrk_p);
  // use PID sub-system
	pidp->usePidSys(pidp->useDedx() | pidp->useTof1() | pidp->useTof2() | pidp->useTofE());
	pidp->identify(pidp->onlyProton()|pidp->onlyPion()|pidp->onlyKaon()|pidp->onlyElectron());
  //pidp->usePidSys(pidp->useDedx() | pidp->useTof1() | pidp->useTof2());
  //pidp->identify(pidp->onlyPionKaonProton());
  pidp->calculate();
  if(pidp->IsPidInfoValid()) {
    prob_pip = pidp->probPion();
    prob_kp  = pidp->probKaon();
		prob_pp  = pidp->probProton();
		prob_ep  = pidp->probElectron();
  }
	return;
}

bool GammaPP::KinematicFit(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> iGam,
					 int iPCharge,
					 int iMCharge,
					 HepLorentzVector & pP4,
					 HepLorentzVector & mP4,
					 HepLorentzVector & gP4,
					 int & ig) {
	float chisq=9999;
	RecMdcKalTrack *pTrk = (*(evtRecTrkCol->begin()+iPCharge))->mdcKalTrack();
	WTrackParameter wvpTrk = WTrackParameter(P_MASS, pTrk->getZHelixP(), pTrk->getZErrorP());
	RecMdcKalTrack *mTrk = (*(evtRecTrkCol->begin()+iMCharge))->mdcKalTrack();
	WTrackParameter wvmTrk = WTrackParameter(P_MASS, mTrk->getZHelixP(), mTrk->getZErrorP());
	pP4 = wvpTrk.p();
	mP4 = wvmTrk.p();
	for(int i = 0; i<iGam.size(); i++){
		RecEmcShower *gTrk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();

		KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setChisqCut(20);
		kmfit->AddTrack(0, wvpTrk);
		kmfit->AddTrack(1, wvmTrk);
		kmfit->AddTrack(2, 0, gTrk);
		kmfit->AddFourMomentum(0, jpsi_cms);
		if(kmfit->Fit() == false) continue;
		if(kmfit->chisq() < chisq){
			chisq = kmfit->chisq();
			pP4 = kmfit->pfit(0);
			mP4 = kmfit->pfit(1);
			gP4 = kmfit->pfit(2);
			ig = iGam[i];
		}
	}
	if(chisq > 20) return false;
  m_4c_chisq = chisq;
	return true;
}

bool GammaPP::ExtraPi0(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> iGam) {
	int nGam = iGam.size();
	float chisq=9999;
	for(int m = 0; m < nGam-1; m++) {
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[m]))->emcShower();
		double eraw1 = g1Trk->energy();
		double phi1 = g1Trk->phi();
		double the1 = g1Trk->theta();
		HepLorentzVector ptrkg1;
		ptrkg1.setPx(eraw1*sin(the1)*cos(phi1));
		ptrkg1.setPy(eraw1*sin(the1)*sin(phi1));
		ptrkg1.setPz(eraw1*cos(the1));
		ptrkg1.setE(eraw1);

		for(int n = m+1; n < nGam; n++) {
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[n]))->emcShower();
			double eraw2 = g2Trk->energy();
			double phi2 = g2Trk->phi();
			double the2 = g2Trk->theta();
			HepLorentzVector ptrkg2;
			ptrkg2.setPx(eraw2*sin(the2)*cos(phi2));
			ptrkg2.setPy(eraw2*sin(the2)*sin(phi2));
			ptrkg2.setPz(eraw2*cos(the2));
			ptrkg2.setE(eraw2);

			HepLorentzVector  ptrkpi0;
			ptrkpi0 = ptrkg1+ptrkg2;
			double m_xmpi0_tem = ptrkpi0.m();
			if(m_xmpi0_tem>0.15||m_xmpi0_tem<0.115)  continue;

			return true;
		}   
	}
	return false;
}

bool GammaPP::ExtraEta(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> iGam) {
	int nGam = iGam.size();
	float chisq=9999;
	for(int m = 0; m < nGam-1; m++) {
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[m]))->emcShower();
		double eraw1 = g1Trk->energy();
		double phi1 = g1Trk->phi();
		double the1 = g1Trk->theta();
		HepLorentzVector ptrkg1;
		ptrkg1.setPx(eraw1*sin(the1)*cos(phi1));
		ptrkg1.setPy(eraw1*sin(the1)*sin(phi1));
		ptrkg1.setPz(eraw1*cos(the1));
		ptrkg1.setE(eraw1);

		for(int n = m+1; n < nGam; n++) {
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[n]))->emcShower();
			double eraw2 = g2Trk->energy();
			double phi2 = g2Trk->phi();
			double the2 = g2Trk->theta();
			HepLorentzVector ptrkg2;
			ptrkg2.setPx(eraw2*sin(the2)*cos(phi2));
			ptrkg2.setPy(eraw2*sin(the2)*sin(phi2));
			ptrkg2.setPz(eraw2*cos(the2));
			ptrkg2.setE(eraw2);

			HepLorentzVector  ptrkpi0;
			ptrkpi0 = ptrkg1+ptrkg2;
			double m_xmpi0_tem = ptrkpi0.m();
			if(m_xmpi0_tem>0.57||m_xmpi0_tem<0.5)  continue;

		}   
	}
	return false;
}


