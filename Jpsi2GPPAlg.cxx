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
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McDecayModeSvc/IMcDecayModeSvc.h"


#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <vector>

class Jpsi2GPP: public Algorithm {
  
public:
  Jpsi2GPP(const std::string&, ISvcLocator*);
  ~Jpsi2GPP(); 
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
  // declare some cuts 

  //  MC truth info
	const int kMaxId; 
	NTuple::Item<int>  m_numParticle;
	NTuple::Array<int> m_pdgid;
	NTuple::Array<int> m_motherindex;

	//NTuple::Item<int> m_gamma;
	//NTuple::Array<float> m_gamma_pxmc;
	//NTuple::Array<float> m_gamma_pymc;
	//NTuple::Array<float> m_gamma_pzmc;
	//NTuple::Array<float> m_gamma_pmc;
	//NTuple::Array<float> m_gamma_emc;
	//NTuple::Array<float> m_gamma_m0mc;
	//NTuple::Array<float> m_gamma_xmc;
	//NTuple::Array<float> m_gamma_ymc;
	//NTuple::Array<float> m_gamma_zmc;
	//NTuple::Array<float> m_gamma_x0mc;
	//NTuple::Array<float> m_gamma_y0mc;
	//NTuple::Array<float> m_gamma_z0mc;

	//NTuple::Item<int> m_gamma_fsr;
	//NTuple::Array<float> m_gamma_fsr_pxmc;
	//NTuple::Array<float> m_gamma_fsr_pymc;
	//NTuple::Array<float> m_gamma_fsr_pzmc;
	//NTuple::Array<float> m_gamma_fsr_pmc;
	//NTuple::Array<float> m_gamma_fsr_emc;
	//NTuple::Array<float> m_gamma_fsr_m0mc;
	//NTuple::Array<float> m_gamma_fsr_xmc;
	//NTuple::Array<float> m_gamma_fsr_ymc;
	//NTuple::Array<float> m_gamma_fsr_zmc;
	//NTuple::Array<float> m_gamma_fsr_x0mc;
	//NTuple::Array<float> m_gamma_fsr_y0mc;
	//NTuple::Array<float> m_gamma_fsr_z0mc;

	//const ecms_4p
	HepLorentzVector ECMS_4P;

  //bool m_isMonteCarlo; 

	McDecayModeSvc* m_svc;
  // define Histograms
  TH1F* h_evtflw; 
	int N1,N2,N3,N4,N5,N6;
  
  // define Trees
  TTree* m_tree;

  // common info 
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;

  
  // neutral tracks
  NTuple::Item<int> m_sh_good;
	//NTuple::Item<float> m_sh_energy;
	//NTuple::Item<float> m_sh_phi;
	//NTuple::Item<float> m_sh_theta;
	//NTuple::Item<float> m_sh_pr;
	//NTuple::Item<float> m_sh_px;
	//NTuple::Item<float> m_sh_py;
	//NTuple::Item<float> m_sh_pz;
	//NTuple::Item<float> m_sh_e3;
	//NTuple::Item<float> m_sh_e5;
	//NTuple::Item<int> m_sh_status;
	//NTuple::Item<int> m_sh_numHits;
	//NTuple::Item<float> m_sh_secondMoment;
	//NTuple::Item<float> m_sh_latMoment;
	//NTuple::Item<float> m_sh_a42Moment;
	//NTuple::Item<float> m_sh_a20Moment;
	//NTuple::Item<float> m_sh_dE;
	//NTuple::Item<float> m_sh_seed;
	//NTuple::Item<float> m_sh_dang;

	NTuple::Item<int> m_ext_pi0;
	NTuple::Item<int> m_ext_eta;
	// charged tracks
	NTuple::Item<float> m_p_px;
	NTuple::Item<float> m_p_py;
	NTuple::Item<float> m_p_pz;
	NTuple::Item<float> m_p_p;
	NTuple::Item<float> m_p_the;
	NTuple::Item<float> m_p_e;

	NTuple::Item<float> m_m_px;
	NTuple::Item<float> m_m_py;
	NTuple::Item<float> m_m_pz;
	NTuple::Item<float> m_m_p;
	NTuple::Item<float> m_m_the;
	NTuple::Item<float> m_m_e;

	NTuple::Item<float> m_the;

	NTuple::Item<float> m_g_px;
	NTuple::Item<float> m_g_py;
	NTuple::Item<float> m_g_pz;
	NTuple::Item<float> m_g_e;
	NTuple::Item<float> m_g_r;
  // functions		
  void book_tree(); 
  void before_execute(); 
  bool buildJpsi2GPP();
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
}; 

//
// module declare
//

DECLARE_ALGORITHM_FACTORY( Jpsi2GPP )
DECLARE_FACTORY_ENTRIES( Jpsi2GPPAlg ) {
  DECLARE_ALGORITHM(Jpsi2GPP);
}

LOAD_FACTORY_ENTRIES( Jpsi2GPPAlg )

//
// constants
//

const double PION_MASS = 0.139570;
const double E_MASS = 0.000511;
const double P_MASS = 0.9382723;
const double PION0_MASS = 0.134977;
const double JPSI_MASS = 3.096916;


//
// member functions
//

Jpsi2GPP::Jpsi2GPP(const std::string& name, ISvcLocator* pSvcLocator) :
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

StatusCode Jpsi2GPP::initialize(){
	//ECMS_4P.set(3.097*sin(0.011),0.,0.,3.097);     
	ECMS_4P.set(3.097*sin(0.011),0.,0.,3.097);     
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

StatusCode Jpsi2GPP::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

	before_execute();
  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader) return StatusCode::FAILURE;

  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();

  
  if(buildJpsi2GPP() == true) m_tuple1->write();

  return StatusCode::SUCCESS; 
}

StatusCode Jpsi2GPP::finalize() {
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

Jpsi2GPP::~Jpsi2GPP() {
}

void Jpsi2GPP::before_execute(){
}

void Jpsi2GPP::book_tree() {
  StatusCode status;
  NTuplePtr nt1(ntupleSvc(), "FILE/Jpsi2GPP");
  if ( nt1 ) m_tuple1 = nt1;
  else {
   m_tuple1 = ntupleSvc()->book ("FILE1/Jpsi2GPP", CLID_ColumnWiseTuple, " N-Tuple example");
    if ( m_tuple1 )    {
			//MC info 
  		status=m_tuple1->addItem("indexmc",m_numParticle,0,100);
			status=m_tuple1->addIndexedItem("pdgid", m_numParticle, m_pdgid);
			status=m_tuple1->addIndexedItem("motheridx", m_numParticle, m_motherindex);

			//truth info
		  //status = m_tuple1->addItem("ngamma_mc",m_gamma,0,100);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_px",m_gamma, m_gamma_pxmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_py",m_gamma, m_gamma_pymc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_pz",m_gamma, m_gamma_pzmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_p",m_gamma, m_gamma_pmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_e",m_gamma, m_gamma_emc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_m0",m_gamma, m_gamma_m0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_x",m_gamma, m_gamma_xmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_y",m_gamma, m_gamma_ymc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_z",m_gamma, m_gamma_zmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_x0",m_gamma, m_gamma_x0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_y0",m_gamma, m_gamma_y0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_z0",m_gamma, m_gamma_z0mc);
		  //
		  //status = m_tuple1->addItem("ngamma_fsr_mc",m_gamma_fsr,0,100);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_px",m_gamma_fsr, m_gamma_fsr_pxmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_py",m_gamma_fsr, m_gamma_fsr_pymc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_pz",m_gamma_fsr, m_gamma_fsr_pzmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_p",m_gamma_fsr, m_gamma_fsr_pmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_e",m_gamma_fsr, m_gamma_fsr_emc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_m0",m_gamma_fsr, m_gamma_fsr_m0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_x",m_gamma_fsr, m_gamma_fsr_xmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_y",m_gamma_fsr, m_gamma_fsr_ymc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_z",m_gamma_fsr, m_gamma_fsr_zmc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_x0",m_gamma_fsr, m_gamma_fsr_x0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_y0",m_gamma_fsr, m_gamma_fsr_y0mc);
		  //status = m_tuple1->addIndexedItem ("mc_gamma_fsr_z0",m_gamma_fsr, m_gamma_fsr_z0mc);
		  
  		//commom info
  		status=m_tuple1->addItem("run",m_run);
  		status=m_tuple1->addItem("event",m_event);
	  
  		//netual tracks
			status=m_tuple1->addItem("ng", m_sh_good);
			status=m_tuple1->addItem("extpi0", m_ext_pi0);
			status=m_tuple1->addItem("exteta", m_ext_eta);
			//status=m_tuple1->addItem("game", m_sh_energy);
			//status=m_tuple1->addItem("gamphi", m_sh_phi);
			//status=m_tuple1->addItem("gamthe", m_sh_theta);
			//status=m_tuple1->addItem("gampx", m_sh_px);
			//status=m_tuple1->addItem("gampy", m_sh_py);
			//status=m_tuple1->addItem("gampz", m_sh_pz);
			//status=m_tuple1->addItem("gameIII", m_sh_e3);
			//status=m_tuple1->addItem("gameV", m_sh_e5);
			//status=m_tuple1->addItem("gamstatus", m_sh_status);
			//status=m_tuple1->addItem("gamhitnum", m_sh_numHits);
			//status=m_tuple1->addItem("gamiindmoment", m_sh_secondMoment);
			//status=m_tuple1->addItem("gamlatmoment", m_sh_latMoment);
			//status=m_tuple1->addItem("gama42", m_sh_a42Moment);
			//status=m_tuple1->addItem("gama20", m_sh_a20Moment);
			//status=m_tuple1->addItem("gamerrorE", m_sh_dE);
			//status=m_tuple1->addItem("gamsed", m_sh_seed);

			//charged tracks
			status=m_tuple1->addItem("pPx", m_p_px);
			status=m_tuple1->addItem("pPy", m_p_py);
			status=m_tuple1->addItem("pPz", m_p_pz);
			status=m_tuple1->addItem("pPp", m_p_p);
			status=m_tuple1->addItem("pPe", m_p_e);
			status=m_tuple1->addItem("pPthe", m_p_the);

			status=m_tuple1->addItem("pMx", m_m_px);
			status=m_tuple1->addItem("pMy", m_m_py);
			status=m_tuple1->addItem("pMz", m_m_pz);
			status=m_tuple1->addItem("pMp", m_m_p);
			status=m_tuple1->addItem("pMe", m_m_e);
			status=m_tuple1->addItem("pMthe", m_m_the);

			status=m_tuple1->addItem("pmthe", m_the);
			status=m_tuple1->addItem("mp_mass", m_mp_mass);

			status=m_tuple1->addItem("gPx", m_g_px);
			status=m_tuple1->addItem("gPy", m_g_py);
			status=m_tuple1->addItem("gPz", m_g_pz);
			status=m_tuple1->addItem("gE", m_g_e);
			status=m_tuple1->addItem("gPr", m_g_r);

		}
    else    {
      return ;
    }
	}
	h_evtflw= new TH1F("hevtflw", "hevtflw", 6, 0, 6);
  return ;
}


bool Jpsi2GPP::buildJpsi2GPP() {
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
 	if (m_run < 0) {
		if(!saveGenInfo(mcParticleCol)) return false;
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
	//storeShowersInfo(evtRecTrkCol, ig);
	storeTracksInfo(pP4, mP4, gP4);
	m_sh_good = iGam.size();

	return true;
}

bool Jpsi2GPP::saveGenInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol) {
  MsgStream log(msgSvc(), name());
	//m_gamma=0;
	//m_gamma_fsr=0;
	if (!mcParticleCol){
		std::cout << "Could not retrieve McParticelCol" << std::endl;
		return false;
	}
  IMcDecayModeSvc* i_svc;
  StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
  if ( sc_DecayModeSvc.isFailure() ){
      log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
      return false;
  }
  m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

	std::vector<int> pdgid, motherindex;
	pdgid.clear();
	motherindex.clear();
  Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
  for (; iter_mc != mcParticleCol->end(); iter_mc++){
    if ((*iter_mc)->primaryParticle()) continue;
  //  if (!(*iter_mc)->decayFromGenerator()) continue;

		if ((*iter_mc)->particleProperty()==443){
			int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
			m_numParticle = pdgid.size();
			for(int i = 0; i < m_numParticle; i++){
				m_pdgid[i]=(pdgid[i]);
				m_motherindex[i]=(motherindex[i]);
			}
		}
		/*
		HepLorentzVector  mctrue_track = (*iter_mc)->initialFourMomentum();
		HepLorentzVector  mctrue_initPost = (*iter_mc)->initialPosition();
		HepLorentzVector  mctrue_end = (*iter_mc)->finalPosition();
		if((*iter_mc)->particleProperty()==22 ){
			m_gamma_pxmc[m_gamma]=mctrue_track.px();
			m_gamma_pymc[m_gamma]=mctrue_track.py();
			m_gamma_pzmc[m_gamma]=mctrue_track.pz();
			m_gamma_pmc[m_gamma]=mctrue_track.rho();
			m_gamma_emc[m_gamma]=mctrue_track.e();
			m_gamma_m0mc[m_gamma]=mctrue_track.m();
			m_gamma_xmc[m_gamma]=mctrue_end.px();
			m_gamma_ymc[m_gamma]=mctrue_end.py();
			m_gamma_zmc[m_gamma]=mctrue_end.pz();
			m_gamma_x0mc[m_gamma]=mctrue_initPost.px();
			m_gamma_y0mc[m_gamma]=mctrue_initPost.py();
			m_gamma_z0mc[m_gamma]=mctrue_initPost.pz();
			m_gamma++;
			}
		if((*iter_mc)->particleProperty()==-22 ){
			m_gamma_fsr_pxmc[m_gamma_fsr]=mctrue_track.px();
			m_gamma_fsr_pymc[m_gamma_fsr]=mctrue_track.py();
			m_gamma_fsr_pzmc[m_gamma_fsr]=mctrue_track.pz();
			m_gamma_fsr_pmc[m_gamma_fsr]=mctrue_track.rho();
			m_gamma_fsr_emc[m_gamma_fsr]=mctrue_track.e();
			m_gamma_fsr_m0mc[m_gamma_fsr]=mctrue_track.m();
			m_gamma_fsr_xmc[m_gamma_fsr]=mctrue_end.px();
			m_gamma_fsr_ymc[m_gamma_fsr]=mctrue_end.py();
			m_gamma_fsr_zmc[m_gamma_fsr]=mctrue_end.pz();
			m_gamma_fsr_x0mc[m_gamma_fsr]=mctrue_initPost.px();
			m_gamma_fsr_y0mc[m_gamma_fsr]=mctrue_initPost.py();
			m_gamma_fsr_z0mc[m_gamma_fsr]=mctrue_initPost.pz();
			m_gamma_fsr++;
			}
			*/
	}
	return true;
}

CLHEP::Hep3Vector Jpsi2GPP::getOrigin() {
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

bool Jpsi2GPP::passVertexSelection(CLHEP::Hep3Vector xorigin,
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

int Jpsi2GPP::selectChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
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

bool Jpsi2GPP::PPid(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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

	return true;
}

bool Jpsi2GPP::PionPid(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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

void Jpsi2GPP::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
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
void Jpsi2GPP::storeShowersInfo(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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
void Jpsi2GPP::storeTracksInfo(HepLorentzVector pP4,
										HepLorentzVector mP4,
										HepLorentzVector gP4){
	m_p_px = pP4.px();
	m_p_py = pP4.py();
	m_p_pz = pP4.pz();
	m_p_p = pP4.rho();
	m_p_the = pP4.pz()/pP4.rho();
	m_p_e = pP4.e();

	m_m_px = mP4.px();
	m_m_py = mP4.py();
	m_m_pz = mP4.pz();
	m_m_p = mP4.rho();
	m_m_the = mP4.pz()/mP4.rho();
	m_m_e = mP4.e();

	m_g_px = gP4.px();
	m_g_py = gP4.py();
	m_g_pz = gP4.pz();
	m_g_e = gP4.rho();
	m_g_r = gP4.perp2();

	m_the = pP4.angle(mP4);
  m_mp_mass = (pP4 + mP4).m();

}

void Jpsi2GPP::calcTrackPID(EvtRecTrackIterator itTrk_p,
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

bool Jpsi2GPP::KinematicFit(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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
		kmfit->AddFourMomentum(0, ECMS_4P);
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
	return true;
}

bool Jpsi2GPP::ExtraPi0(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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

bool Jpsi2GPP::ExtraEta(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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


