#define APP_NAME "MMTool"
#define MyInfo(a,b)    if(m_debugMode<=MSG::INFO   ) std::cout<<"Info    in <MMTool::"<<(a)<<">: "<<(b)<<std::endl;
#define MyWarning(a,b) if(m_debugMode<=MSG::WARNING) std::cout<<"Warning in <MMTool::"<<(a)<<">: "<<(b)<<std::endl;
#define MyError(a,b)   if(m_debugMode<=MSG::ERROR  ) std::cout<<"Error   in <MMTool::"<<(a)<<">: "<<(b)<<std::endl;
#define MyDebug(a,b)   if(m_debugMode<=MSG::DEBUG  ) std::cout<<"Debug   in <MMTool::"<<(a)<<">: "<<(b)<<std::endl;
#define MyAlways(a,b)  if(m_debugMode<=MSG::ALWAYS ) std::cout<<"In         <MMTool::"<<(a)<<">: "<<(b)<<std::endl;

#include<iostream>
#include<MMTool/MMTool.h>

#include"TH2F.h"
#include"TString.h"

namespace MM {

  ////////////////////////////////////////////////////////////////
  MMTool::MMTool(NLepMM mmMode, const MSG::Level& debugMode)
  {
    m_debugMode = debugMode;
    initialize(mmMode);
  }
  ////////////////////////////////////////////////////////////////
  MMTool::~MMTool()
  {
    finalize();
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::initialize(NLepMM mmMode)
  {
    MyDebug("initialize()","Initializing MMTool...");
    initLeps();
    m_mmMode = mmMode;
    m_nDetIsZero=0;

    Bool_t isSupportedMode=kFALSE;
    for(Int_t mode=0; mode<nNLepMM; mode++){
      if(m_mmMode==mode) isSupportedMode=kTRUE;
    }
    if(not isSupportedMode){
      MyError("initialize()","MMTool only supports 2LMM or 3LMM.");
      return;
    }

    if     (m_mmMode==TwoLepMM  ){
      MyDebug("initialize()","MMTool was initialized for two leptons MM (4x4).");
      m_nLepsForMM = 2;
    }
    else if(m_mmMode==ThreeLepMM){
      MyDebug("initialize()","MMTool was initialized for two leptons MM (8x8).");
      m_nLepsForMM = 3;
    }

    return;
  }
  ////////////////////////////////////////////////////////////////
  Bool_t MMTool::prepEffHists(TString fileName)
  {
    //Dummy histograms for testing...
    Int_t nbinsx(0), nbinsy(0);
    Double_t ptmin=0.;
    Double_t ptmax=200.;
    Double_t etamin=0.;
    Double_t etamax=2.5;
    m_realElEff = new TH2F("m_realElEff","",nbinsx,ptmin,ptmax,nbinsy,etamin,etamax);//Efficiency in the (x,y)=(Pt,|Eta|) plane
    m_realMuEff = new TH2F("m_realMuEff","",nbinsx,ptmin,ptmax,nbinsy,etamin,etamax);//Efficiency in the (x,y)=(Pt,|Eta|) plane
    m_fakeElEff = new TH2F("m_fakeElEff","",nbinsx,ptmin,ptmax,nbinsy,etamin,etamax);//Efficiency in the (x,y)=(Pt,|Eta|) plane
    m_fakeMuEff = new TH2F("m_fakeMuEff","",nbinsx,ptmin,ptmax,nbinsy,etamin,etamax);//Efficiency in the (x,y)=(Pt,|Eta|) plane
    for(Int_t binx=1; binx<nbinsx+2; binx++){ //Including overflow bins
      for(Int_t biny=1; biny<nbinsy+2; biny++){ //Including overflow bins
        m_realElEff->SetBinContent(binx,biny,0.686926);
        m_realMuEff->SetBinContent(binx,biny,0.703276);
        m_fakeElEff->SetBinContent(binx,biny,0.117019);
        m_fakeMuEff->SetBinContent(binx,biny,0.129006);
      }
    }
    return true;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::finalize()
  {
    MyDebug("finalize()","Finalizing MMTool...");
    if(m_nDetIsZero!=0){
      MyWarning("finalize()",Form("There are %d evetns ignored for MM due to Determinant==0.",m_nDetIsZero));
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::setLepProp(const Int_t& xthLep, const LepProp& lep)
  {
    if(xthLep<m_nLepsForMM){
      m_leps[xthLep] = lep;
    }else{
      MyError("setLepProp()",Form("Unnecessary leptons are set for MMTool. This might be a potential bug."));
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::setLepProps(const LepProp& lep1, const LepProp& lep2)
  {
    if(m_nLepsForMM==2){
      m_leps[0] = lep1;
      m_leps[1] = lep2;
    }else{
      MyError("setLepProps()",Form("All necessary leptons are not set for MMTool. This might be a potential bug."));
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::setLepProps(const LepProp& lep1, const LepProp& lep2, const LepProp& lep3)
  {
    if(m_nLepsForMM==3){
      m_leps[0] = lep1;
      m_leps[1] = lep2;
      m_leps[2] = lep3;
    }else{
      MyError("setLepProps()",Form("All necessary leptons are not set for MMTool. This might be a potential bug."));
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::initLeps()
  {
    for(Int_t lep=0; lep<3; lep++){
      m_leps[lep].pt=0.;
      m_leps[lep].eta=0.;
      m_leps[lep].flav=LepFlav::undefined;
      m_leps[lep].isSignal=kFALSE;
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  Double_t MMTool::getRealEff(const Int_t& xthLep)
  {
    Double_t eff=0.;
    if(xthLep<m_nLepsForMM){
      if      (m_leps[xthLep].flav==LepFlav::el){
        Int_t tmpBin = m_realElEff->FindBin(m_leps[xthLep].pt/1000.,m_leps[xthLep].eta);
        eff = m_realElEff->GetBinContent(tmpBin);
      }else if(m_leps[xthLep].flav==LepFlav::mu){
        Int_t tmpBin = m_realMuEff->FindBin(m_leps[xthLep].pt/1000.,m_leps[xthLep].eta);
        eff = m_realMuEff->GetBinContent(tmpBin);
      }
    }else{
      MyError("getRealEff()",Form("Unnecessary lepton (id=%d) was accessed. This might be a potential bug.",xthLep));
    }
    return eff;
  }
  ////////////////////////////////////////////////////////////////
  Double_t MMTool::getFakeEff(const Int_t& xthLep)
  {
    Double_t eff=0.;
    if(xthLep<m_nLepsForMM){
      if      (m_leps[xthLep].flav==LepFlav::el){
        Int_t tmpBin = m_fakeElEff->FindBin(m_leps[xthLep].pt/1000.,m_leps[xthLep].eta);
        eff = m_fakeElEff->GetBinContent(tmpBin);
      }else if(m_leps[xthLep].flav==LepFlav::mu){
        Int_t tmpBin = m_fakeMuEff->FindBin(m_leps[xthLep].pt/1000.,m_leps[xthLep].eta);
        eff = m_fakeMuEff->GetBinContent(tmpBin);
      }
    }else{
      MyError("getFakeEff()",Form("Unnecessary lepton (id=%d) was accessed. This might be a potential bug.",xthLep));
    }
    return eff;
  }
  ////////////////////////////////////////////////////////////////
  void MMTool::setEff(Double_t (&r)[3], Double_t (&f)[3], Double_t (&R)[3], Double_t (&F)[3])
  {
    for(Int_t lep=0; lep<m_nLepsForMM; lep++){
      r[lep] = getRealEff(lep);
      f[lep] = getFakeEff(lep);
      R[lep] = 1. - r[lep];
      F[lep] = 1. - f[lep];
    }
    return;
  }
  ////////////////////////////////////////////////////////////////
  Double_t MMTool::getMMWeight()
  {
    Double_t mmWeight=0.;
    MyDebug("getMMWeight()","Calculating MM Weight");

    Double_t r[3]={0.,0.,0.};
    Double_t f[3]={0.,0.,0.};
    Double_t R[3]={0.,0.,0.};
    Double_t F[3]={0.,0.,0.};
    setEff(r,f,R,F);

    if(m_nLepsForMM==2){

      Int_t Nxx[4]={0,0,0,0}; //Nxx[0]=Ntt, Nxx[1]=Ntl, Nxx[2]=Nlt, Nxx[3]=Nll
      if(  m_leps[0].isSignal &&   m_leps[1].isSignal) Nxx[0]=1;
      if(  m_leps[0].isSignal &&  !m_leps[1].isSignal) Nxx[1]=1;
      if( !m_leps[0].isSignal &&   m_leps[1].isSignal) Nxx[2]=1;
      if( !m_leps[0].isSignal &&  !m_leps[1].isSignal) Nxx[3]=1;

      Double_t Nrr =   F[0]*F[1]*Nxx[0] - F[0]*f[1]*Nxx[1] - f[0]*F[1]*Nxx[2] + f[0]*f[1]*Nxx[3];
      Double_t Nrf = - F[0]*R[1]*Nxx[0] + F[0]*r[1]*Nxx[1] + f[0]*R[1]*Nxx[2] - f[0]*r[1]*Nxx[3];
      Double_t Nfr = - R[0]*F[1]*Nxx[0] + R[0]*f[1]*Nxx[1] + r[0]*F[1]*Nxx[2] - r[0]*f[1]*Nxx[3];
      Double_t Nff =   R[0]*R[1]*Nxx[0] - R[0]*r[1]*Nxx[1] - r[0]*R[1]*Nxx[2] + r[0]*r[1]*Nxx[3];

      MyDebug("getMMWeight()",Form("==== Num of evts: Ntt=%d, Ntl=%d, Nlt=%d, Nll=%d",Nxx[0],Nxx[1],Nxx[2],Nxx[3]));
      MyDebug("getMMWeight()",Form("==== Fake Rates : r[0]=%f, f[0]=%f, r[1]=%f, f[1]=%f",r[0],f[0],r[1],f[1]));
      MyDebug("getMMWeight()",Form("==== Fake vector: Nrr=%f, Nrf=%f, Nfr=%f, Nff=%f",Nrr,Nrf,Nfr,Nff));

      if(r[0]-f[0]==0 || r[1]-f[1]==0){
        MyWarning("getMMWeight()","Cannot invert the matrix (i.e. determinant==0)");
        m_nDetIsZero++;
        return 0.;
      }
      Double_t det2leps = 1./(r[0]-f[0])/(r[1]-f[1]);
      mmWeight = det2leps * (r[0]*f[1]*Nrf + f[0]*r[1]*Nfr + f[0]*f[1]*Nff);

    }else if(m_nLepsForMM==3){
      MyError("getMMWeight()","3Lep MM has not been implemented yet.");
    }

    MyDebug("getMMWeight()",Form("Total fake weight: %f",mmWeight));

    return mmWeight;
  }

}//namespace MM
