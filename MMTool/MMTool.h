#ifndef MMTool_MMTool_H
#define MMTool_MMTool_H

#include"AsgTools/MsgLevel.h"
#include<TObject.h>

class TH2F;

namespace MM {

  enum LepFlav{undefined=-1, el=0, mu=1, tau=2, nLepFlavType};
  enum NLepMM{TwoLepMM, ThreeLepMM, nNLepMM};

  struct LepProp{
    LepFlav flav;
    Double_t pt;
    Double_t eta;
    Bool_t isSignal;
  };

  class MMTool
  {
  public:

    MMTool(NLepMM mmMode, const MSG::Level& debugMode=MSG::ERROR);
    ~MMTool();

    Bool_t prepEffHists(TString fileName, Int_t effDimension=0);

    Double_t getMMWeight();

    void setLepProp(const Int_t& xthLep, const LepProp& lep);//xthLep==0->Lead lep, xthLep==1->2nd Lead lep,,,
    void setLepProps(const LepProp& lep1, const LepProp& lep2);
    void setLepProps(const LepProp& lep1, const LepProp& lep2, const LepProp& lep3);

  protected:
    MSG::Level m_debugMode;
    Int_t m_mmMode;
    Int_t m_nLepsForMM;
    Int_t m_nDetIsZero;
    LepProp m_leps[3];

    TH2F* m_realElEff;
    TH2F* m_realMuEff;
    TH2F* m_fakeElEff;
    TH2F* m_fakeMuEff;

    void initialize(NLepMM mmMode);
    void finalize();
    void initLeps();

    void setEff(Double_t (&r)[3], Double_t (&f)[3], Double_t (&R)[3], Double_t (&F)[3]);
    Double_t getRealEff(const Int_t& xthLep);
    Double_t getFakeEff(const Int_t& xthLep);

  };

}// namespace MM

#endif
