## Matrix Method (MM) tool for SUSY EW analyses  

Aiming to provide an easy tool for MM calculation.  
Working with RootCore.  

-- Maintained by TokyoTech group...  
-- Contact: minoru.hirose@cern.ch, daiki.hayakawa@cern.ch  

##Code snipets for usage...  
-- Include the header file.
    #include"MMTool/MMTool.h"

-- Instantiate MMTool and set proper efficiency before the event loop.  
    MM::MMTool m_MMTool = new MM::MMTool(MM::TwoLepMM,m_dbg);  
    m_MMTool->prepEffHists("Path to the file for real/fake efficiency");  
(In this case, MMTool returns the event weight for 2L Matrix Method (4x4 matrix).)

-- In each event, set leptons' property. Then, get MM weight for histograms and fake yields.
    MM::LepProp lepProp[2];
    for(Int_t lep=0; lep<2; lep++){
      lepProp[lep].pt       = (baseline lepton Pt )[xxx]; //in GeV
      lepProp[lep].eta      = (baseline lepton Eta)[xxx];
      lepProp[lep].flav     = (MM::el or MM::mu);
      lepProp[lep].isSignal = (kTRUE if baseline lepton passes the signal criteia);
      m_MMTool->setLepProp(lep,lepProp[lep]);
    }
    Double_t weight = m_MMTool->getMMWeight();

-- Delete MMTool for cleanup.  
    delete m_MMTool;  

Working example can be found on the link below.  
https://github.com/hrsmnr/xAOD_walkthrough  
