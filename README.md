# Lepton_IDSF

N.B. we are using inputs in /eos/cms/store/group/phys_muon/balvarez/TnP_TopEFT_mvaTTHUL_merged/, and the example below is for electron.

1. step1: python3 Step1_makehist.py 

get the histograms of Data and DY MC. There are two kinds of histogram, one is for the events passing the selection "pass", the other one is for the events failing the selection "fail". There are two additional histograms considering the pileup systematic. Finally you should get many root files for different pt and eta bin, and in each root file you should have following histograms:
  KEY: TH1D	TnP_mass_DYpass;1	Pass_etap2p1Top2p4pt60To120
  KEY: TH1D	TnP_mass_DYpass_pileUp;1	Pass_etap2p1Top2p4pt60To120_pileUp
  KEY: TH1D	TnP_mass_DYpass_pileDo;1	Pass_etap2p1Top2p4pt60To120_pileDo
  KEY: TH1D	TnP_mass_DYfail;1	Fail_etap2p1Top2p4pt60To120
  KEY: TH1D	TnP_mass_DYfail_pileUp;1	Fail_etap2p1Top2p4pt60To120_pileUp
  KEY: TH1D	TnP_mass_DYfail_pileDo;1	Fail_etap2p1Top2p4pt60To120_pileDo
  KEY: TH1D	TnP_mass_Mupass;1	Pass_etap2p1Top2p4pt60To120
  KEY: TH1D	TnP_mass_Mufail;1	Fail_etap2p1Top2p4pt60To120

2. step2: following the Step2/README to get the efficiency of data and MC
3. step3: go to directory Step3, run "python3 final_sf.py" to get the final SF. N.B., there are some hard-coded directory names in final_sf.py correspond to the direcotories storing the outputs in step2, always remember to keep them consistent
