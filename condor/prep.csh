# Code to put your input rootfiles referenced in configs onto a condor-friendly storage area
cd /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_10_2_0/src/BStar13TeV/rootfiles/
tar -czvf bstar_presel_rootfiles.tgz TWpreselection*_tau32medium*default.root  TWpreselection*_tau32medium*ttbar.root smooth_QCD_*.root
xrdcp -f bstar_presel_rootfiles.tgz root://cmseos.fnal.gov//store/user/lcorcodi/bstar_presel_rootfiles.tgz
rm bstar_presel_rootfiles.tgz
cd /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_10_2_13/src/2DAlphabet

