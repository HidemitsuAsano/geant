void GetMCdataInfo(     ){

  TFile *_file0 = TFile::Open("ncout.root");
  //TFile *_file0 = TFile::Open("simIMpisigma_nSppim_v101.root");
  //TFile *_file0 = TFile::Open("simIMpisigma_K0n_ns_v20.root");
  std::cout << _file0->Print() << std::endl;

  const int pdg_ep = -11;
  const int pdg_em = 11;
  const int pdg_gamma=22;
  const int pdg_pip=211;
  const int pdg_pim=-211;
  const int pdg_neutron = 2112;
  const int pdg_proton = 2212;

  int nent = ncan_pdg->GetEntries();
  const int electron_bin = ncan_pdg->GetXaxis()->FindBin(pdg_em);
  int nelectron =  ncan_pdg->GetBinContent(electron_bin);
  std::cout << "electron " << (float)nelectron/(float)nent << std::endl;
  
  const int proton_bin = ncan_pdg->GetXaxis()->FindBin(pdg_proton);
  int nproton = ncan_pdg->GetBinContent(proton_bin);
  std::cout << "nproton " << (float)nproton/(float)nent << std::endl;
  
  const int heavy_bin = ncan_pdg->GetXaxis()->FindBin(4000);
  int nheavyion =  ncan_pdg->Integral(heavy_bin,999999);
  std::cout << "nheavy " << (float)nheavyion/(float)nent << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "parent particle study" << std::endl;
  int nent_p = ncan_parentpdg->GetEntries();
  const int gamma_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_gamma);
  int ngamma_p = ncan_parentpdg->GetBinContent(gamma_bin_p);
  std::cout << "gamma "  << (float)ngamma_p/(float)nent_p << std::endl;
  
  const int electron_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_em);
  int nelectron_p =  ncan_parentpdg->GetBinContent(electron_bin_p);
  std::cout << "electron " << (float)nelectron_p/(float)nent_p << std::endl;
   
  const int pim_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_pim);
  int npim_p = ncan_parentpdg->GetBinContent(pim_bin_p);
  std::cout << "pim     " << (float)npim_p/(float)nent_p << std::endl;
  
  const int pip_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_pip);
  int npip_p = ncan_parentpdg->GetBinContent(pip_bin_p);
  std::cout << "pip     " << (float)npip_p/(float)nent_p << std::endl;
  
  const int proton_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_proton);
  int nproton_p = ncan_parentpdg->GetBinContent(proton_bin_p);
  std::cout << "nproton " << (float)nproton_p/(float)nent << std::endl;
  
  const int neutron_bin_p = ncan_parentpdg->GetXaxis()->FindBin(pdg_neutron);     
  int neutron_p = ncan_parentpdg->GetBinContent(neutron_bin_p);
  std::cout << "nneutron " << (float)neutron_p/(float)nent_p << std::endl;
  
  const int heavyion_bin_p = ncan_parentpdg->GetXaxis()->FindBin(4000);
  int nheavyion_p =  ncan_parentpdg->Integral(heavyion_bin_p,999999);
  std::cout << "nheavy " << (float)nheavyion_p/(float)nent_p << std::endl;



}
