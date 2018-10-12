// v11 modified by Kawasaki 
// -> modified by Asano Aug.20th, 2018
// ====================================================================
#include "KnuclAnaManager.hh"
#include "KnuclDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4TwistedTubs.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4TwistedTubs.hh"

#include "G4PhysicalConstants.hh"


//#include "KnuclMaterialManager.hh"
#include "KnuclFieldSetup.hh"
#include "KnuclFieldSetup_uswk.hh"
#include "KnuclFieldSetup_dora.hh"
#include "KnuclFieldSetup_combi.hh"
#include "KnuclCommon.h"
#include "DetectorList.hh"

#define CHECK_OVERLAPS false

// ====================================================================
//
// class description
//
// ====================================================================

/////////////////////////////////////////////////////
KnuclDetectorConstruction::KnuclDetectorConstruction(KnuclAnaManager* ana)
 :  
    fEmFieldSetupUSHIWAKA(0), fEmFieldSetupCDC(0), 
    fEmFieldSetupMap_uswk(0), fEmFieldSetupMap_dora(0),
    fEmFieldSetupMap_combi(0),
    experimentalHall_box(0),  experimentalHall_log(0), experimentalHall_phys(0),
    CDS_tube(0),              CDS_log(0),              CDS_phys(0)
//////////////////////////////////////////////////////
{
  DORAFieldMapFlag    = ana->GetDORAFieldMapFlag();
  USWKFieldMapFlag    = ana->GetUSWKFieldMapFlag();
  FieldInCDC          = ana->GetFieldInCDC();
  FieldInUshiwaka     = ana->GetFieldInUshiwaka();
  neutron_efficiency_mode=false;
  if (ana->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode)
    {
      neutron_efficiency_mode=true;
    }
  materialMgr = 0;
  materialMgr = new KnuclMaterialManager();

  for(int i=0;i<100;i++){
    vnCounterStat.push_back(true);
  }
  ADCthreshold = ana->GetADCThreshold();
  SetGeomParameters(ana->GetConfFileName());
}

//////////////////////////////////////////////////////
KnuclDetectorConstruction::~KnuclDetectorConstruction()
//////////////////////////////////////////////////////
{
  if (fEmFieldSetupUSHIWAKA)  delete fEmFieldSetupUSHIWAKA;
  if (fEmFieldSetupCDC)       delete fEmFieldSetupCDC;
  if (fEmFieldSetupMap_uswk)  delete fEmFieldSetupMap_uswk;
  if (fEmFieldSetupMap_dora)  delete fEmFieldSetupMap_dora;
  if (fEmFieldSetupMap_combi) delete fEmFieldSetupMap_combi;
}

void KnuclDetectorConstruction::SetGeomParameters(std::string filename)
{
  const Int_t MAXCHARS = 512;
  char str[MAXCHARS],str1[MAXCHARS];

  for(int i=0;i<200;i++)
    for(int j=0;j<200;j++)
      for(int k=0;k<20;k++)
	param[i][j][k]=-999;

  std::string BLDCWireMapFileName;
  std::string CDCGeometryFileName;
  std::string GeometryMapFileNameCDS;
  std::string GeometryMapFileNameBL;
  std::string GeometryMapFileNameHall;

  ifstream fin;
  fin.open(filename.data());
  if (!fin.good()){
    cout<<"cant open file: "<<filename<<endl;
    exit(0);
  }
  while (!fin.eof()){
    fin.getline(str, MAXCHARS);
    if( str[0]=='#' ) continue;
    if( sscanf(str,"BLDCWireMap: %s", str1)==1 )
      BLDCWireMapFileName = str1;
    else if( sscanf(str,"CDCGeom: %s", str1)==1 )
      CDCGeometryFileName = str1;
    else if( sscanf(str,"GeomCDS: %s", str1)==1 )
      GeometryMapFileNameCDS = str1;
    else if( sscanf(str,"GeomBL: %s", str1)==1 )
      GeometryMapFileNameBL = str1;
    else if( sscanf(str,"GeomHall: %s", str1)==1 )
      GeometryMapFileNameHall = str1;
  }

  cout << "BLDCWireMap = " << BLDCWireMapFileName     << endl;
  ReadFile(BLDCWireMapFileName);
  cout << "CDCGeometry = " << CDCGeometryFileName     << endl;
  ReadFile(CDCGeometryFileName);
  cout << "GeomCDS     = " << GeometryMapFileNameCDS  << endl;
  ReadFile(GeometryMapFileNameCDS);
  cout << "GeomBL      = " << GeometryMapFileNameBL   << endl;
  ReadFile(GeometryMapFileNameBL);
  cout << "GeomHall    = " << GeometryMapFileNameHall << endl;
  ReadFile(GeometryMapFileNameHall);

}

void KnuclDetectorConstruction::ReadFile(std::string filename)
{
  const int MAXCHAR = 512;
  const int MAXTOKEN = 20;
  const char* DELIMITER = " ";
  
  int cid=-1;
  int seg=-1;
  const int nseg=20;
  double par[20];
  char buf[MAXCHAR];
  ifstream fin;
  
  fin.open(filename.c_str());
  if (!fin.good()){
    std::cerr << " File open fail. [" << filename << "]" << std::endl;
    exit(-1);
  }
  while (!fin.eof()){
    fin.getline(buf, MAXCHAR);
    if( buf[0]=='#' ) continue;
    int n = 0;
    const char* token[MAXTOKEN] = {};
    token[0] = std::strtok(buf, DELIMITER);
    if (token[0] != NULL) {
      if( !strcmp(token[0], "#") ) continue;
      for (n = 1; n < MAXTOKEN; n++){
	token[n] = std::strtok( NULL, DELIMITER);
	if (token[n] == NULL ) break;
      }
    }
    //    if(n!=npar+2) continue;
    for (int i = 0; i < n; i++){
      if(i==0) cid=atoi(token[i]);
      if(i==1) seg=atoi(token[i]);
      if(i >= 2 && i < nseg+2 )
	par[i-2]=atof(token[i]);
    }
    for(int i=0;i<n-2;i++)
      param[cid][seg][i]=par[i];
  }
  fin.close();  
}

G4LogicalVolume* KnuclDetectorConstruction::ConstructShape(Int_t CID){
  //====================================================
  //--- define detector based on CID
  //====================================================
  if(param[CID][0][0]==-999&&param[CID][1][0]==-999){
    std::cout<<CID<<" parameters not defined !!!"<<std::endl;
    return 0;
  }

  DetectorList* dlist=DetectorList::GetInstance();
  Int_t nsegments=dlist->GetNsegs(CID);
  TString name=dlist->GetName(CID);
  char tmpname[100];
  G4Material *medium1= materialMgr->GetMaterial("Air");
  std::string matname=dlist->GetMaterial(CID);
  G4Material *medium2;

  if(dlist->IsChamber(CID))
    medium2 = FindChamberGas(matname);
  else 
    medium2 = materialMgr->GetMaterial(matname);
  if(nsegments==0) medium1=medium2;  
#if 1
  std::cout<<"Construnct "
	   <<std::setw(5)<<CID
	   <<std::setw(10)<<name
	   <<std::setw(20)<<medium1->GetName()
	   <<std::setw(20)<<medium2->GetName()<<std::endl;
#endif
  int seg=0;  
  G4LogicalVolume* assembly_log=0;
  G4LogicalVolume *mother=0;

  G4VisAttributes* att=new G4VisAttributes();
  std::string colname=dlist->GetColor(CID);
  G4Colour tmpcol;
  G4Colour::GetColour(colname.data(), tmpcol);
  att->SetColour(tmpcol);
  //  att->SetForceSolid(true);
  //  att->SetForceWireframe(true);

  if(param[CID][seg][0]!=-999){

    
  //G4CSGSolid* solid=MakeShape(param[CID][seg],name.Data());  //original


    //kawasak    
    G4CSGSolid* solid_tmp=MakeShape(param[CID][seg],name.Data());  

    //--------------------------------------------------------------------------------------------
    //G4double BPC_add_rmin = 0.00*mm;      ///
    //G4double BPC_add_rmax = 84*mm;        ///
    //G4double BPC_add_dz   = 53*mm;        ///
      G4double BPC_add_rmin = 0.00*mm;      ///
      G4double BPC_add_rmax = 145*mm;       ///
      G4double BPC_add_dz   = 52.5*mm;      ///


    G4Tubs *BPC_add       = new G4Tubs("BPC_add",BPC_add_rmin, BPC_add_rmax, BPC_add_dz, 0, twopi);
    G4RotationMatrix* BPC_add_rotCounter = new G4RotationMatrix;
  //G4ThreeVector BPC_add_shift(0, 0, -88*mm);    ///
    G4ThreeVector BPC_add_shift(0, 0, -90*mm);    ///
    G4Transform3D BPC_add_tr(*BPC_add_rotCounter, BPC_add_shift);
    //---------------------------------------------------------------------------------------------

    G4Box* dummy_add= new G4Box("dummy_add",0.0000001, 0.0000001, 0.0000001);


    G4UnionSolid *solid;
    if(CID==CID_BPC)solid  = new G4UnionSolid(name.Data(), solid_tmp, BPC_add, BPC_add_tr);    
    if(CID!=CID_BPC)solid  = new G4UnionSolid(name.Data(), solid_tmp, dummy_add);    
    //kawasak
    


    sprintf(tmpname,"%s_log",name.Data());
    assembly_log=new G4LogicalVolume(solid,medium1,tmpname,0,0,0);
    G4Transform3D assembly_trans=MakeTrans(param[CID][seg],dlist->IsChamber(CID));
    if(nsegments==0&&!dlist->GetFlag(CID,1))   assembly_log->SetVisAttributes(att);
    else           assembly_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    if((int)param[CID][seg][10]>=0){
      std::string mothername=dlist->GetName((int)param[CID][seg][10]);
      sprintf(tmpname,"%s_log",mothername.data());
      mother=G4LogicalVolumeStore::GetInstance()->GetVolume(tmpname);
      if(!mother) mother=experimentalHall_log;  
      new G4PVPlacement(assembly_trans,assembly_log,name.Data(),mother,false,0,CHECK_OVERLAPS);
    }else if((int)param[CID][seg][10]==-1){
      experimentalHall_phys=new G4PVPlacement(assembly_trans,assembly_log,name.Data(),0,false,0,CHECK_OVERLAPS);
      return assembly_log;
    }
  }
  seg=1;
  if(!mother){
    if((int)param[CID][seg][10]>=0){
      std::string mothername=dlist->GetName((int)param[CID][seg][10]);
      sprintf(tmpname,"%s_log",mothername.data());
      mother=G4LogicalVolumeStore::GetInstance()->GetVolume(tmpname);
    }
  }
  if(!assembly_log) assembly_log=mother;
  if(!assembly_log) return 0;
  
  if(dlist->IsChamber(CID)){
    CreateChamberCells(CID,assembly_log,false);
  }else{
    if(nsegments==0) return assembly_log;
    G4LogicalVolume *segment_log=0;
    if(dlist->GetFlag(CID,0)){
      if(param[CID][seg][0]!=-999){
	G4CSGSolid* solid=MakeShape(param[CID][seg],name.Data());
	sprintf(tmpname,"s%s_log",name.Data());
	segment_log=new G4LogicalVolume(solid,medium2,tmpname,0,0,0);
	if(dlist->IsCounter(CID))
	  segment_log->SetSensitiveDetector(counterSD);
	if(!dlist->GetFlag(CID,1))   assembly_log->SetVisAttributes(att);
	else                        assembly_log->SetVisAttributes(G4VisAttributes::GetInvisible());
      }
    }
    for (seg=1; seg<=nsegments; seg++){
      if(param[CID][seg][0]==-999) continue;
      int volno=seg-1;
      if(!dlist->GetFlag(CID,0)){
	G4CSGSolid* solid=MakeShape(param[CID][seg],name.Data());
	sprintf(tmpname,"s%s_log%02d",name.Data(),seg);
	segment_log=new G4LogicalVolume(solid,medium2,tmpname,0,0,0);
	if(dlist->IsCounter(CID))
	  segment_log->SetSensitiveDetector(counterSD);
	if(!dlist->GetFlag(CID,1))   assembly_log->SetVisAttributes(att);
	else                         assembly_log->SetVisAttributes(G4VisAttributes::GetInvisible());
	//      volno=0;
      }
      G4Transform3D segment_trans = MakeTrans(param[CID][seg]);
      sprintf(tmpname,"s%s_phys%02d",name.Data(),seg);      
      new G4PVPlacement(segment_trans,segment_log,tmpname,assembly_log,false,volno,CHECK_OVERLAPS);
    }
  }
  return assembly_log;
}

G4CSGSolid* KnuclDetectorConstruction::MakeShape(double *par,const char* name){
  Double_t z    = par[9]*cm/2.0;
  G4CSGSolid *shape=0;
  char tmpname[100];
  if(z>0){
    Double_t rmin = par[6]*cm;
    Double_t rmax = par[7]*cm;
    Double_t phi  = par[8]*degree;    
    sprintf(tmpname,"%s_tube",name);
    shape= new G4Tubs(tmpname,rmin, rmax, z,0,phi);
  }else{
    Double_t x = par[6]*cm/2.0;
    Double_t y = par[7]*cm/2.0;
    z    = par[8]*cm/2.0;
    sprintf(tmpname,"%s_box",name);
    if(z>0)  shape= new G4Box(tmpname,x,y,z);    
  }
  return shape;
}

G4Transform3D KnuclDetectorConstruction::MakeTrans(double *par, bool ZFIRST){
  G4ThreeVector pos(par[0]*cm,par[1]*cm,par[2]*cm);
  if(par[9]!=0&&par[8]<90){
    Double_t pos_r   = par[0]*cm;
    Double_t pos_phi = par[1]*degree;
    Double_t pos_z   = par[2]*cm;
    pos.set(pos_r,0,pos_z);
    pos.rotateZ(pos_phi);
  }
  G4RotationMatrix rot;
  if(ZFIRST){
    rot.rotateZ(par[5]*degree);
    rot.rotateX(par[3]*degree);
    rot.rotateY(par[4]*degree);
  }else{
    rot.rotateX(par[3]*degree);
    rot.rotateY(par[4]*degree);
    rot.rotateZ(par[5]*degree);
  }
  G4Transform3D trans(rot,pos); 
  return trans;
}

G4Material* KnuclDetectorConstruction::FindChamberGas(const TString& chamberType)
{
  if ( chamberType == "CDCGas" ){
    return materialMgr-> GetMaterial("ArEthan_50_50");
  } else
  if ( chamberType == "BLDCGas" ){
    return materialMgr-> GetMaterial("ArIsobutane_80_20");
  } else {
    G4cout << "!!! Selected Chamber type " << chamberType
	   << " is not supported yet. " << G4endl; 
    return materialMgr-> GetMaterial("interGalactic");
  }
}

void KnuclDetectorConstruction::SetMagneticFields(){
  // ==============================================================
  // magnetic field
  // ==============================================================

  G4double fCDSField    = FieldInCDC   *tesla;  
  G4ThreeVector fCDC   (0.0, 0.0,          fCDSField);  
  G4double fUshiwakaField =  FieldInUshiwaka*tesla;  
  G4ThreeVector fUSHIWAKA(0.0, fUshiwakaField, 0.0      );

  G4cout<<"DORAFieldMapFlag = "<<DORAFieldMapFlag<<G4endl;
  G4cout<<"USWKFieldMapFlag = "<<USWKFieldMapFlag<<G4endl;

  if(DORAFieldMapFlag == 0 && USWKFieldMapFlag == 0){
    fEmFieldSetupUSHIWAKA = new KnuclFieldSetup(fUSHIWAKA) ;
    fEmFieldSetupCDC      = new KnuclFieldSetup(fCDC) ;
  }else if(DORAFieldMapFlag == 0 && USWKFieldMapFlag == 1){
    fEmFieldSetupMap_uswk = new KnuclFieldSetup_uswk("ushiwaka_field_map.txt", 0.*cm, 0.*cm, 250.*cm) ;
    fEmFieldSetupCDC      = new KnuclFieldSetup(fCDC) ;
  }else if(DORAFieldMapFlag == 1 && USWKFieldMapFlag == 0){
    fEmFieldSetupUSHIWAKA = new KnuclFieldSetup(fUSHIWAKA) ;
    //normalized to 0.7 Tesla at (0., 0., 0)
    fEmFieldSetupMap_dora = new KnuclFieldSetup_dora("Magnet_Field_Doraemon_3.txt", 0.*cm, 0.*cm, 0.*cm, FieldInCDC) ;
  }else{
    //normalized to 0.7 Tesla at (0., 0., 0)
    fEmFieldSetupMap_combi = new KnuclFieldSetup_combi("ushiwaka_field_map.txt", 0.*cm, 0.*cm, 250.*cm, "Magnet_Field_Doraemon_3.txt", 0.*cm, 0.*cm, 0.*cm, FieldInCDC) ;
  }
}

G4VPhysicalVolume*  KnuclDetectorConstruction::SetupForNCEfficencyStudy()
{

  // 
  // Fetch material information for Scintillator
  //
  G4Material* NE213     = materialMgr-> GetMaterial("NE213");
   
  //
  G4double LS_inner_diameter =  0.0    *cm;
  G4double LS_outer_diameter = 12.7/2.0*cm;
  G4double LS_length         = 12.7/2.0*cm;
  G4double LS_phi_start      =  0.0;
  G4double LS_phi_end        = twopi;
  G4Tubs *LS_tubs = new G4Tubs("LS_tubs", LS_inner_diameter, LS_outer_diameter, 
			 LS_length,LS_phi_start, LS_phi_end);

  G4LogicalVolume *LS_log = new G4LogicalVolume(LS_tubs, NE213, "LS_log", 0, 0, 0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,param[CID_NC][0][2]*cm),
		    LS_log, "LS", experimentalHall_log, false, 0, CHECK_OVERLAPS);
  
  LS_log->SetSensitiveDetector(counterSD);

  //--- Visualization ---//                                                                         
  G4VisAttributes *LS_att = new G4VisAttributes(Red);
  //    LS_att->SetForceWireframe(true);
  LS_log->SetVisAttributes(LS_att);

#if 0
  ConstructShape(CID_NC);
#endif
    G4cout << "neutron_efficiency_mode: KnuclDetectorConstruction completed" << G4endl;
    
    return experimentalHall_phys;  
}

void KnuclDetectorConstruction::CreateCDCCells(G4LogicalVolume *cdc_log,bool WIRE)
{
  G4Material* W                   = materialMgr-> GetMaterial("Tungsten");
  G4Material* Al                  = materialMgr-> GetMaterial("Aluminum");
  G4Material* ChamberGas = FindChamberGas("CDCGas");  

  int cid=CID_CDC;
  DetectorList* dlist=DetectorList::GetInstance();
  //  const char* name=dlist->GetName(cid).data();

  if(param[cid][0][0]==-999||param[cid][1][0]==-999){
    std::cout<<"CDC parameters not defined !!!"<<std::endl;
    exit(0);
  }
  G4ThreeVector xyzCounter(0,0,0);
  G4double cdc_z     = param[cid][0][9]*cm/2.; 

  const G4double cdc_wire_dist = 0.45*cm;
  const G4double cdc_off = cdc_wire_dist;

  char name_sol[64];
  char name_log[64];
  char name_phy[64];

  G4VisAttributes* att=new G4VisAttributes(!dlist->GetFlag(cid,1));
  std::string colname=dlist->GetColor(cid);
  G4Colour tmpcol;
  G4Colour::GetColour(colname.data(), tmpcol);
  att->SetColour(tmpcol);
  //  att->SetForceSolid(true);
  att->SetForceWireframe(true);

  // CDC cell
  {
    for (G4int i=0; i<15; i++){
      G4double cdc_radius     = param[cid][i+1][1]*cm; 
      G4double cdc_phi0       = param[cid][i+1][2]*degree; 
      G4double cdc_dphi       = param[cid][i+1][3]*degree; 
      G4double cdc_cell_twist = -param[cid][i+1][4]*degree; 
      G4int    cdc_nwire      = (int)param[cid][i+1][0];
      G4double rmin = cdc_radius-cdc_off;
      G4double rmax = cdc_radius+cdc_off;

      //      std::cout<<i+1<<"  "<<cdc_nwire<<"  "<<cdc_radius<<"  "<<cdc_phi0/degree<<"  "<<cdc_phi0/degree<<"  "<<cdc_cell_twist/degree<<std::endl;      
      sprintf(name_sol,"CDC_sol_%2.2d",i);
      sprintf(name_log,"CDC_log_%2.2d",i);
      
      G4TwistedTubs *tubs = new G4TwistedTubs(name_sol, cdc_cell_twist, rmin, rmax,
					      cdc_z, cdc_dphi);
      G4LogicalVolume *log = new G4LogicalVolume(tubs,ChamberGas,name_log,0,0,0);
      log->SetSensitiveDetector(chamberSD);
      log->SetVisAttributes(att);      
      for (G4int j=0; j<cdc_nwire; j++){
	sprintf(name_phy,"CDC_phys%2.2d%03d",i+1,j+1);
	G4RotationMatrix* rotCounter = new G4RotationMatrix;
	G4double rotangle=cdc_dphi*j+cdc_phi0-cdc_cell_twist/2.;
	rotCounter->rotateZ(rotangle);	
	G4Transform3D posCounter(*rotCounter, xyzCounter);
	new G4PVPlacement(posCounter,log, name_phy, cdc_log, false, j, CHECK_OVERLAPS);  
      }
    }
  }  
    
  if(!WIRE) return;
  //*****************//
  //*** CDC wires ***//
  //*****************//
  G4double r0 = 17.5*cm;

  G4double  radius[67];
  G4int     nwires[67];
  G4double  tilt[67];
  for (int i=0; i<67; i++){
    if(i==0){
      radius[i] = r0;
    }else{
      radius[i] = radius[i-1]+cdc_wire_dist;
    }
    if(i==1 || i==12 || i==13 || i==21 || i==22 || i==30 || i==31 ||
       i==39 || i==40 || i==48 || i==49 || i==57 || i==58 || i==66){
      radius[i] += 0.2*cm;
    }
    
    if (i<12           ) {
       nwires[i] =  72;
       tilt[i] = 0;
    }
    else
      if (i>=12 && i<21  ) {
       nwires[i] =  90;
       //tilt[i] = i==12 ? 0 : twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]));
       tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
      }
    else
      if (i>=21 && i<30  ) { 
       nwires[i] = 100;
       //tilt[i] = i==21 ? 0 : -(twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i])));
       tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
      }
    else
      if (i>=30 && i<39  ) { 
       nwires[i] = 120;
       tilt[i] = 0;
      }
    else
      if (i>=39 && i<48  ) { 
       nwires[i] = 150;
       //tilt[i] = i==39 ? 0 : twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]));
       tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
      }
    else
      if (i>=48 && i<57  ) { 
       nwires[i] = 160;
       //tilt[i] = i==48 ? 0 : -(twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i])));
       tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
      }
    else
      if (i>=57 && i<=66 ) {
       nwires[i] = 180;
       tilt[i] = 0;
      }
    else {
      nwires[i] =   0;
    }
  }
  //for (int i=0; i<67; i++){
  //G4cout << radius[i] << " " << nwires[i] << " " << tilt[i]/deg << G4endl;
  //}   
 
  G4int total_sense_wires = 0;
  G4int total_field_wires = 0;

  G4int   idd = 0;
  G4int   IsSenseWire = 0;

  //--- Visualization ---//
  G4VisAttributes *CDCSenseWire_att = new G4VisAttributes(Red);
  //CDCSenseWire_att->SetForceWireframe(true);
  CDCSenseWire_att->SetForceSolid(true);
  G4VisAttributes *CDCFieldWire_att = new G4VisAttributes(Green);
  //CDCFieldWire_att->SetForceWireframe(true);
  CDCFieldWire_att->SetForceSolid(true);

  for (G4int i=0; i<67; i++) {
    if ( i== 0 || i== 1 || i== 3 || i== 5 || i== 7 || i== 9|| i==11 || i==12 || i==13 || i==15 || i==17 ||
         i==19 || i==21 || i==22 || i==24 || i==26 || i==28|| i==30 || i==31 || i==33 || i==35 || i==37 ||
         i==39 || i==40 || i==42 || i==44 || i==46 || i==48|| i==49 || i==51 || i==53 || i==55 || i==57 ||
         i==58 || i==60 || i==62 || i==64 || i==66 ) {
      idd = 0;
    } else {
      idd = 1;
    }

    if (i== 3 || i== 6 || i== 9 || i==15 || i==18 || i==24 || i==27 || i==33 || i==36 || 
        i==42 || i==45 || i==51 || i==54 || i==60 || i==63             )
      IsSenseWire = 1; 
    else 
      IsSenseWire = 0;

    G4Tubs *sense_tub = new G4Tubs("CDC_sense_Wire", 0.0, 0.03*mm/2., cdc_z/cos(fabs(tilt[i])), 0.0, twopi);
    G4Tubs *field_tub = new G4Tubs("CDC_field_Wire", 0.0, 0.10*mm/2., cdc_z/cos(fabs(tilt[i])), 0.0, twopi);

    G4LogicalVolume* sense_log =  new G4LogicalVolume(sense_tub, W,  "CDCSWire_log", 0,0,0);
    G4LogicalVolume* field_log =  new G4LogicalVolume(field_tub, Al, "CDCFWire_log", 0,0,0);
    
    //    sense_log->SetVisAttributes(CDCSenseWire_att);
    //    field_log->SetVisAttributes(CDCFieldWire_att);
    sense_log->SetVisAttributes(G4VisAttributes::Invisible);
    field_log->SetVisAttributes(G4VisAttributes::Invisible);

    G4double unit_ang = twopi/(G4double)nwires[i];
    G4double Rwire    = radius[i];

    for (int j=0; j<nwires[i]; j++){
      G4double xxx = 0.0;
      G4double yyy = 0.0;
      G4double ang = 0.0;
      if (idd==0) {
	ang = unit_ang*(G4double)j;
        xxx = Rwire*cos(ang);
        yyy = Rwire*sin(ang); 
      } else 
      if (idd==1) {
	ang = unit_ang*(G4double)j + unit_ang/2.0;
        xxx = Rwire*cos(ang);
        yyy = Rwire*sin(ang); 
      }
      G4ThreeVector wire_pos(xxx, yyy, 0.0);
      G4RotationMatrix* wire_rot = new G4RotationMatrix;
      wire_rot->rotateX(-tilt[i]);
      wire_rot->rotateZ(ang);
      G4Transform3D wire_geom(*wire_rot, wire_pos);
      //  /**
      if (IsSenseWire == 1 ) {
        sprintf(name_phy,"CDCSWire_phy%06d", total_sense_wires);
	new G4PVPlacement(wire_geom, sense_log, name_phy, cdc_log, false, 0, CHECK_OVERLAPS ); 	
        total_sense_wires++;
      }
      else 
      if (IsSenseWire == 0 ) 
      {
        sprintf(name_phy,"CDCFWire_phy%06d", total_field_wires);
	new G4PVPlacement(wire_geom, field_log, name_phy, cdc_log, false, 0, CHECK_OVERLAPS );
	total_field_wires++;
      } else {
        G4cout << "Wrong IsSenseWire assigned " << IsSenseWire << G4endl;
      }
      
    }
  }  
  //G4cout << total_sense_wires << " " << total_field_wires << G4endl;
  //G4cout << "Number of wires generated = " << total_wires << G4endl;
  return;
}

void KnuclDetectorConstruction::CreateBLCCells(int cid,G4LogicalVolume* blc_log, bool WIRE)
{
  if(param[cid][0][0]==-999||param[cid][1][0]==-999){
    std::cout<<cid<<" parameters not defined !!!"<<std::endl;
    return;
  }
  DetectorList* dlist=DetectorList::GetInstance();
  const char* name=dlist->GetName(cid).data();
  G4Material* ChamberGas = FindChamberGas("BLDCGas");
  
  const G4double BLC_cell_size = param[cid][1][4]*cm;
  G4double BLC_seg_x = fabs(BLC_cell_size)/2;
  G4double BLC_seg_y = param[cid][1][5]*cm/2;
  G4double BLC_seg_z = 0.000001/2*mm;

  char box_name[20];
  char log_name[20];
  char phys_name[20];
  sprintf(box_name, "s%s_box", name);
  G4Box *box = new G4Box(box_name,BLC_seg_x,BLC_seg_y,BLC_seg_z);
  sprintf(log_name, "s%s_log", name);
  G4LogicalVolume* log = new G4LogicalVolume(box, ChamberGas, log_name, 0,0,0);
  log->SetSensitiveDetector(chamberSD);  

  G4VisAttributes* att=new G4VisAttributes();
  std::string colname=dlist->GetColor(cid);
  G4Colour tmpcol;
  G4Colour::GetColour(colname.data(), tmpcol);
  att->SetColour(tmpcol);
  //  att->SetForceSolid(true);
  att->SetForceWireframe(true);
  log->SetVisAttributes(att);      

  for (G4int i=1; i<=dlist->GetNlayers(cid); i++){
    for (G4int j=1; j<=dlist->GetNwires(cid); j++){
      sprintf(phys_name, "s%s_phys%d%02d", name, i, j);
      G4double xyflag=param[cid][i][2];
      G4double tmpx = param[cid][i][3]*cm+BLC_cell_size*(j-1);
      G4double tmpz = param[cid][i][1]*cm;
      G4ThreeVector xyzCounter;
      if(xyflag==0) xyzCounter = G4ThreeVector(tmpx, 0, tmpz);
      else          xyzCounter = G4ThreeVector(0, tmpx, tmpz);
      G4RotationMatrix* rotCounter = new G4RotationMatrix;
      rotCounter->rotateZ(xyflag*90*deg);
      G4Transform3D posCounter(*rotCounter, xyzCounter);
      new G4PVPlacement(posCounter, log, phys_name, blc_log, false, j-1, CHECK_OVERLAPS);
    }
  }  
  if(!WIRE) return;
}

void KnuclDetectorConstruction::CreateFDCCells(int cid,G4LogicalVolume* blc_log, bool WIRE)
{
  if(param[cid][0][0]==-999||param[cid][1][0]==-999){
    std::cout<<cid<<" parameters not defined !!!"<<std::endl;
    return;
  }
  DetectorList* dlist=DetectorList::GetInstance();
  const char* name=dlist->GetName(cid).data();
  G4Material* ChamberGas = FindChamberGas("BLDCGas");
  
  const G4double BLC_cell_size = param[cid][1][4]*cm;
  G4double BLC_seg_x = fabs(BLC_cell_size)/2;
  G4double BLC_seg_y = param[cid][1][5]*cm/2;
  G4double BLC_seg_z = 0.000001/2*mm;

  char box_name[20];
  char log_name[20];
  char phys_name[20];  
  sprintf(box_name, "s%s_box", name);
  G4Box *box = new G4Box(box_name,BLC_seg_x,BLC_seg_y,BLC_seg_z);
  sprintf(log_name, "s%s_log", name);
  G4LogicalVolume* log = new G4LogicalVolume(box, ChamberGas, log_name, 0,0,0);
  log->SetSensitiveDetector(chamberSD);  

  G4VisAttributes* att=new G4VisAttributes();
  std::string colname=dlist->GetColor(cid);
  G4Colour tmpcol;
  G4Colour::GetColour(colname.data(), tmpcol);
  att->SetColour(tmpcol);
  //  att->SetForceSolid(true);
  att->SetForceWireframe(true);
  log->SetVisAttributes(att);      

  for (G4int i=1; i<=dlist->GetNlayers(cid); i++){
    for (G4int j=1; j<=dlist->GetNwires(cid); j++){
      sprintf(phys_name, "s%s_phys%d%02d", name, i, j);
      G4double rotangle=param[cid][i][7]*degree;
      G4double xyflag=param[cid][i][2];
      G4double tmpx = param[cid][i][3]*cm+BLC_cell_size*(j-1);
      if(i==3||i==4) tmpx = param[cid][i][3]*cm+BLC_cell_size*(j-1);
      else tmpx = (param[cid][i][3]*cm+BLC_cell_size*(j-1))/cos(rotangle);
      G4double tmpz = param[cid][i][1]*cm;
      G4ThreeVector xyzCounter;
      if(xyflag==0) xyzCounter = G4ThreeVector(tmpx, 0, tmpz);
      else          xyzCounter = G4ThreeVector(0, tmpx, tmpz);
      G4RotationMatrix* rotCounter = new G4RotationMatrix;
      rotCounter->rotateZ(rotangle);
      G4Transform3D posCounter(*rotCounter, xyzCounter);
      new G4PVPlacement(posCounter, log, phys_name, blc_log, false, j-1, CHECK_OVERLAPS);
    }
  }  
  if(!WIRE) return;
}

void KnuclDetectorConstruction::CreateBPCCells(int cid,G4LogicalVolume* blc_log, bool WIRE)
{
  if(param[cid][0][0]==-999||param[cid][1][0]==-999){
    std::cout<<cid<<" parameters not defined !!!"<<std::endl;
    return;
  }
  DetectorList* dlist=DetectorList::GetInstance();
  const char* name=dlist->GetName(cid).data();
  G4Material* ChamberGas = FindChamberGas("BLDCGas");
  
  const G4double BLC_cell_size = param[cid][1][4]*cm;
  G4double BLC_seg_x = fabs(BLC_cell_size)/2;
//G4double BLC_r     = (param[cid][0][7]-0.2)*cm;        //original
  G4double BLC_r     = (param[cid][0][7]-0.2 -2.4)*cm;   //kawasak

  G4double BLC_seg_z = 0.000001/2*mm;
  char box_name[20];
  char log_name[20];
  char phys_name[20];

  G4VisAttributes* att=new G4VisAttributes();
  std::string colname=dlist->GetColor(cid);
  G4Colour tmpcol;
  G4Colour::GetColour(colname.data(), tmpcol);
  att->SetColour(tmpcol);
  //  att->SetForceSolid(true);
  att->SetForceWireframe(true);

  G4Box *box;
  G4LogicalVolume *log;
  for (G4int i=1; i<=dlist->GetNlayers(cid); i++){
    for (G4int j=1; j<=dlist->GetNwires(cid); j++){
      sprintf(box_name, "s%s_%d_%d_box", name,i,j);
      sprintf(log_name, "s%s_%d_%d_log", name,i,j);
      G4double tmpz = param[cid][i][1]*cm;
      G4double xyflag=param[cid][i][2];
      G4double tmpx = param[cid][i][3]*cm+BLC_cell_size*(j-1);
      G4double tmpx2 = tmpx>0 ? param[cid][i][3]*cm+BLC_cell_size*(j-1+0.5) :
	param[cid][i][3]*cm+BLC_cell_size*(j-1-0.5);
      G4double BLC_seg_y = BLC_r*sqrt(1-tmpx2*tmpx2/BLC_r/BLC_r);
      //std::cerr<<BLC_r<<" "<<i<<" "<<j<<" "<<tmpx<<" "<<tmpx2<<" "<<BLC_seg_y<<std::endl;
      box = new G4Box( box_name, BLC_seg_x,BLC_seg_y,BLC_seg_z);
      log = new G4LogicalVolume( box, ChamberGas, log_name, 0,0,0);
      log->SetSensitiveDetector(chamberSD);  
      log->SetVisAttributes(att);      

      sprintf(phys_name, "s%s_phys%d%02d", name, i, j);
      G4ThreeVector xyzCounter;
      if(xyflag==0) xyzCounter = G4ThreeVector(tmpx, 0, tmpz);
      else          xyzCounter = G4ThreeVector(0, tmpx, tmpz);
      G4RotationMatrix* rotCounter = new G4RotationMatrix;
      rotCounter->rotateZ(xyflag*90*deg);
      G4Transform3D posCounter(*rotCounter, xyzCounter);
      new G4PVPlacement(posCounter, log, phys_name, blc_log, false, j-1, CHECK_OVERLAPS);
    }
  } 

  
  //LAYER FLAME
  G4Material* G10 = materialMgr-> GetMaterial("CFRP");     ///



  //for cathode layer

  //G4double rmin1 = 58.5*mm;  ///
  //G4double rmax1 = 78.5*mm;  ///
  //G4double dz1   = 0.50*mm;  ///

    G4double rmin1 = 100*mm;    ///
    G4double rmax1 = 118.5*mm;  ///
    G4double dz1   = 0.50*mm;   ///


  G4Tubs *lay1  = new G4Tubs("lay1",rmin1, rmax1, dz1, 0, twopi);
  G4LogicalVolume *ca_log = new G4LogicalVolume(lay1, G10, "ca_log", 0, 0, 0); 



  //for anode & spacer (cathode layer -> width x 3)

  G4double rmin2 = rmin1; 
  G4double rmax2 = rmax1;  
  G4double dz2   = dz1*3;  

  G4Tubs *lay2  = new G4Tubs("lay2",rmin2, rmax2, dz2, 0, twopi);


  //for spacer

  //G4double box1_dx = 100*mm;     ///
  //G4double box1_dy = 13.777*mm;  ///
  //G4double box1_dz = 1.5*mm;     ///

    G4double box1_dx = 200*mm;     ///
    G4double box1_dy = 34.26*mm;  ///
    G4double box1_dz = 1.5*mm;     ///



  G4Box *box1  = new G4Box("box1", box1_dx, box1_dy, box1_dz);


  G4RotationMatrix* box1_rotCounter = new G4RotationMatrix;

//G4ThreeVector shift1(0, 53.767*mm, 0);    ///
  G4ThreeVector shift1(0, 84.24*mm, 0);    ///

  G4Transform3D tr1(*box1_rotCounter, shift1);

  G4IntersectionSolid *int1 = new  G4IntersectionSolid("int1", lay2, box1 ,tr1);


//G4ThreeVector shift2(0, -53.767*mm, 0);    ///
  G4ThreeVector shift2(0, -84.24*mm, 0);    ///
  G4Transform3D tr2(*box1_rotCounter, shift2);

  G4IntersectionSolid *int2 = new  G4IntersectionSolid("int2", lay2, box1 ,tr2);


  G4UnionSolid *sp_int  = new G4UnionSolid("sp_int", int1, int2);
  G4UnionSolid *sp      = new G4UnionSolid("sp_tmp", lay1, sp_int);




  //for anode
  //G4double box2_dx = 100*mm; ///
  //G4double box2_dy = 55*mm;  ///
  //G4double box2_dz = 1*mm;   ///

    G4double box2_dx = 200*mm; ///
    G4double box2_dy = 80*mm;  ///
    G4double box2_dz = 1*mm;   ///


  G4Box *box2  = new G4Box("box2", box2_dx, box2_dy, box2_dz);


  G4RotationMatrix* box2_rotCounter = new G4RotationMatrix;

  G4ThreeVector shift3(0,-15*mm, 0.5*mm);    ///
  G4Transform3D tr3(*box2_rotCounter, shift3);

  G4SubtractionSolid *an = new G4SubtractionSolid("an", lay2, box2, tr3);


  G4LogicalVolume *an_log = new G4LogicalVolume(an, G10, "an_log", 0, 0, 0);
  G4LogicalVolume *sp_log = new G4LogicalVolume(sp, G10, "sp_log", 0, 0, 0);


  for (G4int i=1; i<=dlist->GetNlayers(cid); i++){

  G4RotationMatrix* flame_rotCounter = new G4RotationMatrix;

  G4double flame_tmpz;
 

 //flame_tmpz  = (param[cid][i][1] )*cm;     ///
   flame_tmpz  = 27.5*mm;                    ///

  G4ThreeVector ca_xyzCounter(0., 0., flame_tmpz -8*(i-1)*mm );
  G4Transform3D ca_posCounter(*flame_rotCounter, ca_xyzCounter);
  new G4PVPlacement(ca_posCounter, ca_log, "ca_phy", blc_log, false, i-1, CHECK_OVERLAPS);

  G4ThreeVector ca2_xyzCounter(0., 0., flame_tmpz -8*(i-1)*mm -4*mm);
  G4Transform3D ca2_posCounter(*flame_rotCounter, ca2_xyzCounter);
  new G4PVPlacement(ca2_posCounter, ca_log, "ca_phy", blc_log, false, i-1+8, CHECK_OVERLAPS);



  G4ThreeVector an_xyzCounter(0., 0., flame_tmpz -8*(i-1)*mm -2*mm);
  G4Transform3D an_posCounter(*flame_rotCounter, an_xyzCounter);
  new G4PVPlacement(an_posCounter, an_log, "an_phy", blc_log, false, i-1, CHECK_OVERLAPS);



  G4ThreeVector sp_xyzCounter(0., 0., flame_tmpz -8*(i-1)*mm -6*mm);
  G4Transform3D sp_posCounter(*flame_rotCounter, sp_xyzCounter);
  new G4PVPlacement(sp_posCounter, sp_log, "sp_phy", blc_log, false, i-1, CHECK_OVERLAPS);
  }
  
  
  
  //ASD
  G4Material* GFRP = materialMgr-> GetMaterial("CFRP");   ///

  G4double asd_dx  = 4.05*mm;    ///
  G4double asd_dy  = 22.5*mm;    ///
  G4double asd_dz  = 37.75*mm;   ///


  G4Box *asd   = new G4Box("asd", asd_dx, asd_dy, asd_dz);
  G4LogicalVolume *asd_log = new G4LogicalVolume( asd, GFRP, "asd_log", 0,0,0);


  G4RotationMatrix* asd_rotCounter = new G4RotationMatrix;
//G4double asd_dphi = 45*deg;    ///
  G4double asd_dphi = 22.5*deg;    ///


  G4ThreeVector asd_xyzCounter;
//asd_xyzCounter = G4ThreeVector((73-4.05)*mm, 0, -101.05*mm); ///
  asd_xyzCounter = G4ThreeVector((134-4.05)*mm, 0, -103.65*mm); ///


//char asd_phy_name[8]; 
  char asd_phy_name[16]; 


//for(int i=0; i<8; i++){
  for(int i=0; i<16; i++){
  sprintf(asd_phy_name, "asd_phy_%s_%d", name,i+1);

  asd_rotCounter->rotateZ(asd_dphi*i +asd_dphi/2);    
  asd_xyzCounter.rotateZ(asd_dphi*i +asd_dphi/2);
  
  G4Transform3D asd_posCounter(*asd_rotCounter, asd_xyzCounter);
  new G4PVPlacement(asd_posCounter, asd_log, asd_phy_name, blc_log, false, i, CHECK_OVERLAPS);
  }
  
  

  
  //WIRE

  G4Material* Tun                   = materialMgr-> GetMaterial("Tungsten"); 
  G4Material* CuBe                  = materialMgr-> GetMaterial("Tungsten");      ///

  G4double sense_phi =0.0125*mm;  ///
  G4double field_phi =0.05*mm;    ///


  char sense_tub_name[200];
  char sense_log_name[200];
  char sense_phy_name[200];

  char field_tub_name[200];
  char field_log_name[200];
  char field_phy_name[200];


  for (G4int i=1; i<=dlist->GetNlayers(cid); i++){
  for (G4int j=1; j<=dlist->GetNwires(cid);  j++){


      sprintf(sense_tub_name, "sense_tub_%s_%d_%d", name,i,j);
      sprintf(sense_log_name, "sense_log_%s_%d_%d", name,i,j);
      sprintf(sense_phy_name, "sense_phy_%s_%d_%d", name,i,j);

      sprintf(field_tub_name, "field_tub_%s_%d_%d", name,i,j);
      sprintf(field_log_name, "field_log_%s_%d_%d", name,i,j);
      sprintf(field_phy_name, "field_phy_%s_%d_%d", name,i,j);


      G4double tmpz  = param[cid][i][1]*cm;
      G4double xyflag= param[cid][i][2];


      G4double sense_tmpx  = param[cid][i][3]*cm+BLC_cell_size*(j-1);
      G4double sense_tmpx2 = sense_tmpx>0 ? param[cid][i][3]*cm+BLC_cell_size*(j-1+0.5) :
	param[cid][i][3]*cm+BLC_cell_size*(j-1-0.5);
      G4double sense_BLC_seg_y = BLC_r*sqrt(1-sense_tmpx2*sense_tmpx2/BLC_r/BLC_r);


      G4double field_tmpx  = param[cid][i][3]*cm+BLC_cell_size*(j-1) +BLC_cell_size*0.5;  ///
      G4double field_tmpx2 = field_tmpx>0 ? param[cid][i][3]*cm+BLC_cell_size*(j-1+0.5) :
	param[cid][i][3]*cm+BLC_cell_size*(j-1-0.5);
      G4double field_BLC_seg_y = BLC_r*sqrt(1-field_tmpx2*field_tmpx2/BLC_r/BLC_r);


      G4Tubs *sense_tub = new G4Tubs(sense_tub_name, 0, sense_phi/2, sense_BLC_seg_y, 0, twopi);
      G4Tubs *field_tub = new G4Tubs(field_tub_name, 0, field_phi/2, field_BLC_seg_y, 0, twopi);

      G4LogicalVolume* sense_log =  new G4LogicalVolume(sense_tub, Tun,  sense_log_name, 0,0,0);
      G4LogicalVolume* field_log =  new G4LogicalVolume(field_tub, CuBe, field_log_name, 0,0,0);



      G4RotationMatrix* rotCounter = new G4RotationMatrix;
      rotCounter->rotateX(90*deg);
      rotCounter->rotateZ(xyflag*90*deg);


      G4ThreeVector sense_xyzCounter;
      if(xyflag==0) sense_xyzCounter = G4ThreeVector(sense_tmpx, 0, tmpz);
      else          sense_xyzCounter = G4ThreeVector(0, sense_tmpx, tmpz);
      G4Transform3D sense_posCounter(*rotCounter, sense_xyzCounter);


      G4ThreeVector field_xyzCounter;
      if(xyflag==0) field_xyzCounter = G4ThreeVector(field_tmpx, 0, tmpz);
      else          field_xyzCounter = G4ThreeVector(0, field_tmpx, tmpz);
      G4Transform3D field_posCounter(*rotCounter, field_xyzCounter);


      //new G4PVPlacement(sense_posCounter, sense_log, sense_phy_name, blc_log, false, (i-1)*16+(j-1), CHECK_OVERLAPS);
      //new G4PVPlacement(field_posCounter, field_log, field_phy_name, blc_log, false, (i-1)*16+(j-1), CHECK_OVERLAPS);
        new G4PVPlacement(sense_posCounter, sense_log, sense_phy_name, blc_log, false, (i-1)*32+(j-1), CHECK_OVERLAPS);
        new G4PVPlacement(field_posCounter, field_log, field_phy_name, blc_log, false, (i-1)*32+(j-1), CHECK_OVERLAPS);


  }
  }
  
 
  
  //HOUSING
  
  
  G4Material* Al                  = materialMgr-> GetMaterial("Aluminum");

  //end_cup
  //G4double end_cup_rmin = 58.5*mm;  ///
  //G4double end_cup_rmax = 84.0*mm;  ///
  //G4double end_cup_dz   = 1.00*mm;  ///
    G4double end_cup_rmin = 100*mm;  ///
    G4double end_cup_rmax = 145*mm;  ///
    G4double end_cup_dz   = 1.5*mm;  ///



  G4ThreeVector end_cup_xyzCounter;
//end_cup_xyzCounter = G4ThreeVector(0, 0, 34.0*mm); ///
  end_cup_xyzCounter = G4ThreeVector(0, 0, 34.5*mm); ///
  G4RotationMatrix* end_cup_rotCounter = new G4RotationMatrix;
  G4Transform3D end_cup_posCounter(*end_cup_rotCounter, end_cup_xyzCounter);


  G4Tubs *end_cup  = new G4Tubs("end_cup", end_cup_rmin, end_cup_rmax, end_cup_dz, 0, twopi);
  G4LogicalVolume *end_cup_log = new G4LogicalVolume(end_cup, Al, "end_cup_log", 0, 0, 0);  ///
  new G4PVPlacement(end_cup_posCounter, end_cup_log, "end_cup_phy", blc_log, false, 0, CHECK_OVERLAPS);



  //side_cover
  //G4double side_cover_rmin = 83.0*mm;  ///
  //G4double side_cover_rmax = 84.0*mm;  ///
  //G4double side_cover_dz   = 38.1*mm;  ///
    G4double side_cover_rmin = 144*mm;  ///
    G4double side_cover_rmax = 145*mm;  ///
    G4double side_cover_dz   = 38.1*mm;  ///


  G4ThreeVector side_cover_xyzCounter;
  side_cover_xyzCounter = G4ThreeVector(0, 0, 5.1*mm); ///
  G4RotationMatrix* side_cover_rotCounter = new G4RotationMatrix;
  G4Transform3D side_cover_posCounter(*side_cover_rotCounter, side_cover_xyzCounter);


  G4Tubs *side_cover  = new G4Tubs("side_cover", side_cover_rmin, side_cover_rmax, side_cover_dz, 0, twopi);
  G4LogicalVolume *side_cover_log = new G4LogicalVolume(side_cover, Al, "side_cover_log", 0, 0, 0);  ///
  new G4PVPlacement(side_cover_posCounter, side_cover_log, "side_cover_phy", blc_log, false, 0, CHECK_OVERLAPS);

  

  //base
  //G4double base_rmin = 58.5*mm;  ///
  //G4double base_rmax = 84.0*mm;  ///
  //G4double base_dz   = 5.75*mm;  ///
    G4double base_rmin = 100*mm;   ///
    G4double base_rmax = 145*mm;   ///
    G4double base_dz   = 6.5*mm;   ///


  G4ThreeVector base_xyzCounter;
//base_xyzCounter = G4ThreeVector(0, 0, -48.95*mm); ///
  base_xyzCounter = G4ThreeVector(0, 0, -49.7*mm);  ///
  G4RotationMatrix* base_rotCounter = new G4RotationMatrix;
  G4Transform3D base_posCounter(*base_rotCounter, base_xyzCounter);


  G4Tubs *base  = new G4Tubs("base", base_rmin, base_rmax, base_dz, 0, twopi);
  G4LogicalVolume *base_log = new G4LogicalVolume(base, Al, "base_log", 0, 0, 0);  ///
  new G4PVPlacement(base_posCounter, base_log, "base_phy", blc_log, false, 0, CHECK_OVERLAPS);

  
  //fix_ring
  //G4double fix_ring_rmin = 58.5*mm;    ///
  //G4double fix_ring_rmax = 84.0*mm;    ///
  //G4double fix_ring_dz   = 1.00*mm;    ///
  //G4double fix_ring_rmin = 119.5*mm;   ///
    G4double fix_ring_rmin = 132.25*mm;  ///

    G4double fix_ring_rmax = 145*mm;     ///
    G4double fix_ring_dz   = 1.00*mm;    ///

  G4ThreeVector fix_ring_xyzCounter;
//fix_ring_xyzCounter = G4ThreeVector(0, 0, -136.263*mm); ///
  fix_ring_xyzCounter = G4ThreeVector(0, 0, -138.23*mm);  ///
  G4RotationMatrix* fix_ring_rotCounter = new G4RotationMatrix;
  G4Transform3D fix_ring_posCounter(*fix_ring_rotCounter, fix_ring_xyzCounter);


  G4Tubs *fix_ring  = new G4Tubs("fix_ring", fix_ring_rmin, fix_ring_rmax, fix_ring_dz, 0, twopi);
  G4LogicalVolume *fix_ring_log = new G4LogicalVolume(fix_ring, Al, "fix_ring_log", 0, 0, 0);  ///
  new G4PVPlacement(fix_ring_posCounter, fix_ring_log, "fix_ring_phy", blc_log, false, 0, CHECK_OVERLAPS);
  


  //pole
  //G4double pole_rmin = 0*mm;  
  //G4double pole_rmax = 3.175*mm;    ///
  //G4double pole_dz   = 42.7815*mm;  ///
  //G4double pole_dphi = 45*deg;      ///
    G4double pole_rmin = 0*mm;  
    G4double pole_rmax = 3.175*mm;    ///
    G4double pole_dz   = 43.015*mm;  ///
    G4double pole_dphi = 22.5*deg;      ///


  G4ThreeVector pole_xyzCounter;
  //pole_xyzCounter = G4ThreeVector(81*mm, 0, -97.482*mm); ///
  pole_xyzCounter = G4ThreeVector(141.5*mm, 0, -99.215*mm); ///
  G4RotationMatrix* pole_rotCounter = new G4RotationMatrix;


  G4Tubs *pole  = new G4Tubs("pole", pole_rmin, pole_rmax, pole_dz, 0, twopi);
  G4LogicalVolume *pole_log = new G4LogicalVolume(pole, Al, "pole_log", 0, 0, 0);  ///

//for(int i=0; i<8; i++){
  for(int i=0; i<16; i++){

 //pole_rotCounter->rotateZ(pole_dphi*i +6*deg);
 //pole_xyzCounter.rotateZ(pole_dphi*i +6*deg);
   pole_rotCounter->rotateZ(pole_dphi*i);
   pole_xyzCounter.rotateZ(pole_dphi*i );


  G4Transform3D pole_posCounter(*pole_rotCounter, pole_xyzCounter);
  new G4PVPlacement(pole_posCounter, pole_log, "pole_phy", blc_log, false, 0, CHECK_OVERLAPS);
  }
  

  std::cout<<"param[cid][1][4]: "<<param[cid][1][4]<<std::endl;
  std::cout<<"param[cid][0][7]: "<<param[cid][0][7]<<std::endl;
 
  

  if(!WIRE) return;
}

void KnuclDetectorConstruction::CreateChamberCells(int cid,G4LogicalVolume* log, bool WIRE)
{
  if(cid==CID_CDC) CreateCDCCells(log,WIRE);
  else if(cid==CID_BLC1a){
    CreateBLCCells(cid,log,WIRE);
    CreateBLCCells(cid+1,log,WIRE); // CID_BLC1b
  }
  else if(cid==CID_BLC2a){
    CreateBLCCells(cid,log,WIRE);
    CreateBLCCells(cid+1,log,WIRE); // CID_BLC2b
  }
  else if(cid==CID_BPC)   CreateBPCCells(cid,log,WIRE);
  else if(cid==CID_FDC1)  CreateFDCCells(cid,log,WIRE);
  else if(cid==CID_FDC2)  CreateFDCCells(cid,log,WIRE);
  else{
    std::cout<< "[Warning] no definition for cid = "<<cid<<std::endl;
  }
}
//////////////////////////////////////////////////////
G4VPhysicalVolume* KnuclDetectorConstruction::Construct()
//////////////////////////////////////////////////////
{

  // Define magnetic field
  SetMagneticFields(); 

  // ==============================================================
  // Definition of sensitive detectors
  // ==============================================================

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  counterSD   = new KnuclCounterSD  ("counterSD", this, ADCthreshold);
  chamberSD   = new KnuclChamberSD  ("chamberSD", this, param[0]);
  SDman->AddNewDetector(counterSD);
  SDman->AddNewDetector(chamberSD);

  // ==============================================================
  // experimental hall (world volume)
  //   --- beam line along z axis ---
  // ==============================================================

  experimentalHall_log=ConstructShape(CID_Hall);

  // Special geometry for NC efficiency studies
  if(neutron_efficiency_mode) return SetupForNCEfficencyStudy();

  G4LogicalVolume* cds=ConstructShape(CID_Doraemon);
  if(DORAFieldMapFlag == 0){
    G4FieldManager* CDCFieldMgr = fEmFieldSetupCDC->GetFieldManager();
    cds->SetFieldManager(CDCFieldMgr,0);
  }

  G4LogicalVolume* uswk=ConstructShape(CID_USWK);
  if(USWKFieldMapFlag == 0){
    G4FieldManager* UshiwakaFieldMgr = fEmFieldSetupUSHIWAKA->GetFieldManager(); 
    uswk->SetFieldManager(UshiwakaFieldMgr,0);
  }

  ConstructBeamDump();

  if (GetCounterStatus(CID_NC)  ) ConstructShape(CID_NC);
  if (GetCounterStatus(CID_CVC) ) ConstructShape(CID_CVC);
  if (GetCounterStatus(CID_PC)  ) ConstructShape(CID_PC);
  
  if (GetCounterStatus(CID_CDC) ) ConstructShape(CID_CDC);
  if (GetCounterStatus(CID_CDH) ) ConstructShape(CID_CDH);
//if (GetCounterStatus(CID_IH)  ) ConstructShape(CID_IH);   kawasak

  if (GetCounterStatus(CID_BVC) ) ConstructShape(CID_BVC);
  if (GetCounterStatus(CID_FDC1)) ConstructShape(CID_FDC1);
  
  if (GetCounterStatus(CID_BLC2)) ConstructShape(CID_BLC2a);
  if (GetCounterStatus(CID_T0)  ) ConstructShape(CID_T0);
  if (GetCounterStatus(CID_AC)  ) ConstructShape(CID_AC);
  if (GetCounterStatus(CID_BPC) ) ConstructShape(CID_BPC);
  if (GetCounterStatus(CID_BPD) ) ConstructShape(CID_BPD);
  if (GetCounterStatus(CID_DEF) ) ConstructShape(CID_DEF);

  ConstructTarget();

  G4cout << "KnuclDetectorConstruction completed" << G4endl;

  return experimentalHall_phys;
}

void KnuclDetectorConstruction::ConstructBeamDump()
{
  ConstructShape(CID_Floor);
  ConstructShape(CID_BeamDump);
  ConstructShape(CID_SideDump);
  ConstructShape(CID_NShield);
  ConstructShape(CID_SideCon);
  ConstructShape(CID_DoorCon);
}
void KnuclDetectorConstruction::ConstructTarget()
{
  ConstructShape(CID_CDCCFRP);
  ConstructShape(CID_CDCMylar);
  ConstructShape(CID_TarChm);
  ConstructShape(CID_RadS);
  ConstructShape(CID_TarCFRP);
  ConstructShape(CID_TarCap);
  ConstructShape(CID_TarCell);
  ConstructShape(CID_TarRing);
  //  ConstructShape(CID_Target);
  ConstructShape(CID_Fiducial);
  ConstructShape(CID_CellBe);
  ConstructShape(CID_CellAlBe);
  ConstructShape(CID_CellRing);
  ConstructShape(CID_BShield);
  ConstructShape(CID_BFrange);
}
