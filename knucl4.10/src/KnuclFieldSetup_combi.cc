#include "KnuclFieldSetup_combi.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4SystemOfUnits.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

KnuclFieldSetup_combi::  KnuclFieldSetup_combi(const char* uswk_file, 
	      double uswk_xOffset, double uswk_yOffset, double uswk_zOffset,
	      const char* dora_file, 
	      double dora_xOffset, double dora_yOffset, double dora_zOffset,
	      double dora_norm)
  : fChordFinder(0), fEquation(0), fMagneticField(0), fStepper(0)
{    
  fMinStep     = 0.25*mm ;
  fMagneticField = new combi_field(uswk_file, 
				   uswk_xOffset, 
				   uswk_yOffset, 
				   uswk_zOffset,
				   dora_file, 
				   dora_xOffset, 
				   dora_yOffset, 
				   dora_zOffset,
				   dora_norm); 
  
  G4FieldManager   *pFieldMgr;
  pFieldMgr=G4TransportationManager::GetTransportationManager()->GetFieldManager();
  pFieldMgr->SetDeltaOneStep(fMinStep);
  G4ChordFinder *pChordFinder = new G4ChordFinder(fMagneticField);
  pFieldMgr->SetChordFinder( pChordFinder );
  pFieldMgr->SetDetectorField(fMagneticField);

}

KnuclFieldSetup_combi::~KnuclFieldSetup_combi()
{
  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
  if(fEquation)      delete fEquation;
  G4cout<<"KnuclFieldSetup_combi::~KnuclFieldSetup_combi()"<<G4endl;
  //  getchar();
}

combi_field::combi_field(const char* uswk_file, 
				   double uswk_xOffset, 
				   double uswk_yOffset, 
				   double uswk_zOffset,
				   const char* dora_file, 
				   double dora_xOffset, 
				   double dora_yOffset, 
				   double dora_zOffset,
				   double dora_norm)
{
  puswk = new sub_field_uswk(uswk_file, uswk_xOffset, uswk_yOffset, uswk_zOffset);
  pdora = new sub_field_dora(dora_file, dora_xOffset, dora_yOffset, dora_zOffset, dora_norm);

  {
    G4cout << "sample value from combi field map: z=250cm " <<G4endl;
    const double test_point[4]={0*cm, 0*cm, 250*cm, 0};
    double test_field[3]={-1, -1, -1};
    GetFieldValue(test_point, test_field);
    G4cout << test_field[0]/gauss << " " << test_field[1]/gauss << " " << test_field[2]/gauss << " gauss "
	   << "\n-----------------------------------------------------------" << endl;
  }

  {
    G4cout << "sample value from combi field map: z=-50cm " <<G4endl;
    const double test_point[4]={0*cm, 0*cm, -50*cm, 0};
    double test_field[3]={-1, -1, -1};
    GetFieldValue(test_point, test_field);
    G4cout << test_field[0]/gauss << " " << test_field[1]/gauss << " " << test_field[2]/gauss << " gauss "
	   << "\n-----------------------------------------------------------" << endl;
  }

  {
    G4cout << "sample value from combi field map: z=40cm " <<G4endl;
    const double test_point[4]={0*cm, 0*cm, 40*cm, 0};
    double test_field[3]={-1, -1, -1};
    GetFieldValue(test_point, test_field);
    G4cout << test_field[0]/gauss << " " << test_field[1]/gauss << " " << test_field[2]/gauss << " gauss "
	   << "\n-----------------------------------------------------------" << endl;
  }


}

void combi_field::GetFieldValue(const double Point[4], double *Bfield) const
{

  if(Point[2] > 100.*cm)
    puswk->GetFieldValue(Point, Bfield);
  else
    pdora->GetFieldValue(Point, Bfield);  

}

sub_field_uswk::sub_field_uswk( const char* filename, double xOffset, double yOffset, double zOffset) 
  :fXoffset(xOffset),fYoffset(yOffset),fZoffset(zOffset),invertX(false),invertY(false),invertZ(false)
{    
 
  double lenUnit= cm;
  double fieldUnit= gauss; 
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 
  ifstream file( filename ); // Open the file for reading.
  
  //char buffer[256];
  // Ignore first blank line
  //  file.getline(buffer,256);
  
  // Read table dimensions 
  file >> nx >> ny >> nz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: " 
	 << nx << " " << ny << " " << nz << " ] "
	 << endl;

  // Set up storage space for table
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }
  
  // Ignore other header information    
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  //  do {
  //    file.getline(buffer,256);
  //  } while ( buffer[1]!='0');
  
  // Read in the data
  G4cout << " Read values of field from file " << filename << endl; 
  //  double xval,yval,zval,bx,by,bz;
  float xval,yval,zval,bx,by,bz;
  //  double permeability; // Not used in this example.
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
	//        file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;

        file >> xval >> yval >> zval >> bx >> by >> bz;

	// G4cout << " xval = " << xval << " yval = " << yval << " zval = " << zval << endl;
	// G4cout << " Bx   = " << bx << " By = " << by << " Bz = " << bz << endl;
	// if(xval==43 && yval==23 && zval==98)
	//   getchar();

        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
      }
    }
  }
  file.close();

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by X " << xOffset/cm << " cm " << endl
	 << "\n ---> The field will be offset by Y " << yOffset/cm << " cm " << endl
	 << "\n ---> The field will be offset by Z " << zOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm "
	 << "\n-----------------------------------------------------------" << endl;

  const double test_point[4]={1*cm, 0*cm, 250*cm, 0};
  double test_field[3]={-1, -1, -1};
  GetFieldValue(test_point, test_field);
  G4cout << "sample value from field map " 
    	 << test_field[0]/gauss << " " << test_field[1]/gauss << " " << test_field[2]/gauss << " gauss "
	 << "\n-----------------------------------------------------------" << endl;

}

//#define DEBUG_INTERPOLATING_FIELD 1

void sub_field_uswk::GetFieldValue(const double point[4],
				      double *Bfield ) const
{

  double x = point[0] - fXoffset;
  double y = point[1] - fYoffset;
  double z = point[2] - fZoffset;

  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    double valx0z0, mulx0z0, valx1z0, mulx1z0;
    double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}

sub_field_dora::sub_field_dora( const char* filename, double xOffset, double yOffset, double zOffset, double norm) 
  :fXoffset(xOffset),fYoffset(yOffset),fZoffset(zOffset),invertX(false),invertY(false),invertZ(false), norm_factor(1.0)
{    
 
  double lenUnit= cm;
  double fieldUnit= gauss; 
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 
  ifstream file( filename ); // Open the file for reading.
  
  //char buffer[256];
  // Ignore first blank line
  //  file.getline(buffer,256);
  
  // Read table dimensions 
  file >> nx >> ny >> nz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: " 
	 << nx << " " << ny << " " << nz << " ] "
	 << endl;

  // Set up storage space for table
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }
  
  // Ignore other header information    
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  //  do {
  //    file.getline(buffer,256);
  //  } while ( buffer[1]!='0');
  
  // Read in the data
  G4cout << " Read values of field from file " << filename << endl; 
  //  double xval,yval,zval,bx,by,bz;
  float xval,yval,zval,bx,by,bz;
  //  double permeability; // Not used in this example.
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
	//        file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;

        file >> xval >> yval >> zval >> bx >> by >> bz;

	// G4cout << " xval = " << xval << " yval = " << yval << " zval = " << zval << endl;
	// G4cout << " Bx   = " << bx << " By = " << by << " Bz = " << bz << endl;
	// if(xval==43 && yval==23 && zval==98)
	//   getchar();

        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
      }
    }
  }
  file.close();

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by X " << xOffset/cm << " cm " << endl
	 << "\n ---> The field will be offset by Y " << yOffset/cm << " cm " << endl
	 << "\n ---> The field will be offset by Z " << zOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm "
	 << "\n-----------------------------------------------------------" << endl;

  const double test_point[4]={0.*cm, 0.*cm, 0.*cm, 0.};
  double test_field[3]={-1., -1., -1.};
  GetFieldValue(test_point, test_field);
  G4cout << "calculate normalization factor " <<G4endl;
  G4cout << "before normalization " <<G4endl;
  G4cout << test_field[0]/tesla << " " << test_field[1]/tesla << " " << test_field[2]/tesla << " gauss "<< endl;

  norm_factor = norm/(test_field[2]/tesla);
  G4cout << "norm_factor = " << norm_factor <<G4endl;

  G4cout << "after normalization " <<G4endl;
  GetFieldValue(test_point, test_field);
  G4cout << test_field[0]/tesla << " " << test_field[1]/tesla << " " << test_field[2]/tesla << " tesla "<< endl;

}

//#define DEBUG_INTERPOLATING_FIELD 1

void sub_field_dora::GetFieldValue(const double point[4],
				      double *Bfield ) const
{

  double x = point[0] - fXoffset;
  double y = point[1] - fYoffset;
  double z = point[2] - fZoffset;

  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    double valx0z0, mulx0z0, valx1z0, mulx1z0;
    double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

    Bfield[0] = norm_factor * Bfield[0];
    Bfield[1] = norm_factor * Bfield[1];
    Bfield[2] = norm_factor * Bfield[2];

  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}

