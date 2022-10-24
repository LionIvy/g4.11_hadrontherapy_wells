//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file TubesParameterisation.cc
/// \brief Implementation of the TubesParameterisation class

#include "TubesParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TubesParameterisation::TubesParameterisation(
        G4int    StepNumber,
        G4double ZcenterPos[],          //  Z of center of first
        G4double Height[],        //  Z spacing of centers
        G4double Rmin[],
        G4double Rmax[])
 : G4VPVParameterisation()
{

    fStepNumber= StepNumber;
    fZcenterPos= ZcenterPos;          //  Z of center of first
    fHeight = Height;        //  Z spacing of centers
    fRmin = Rmin;
    fRmax = Rmax;
    //G4double* N=Rmax;
//   fNoChambers =  noChambers;
//   fStartZ     =  startZ;
//   fHalfWidth  =  0.5*widthChamber;
//   fSpacing    =  spacingZ;
//   fRmaxFirst = 0.5 * lengthInitial;
//   if( noChambers > 0 ){
//      fRmaxIncr =  0.5 * (lengthFinal-lengthInitial)/(noChambers-1);
//      if (spacingZ < widthChamber) {
//         G4Exception("TubesParameterisation::TubesParameterisation()",
//                     "InvalidSetup", FatalException,
//                     "Width>Spacing");
//      }
//   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TubesParameterisation::~TubesParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TubesParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // Note: copyNo will start with zero!
  G4double Zposition = fZcenterPos[copyNo];
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TubesParameterisation::ComputeDimensions
(G4Tubs& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
{
  // Note: copyNo will start with zero!
  //G4double rmax = fRmaxFirst + copyNo * fRmaxIncr;
  trackerChamber.SetInnerRadius(fRmin[copyNo]);
  trackerChamber.SetOuterRadius(fRmax[copyNo]);
//  G4double ZcenterPos,          //  Z of center of first
//  G4double Height,        //  Z spacing of centers
//  G4double Rmax
  trackerChamber.SetZHalfLength(0.5*fHeight[copyNo]);
  trackerChamber.SetStartPhiAngle(0.*deg);
  trackerChamber.SetDeltaPhiAngle(360.*deg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
