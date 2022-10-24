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
/// \file BoxesParameterisation.cc
/// \brief Implementation of the BoxesParameterisation class

#include "BoxesParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BoxesParameterisation::BoxesParameterisation(
        G4int    ObjNumber,
        G4double YcenterPos[],
        G4double ZcenterPos[],          //  Z of center of first
        G4double XWidth,
        G4double YWidth[],
        G4double Height[])
 : G4VPVParameterisation()
{

    fObjNumber= ObjNumber;

    fYcenterPos= YcenterPos;
    fZcenterPos= ZcenterPos;          //  Z of center of first
    fXWidth=XWidth;
    fYWidth=YWidth;
    fHeight = Height;        //  Z spacing of centers


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BoxesParameterisation::~BoxesParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BoxesParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // Note: copyNo will start with zero!

  G4double Yposition = fYcenterPos[copyNo];
  G4double Zposition = fZcenterPos[copyNo];
  G4ThreeVector origin(Yposition,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BoxesParameterisation::ComputeDimensions
(G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
{
  // Note: copyNo will start with zero!
  //G4double rmax = fRmaxFirst + copyNo * fRmaxIncr;

//  G4double ZcenterPos,          //  Z of center of first
//  G4double Height,        //  Z spacing of centers
//  G4double Rmax

   G4double fXHWidth=0.5*fXWidth;
   G4double fYHWidth=0.5*fYWidth[copyNo];
   G4double fHHeight=0.5*fHeight[copyNo];
  trackerChamber.SetYHalfLength(fXHWidth);
  trackerChamber.SetXHalfLength(fYHWidth);
  trackerChamber.SetZHalfLength(fHHeight);

//  trackerChamber.SetZHalfLength(0.5*fHeight[copyNo]);
//  trackerChamber.SetStartPhiAngle(0.*deg);
//  trackerChamber.SetDeltaPhiAngle(360.*deg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
