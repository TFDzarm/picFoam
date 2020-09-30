/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 picFoam
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "picParcel.H"
#include "PICCloud.H"

#include "VariableHardSphere.H"
#include "HardSphere.H"

#include "IsotropicScattering.H"
#include "NoBinaryCollision.H"

#include "uniformBackground.H"
#include "noBackground.H"

#include "NanbuCorrection.H"
#include "SentokuKempCorrection.H"
#include "SentokuKempNonRelCorrection.H"
#include "NoCorrection.H"

#include "picParcel.H"
#include "PICCloud.H"
#include "NoCoulombCollision.H"
#include "PerezRelativisticCollision.H"
#include "NanbuCollision.H"
#include "TakizukaAbePairing.H"
#include "NanbuYonemuraPairing.H"

#include "NoElectronNeutralCollision.H"
#include "NanbuElectronNeutralCollision.H"
#include "RelativisticElectronNeutralCollision.H"

#include "RajuElastic.H"
#include "RajuExcitation.H"
#include "RajuIonization.H"
#include "KimCarlsonCrossSections.H"
#include "StraubIonization.H"
#include "WetzelIonization.H"
#include "BerkeleyArgonElastic.H"
#include "BerkeleyArgonExcitation.H"
#include "BerkeleyArgonIonization.H"
#include "BerkeleyHeliumElastic.H"
#include "BerkeleyHeliumExcitation.H"
#include "BerkeleyHeliumIonization.H"
#include "NoCrossSection.H"
#include "BrusaCrossSection.H"
#include "FixedValueCrossSection.H"

#include "NoIonization.H"
#include "PerezIonizationModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef PICCloud<picParcel> CloudType;

    //Background gas model
    makeBackgroundGasModel(CloudType);
    makeBackgroundGasModelType(noBackground,CloudType);
    makeBackgroundGasModelType(uniformBackground,CloudType);

    //Binary collisions (neutral and neutral ion)
    makeBinaryCollisionModel(CloudType);
    makeBinaryCollisionModelType(NoBinaryCollision, CloudType);
    makeBinaryCollisionModelType(IsotropicScattering, CloudType);

    // Cross section models
    makeTotalCrossSectionModel(CloudType);
    makeTotalCrossSectionModelType(HardSphere, CloudType);
    makeTotalCrossSectionModelType(VariableHardSphere, CloudType);


    //Pairing algorithm for coulomb collisons
    makePairingAlgorithm(PICCloud<picParcel>);
    makePairingAlgorithmType(TakizukaAbePairing, CloudType);
    makePairingAlgorithmType(NanbuYonemuraPairing, CloudType);

    //Coulomb collision models
    makeCoulombCollisionModel(PICCloud<picParcel>);
    makeCoulombCollisionModelType(NoCoulombCollision, CloudType);
    makeCoulombCollisionModelType(PerezRelativisticCollision, CloudType);
    makeCoulombCollisionModelType(NanbuCollision, CloudType);

    //Coulomb collision ionization model
    makeIonizationModel(PICCloud<picParcel>)
    makeIonizationModelType(NoIonization, CloudType)
    makeIonizationModelType(PerezIonizationModel, CloudType)


    //Correction models for unequal weights
    makeWeightCorrectionModel(CloudType);
    makeWeightCorrectionModelType(NanbuWeightCorrectionModel,CloudType);
    makeWeightCorrectionModelType(SentokuKempCorrection,CloudType);
    makeWeightCorrectionModelType(SentokuKempNonRelCorrection,CloudType);
    makeWeightCorrectionModelType(NoWeightCorrectionModel,CloudType);


    //Electron neutral collison models
    makeElectronNeutralCollisionModel(PICCloud<picParcel>);
    makeElectronNeutralCollisionModelType(NoElectronNeutralCollision, CloudType);
    makeElectronNeutralCollisionModelType(NanbuElectronNeutralCollision, CloudType);
    makeElectronNeutralCollisionModelType(RelativisticElectronNeutralCollision, CloudType);

    //CrossSection Models for electron neutral collision
    typedef CrossSectionModel<CloudType, crossSectionType::ElectronElasticCS> elasticCS;
    makeCSModel(elasticCS);

    typedef CrossSectionModel<CloudType, crossSectionType::ElectronExciationCS> exciationCS;
    makeCSModel(exciationCS);

    typedef CrossSectionModel<CloudType, crossSectionType::ElectronIonizationCS> ionizationCS;
    makeCSModel(ionizationCS);

    //Raju Argon Elastic, Excitation, Ionization
    makeCSType(RajuElasticCS,CloudType, crossSectionType::ElectronElasticCS);
    makeCSType(RajuExcitationCS,CloudType, crossSectionType::ElectronExciationCS);
    makeCSType(RajuIonizationCS,CloudType, crossSectionType::ElectronIonizationCS);

    //Berkeley model of XPDP1 https://ptsg.egr.msu.edu/
    //Argon
    makeCSType(BerkeleyArgonElasticCS,CloudType, crossSectionType::ElectronElasticCS);
    makeCSType(BerkeleyArgonExcitationCS,CloudType, crossSectionType::ElectronExciationCS);
    makeCSType(BerkeleyArgonIonizationCS,CloudType, crossSectionType::ElectronIonizationCS);

    //Helium
    makeCSType(BerkeleyHeliumElasticCS,CloudType, crossSectionType::ElectronElasticCS);
    makeCSType(BerkeleyHeliumExcitationCS,CloudType, crossSectionType::ElectronExciationCS);
    makeCSType(BerkeleyHeliumIonizationCS,CloudType, crossSectionType::ElectronIonizationCS);

    //Fixed cross section model
    makeCSTypeMulti(FixedValueCrossSection, elasticCS, CloudType, crossSectionType::ElectronElasticCS);
    makeCSTypeMulti(FixedValueCrossSection, exciationCS, CloudType, crossSectionType::ElectronExciationCS);
    makeCSTypeMulti(FixedValueCrossSection, ionizationCS, CloudType, crossSectionType::ElectronIonizationCS);

    //Argon ionization
    makeCSType(KimCarlsonCrossSection,CloudType, crossSectionType::ElectronIonizationCS);
    makeCSType(StraubIonizationCS,CloudType, crossSectionType::ElectronIonizationCS);
    makeCSType(WetzelIonizationCs,CloudType, crossSectionType::ElectronIonizationCS);

    //Brusa Ne,Ar,Kr,Xe
    makeCSTypeMulti(BrusaCrossSection, elasticCS, CloudType, crossSectionType::ElectronElasticCS);
    makeCSTypeMulti(BrusaCrossSection, exciationCS, CloudType, crossSectionType::ElectronExciationCS);
    makeCSTypeMulti(BrusaCrossSection, ionizationCS, CloudType, crossSectionType::ElectronIonizationCS);

    //Zero cross section
    makeCSTypeMulti(NoCrossSection,elasticCS,CloudType, crossSectionType::ElectronElasticCS);
    makeCSTypeMulti(NoCrossSection,exciationCS,CloudType, crossSectionType::ElectronExciationCS);
    makeCSTypeMulti(NoCrossSection,ionizationCS,CloudType, crossSectionType::ElectronIonizationCS);

}


// ************************************************************************* //
