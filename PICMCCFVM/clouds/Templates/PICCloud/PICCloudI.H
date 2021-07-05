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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::word& Foam::PICCloud<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::PICCloud<ParcelType>::mesh() const
{
    return mesh_;
}

template<class ParcelType>
inline const Foam::IOdictionary&
Foam::PICCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::PICCloud<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline const Foam::List<Foam::label>&
Foam::PICCloud<ParcelType>::fieldCalculation() const
{
    return fieldCalculation_;
}

template<class ParcelType>
inline const Foam::List<Foam::label>&
Foam::PICCloud<ParcelType>::chargedSpecies() const
{
    return chargedSpecies_;
}

template<class ParcelType>
inline const Foam::List<Foam::label>&
Foam::PICCloud<ParcelType>::neutralSpecies() const
{
    return neutralSpecies_;
}

template<class ParcelType>
inline const Foam::List<Foam::label>&
Foam::PICCloud<ParcelType>::ionSpecies() const
{
    return ionSpecies_;
}

template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::nParticleEqWeight() const
{
    return nParticleEqWeight_;
}

template<class ParcelType>
inline const Foam::List<Foam::List<Foam::DynamicList<ParcelType*>>>&
Foam::PICCloud<ParcelType>::sortedCellOccupancy() const
{
    return sortedCellOccupancy_;
}

template<class ParcelType>
inline Foam::List<Foam::List<Foam::DynamicList<ParcelType*>>>&
Foam::PICCloud<ParcelType>::sortedCellOccupancy()
{
    return sortedCellOccupancy_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::PICCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::PICCloud<ParcelType>::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


template<class ParcelType>
inline Foam::Random& Foam::PICCloud<ParcelType>::rndGen()
{
    return rndGen_;
}

template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::qBF()
{
    return q_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::PICCloud<ParcelType>::fDBF()
{
    return fD_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::rhoNBF()
{
    return rhoN_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::rhoMBF()
{
    return rhoM_.boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::linearKEBF()
{
    return linearKE_.boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::rhoNSpeciesBF(label fieldIndex)
{
    return rhoNSpecies_[fieldIndex].boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::rhoMSpeciesBF(label fieldIndex)
{
    return rhoMSpecies_[fieldIndex].boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::PICCloud<ParcelType>::linearKESpeciesBF(label fieldIndex)
{
    return linearKESpecies_[fieldIndex].boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::PICCloud<ParcelType>::momentumSpeciesBF(label fieldIndex)
{
    return momentumSpecies_[fieldIndex].boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::PICCloud<ParcelType>::momentumBF()
{
    return momentum_.boundaryFieldRef();
}

template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::PICCloud<ParcelType>::jBF()
{
    return j_.boundaryFieldRef();
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::PICCloud<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}

template<class ParcelType>
inline const Foam::ElectronNeutralCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::electronNeutralCollision() const
{
    return electronNeutralCollisionModel_;
}


template<class ParcelType>
inline Foam::ElectronNeutralCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::electronNeutralCollision()
{
    return electronNeutralCollisionModel_();
}

template<class ParcelType>
inline const Foam::IonNeutralCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::ionNeutralCollision() const
{
    return ionNeutralCollisionModel_;
}


template<class ParcelType>
inline Foam::IonNeutralCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::ionNeutralCollision()
{
    return ionNeutralCollisionModel_();
}

template<class ParcelType>
inline Foam::BackgroundGasModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::backgroundGas()
{
    return backgroundGas_();
}

template<class ParcelType>
inline const Foam::BackgroundGasModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::backgroundGas() const
{
    return backgroundGas_;
}


template<class ParcelType>
inline const Foam::WallReflectionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::wallReflection() const
{
    return wallReflectionModel_;
}


template<class ParcelType>
inline Foam::WallReflectionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::wallReflection()
{
    return wallReflectionModel_();
}
template<class ParcelType>
inline const Foam::MaxwellSolver<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::maxwellSolver() const
{
    return maxwellSolver_;
}

template<class ParcelType>
inline Foam::MaxwellSolver<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::maxwellSolver()
{
    return maxwellSolver_();
}

template<class ParcelType>
inline const Foam::ParticlePusher<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::particlePusher() const {
    return particlePusher_;
}

template<class ParcelType>
inline Foam::ParticlePusher<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::particlePusher(){
        return particlePusher_();
}

template<class ParcelType>
inline const Foam::ParticleMerging<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::particleMerging() const {
    return particleMerging_;
}

template<class ParcelType>
inline Foam::ParticleMerging<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::particleMerging(){
        return particleMerging_();
}

template<class ParcelType>
inline const Foam::CoulombCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::coulombCollision() const {
    return coulombCollisionModel_;
}

template<class ParcelType>
inline Foam::CoulombCollisionModel<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::coulombCollision() {
    return coulombCollisionModel_();
}

template<class ParcelType>
inline const Foam::BoundaryEventModelList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::boundaryEventModels() const {
    return boundaryEvents_;
}

template<class ParcelType>
inline Foam::BoundaryEventModelList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::boundaryEventModels() {
    return boundaryEvents_;
}

template<class ParcelType>
inline const Foam::BoundaryModelList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::boundaryModels() const {
    return boundaryModels_;
}

template<class ParcelType>
inline Foam::BoundaryModelList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::boundaryModels() {
    return boundaryModels_;
}

template<class ParcelType>
inline const Foam::DiagnosticsList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::Diagnostics() const {
    return parcelDiagnostics_;
}

template<class ParcelType>
inline Foam::DiagnosticsList<Foam::PICCloud<ParcelType>>&
Foam::PICCloud<ParcelType>::Diagnostics() {
    return parcelDiagnostics_;
}

//----------------------------------------------------------------------//

template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass()*p.nParticle();
    }

    return sysMass;
}


template<class ParcelType>
inline Foam::vector Foam::PICCloud<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += p.nParticle()*cP.mass()*p.U();
    }

    return linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::PICCloud<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*p.nParticle()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::PICCloud<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei()*p.nParticle();
    }

    return internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::PICCloud<ParcelType>::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
    return tfld();
}


template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::PICCloud<ParcelType>::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(3.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline Foam::scalar
Foam::PICCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::PICCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(2.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}

template<class ParcelType>
inline Foam::scalar
Foam::PICCloud<ParcelType>::maxwellianComponent(
        scalar temperature,
        scalar mass
) const
{
    return sqrt(physicoChemical::k.value()*temperature/mass);
}

template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::lorentzFactor( Foam::scalar velocity) const
{
    return 1.0/sqrt(1.0-sqr(velocity/universal::c.value()));
}

template<class ParcelType>
inline Foam::scalar Foam::PICCloud<ParcelType>::lorentzFactor( Foam::scalar momentum, Foam::scalar mass) const
{
    return sqrt(1.0+sqr(momentum/(mass*universal::c.value())));
}

template<class ParcelType>
inline const Foam::volScalarField& Foam::PICCloud<ParcelType>::q() const
{
    return q_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::PICCloud<ParcelType>::fD() const
{
    return fD_;
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::rhoN() const
{
    return rhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::PICCloud<ParcelType>::rhoM() const
{
    return rhoM_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::picRhoN() const
{
    return picRhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::linearKE() const
{
    return linearKE_;
}

template<class ParcelType>
inline const Foam::volVectorField& Foam::PICCloud<ParcelType>::momentum() const
{
    return momentum_;
}

template<class ParcelType>
inline Foam::volScalarField& Foam::PICCloud<ParcelType>::rhoCharge()
{
    return chargeDensity_().field();
}

template<class ParcelType>
inline const Foam::volScalarField& Foam::PICCloud<ParcelType>::rhoCharge() const
{
    return chargeDensity_().field();
}

template<class ParcelType>
inline Foam::volVectorField&
Foam::PICCloud<ParcelType>::electricField()
{
    return eField_().field();
}

template<class ParcelType>
inline const Foam::volVectorField&
Foam::PICCloud<ParcelType>::electricField() const
{
    return eField_().field();
}

template<class ParcelType>
inline Foam::FieldWeigthing& Foam::PICCloud<ParcelType>::eFieldWeighting()
{
    return eField_();
}

template<class ParcelType>
inline const Foam::FieldWeigthing& Foam::PICCloud<ParcelType>::eFieldWeighting() const
{
    return eField_;
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::spaceChargeDensity() const
{
    return spaceChargeDensity_;
}

template<class ParcelType>
inline const bool&
Foam::PICCloud<ParcelType>::syncVelocityAtBoundary() const
{
    return syncVelocityAtBoundary_;
}

template<class ParcelType>
inline const bool&
Foam::PICCloud<ParcelType>::isInitializing() const
{
    return isInitializing_;
}

template<class ParcelType>
inline Foam::dictionary&
Foam::PICCloud<ParcelType>::picInitialiseDict()
{
    return picInitialiseDict_;
}

template<class ParcelType>
inline const Foam::dictionary&
Foam::PICCloud<ParcelType>::picInitialiseDict() const
{
    return picInitialiseDict_;
}

template<class ParcelType>
inline Foam::volScalarField&
Foam::PICCloud<ParcelType>::elpotentialField()
{
    return phiE_;
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::elpotentialField() const
{
    return phiE_;
}

template<class ParcelType>
inline Foam::volVectorField&
Foam::PICCloud<ParcelType>::magneticField()
{
    return B_;
}

template<class ParcelType>
inline const Foam::volVectorField&
Foam::PICCloud<ParcelType>::magneticField() const
{
    return B_;
}

template<class ParcelType>
inline Foam::volScalarField&
Foam::PICCloud<ParcelType>::jouleHeat()
{
    return jE_;
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::PICCloud<ParcelType>::jouleHeat() const
{
    return jE_;
}

template<class ParcelType>
inline const Foam::scalarField&
Foam::PICCloud<ParcelType>::cellLengthScale() const
{
    return cellLengthScale_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::PICCloud<ParcelType>::avgCellLengthScale() const
{
    return avgcellLengthScale_;
}

template<class ParcelType>
inline bool& Foam::PICCloud<ParcelType>::printVelocityWarning()
{
    return printVelocityWarning_;
}

template<class ParcelType>
inline const bool& Foam::PICCloud<ParcelType>::printVelocityWarning() const
{
    return printVelocityWarning_;
}


template<class ParcelType>
inline const Foam::label& Foam::PICCloud<ParcelType>::electronTypeId() const
{
    return electronTypeId_;
}

template<class ParcelType>
inline const Foam::label& Foam::PICCloud<ParcelType>::neutralTypeId(Foam::label i) const
{
    return neutralTypeList_[i];
}

template<class ParcelType>
inline const Foam::label& Foam::PICCloud<ParcelType>::ionTypeId(Foam::label i) const
{
    return ionTypeList_[i];
}

template<class ParcelType>
inline const Foam::List<Foam::label>& Foam::PICCloud<ParcelType>::ionTypeIdList() const
{
    return ionTypeList_;
}

template<class ParcelType>
inline void Foam::PICCloud<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}

template<class ParcelType>
inline const Foam::scalar& Foam::PICCloud<ParcelType>::warnCellTrajectory() const
{
    return warnCellTrajectory_;
}

// ************************************************************************* //
