/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "FloatingPotential.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"
#include "fixedGradientFvPatchField.H"
#include "fixedValueFvPatchField.H"
#include "processorPolyPatch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FloatingPotential<CloudType>::FloatingPotential
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud,associatedPatches),
    Qconv_(0.0),
    sigma_(0.0),
    patchId_(associatedPatches[0]),
    patchArea_(),
    boundaryCondition_(nullptr)
{
    volScalarField& phiE = cloud.elpotentialField();

    if(associatedPatches.size() > 1)
    {
        FatalErrorInFunction
                << "model " << typeName << " does not support mutiple patches" << nl
                << abort(FatalError);
    }
    const polyPatch& patch = cloud.mesh().boundaryMesh()[patchId_];

    if(!isA<circuitBoundaryFvPatchField>(phiE.boundaryField()[patchId_]))
        FatalErrorInFunction << "Expected circuitBoundary boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

    boundaryCondition_ = &(refCast<circuitBoundaryFvPatchField>(phiE.boundaryFieldRef()[patchId_]));

    const scalarField magSf(mag(patch.faceAreas()));
    patchArea_ = sum(magSf);
    reduce(patchArea_, sumOp<scalar>());
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FloatingPotential<CloudType>::~FloatingPotential()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FloatingPotential<CloudType>::preUpdate_Boundary()
{
    //Surface update
    reduce(Qconv_, sumOp<scalar>());
    sigma_ += Qconv_/patchArea_;

    boundaryCondition_->circuitGradient() = sigma_/constant::electromagnetic::epsilon0.value();

    Qconv_ = 0.0;
}

template<class CloudType>
void Foam::FloatingPotential<CloudType>::postUpdate_Boundary()
{
    Info  << "[" << typeName << "] sigma: " << sigma_ << " As/m^2" << endl;
}

template<class CloudType>
void Foam::FloatingPotential<CloudType>::injection() {}



template<class CloudType>
bool Foam::FloatingPotential<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    const scalar charge = p.nParticle()*p.charge();

    td.keepParticle = false;
    Qconv_ += charge;
    return true;
}

template<class CloudType>
void Foam::FloatingPotential<CloudType>::particleEjection(typename CloudType::parcelType& p, label patchId)
{
    const scalar Q = p.charge()*p.nParticle();
    Qconv_ -= Q;
}

// ************************************************************************* //
