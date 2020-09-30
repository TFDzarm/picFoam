/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "fvCFD.H"
#include "ElectroStaticSolver.H"
#include "electromagneticConstants.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectroStaticSolver<CloudType>::ElectroStaticSolver
(
    const dictionary& dict,
    CloudType& cloud
)
:
    MaxwellSolver<CloudType>(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectroStaticSolver<CloudType>::~ElectroStaticSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ElectroStaticSolver<CloudType>::solveFields()
{
    CloudType& cloud(this->owner());

    volVectorField& E(cloud.electricField());
    volScalarField& phiE(cloud.elpotentialField());
    volScalarField& rhoCharge(cloud.rhoCharge());

    solve
    (
        fvm::laplacian(phiE) + rhoCharge/constant::electromagnetic::epsilon0
    );
    phiE.correctBoundaryConditions();

    E = -fvc::grad(phiE);
    E.correctBoundaryConditions();
    cloud.eFieldWeighting().update();//FieldWeighting
}


// ************************************************************************* //
