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

#include "circuitBoundaryFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::circuitBoundaryFvPatchField::circuitBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    circuitGradient_(0.0),
    area_(0.0)
{
    const scalarField magSf(mag(patch().Sf()));
    area_ = sum(magSf);
    reduce(area_, sumOp<scalar>());
}


Foam::circuitBoundaryFvPatchField::circuitBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, dict),
    circuitGradient_(0.0),
    area_(0.0)
{
    const scalarField magSf(mag(patch().Sf()));
    area_ = sum(magSf);
    reduce(area_, sumOp<scalar>());

    scalar Q = dict.lookupOrDefault("Q",0.0);
    if(Q != 0.0)
    {
        circuitGradient_ = Q/area_/constant::electromagnetic::epsilon0.value();//Set the inital value read from boundary definition
    }


    this->evaluate();
}


Foam::circuitBoundaryFvPatchField::circuitBoundaryFvPatchField
(
    const circuitBoundaryFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    circuitGradient_(ptf.circuitGradient_),
    area_(ptf.area_)

{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}

Foam::circuitBoundaryFvPatchField::circuitBoundaryFvPatchField
(
    const circuitBoundaryFvPatchField& ptf
)
:
    fvPatchField<scalar>(ptf),
    circuitGradient_(ptf.circuitGradient_),
    area_(ptf.area_)
{}


Foam::circuitBoundaryFvPatchField::circuitBoundaryFvPatchField
(
    const circuitBoundaryFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    circuitGradient_(ptf.circuitGradient_),
    area_(ptf.area_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circuitBoundaryFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}

void Foam::circuitBoundaryFvPatchField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
}

void Foam::circuitBoundaryFvPatchField::evaluate(const Pstream::commsTypes)//Called every timestep
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<scalar>::operator=
    (
         this->patchInternalField() + circuitGradient_/this->patch().deltaCoeffs()//Update the field on the patch
    );

    fvPatchField<scalar>::evaluate();
}



Foam::tmp<Foam::Field<Foam::scalar>>
Foam::circuitBoundaryFvPatchField::valueInternalCoeffs//Used by the solver
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>(new Field<scalar>(this->size(), pTraits<scalar>::one));
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::circuitBoundaryFvPatchField::valueBoundaryCoeffs//Used by the solver
(
    const tmp<scalarField>&
) const
{
    return circuitGradient()/this->patch().deltaCoeffs();
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::circuitBoundaryFvPatchField::gradientInternalCoeffs() const//Used by the solver
{
    return tmp<Field<scalar>>
    (
        new Field<scalar>(this->size(), Zero)
    );
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::circuitBoundaryFvPatchField::gradientBoundaryCoeffs() const//Used by the solver
{
    return tmp<Field<scalar>>
    (
        new Field<scalar>(this->size(), circuitGradient())
    );
}

void Foam::circuitBoundaryFvPatchField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);//writes patch and type
    writeEntry(os,"value",*this);
    writeEntry(os,"Q",circuitGradient_*area_*constant::electromagnetic::epsilon0.value());
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        circuitBoundaryFvPatchField
    )
}

// ************************************************************************* //
