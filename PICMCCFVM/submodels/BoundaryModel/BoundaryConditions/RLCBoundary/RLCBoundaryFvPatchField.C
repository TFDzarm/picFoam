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

#include "RLCBoundaryFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::RLCBoundaryFvPatchField::RLCBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    coeff1_(p.size(), Zero),
    coeff2_(p.size(), Zero),
    Q_(0.0),
    Q1_(0.0),
    Q2_(0.0),
    Q3_(0.0),
    sigma_(0.0)
{}


Foam::RLCBoundaryFvPatchField::RLCBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, dict),
    coeff1_(p.size(), Zero),
    coeff2_(p.size(), Zero),
    Q_(dict.lookupOrDefault<scalar>("Q",0.0)),//Read inital values from the boundary definition
    Q1_(dict.lookupOrDefault<scalar>("Q1",0.0)),
    Q2_(dict.lookupOrDefault<scalar>("Q2",0.0)),
    Q3_(dict.lookupOrDefault<scalar>("Q3",0.0)),
    sigma_(dict.lookupOrDefault<scalar>("sigma",0.0))
{
    this->evaluate();
}


Foam::RLCBoundaryFvPatchField::RLCBoundaryFvPatchField
(
    const RLCBoundaryFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    coeff1_(mapper(ptf.coeff1_)),
    coeff2_(mapper(ptf.coeff2_)),
    Q_(ptf.Q_),
    Q1_(ptf.Q1_),
    Q2_(ptf.Q2_),
    Q3_(ptf.Q3_),
    sigma_(ptf.sigma_)
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

Foam::RLCBoundaryFvPatchField::RLCBoundaryFvPatchField
(
    const RLCBoundaryFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    coeff1_(ptf.coeff1_),
    coeff2_(ptf.coeff2_),
    Q_(ptf.Q_),
    Q1_(ptf.Q1_),
    Q2_(ptf.Q2_),
    Q3_(ptf.Q3_),
    sigma_(ptf.sigma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RLCBoundaryFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
    m(coeff1_, coeff1_);
    m(coeff2_, coeff2_);
}

void Foam::RLCBoundaryFvPatchField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
    const RLCBoundaryFvPatchField& fgptf =
            refCast<const RLCBoundaryFvPatchField>(ptf);

    coeff1_.rmap(fgptf.coeff1_, addr);
    coeff2_.rmap(fgptf.coeff2_, addr);
}

void Foam::RLCBoundaryFvPatchField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<scalar>::operator=((1.0-this->coeff1_/this->patch().deltaCoeffs())*this->patchInternalField());//Update the field values in the patch (coeffs is set by the PICMCC lib)
    fvPatchField<scalar>::evaluate();
}



Foam::tmp<Foam::Field<Foam::scalar>>
Foam::RLCBoundaryFvPatchField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;
    return tmp<Field<scalar>>(new Field<scalar>(this->size(), Zero));
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::RLCBoundaryFvPatchField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;
    return tmp<Field<scalar>>(new Field<scalar>(this->size(), Zero));
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::RLCBoundaryFvPatchField::gradientInternalCoeffs() const
{
    return -this->coeff1_;// - 1/(alpha0*epsilon0*A+dx)
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::RLCBoundaryFvPatchField::gradientBoundaryCoeffs() const
{
    return this->coeff1_*this->coeff2_;// 1/(alpha0*epsilon0*A+dx) * (...)
}

void Foam::RLCBoundaryFvPatchField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);//writes patch and type
    writeEntry(os,"value",*this);
    writeEntry(os,"Q",this->Q_);
    writeEntry(os,"Q1",this->Q1_);
    writeEntry(os,"Q2",this->Q2_);
    writeEntry(os,"Q3",this->Q3_);
    writeEntry(os,"sigma",this->sigma_);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        RLCBoundaryFvPatchField
    )
}

// ************************************************************************* //
