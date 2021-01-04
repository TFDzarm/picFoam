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

#include "PICParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const std::size_t Foam::PICParcel<ParcelType>::sizeofFields_
(
    sizeof(PICParcel<ParcelType>) - sizeof(ParcelType)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PICParcel<ParcelType>::PICParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    U_(Zero),
    typeId_(-1),
    chargeModifier_(1),
    charge_(0.0),
    nParticle_(1.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            typeId_ = readLabel(is);
            chargeModifier_ = readLabel(is);
            charge_ = readScalar(is);
            nParticle_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&U_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "PICParcel<ParcelType>::PICParcel"
        "(const Cloud<ParcelType>&, Istream&, bool)"
    );
}

//Read fields on simulation start
template<class ParcelType>
void Foam::PICParcel<ParcelType>::readFields(Cloud<PICParcel<ParcelType>>& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, U);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, typeId);

    IOField<label> chargeModifier(c.fieldIOobject("chargeModifier", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, chargeModifier);

    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, nParticle);

    label i = 0;
    forAllIter(typename Cloud<PICParcel<ParcelType>>, c, iter)
    {
        PICParcel<ParcelType>& p = iter();

        p.U_ = U[i];
        p.typeId_ = typeId[i];
        p.chargeModifier_ = chargeModifier[i];

        p.nParticle_ = nParticle[i];

        i++;
    }
}

//Called on writeTime
template<class ParcelType>
void Foam::PICParcel<ParcelType>::writeFields
(
    const Cloud<PICParcel<ParcelType>>& c
)
{
    ParcelType::writeFields(c);

    const PICCloud<PICParcel<ParcelType>>& picCloud(dynamicCast<const PICCloud<PICParcel<ParcelType>>>(c));
    label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);

    IOField<label> chargeModifier(c.fieldIOobject("chargeModifier", IOobject::NO_READ), np);
    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<PICParcel<ParcelType>>, c, iter)
    {
        const PICParcel<ParcelType>& p = iter();

        U[i] = p.U();
        typeId[i] = p.typeId();

        chargeModifier[i] = p.chargeModifier();
        nParticle[i] = p.nParticle();

        i++;
    }

    bool writeFields = np > 0;
    if(picCloud.isInitializing() && !writeFields)//We cannot overwrite existing fields with empty ones
    {
        Info << "WARNING no particles initialized make sure to remove old lagrangian fields" << endl;
    }

    U.write(writeFields);
    typeId.write(writeFields);
    chargeModifier.write(writeFields);
    nParticle.write(writeFields);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//This is used in parallel COM
template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PICParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType& >(p)
            << token::SPACE << p.U()
            << token::SPACE << p.typeId()
            << token::SPACE << p.chargeModifier()
            << token::SPACE << p.charge()
            << token::SPACE << p.nParticle();
    }
    else
    {
        os  << static_cast<const ParcelType& >(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            PICParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const PICParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
