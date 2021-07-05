/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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

Class
    Foam::SmirnovChargeEx

Description
    Simple ChargeExchange cross sections
    Boris M. Smirnov - Theory of Gas Discharge Plasma p.63

\*---------------------------------------------------------------------------*/

#ifndef SmirnovChargeEx_H
#define SmirnovChargeEx_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class SmirnovChargeEx Declaration
\*---------------------------------------------------------------------------*/

//sigma = C ln^2(v0/v)
// => sigma = (sqrt(C)*ln(v0)-sqrt(c)ln(v))^2
// parameter1 = sqrt(C)*ln(v0)
// parameter2 = sqrt(C)
//e.g. Helium in Angstrom^2
// parameter1: 5.5
// parameter2: 0.28
template<class CloudType, crossSectionType Type>
class SmirnovChargeExCS
:
    public CrossSectionModel<CloudType,Type>
{
    //
    scalar coeff1_;
    scalar coeff2_;
public:

    //- Runtime type information
    TypeName("Smirnov");


    // Constructors

        //- Construct from dictionary
        SmirnovChargeExCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~SmirnovChargeExCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new SmirnovChargeExCS<CloudType,Type>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SmirnovChargeEx.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
