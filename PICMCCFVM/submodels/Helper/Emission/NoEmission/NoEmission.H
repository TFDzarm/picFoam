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

Class
    Foam::NoEmission

Description
    No emission.

SourceFiles
    NoEmission.C

\*---------------------------------------------------------------------------*/

#ifndef NoEmission_H
#define NoEmission_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NoEmission Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoEmission: public EmissionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("NoEmission");


    NoEmission
    (
        const dictionary& dict,
        CloudType& cloud
    );

    virtual autoPtr<EmissionModel<CloudType>> clone() const
        {
            return autoPtr<EmissionModel<CloudType>>
            (
                new NoEmission<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~NoEmission();


    // Member Functions

    //- Do nothing
    virtual void emission();

    //- Do nothing
    virtual void collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td){}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoEmission.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
