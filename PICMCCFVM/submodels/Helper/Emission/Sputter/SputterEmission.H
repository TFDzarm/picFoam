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
    Foam::SputterEmission

Description
    Sputter from boundary with given probability.
    Like the Sputter BoundaryEvent model.

SourceFiles
    SputterEmission.C

\*---------------------------------------------------------------------------*/

#ifndef SputterEmission_H
#define SputterEmission_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SputterEmission Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SputterEmission: public EmissionModel<CloudType>
{
    //- List of all species and which species(typeIds) they should sputter
    labelListList species_;

    //- Energy of sputterd parcels
    scalarListList energies_;

    //- Probability for sputtering
    scalarListList probability_;
public:

    //- Runtime type information
    TypeName("Sputter");


    SputterEmission
    (
        const dictionary& dict,
        CloudType& cloud
    );

    virtual autoPtr<EmissionModel<CloudType>> clone() const
        {
            return autoPtr<EmissionModel<CloudType>>
            (
                new SputterEmission<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~SputterEmission();


    // Member Functions
    //- Emission of particles from boundary
    virtual void emission();

    //- Collisional emission e.g. sputter
    virtual void collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SputterEmission.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
