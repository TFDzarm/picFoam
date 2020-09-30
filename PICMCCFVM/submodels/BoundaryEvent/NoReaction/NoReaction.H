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
    Foam::NoReaction

Description
    Empty model does nothing.

SourceFiles
    NoReaction.C

\*---------------------------------------------------------------------------*/

#ifndef NoReaction_H
#define NoReaction_H

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NoReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoReaction :
        public BoundaryEvent<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    NoReaction(const NoReaction<CloudType>& im);

    //- Construct from dictionary
        NoReaction
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~NoReaction();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new NoReaction<CloudType>(*this)
        );
    }

    // Member Functions

        //- Do nothing
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
