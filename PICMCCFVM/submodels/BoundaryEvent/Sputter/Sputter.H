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
    Foam::SputterEvent

Description
    Sputter event: Create particles with a given probability when the associated patch is hit.
    Used for e.g. wall boundaries. Special boundary models use the EmissionModel class, but can also use this.

SourceFiles
    Sputter.C

\*---------------------------------------------------------------------------*/

#ifndef SputterEvent_H
#define SputterEvent_H

#include "BoundaryEvent.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SputterEvent Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SputterEvent :
        public BoundaryEvent<CloudType>
{
    //Private data

        //- List for all typeId containing a list of typeIds to sputter
        List<List<List<label>>> sputterSpecies_;

        //- List for all typeId containing a list of inital velocity
        List<List<List<scalar>>> initalEnergies_;

        //- List contatining the probability for the sputtering to occur
        List<List<List<scalar>>> sputterProbability_;

        //- Patch normals
        List<List<vector>> patchNormal_;

        //- Tris that belong to a face on a patch
        List<List<List<triFace>>> patchTriFaces_;

        //- Should we print info?
        bool printInfo_;

        //- Current timestep sputterd parcels
        List<Field<scalar>> sputterdParcels_;

        //- Total sum of sputterd parcels
        List<Field<scalar>> totalSputterdParcels_;

        //- Running average of sputterd parcels
        List<Field<scalar>> runningAverage_;

        //- Time when averaging began
        scalar averageStart_;
public:

    //- Runtime type information
    TypeName("Sputter");


    // Constructors

    SputterEvent(const SputterEvent<CloudType>& im);

    //- Construct from dictionary
        SputterEvent
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~SputterEvent();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new SputterEvent<CloudType>(*this)
        );
    }

    // Member Functions
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Sputter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
