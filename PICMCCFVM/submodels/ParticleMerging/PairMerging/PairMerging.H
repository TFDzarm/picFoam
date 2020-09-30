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
    Foam::PairMerging

Description
    Merge two particle time that have smiliar velocities.

    Model inspired by:

    M. M. Hopkins et al.
    Challenges to Simulating Vacuum Arc Discharge
    Report: SAND2013-5329C
            478678
    https://www.osti.gov/biblio/1115246
    Slide 15-16

\*---------------------------------------------------------------------------*/

#ifndef PairMerging_H
#define PairMerging_H

#include "ParticleMerging.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class NoInflow Declaration
\*---------------------------------------------------------------------------*/

//Inspired by Challenges to Simulating Vacuum Arc Discharge (M.Hopkins 2013)
template<class CloudType>
class PairMerging
:
    public ParticleMerging<CloudType>
{
    //- Try to reduce down to ParcelNumber*reductionFactor_
    scalar reductionFactor_;

    //- Maximum allowed numbers of parcel before merging kicks in
    label cellParcelMax_;

    //- Maximum number of checks
    label maxChecks_;

    //- Number of neighbors in velocity sorted list
    label nCheckNeighbors_;

    //- Only merge if below this angle difference
    scalar angleDifference_;

    //- Only merge particles with smiliar velocity (maximum deviation given by this ratio)
    scalar velocityRatio_;

    //- Save error that is made
    scalar cumulativeKinEnergyError_;

    //- Species to be merged
    List<label> mergeSpecies_;

    //- Maximum nParticles deviation
    List<scalar> nParticleDev_;
public:

    //- Runtime type information
    TypeName("PairMerging");


    // Constructors

        //- Construct from dictionary
        PairMerging
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PairMerging();


    // Member Functions
        //- Move the parcel
        virtual bool active() const {
            return true;
        }

        //- Are conditions for merging met?
        virtual bool shouldMerge() const;

        //- Merge particle
        virtual void mergeParticles();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PairMerging.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
