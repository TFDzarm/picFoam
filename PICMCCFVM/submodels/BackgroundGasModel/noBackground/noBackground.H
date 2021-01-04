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
    Foam::noBackground

Description
    No background, do nothing!

SourceFiles
    noBackground.C

\*---------------------------------------------------------------------------*/

#ifndef noBackground_H
#define noBackground_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BackgroundGasModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noBackground Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class noBackground
:
        public BackgroundGasModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        noBackground
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~noBackground();


    // Member Functions
        virtual bool active() const;

        virtual void update() {}
        virtual void collisionUpdate(CollisionEvent type, label celli, const vector& preU, const vector& postU) {}

        virtual label species() const;
        virtual scalar nParticle() const;

        virtual tmp<scalarField> numberDensity() const;
        virtual tmp<scalarField> temperature() const;

        virtual vector sampleVelocity(label celli);

        virtual void initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "noBackground.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
