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
    Foam::uniformBackground

Description
    Uniform background model.
    Density and Temperature are constant in the whole domain.
    Velocities are sampled from the Maxwell-Boltzmann distribution

SourceFiles
    uniformBackground.C

\*---------------------------------------------------------------------------*/

#ifndef uniformBackground_H
#define uniformBackground_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BackgroundGasModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class uniformBackground Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class uniformBackground
:
        public BackgroundGasModel<CloudType>
{
    // Private data

        //Background species
        label species_;

        //Uniform background fields
        scalarField numberDensity_;
        scalarField temperature_;

        //Equivalent particle
        scalar nParticle_;

        //Dirft velocity
        vector velocity_;
public:

    //- Runtime type information
    TypeName("uniformBackground");

    // Constructors

        //- Construct from components
        uniformBackground
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~uniformBackground();


    // Member Functions
        virtual bool active() const;

        virtual void update() {}
        virtual void collisionUpdate(CollisionEvent type, label celli, const vector& preU, const vector& postU) {}

        virtual label species() const;
        virtual scalar nParticle() const;

        virtual tmp<scalarField> numberDensity() const ;
        virtual tmp<scalarField> temperature() const;

        virtual vector sampleVelocity(label celli);

        virtual void initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "uniformBackground.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
