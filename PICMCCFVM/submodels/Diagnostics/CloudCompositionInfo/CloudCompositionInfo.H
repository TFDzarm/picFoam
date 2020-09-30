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
    Foam::CloudCompositionInfo

Description
    Print statistics on cloud composition.

SourceFiles
    CloudCompositionInfo.C

\*---------------------------------------------------------------------------*/

#ifndef CloudCompositionInfo_H
#define CloudCompositionInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CloudCompositionInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudCompositionInfo :
        public DiagnosticInfo<CloudType>
{
//Private data

    //- Current number of parcels for all species
    List<label> nParcelsTypes_;

    //- Current number of particles for all species
    List<scalar> nParticleTypes_;

    //- Total volume of the domain
    scalar totalVolume_;

    //- Should we also print the number density?
    Switch printNumberDensity_;
public:

    //- Runtime type information
    TypeName("Composition");


    // Constructors

    CloudCompositionInfo(const CloudCompositionInfo<CloudType>& im);

    //- Construct from dictionary
        CloudCompositionInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~CloudCompositionInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new CloudCompositionInfo<CloudType>(*this)
        );
    }

    // Member Functions

        //Called by DiagnosticsList so we loop only once through all particles
        virtual void gatherDiagnostic(const typename CloudType::parcelType& p);

        //- Print info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudCompositionInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
