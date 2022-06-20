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
    Foam::ChargeInfo

Description
    Calculates the total charge in the system.

SourceFiles
    ChargeInfo.C

\*---------------------------------------------------------------------------*/

#ifndef ChargeInfo_H
#define ChargeInfo_H

#include "DiagnosticInfo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ChargeInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ChargeInfo :
        public DiagnosticInfo<CloudType>
{
public:

    //- Runtime type information
    TypeName("Charge");


    // Constructors

    ChargeInfo(const ChargeInfo<CloudType>& im);

    //- Construct from dictionary
        ChargeInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~ChargeInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new ChargeInfo<CloudType>(*this)
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
    #include "ChargeInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
