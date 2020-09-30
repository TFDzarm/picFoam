/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::WallInteractionModel

Description
    Templated wall interaction model class

SourceFiles
    WallInteractionModel.C
    WallInteractionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PrintParcelInfo_H
#define PrintParcelInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class WallInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PrintParcelInfo :
        public DiagnosticInfo<CloudType>
{
    Switch printPosition_;
    Switch printVelocity_;
    Switch printSums_;

    vector sumPosition_;
    vector sumVelocity_;
public:

    //- Runtime type information
    TypeName("ParcelInfo");


    // Constructors

    PrintParcelInfo(const PrintParcelInfo<CloudType>& im);

    //- Construct from dictionary
        PrintParcelInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PrintParcelInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new PrintParcelInfo<CloudType>(*this)
        );
    }
    // Member Functions
    virtual void gatherDiagnostic(const typename CloudType::parcelType& p);

    //- Print info
    virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PrintParcelInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
