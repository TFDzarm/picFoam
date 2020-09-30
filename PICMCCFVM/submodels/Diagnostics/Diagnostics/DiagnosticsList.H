/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::DiagnosticsList

Description
    List of DiagnosticInfo models

SourceFiles
    DiagnosticsList.C

\*---------------------------------------------------------------------------*/

#ifndef DiagnosticsList_H
#define DiagnosticsList_H

#include "PtrList.H"
#include "DiagnosticInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DiagnosticsList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DiagnosticsList
:
    public PtrList<DiagnosticInfo<CloudType>>
{
protected:
    // reference to the owner cloud class
    CloudType& owner_;
public:

    // Constructors

        //- Construct null from owner
        DiagnosticsList(CloudType& owner);

        //- Construct from dictionary and cloud owner
        DiagnosticsList(const dictionary& dict, CloudType& owner);
        void setupModels(const dictionary& dict, CloudType& owner);

        //- Construct copy
        DiagnosticsList(const DiagnosticsList<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<DiagnosticsList<CloudType>> clone() const
        {
            return autoPtr<DiagnosticsList<CloudType>>
            (
                new DiagnosticsList<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DiagnosticsList();

    //- Return the owner cloud object
    const CloudType& owner() const;

    //- Return non-const access to the owner cloud object
    CloudType& owner();



    // Member Functions

        //- Loop through all particles, call gatherDiagnostic() and info() on all models
        void info();

   //void write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagnosticsList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


