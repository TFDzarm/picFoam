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
    Foam::DiagnosticInfo

Description
    Templated DiagnosticInfo class

SourceFiles
    DiagnosticInfo.C
    DiagnosticInfoNew.C

\*---------------------------------------------------------------------------*/

#ifndef DiagnosticInfo_H
#define DiagnosticInfo_H

#include "timeControl.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DiagnosticInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DiagnosticInfo
{
protected:
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Execution control
        Foam::timeControl diagnosticControl_;

        //Foam::timeControl writeControl_;
public:

    //- Runtime type information
    TypeName("DiagnosticInfo")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DiagnosticInfo,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    )


    // Constructors
        //- Construct from components
        DiagnosticInfo
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct and return a clone
        virtual autoPtr<DiagnosticInfo<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DiagnosticInfo();


    //- Selector
    static autoPtr<DiagnosticInfo<CloudType>> New
    (
        const word& info,
        const dictionary& dict,
        CloudType& owner

    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return non-const access to the owner cloud object
        CloudType& owner();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //Called by DiagnosticsList so we loop only once through all particles
        virtual void gatherDiagnostic(const typename CloudType::parcelType& p) = 0;

        //- Should we execute, checks diagnosticControl_
        bool shouldExecute();

        //- Print info
        virtual void info() = 0;

        //bool shouldWrite();
        //virtual void write() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDiagnosticInfo(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(DiagnosticInfo<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        DiagnosticInfo<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeDiagnosticInfoType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    DiagnosticInfo<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagnosticInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
