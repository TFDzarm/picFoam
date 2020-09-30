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
    Foam::CrossSectionList

Description
    List of cross section models

SourceFiles
    CrossSectionList.C

\*---------------------------------------------------------------------------*/

#ifndef CrossSectionModelList_H
#define CrossSectionModelList_H

#include "PtrList.H"
#include "CrossSectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CrossSectionList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class CrossSectionList
:
    public PtrList<CrossSectionModel<CloudType,Type>>
{
public:
    // Constructors

        //- Construct null from owner
        CrossSectionList(CloudType& owner);

        //- Construct from dictionary and cloud owner
        CrossSectionList(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CrossSectionList(const CrossSectionList<CloudType,Type>& im);

        //- Construct and return a clone
        virtual autoPtr<CrossSectionList<CloudType,Type>> clone() const
        {
            return autoPtr<CrossSectionList<CloudType,Type>>
            (
                new CrossSectionList<CloudType,Type>(*this)
            );
        }


    //- Destructor
    virtual ~CrossSectionList();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CrossSectionModelList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


