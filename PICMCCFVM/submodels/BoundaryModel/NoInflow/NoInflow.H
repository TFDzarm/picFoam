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
    Foam::NoInflow

Description
    Not inserting any particles

    Modified OpenFOAM DSMC implementaion to fit in this framework.
\*---------------------------------------------------------------------------*/

#ifndef NoInflow_H
#define NoInflow_H

#include "BoundaryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class NoInflow Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoInflow
:
    public BoundaryModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoInflow
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );

        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new NoInflow<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~NoInflow();


    // Member Functions

        //- Introduce particles (none in this case)
        virtual void injection();

        virtual bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoInflow.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
