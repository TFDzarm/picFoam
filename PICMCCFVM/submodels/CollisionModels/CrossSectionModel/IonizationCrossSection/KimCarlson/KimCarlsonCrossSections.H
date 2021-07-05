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
    Foam::KimCarlsonCrossSection

Description
    Use the technique described in

    F. Pérez et al.
    Improved modeling of relativistic collisions and collisional ionization in particle-in-cell codes
    https://doi.org/10.1063/1.4742167

    to calculate a mean cross section from the binding energies of the species.

\*---------------------------------------------------------------------------*/

#ifndef KimCarlsonSections_H
#define KimCarlsonSections_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class HardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class KimCarlsonCrossSection
:
    public CrossSectionModel<CloudType,Type>
{
//Private data

    //- Calculate the binding energy
    scalar bindingEnergyCarlson(label Zstar, label k);

    //- Atomic number of the species
    label atomicNumber_;

    //- Binding energies provided in the dictionary
    scalarList bindingEneries_;

    //- Calculated ionization energies used for interpolation
    scalarList ionizationEneries_;

    //- Calculated ionization cross section used for interpolation
    scalarList ionizationk0CrossSection_;

public:

    //- Runtime type information
    TypeName("KimCarlson");


    // Constructors

        //- Construct from dictionary
        KimCarlsonCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~KimCarlsonCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new KimCarlsonCrossSection<CloudType,Type>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KimCarlsonCrossSections.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
