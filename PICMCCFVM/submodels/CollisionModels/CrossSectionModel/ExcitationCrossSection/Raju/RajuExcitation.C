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

\*---------------------------------------------------------------------------*/

#include "RajuExcitation.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::RajuExcitationCS<CloudType,Type>::RajuExcitationCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId),
    excitationCS_()

{
    Info << "       |= WARNING the CrossSectionModel " << typeName << " is only valid for species of type argon!" << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Total and Elastic CrossSection
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //Create interpolation table from datasets
    label tableSize = 100;
    scalar tableStepSize = (tableSize-1)/::log(1000/0.08);

    DynamicList<scalar> excitationCS;

    excitationCS.append(raju_Qex.first());

    for(label i = 1; i < tableSize-1; i++)
    {
        scalar E = ::exp( scalar(i) / tableStepSize ) * 0.08;

        label index = 0;
        for(label j = 0; j < raju_Energies.size();j++)
        {
            if(raju_Energies[j] >= E)
                break;
            index++;
        }
        scalar Qex = raju_Qex[index-1] + (E-raju_Energies[index-1])*(raju_Qex[index] - raju_Qex[index-1])/(raju_Energies[index]-raju_Energies[index-1]);
        excitationCS.append(Qex);
    }
    //Save interpolation table
    excitationCS_.transfer(excitationCS);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::RajuExcitationCS<CloudType,Type>::~RajuExcitationCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::RajuExcitationCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    if(eVEnergy == 0.0)
        return 0.0;

    scalar Qex=0.0;
    scalar s = 10.49453212*::log(eVEnergy/0.08);
    if( s < 0.0 ) return 0.0;
    if( s < 99.0 )//interpolate
    {
       label l = label(s);
       scalar d = s-scalar(l);
       Qex = (excitationCS_[l+1]-excitationCS_[l])*d + excitationCS_[l];
    }
    else//extrapolate
    {
        Qex = excitationCS_[99];
    }
    if(Qex <= 0.0)
       return 0.0;

    Qex *= 1E-20;
    return Qex;
}

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::RajuExcitationCS<CloudType,Type>::threshold() const
{
    return 12.0;
}

// ************************************************************************* //
