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

#include "StraubIonization.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::StraubIonizationCS<CloudType,Type>::StraubIonizationCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId),
    ionizationCS_()
{
    Info << "WARNING the CrossSectionModel " << typeName << " is only valid for species of type argon!" << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Ionization CrossSection
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //Create a lookup table for interpolation
    //DynamicList<scalar> ionEnergies;
    DynamicList<scalar> ionizationCS;

    label tableSize = 100;
    scalar tableStepSize = (tableSize-1)/::log(1000/17);
    label minEi = straub_Energies.size()-straub_Qi1.size();

    ionizationCS.append(straub_Qi1.first());

    for(label i = 1; i < tableSize-1; i++)
    {
        scalar E = ::exp( scalar(i) / tableStepSize ) * 17;

        label index = 0;
        for(label j = 0; j < straub_Energies.size();j++)
        {
            if(straub_Energies[j] >= E)
                break;
            index++;
        }
        scalar Qi = straub_Qi1[(index-minEi)-1] + (E-straub_Energies[index-1])*(straub_Qi1[(index-minEi)] - straub_Qi1[(index-minEi)-1])/(straub_Energies[index]-straub_Energies[index-1]);
        ionizationCS.append(Qi);
    }
    ionizationCS.append(straub_Qi1.last());
    //Save the table
    ionizationCS_.transfer(ionizationCS);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::StraubIonizationCS<CloudType,Type>::~StraubIonizationCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::StraubIonizationCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qi=0.0;
    scalar s = 24.29720974*::log(eVEnergy/17);
    if( s < 0.0 ) return 0.0;
    if( s < 99.0 )//interpolate
    {
       label l = label(s);
       scalar d = s-scalar(l);
       Qi = (ionizationCS_[l+1]-ionizationCS_[l])*d + ionizationCS_[l];
    }
    else
    {
        //scalar d = s - 99.0;
        //Qi = (ionizationCS_[99]-ionizationCS_[98])*d + ionizationCS_[99];
        Qi = ionizationCS_[99];//return the last entry
    }
    if(Qi <= 0.0)
       return 0.0;

    Qi *= 1E-20;
    return Qi;
}

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::StraubIonizationCS<CloudType,Type>::threshold() const
{
    return 17.0;
}

// ************************************************************************* //
