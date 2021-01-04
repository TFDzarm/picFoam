/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "WetzelIonization.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::WetzelIonizationCs<CloudType,Type>::WetzelIonizationCs
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
    scalar tableStepSize = (tableSize-1)/::log(200/15);
    label minEi = wetzel_Energies.size()-wetzel_Qi1.size();

    ionizationCS.append(wetzel_Qi1.first());

    for(label i = 1; i < tableSize-1; i++)
    {
        scalar E = ::exp( scalar(i) / tableStepSize ) * 15;
        //ionEnergies.append(E);

        label index = 0;
        for(label j = 0; j < wetzel_Energies.size();j++)
        {
            if(wetzel_Energies[j] >= E)
                break;
            index++;
        }
        scalar Qi = wetzel_Qi1[(index-minEi)-1] + (E-wetzel_Energies[index-1])*(wetzel_Qi1[(index-minEi)] - wetzel_Qi1[(index-minEi)-1])/(wetzel_Energies[index]-wetzel_Energies[index-1]);
        ionizationCS.append(Qi);
    }
    ionizationCS.append(wetzel_Qi1.last());
    //Save the table
    ionizationCS_.transfer(ionizationCS);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::WetzelIonizationCs<CloudType,Type>::~WetzelIonizationCs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::WetzelIonizationCs<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qi=0.0;
    scalar s = 38.21999573*::log(eVEnergy/15.0);
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
Foam::scalar Foam::WetzelIonizationCs<CloudType,Type>::threshold() const
{
    return 15.0;
}

// ************************************************************************* //
