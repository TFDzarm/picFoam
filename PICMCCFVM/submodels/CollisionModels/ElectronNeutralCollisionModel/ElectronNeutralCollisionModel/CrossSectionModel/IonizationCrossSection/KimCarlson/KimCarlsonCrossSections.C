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

#include "KimCarlsonCrossSections.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::KimCarlsonCrossSection<CloudType,Type>::KimCarlsonCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(dict, cloud, typeName ,associatedTypeId),
    atomicNumber_(readLabel(this->coeffDict().lookup("atomicNumber"))),
    bindingEneries_(this->coeffDict().lookup("bindingEnergy")),
    ionizationEneries_(this->coeffDict().lookup("ionizationEnergy")),
    ionizationk0CrossSection_()
{
    if( bindingEneries_.size() != atomicNumber_ )
        FatalErrorInFunction << "Expected more data points in list bindingEnergy"  << abort(FatalError);
    if(ionizationEneries_.size() != atomicNumber_)
        FatalErrorInFunction << "Expected more data points in list ionizationEneries"  << abort(FatalError);

    namespace cm = constant::mathematical;

    scalar r_e = 2.8179403227e-15;
    label tableSize = 100;
    //1eV-10MeV
    scalar tableStepSize = (tableSize-1)/::log(1e7/1);

    scalar E, Bk, eps, eps_dash, u_dash, b_dash,beta_eps2,beta_b2,beta_u2,sigma0,ediv2,A1,A2,A3,sigmak;
    label N;
    for(label Zstar = 0; Zstar < 1; Zstar++)//Only calculate first ionization cross section
    {
        ionizationk0CrossSection_ = List<scalar>(tableSize,0.0);

        for(label i = 0; i < tableSize; i++)
        {
            E = ::exp( scalar(i)/ tableStepSize);
            N = 1;
            for(label k = 0; k < atomicNumber_-Zstar; k++)
            {
                Bk = bindingEnergyCarlson(Zstar,k);
                if(k < atomicNumber_-Zstar-1) {
                    if(Bk == bindingEnergyCarlson(Zstar,k+1)){
                        N++;
                        continue;
                    }
                }
                eps = E/Bk;

                if(eps > 1.0)
                {
                    b_dash = Bk/510997.37;
                    eps_dash = E/510997.37;
                    u_dash = b_dash;//mean kinetic energy of orbital k: unknown use Uk=Bk as in Smilie

                    beta_eps2 = 1.0-1.0/((1.0+eps_dash)*(1.0+eps_dash));
                    beta_b2 = 1.0-1.0/((1.0+b_dash)*(1.0+b_dash));
                    beta_u2 = 1.0-1.0/((1.0+u_dash)*(1.0+u_dash));

                    sigma0 = 2.0*cm::pi*r_e*r_e*N/(b_dash*(beta_eps2+beta_b2+beta_u2));

                    ediv2 = (1.0+eps_dash*0.5); ediv2 *= ediv2;

                    A1 = 1.0/(1.0+eps)*(1.0+2*eps_dash)/ediv2;
                    A2 = (eps-1.0)/2.0 * b_dash*b_dash/ediv2;
                    A3 = ::log(beta_eps2/(1.0-beta_eps2))-beta_eps2-::log(2.0*b_dash);

                    sigmak = sigma0*(0.5*A3*(1.0-1.0/(eps*eps))+1.0-1.0/eps+A2-A1*::log(eps));

                    ionizationk0CrossSection_[i] += sigmak;
                }
                N=1;
            }
        }
    }

}

template<class CloudType, Foam::crossSectionType Type>
scalar Foam::KimCarlsonCrossSection<CloudType,Type>::bindingEnergyCarlson(label Zstar, label k)
{
    return ionizationEneries_[Zstar] - bindingEneries_[atomicNumber_-Zstar-1] + bindingEneries_[k];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::KimCarlsonCrossSection<CloudType,Type>::~KimCarlsonCrossSection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::KimCarlsonCrossSection<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    if(eVEnergy <= 0.0)
        return 0.0;

    scalar s,d;
    scalar crossSection = 0.0;
    label l;

    s = 6.142165*::log(eVEnergy);
    if( s < 0.0 )
        return crossSection;

    if( s < 99.0 )//interpolate
    {
       l = label(s);
       d = s-scalar(l);
       crossSection = (ionizationk0CrossSection_[l+1]-ionizationk0CrossSection_[l])*d + ionizationk0CrossSection_[l];
    }
    else//extrapolate
    {
        d = s - 99.0;
        crossSection = (ionizationk0CrossSection_[99]-ionizationk0CrossSection_[98])*d + ionizationk0CrossSection_[99];
    }
    if(crossSection < 0.0)//Cannot ionize
       crossSection = 0.0;

    return crossSection;
}

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::KimCarlsonCrossSection<CloudType,Type>::threshold() const
{
    return ionizationEneries_[0];
}

// ************************************************************************* //
