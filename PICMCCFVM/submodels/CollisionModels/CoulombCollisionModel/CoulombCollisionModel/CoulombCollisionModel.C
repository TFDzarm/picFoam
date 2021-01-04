/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "CoulombCollisionModel.H"
#include "PairingAlgorithm.H"
#include "IonizationModel.H"
#include "constants.H"
#include "WeightCorrectionModel.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CoulombCollisionModel<CloudType>::CoulombCollisionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    pairingAlgorithm_(),
    ionizationModel_(),
    allowIntraCollision_(true),
    calculateDebyeLength_(true),
    weightCorrection_(),
    coulombLog_()
{}


template<class CloudType>
Foam::CoulombCollisionModel<CloudType>::CoulombCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    pairingAlgorithm_//Choose and construct the model
    (
        PairingAlgorithm<CloudType>::New
        (
            coeffDict_,
            owner,
            *this
        )
    ),
    ionizationModel_//Choose and construct the model
    (
        IonizationModel<CloudType>::New
        (
            coeffDict_,
            owner
        )
    ),
    allowIntraCollision_(readBool(coeffDict_.lookup("allowIntraCollision"))),
    calculateDebyeLength_(true),
    weightCorrection_//Choose and construct the model
    (
        WeightCorrectionModel<CloudType>::New
        (
            coeffDict_,
            owner
        )
    ),
    coulombLog_()
{
    label nSpecies = owner.typeIdList().size();

    //Simply save the coulomb logarithm twice ignoring the order
    coulombLog_.setSize(nSpecies,List<scalar>(nSpecies,0.0));

    //Return empty dict if entry is not found
    const dictionary& coulombLogdict = coeffDict_.subOrEmptyDict("Collisions");
    label fixedCoulombLog = 0;

    //Read user provided coulomb logarithm
    forAllConstIter(IDLList<entry>, coulombLogdict, iter)
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();
            wordList cP = subDict.lookup("collisionPartner");
            if(cP.size() != 2)
                FatalErrorInFunction << "Group of collisionPartners contains more than two species" << nl << cP << abort(FatalError);
            label p1 = findIndex(owner.typeIdList(),cP[0]);
            label p2 = findIndex(owner.typeIdList(),cP[1]);
            scalar cLog = readScalar(subDict.lookup("coulombLog"));
            if(coulombLog_[p1][p2] != 0 && coulombLog_[p1][p2] != cLog)
            {
                FatalErrorInFunction << "Ambiguous definition for the coulombLog of species" << cP[0] << " and species " << cP[1] << " (" << coulombLog_[p1][p2] << " and " << cLog << ")" << abort(FatalError);
            }
            else if(coulombLog_[p1][p2] == cLog)
                continue;

            //Save the coulomb logarithm
            coulombLog_[p1][p2] = cLog;coulombLog_[p2][p1] = cLog;
            if(coulombLog_[p1][p2] > 0.0)
            {
                fixedCoulombLog++;
                Info << "    CoulombCollisionModel: Using fixed coulombLog=" << coulombLog_[p1][p2] << " for the collisions of species " << cP[0] << " and species " << cP[1] << endl;
            }
        }
    }
    //If the user provided a logarithm for all possible collisions no calculations need to be done
    if(fixedCoulombLog == nSpecies*(nSpecies+1)/2)
        calculateDebyeLength_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CoulombCollisionModel<CloudType>::~CoulombCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::CoulombCollisionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType&
Foam::CoulombCollisionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary&
Foam::CoulombCollisionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::CoulombCollisionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType>
void Foam::CoulombCollisionModel<CloudType>::handleCollisions()
{
    Info << "[Coulomb Collisions]" << nl;
    pairingAlgorithm_().pairANDcollide();
}


template<class CloudType>
scalar Foam::CoulombCollisionModel<CloudType>::debyeLength(label celli)
{
    const CloudType& cloud(this->owner());
    scalar debyeLength = 0.0;
    scalar debye(0.0), n_max(0.0);

    //Go through all charged species
    forAll(cloud.chargedSpecies(),si)
    {
        scalar T(0.0), n(0.0);
        label speci = cloud.chargedSpecies()[si];
        scalar mass = cloud.constProps(speci).mass();

        n = 0.0;

        scalar charge = 0.0;
        scalar vSqr(0.0);
        //vector v = Zero;

        //Calculate average properties
        forAll(cloud.sortedCellOccupancy()[celli][speci],parti)
        {
            typename CloudType::parcelType* p = cloud.sortedCellOccupancy()[celli][speci][parti];
            vSqr += (p->U() & p->U())*p->nParticle();
            charge += p->charge()*p->nParticle();
            //v += p->U()*p->nParticle();
            n += p->nParticle();
        }

        if(n <= 0.0)
            continue;
        //v /= n;
        charge /= n;//avg charge
        vSqr /= n;//avg vSqr
        n /= cloud.mesh().cellVolumes()[celli];//number density

        if(n > n_max)
            n_max = n;

        //Temperatur according to Nanbu assumes Maxwaillian distribution (we ignore drift velocity this prevent error at low parcel counts)
        T = mass/(3.0*constant::physicoChemical::k.value())*(vSqr/*-(v&v)*/);

        if(T <= 0.0)
            continue;

        debye += n*charge*charge/constant::electromagnetic::epsilon0.value() * 1.0/(T*constant::physicoChemical::k.value());

    }
    if(debye <= 0.0)
        return debyeLength;

    //mean interatomic distance (Wigner-Seitz radius)
    scalar r_min = ::pow(4.0*constant::mathematical::pi*n_max/3.0,-1.0/3.0); 

    debyeLength = max(r_min,::sqrt(1/debye));

    return debyeLength;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoulombCollisionModelNew.C"

// ************************************************************************* //
