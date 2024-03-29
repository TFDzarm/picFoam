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
    coulombLog_(),
    average_coulombLog_(0.0),
    debyeLength_accountForDrift_(true),
    debyeLength_Species_(),
    average_b_(0.0),
    calculate_bmin_(true)
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
    coulombLog_(),
    average_coulombLog_(0.0),
    debyeLength_accountForDrift_(true),
    debyeLength_Species_(),
    average_b_(0.0),
    calculate_bmin_(true)
{

    //Check if a DebyeLength subDict was provided and if it contains a species list
    const dictionary& debyeLengthdict = coeffDict_.subOrEmptyDict("DebyeLength");
    wordList debyeWordList = debyeLengthdict.lookupOrDefault("species",wordList::null());
    if(debyeWordList.empty())//If not consider all species for the calculation of the Debye length
    {
        Info << "|->    Consider all charged species in the calculation of the debye length." << endl;
        List<label> copy = owner.chargedSpecies().clone();
        debyeLength_Species_.transfer(copy);
    }
    else//Else only consider species given in the list
    {
        Info << "|->    Only consider species in list " << debyeWordList << " in the calculation of the Debye length." << endl;

        DynamicList<label> debyeLength_Species;
        forAll(debyeWordList,i)
        {
            word species = debyeWordList[i];
            label typeId = findIndex(owner.typeIdList(),species);

            if(typeId < 0)
                FatalErrorInFunction << "Unable to find the species " << species << nl << abort(FatalError);

            if(owner.constProps(typeId).charge() == 0.0)
                FatalErrorInFunction << "Species " << species << " considered for calculation of the Debye length is not a particle with charge" << nl << abort(FatalError);

            debyeLength_Species.append(typeId);
        }
        debyeLength_Species_.transfer(debyeLength_Species);
    }

    //Check if we consider bmin in the Coulomb log calculation else use default value true
    calculate_bmin_.readIfPresent("calculate_bmin",debyeLengthdict);
    if(calculate_bmin_)
        Info << "|->    Consider bmin in Coulomb logarithm calculation" << endl;
    else
        Info << "|->    Ignore bmin in Coulomb logarithm calculation" << endl;

    //Check if we calculate the species temperatures while accounting for the drift velocity
    debyeLength_accountForDrift_.readIfPresent("accountForDrift",debyeLengthdict);
    if(!debyeLength_accountForDrift_)
        Info << "|->    Do not incorporate dirft velocities in Debye length calculation" << endl;
    else
        Info << "|->    Incorporate the drift velocities in Debye length calculation" << endl;

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

            if(p1 < 0 || p2 < 0)
                FatalErrorInFunction << "Unable to find the typeId of one or both species in pair (" << cP[0] << " " << cP[1] << ")" << nl << abort(FatalError);

            if(owner.constProps(p1).charge() == 0.0 || owner.constProps(p2).charge() == 0.0)
                FatalErrorInFunction << "Collision between species " << cP[0] << " and " << cP[1] << " is not a Coulomb collision" << nl << abort(FatalError);

            scalar cLog = readScalar(subDict.lookup("coulombLog"));
            if(coulombLog_[p1][p2] != 0 && coulombLog_[p1][p2] != cLog)
            {
                FatalErrorInFunction << "Ambiguous definition for the coulombLog of species " << cP[0] << " and species " << cP[1] << " (" << coulombLog_[p1][p2] << " and " << cLog << ")" << abort(FatalError);
            }
            else if(coulombLog_[p1][p2] == cLog)
                continue;

            //Save the coulomb logarithm
            coulombLog_[p1][p2] = cLog;coulombLog_[p2][p1] = cLog;
            if(coulombLog_[p1][p2] > 0.0)
            {
                fixedCoulombLog++;
                Info << "|->    Using fixed coulombLog = " << coulombLog_[p1][p2] << " for the collisions of species " << cP[0] << " and species " << cP[1] << endl;
            }
        }
    }

    //Use all charged species instead of debyeLength_Species_(may be not all charged species),
    //so we still calculate the Debye length for species where no Coulomb logarithm is provided
    label nChargedSpecies = owner.chargedSpecies().size();

    //If the user provided a logarithm for all possible collisions no calculations need to be done
    if(fixedCoulombLog == nChargedSpecies*(nChargedSpecies+1)/2)
    {
        calculateDebyeLength_ = false;
        Info << "|->    Provided a coulomb logarithm for every possible collision of charged species. Skipping calculation of Debye length" << endl;
    }
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

    //Go through all considered species
    forAll(debyeLength_Species_,si)
    {
        scalar T(0.0), n(0.0);
        label speci = debyeLength_Species_[si];
        scalar mass = cloud.constProps(speci).mass();

        n = 0.0;

        scalar charge = 0.0;
        scalar vSqr(0.0);
        vector v = Zero;
        scalar vDrift2 = 0.0;

        //Calculate average properties
        forAll(cloud.sortedCellOccupancy()[celli][speci],parti)
        {
            typename CloudType::parcelType* p = cloud.sortedCellOccupancy()[celli][speci][parti];
            vSqr += (p->U() & p->U())*p->nParticle();
            charge += p->charge()*p->nParticle();
            v += p->U()*p->nParticle();
            n += p->nParticle();
        }

        if(n <= 0.0)
            continue;
        v /= n;
        charge /= n;//avg charge
        vSqr /= n;//avg vSqr
        n /= cloud.mesh().cellVolumes()[celli];//number density

        if(n > n_max)
            n_max = n;

        if(debyeLength_accountForDrift_)//Do we add the drift velocity?
            vDrift2 = (v&v);

        //Temperatur according to Nanbu assumes Maxwaillian distribution
        T = mass/(3.0*constant::physicoChemical::k.value())*(vSqr-vDrift2);

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
