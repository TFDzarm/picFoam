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

#include "TakizukaAbePairing.H"
#include "IonizationModel.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TakizukaAbePairing<CloudType>::TakizukaAbePairing
(
    const dictionary& dict,
    CloudType& cloud,
    CoulombCollisionModel<CloudType>& collisionModel
)
:
    PairingAlgorithm<CloudType>(dict,cloud,collisionModel,typeName)
{
    //nParticleEqWeight is zero if the species have different weights
    if(cloud.nParticleEqWeight() == 0.0)
       FatalErrorInFunction << "Model TakizukaAbe is only valid for equal weighted particles" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TakizukaAbePairing<CloudType>::~TakizukaAbePairing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TakizukaAbePairing<CloudType>::pairANDcollide()
{
    namespace cm = constant::mathematical;
    namespace ce = constant::electromagnetic;
    namespace cu = constant::universal;

    CloudType& cloud(this->cloud());
    IonizationModel<CloudType>& ionizationModel(this->coulombCollisionModel().ionizationModel());
    const polyMesh& mesh(cloud.mesh());

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());
    const List<label>& chargedSpecies(cloud.chargedSpecies());

    const scalar& nParticle(cloud.nParticleEqWeight());

    label coulombCollisions = 0;
    label ionizations = 0;
    scalar average_DebyeLength = 0.0;
    forAll(sortedCellOccupancy, celli)// Go through all cells
    {
        scalar debyeLength = 0.0;
        if(this->coulombCollisionModel().calculateDebyeLength()) {
            debyeLength = this->coulombCollisionModel().debyeLength(celli);
            average_DebyeLength += debyeLength;
        }

        forAll(chargedSpecies,i)//Species 1
        {
            label alpha = chargedSpecies[i];
            forAllReverse(chargedSpecies,j)//Reverse species 2...
            {
                if(j < i)//... so that we can break and the same species do not collide twice
                    break;

                //Collision partner species
                label beta = chargedSpecies[j];

                //Number of particles in this cell
                label nAlpha= sortedCellOccupancy[celli][alpha].size();
                label nBeta = sortedCellOccupancy[celli][beta].size();

                if(alpha == beta && nAlpha > 1)//Collision of species alpha with iself
                {
                    if(!this->coulombCollisionModel().allowIntraCollision())
                        continue;

                    //Shuffle colliding parcel ids
                    List<label> collisionSelect(nAlpha);
                    forAll(collisionSelect,pos)
                        collisionSelect[pos] = pos;
                    shuffle(collisionSelect);

                    //This are constant in the loop
                    scalar n = nAlpha*nParticle/mesh.cellVolumes()[celli];

                    if(nAlpha % 2 != 0)//Odd number of parcels
                    {
                        //Method by Takizuka and Abe: Collide Parcel 1 twice (Case 1b)
                        typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[0]];
                        typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[1]];

                        scalar n12 = (nAlpha+3.0)*nParticle/mesh.cellVolumes()[celli];//pairs (nAlpha+3.0)/2.0? -> should be (nAlpha+3.0) so we account for 2 (eg 7: 1-2,3-4,5-6,[5-7],[6-7] == 5 collisions (7+3)/2 paare) extra collisons else if we divide by two  n*n/n21 approx 2*n
                        //P with Q
                        this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n12);

                        typename CloudType::parcelType* parcelR = sortedCellOccupancy[celli][alpha][collisionSelect[2]];

                        //Q with R
                        this->coulombCollisionModel().collide(parcelQ,parcelR,debyeLength,n,n,n12);

                        //P with R
                        this->coulombCollisionModel().collide(parcelP,parcelR,debyeLength,n,n,n12);

                        for (label x=3; x<(collisionSelect.size()-1); x+=2)
                        {
                            parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];
                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n12);
                        }
                        coulombCollisions += (2+collisionSelect.size()/2);

                    }
                    else//Even number of parcels
                    {
                        for (label x=0; x<(collisionSelect.size()-1); x+=2)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];

                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n);//n12=n because n real particle collide N/2 + N/2
                        }
                        coulombCollisions += collisionSelect.size()/2;
                    }

                }
                else if(alpha != beta && min(nAlpha,nBeta) != 0)//Collision of species alpha with beta
                {
                    namespace cu = constant::universal;

                    if(nAlpha == nBeta)
                    {

                        DynamicList<typename CloudType::parcelType*> pAlpha = sortedCellOccupancy[celli][alpha];
                        scalar nP = nAlpha*nParticle/mesh.cellVolumes()[celli];
                        scalar nQ = nBeta*nParticle/mesh.cellVolumes()[celli];
                        ionizationModel.initIonization(alpha,beta,nP,nQ);

                        //Shuffle colliding parcel ids
                        List<label> collisionSelect(nAlpha);
                        forAll(collisionSelect,pos)
                            collisionSelect[pos] = pos;
                        shuffle(collisionSelect);

                        for(label indx = 0; indx < nAlpha; indx++)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[indx]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][beta][collisionSelect[nAlpha-indx-1]];

                            //prepare ionization nei/ne
                            ionizationModel.prepareNumberDensities(*parcelP,*parcelQ);
                        }

                        forAll(pAlpha, indx)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[indx]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][beta][collisionSelect[nAlpha-indx-1]];

                            ionizationModel.prepareIonization(*parcelP,*parcelQ);
                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,nP,nQ,nP);// s= nP*nQ/nP * dt => s= nQ *dt = nP * dt
                            if(ionizationModel.ionize())
                                ionizations++;

                        }
                        coulombCollisions+=nAlpha;
                    }
                    else
                    {
                        //Method by Takizuka and Abe: Collide Parcel 1 twice (Case 2b)
                        label minSpecies = (nAlpha < nBeta)? alpha: beta;
                        label maxSpecies = (nAlpha > nBeta)? alpha: beta;
                        label nMax = max(nAlpha,nBeta);
                        label nMin = min(nAlpha,nBeta);

                        scalar nMinSpecies = nMin*nParticle/mesh.cellVolumes()[celli];
                        scalar nMaxSpecies = nMax*nParticle/mesh.cellVolumes()[celli];


                        //Shuffle colliding parcel ids
                        List<label> collisionSelectMax(nMax);
                        List<label> collisionSelectMin(nMin);
                        forAll(collisionSelectMax,pos) {
                            collisionSelectMax[pos] = pos;
                            if(pos < nMin)
                                collisionSelectMin[pos] = pos;
                        }
                        shuffle(collisionSelectMin);
                        shuffle(collisionSelectMax);


                        label delta = nMax/nMin;
                        scalar delta_r = (static_cast<scalar>(nMax%nMin))/nMin;

                        label nMinGroup1 = delta_r*nMin;
                        label nMaxGroup1 = (delta+1)*delta_r*nMin;

                        label nMinGroup2 = (1.0-delta_r)*nMin;
                        label nMaxGroup2 = delta*(1.0-delta_r)*nMin;

                        if(nMinGroup1 > 0)
                        {
                            label g1select1 = 0;//minSpecies
                            label g1select2 = 0;

                            ionizationModel.initIonization(minSpecies,maxSpecies,nMinSpecies,nMaxSpecies);
                            for(label g1 = 0; g1 < nMaxGroup1;g1++)
                            {
                                if(g1select1 >= nMinGroup1)
                                    g1select1 = 0;

                                typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[g1select1]];
                                typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[g1select2]];
                                //prepare ionization nei/ne
                                ionizationModel.prepareNumberDensities(*parcelP,*parcelQ);

                                g1select1++;
                                g1select2++;
                            }
                            g1select1 = 0;
                            g1select2 = 0;

                            for(label g1 = 0; g1 < nMaxGroup1;g1++)
                            {
                                if(g1select1 >= nMinGroup1)
                                    g1select1 = 0;

                                typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[g1select1]];//minSpecies
                                typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[g1select2]];



                                ionizationModel.prepareIonization(*parcelP,*parcelQ);//prepare: determines the ion species
                                this->coulombCollisionModel().collide(parcelQ,parcelP,debyeLength,nMaxSpecies,nMinSpecies,nMaxSpecies);//Advance max species (pass in Q first)!!!
                                if(ionizationModel.ionize())
                                    ionizations++;

                                g1select1++;
                                g1select2++;
                            }
                        }

                        label g2select1 = nMinGroup1;//minSpecies
                        label g2select2 = nMaxGroup1;

                        ionizationModel.initIonization(minSpecies,maxSpecies,nMinSpecies,nMaxSpecies);
                        for(label g2 = 0; g2 < nMaxGroup2;g2++)
                        {
                            if(g2select1 >= (nMinGroup2+nMinGroup1))
                                g2select1 = nMinGroup1;

                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[g2select1]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[g2select2]];
                            //prepare ionization nei/ne
                            ionizationModel.prepareNumberDensities(*parcelP,*parcelQ);

                            g2select1++;
                            g2select2++;
                        }
                        g2select1 = nMinGroup1;
                        g2select2 = nMaxGroup1;

                        for(label g2 = 0; g2 < nMaxGroup2;g2++)
                        {
                            if(g2select1 >= (nMinGroup2+nMinGroup1))
                                g2select1 = nMinGroup1;

                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[g2select1]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[g2select2]];



                            ionizationModel.prepareIonization(*parcelP,*parcelQ);
                            this->coulombCollisionModel().collide(parcelQ,parcelP,debyeLength,nMaxSpecies,nMinSpecies,nMaxSpecies);//Advance max species (pass in Q first)!!!
                            if(ionizationModel.ionize())
                               ionizations++;

                            g2select1++;
                            g2select2++;
                        }
                        coulombCollisions += (nMaxGroup1 + nMaxGroup2);
                    }

                }
            }
        }
    }

    reduce(coulombCollisions, sumOp<label>());
    reduce(ionizations, sumOp<label>());
    reduce(average_DebyeLength, sumOp<scalar>());
    reduce(this->coulombCollisionModel().average_impactParameter(), sumOp<scalar>());
    reduce(this->coulombCollisionModel().average_coulombLog(), sumOp<scalar>());

    if (coulombCollisions)
    {
        Info << "    Coulomb Collisions              = " << coulombCollisions << nl
             << "    Ionizations                     = " << ionizations << nl
             << "    Average Debye length            = " << average_DebyeLength/this->nCellsGlobal() << nl
             << "    Average b0                      = " << this->coulombCollisionModel().average_impactParameter()/coulombCollisions << nl
             << "    Average Coulomb logarithm       = " << this->coulombCollisionModel().average_coulombLog()/coulombCollisions << nl << endl;
    }
    else
    {
        Info<< "    No Coulomb collisions" << nl << endl;
    }

    //Reset
    this->coulombCollisionModel().average_coulombLog() = 0.0;
    this->coulombCollisionModel().average_impactParameter() = 0.0;
}


// ************************************************************************* //
