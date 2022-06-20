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

#include "NanbuYonemuraPairing.H"
#include "IonizationModel.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuYonemuraPairing<CloudType>::NanbuYonemuraPairing
(
    const dictionary& dict,
    CloudType& cloud,
    CoulombCollisionModel<CloudType>& collisionModel
)
:
    PairingAlgorithm<CloudType>(dict,cloud,collisionModel,typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuYonemuraPairing<CloudType>::~NanbuYonemuraPairing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NanbuYonemuraPairing<CloudType>::pairANDcollide()
{
    namespace cm = constant::mathematical;
    namespace ce = constant::electromagnetic;
    namespace cu = constant::universal;

    CloudType& cloud(this->cloud());
    IonizationModel<CloudType>& ionizationModel(this->coulombCollisionModel().ionizationModel());
    const polyMesh& mesh(cloud.mesh());

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());
    const List<label>& chargedSpecies(cloud.chargedSpecies());

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
                label nParcelAlpha= sortedCellOccupancy[celli][alpha].size();
                label nParcelBeta = sortedCellOccupancy[celli][beta].size();

                //Count weights for nAlpha number density of species alpha
                scalar nAlpha = 0.0;
                forAll(sortedCellOccupancy[celli][alpha],iAlpha)
                {
                    typename CloudType::parcelType* p = sortedCellOccupancy[celli][alpha][iAlpha];
                    nAlpha += p->nParticle();
                }

                //Count weights for nBeta number density of species beta
                scalar nBeta = 0.0;
                if(alpha != beta)
                {
                    forAll(sortedCellOccupancy[celli][beta],iBeta)
                    {
                        typename CloudType::parcelType* p = sortedCellOccupancy[celli][beta][iBeta];
                        nBeta += p->nParticle();
                    }
                }


                if(alpha == beta && nParcelAlpha > 1)//Collision of species alpha with itself
                {
                    if(!this->coulombCollisionModel().allowIntraCollision())
                        continue;

                    //Shuffle colliding parcel ids
                    List<label> collisionSelect(nParcelAlpha);
                    forAll(collisionSelect,pos)
                        collisionSelect[pos] = pos;
                    shuffle(collisionSelect);

                    //Number density n of species alpha
                    scalar n = nAlpha/mesh.cellVolumes()[celli];

                    if(nParcelAlpha % 2 != 0)//Odd number of parcels
                    {
                        scalar nAlphaAlpha = 0.0;

                        //Methode by Nanbu and Yonemura collide last one with first one Fig3b (one extra collison)
                        typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[collisionSelect.size()-1]];
                        typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[0]];

                        nAlphaAlpha += 2.0*min(parcelP->nParticle(),parcelQ->nParticle());

                        for (label x=0; x<(collisionSelect.size()-2); x+=2)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];

                            nAlphaAlpha += 2.0*min(parcelP->nParticle(),parcelQ->nParticle());
                        }
                        scalar n12 = nAlphaAlpha/mesh.cellVolumes()[celli];

                        for (label x=0; x<(collisionSelect.size()-2); x+=2)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];

                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n12);
                        }
                        //The extra collision
                        this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n12);

                        coulombCollisions += (1+collisionSelect.size()/2);

                    }
                    else//Even number of parcels
                    {
                        //Count number of real particle that collide
                        scalar nAlphaAlpha = 0.0;
                        for (label x=0; x<(collisionSelect.size()-1); x+=2)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];
                            nAlphaAlpha += 2.0*min(parcelP->nParticle(),parcelQ->nParticle());//2* because if w1 = 3 and w2=5 -> min 3 real ones of p1 collide with 3 real ones of p2
                        }
                        scalar n12 = nAlphaAlpha/mesh.cellVolumes()[celli];

                        for (label x=0; x<(collisionSelect.size()-1); x+=2)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[x]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][alpha][collisionSelect[x+1]];

                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,n,n,n12);
                        }
                        coulombCollisions += collisionSelect.size()/2;
                    }

                }
                else if(alpha != beta && min(nParcelAlpha,nParcelBeta) != 0)//Collision of species alpha with beta
                {

                    if(nParcelAlpha == nParcelBeta)
                    {
                        //Number densities
                        scalar nP = nAlpha/mesh.cellVolumes()[celli];
                        scalar nQ = nBeta/mesh.cellVolumes()[celli];

                        List<label> collisionSelect(nParcelAlpha);
                        forAll(collisionSelect,pos)
                            collisionSelect[pos] = pos;
                        shuffle(collisionSelect);

                        ionizationModel.initIonization(alpha,beta,nP,nQ);
                        //Count hybrid number density of real particle that collide
                        scalar nAlphaBeta = 0.0;
                        for(label indx = 0; indx < nParcelAlpha; indx++)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[indx]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][beta][collisionSelect[nParcelAlpha-indx-1]];
                            nAlphaBeta += min(parcelP->nParticle(),parcelQ->nParticle());
                            //prepare ionization nei
                            ionizationModel.prepareNumberDensities(*parcelP,*parcelQ);
                        }
                        scalar nPQ = nAlphaBeta/mesh.cellVolumes()[celli];
                        for(label indx = 0; indx < nParcelAlpha; indx++)
                        {
                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][alpha][collisionSelect[indx]];
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][beta][collisionSelect[nParcelAlpha-indx-1]];

                            ionizationModel.prepareIonization(*parcelP,*parcelQ);

                            this->coulombCollisionModel().collide(parcelP,parcelQ,debyeLength,nP,nQ,nPQ);


                            if(ionizationModel.ionize())
                                ionizations++;

                        }
                        coulombCollisions+=nParcelAlpha;
                    }
                    else
                    {
                        label minSpecies = (nParcelAlpha < nParcelBeta)? alpha: beta;
                        label maxSpecies = (nParcelAlpha > nParcelBeta)? alpha: beta;
                        label nMax = max(nParcelAlpha,nParcelBeta);
                        label nMin = min(nParcelAlpha,nParcelBeta);

                        List<label> collisionSelectMax(nMax);
                        List<label> collisionSelectMin(nMin);
                        forAll(collisionSelectMax,pos) {
                            collisionSelectMax[pos] = pos;
                            if(pos < nMin)
                                collisionSelectMin[pos] = pos;
                        }
                        shuffle(collisionSelectMin);
                        shuffle(collisionSelectMax);

                        scalar nMinSpecies = 0.0;
                        scalar nMaxSpecies = 0.0;
                        if(maxSpecies == alpha) {
                            nMaxSpecies = nAlpha/mesh.cellVolumes()[celli];
                            nMinSpecies = nBeta/mesh.cellVolumes()[celli];
                        }
                        else {
                            nMaxSpecies = nBeta/mesh.cellVolumes()[celli];
                            nMinSpecies = nAlpha/mesh.cellVolumes()[celli];
                        }


                        ionizationModel.initIonization(minSpecies,maxSpecies,nMinSpecies,nMaxSpecies);

                        scalar nAlphaBeta = 0.0;
                        label k = 0;
                        for(label l = 0; l < nMax; l++)
                        {
                            if(k >= nMin)
                                k = 0;

                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[k]];//minSpecies
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[l]];

                            nAlphaBeta+= min(parcelP->nParticle(),parcelQ->nParticle());
                            //prepare ionization nei
                            ionizationModel.prepareNumberDensities(*parcelP,*parcelQ);
                            k++;
                        }

                        scalar n12 = nAlphaBeta/mesh.cellVolumes()[celli];

                        k = 0;
                        for(label l = 0; l < nMax; l++)
                        {
                            if(k >= nMin)
                                k = 0;

                            typename CloudType::parcelType* parcelP = sortedCellOccupancy[celli][minSpecies][collisionSelectMin[k]];//minSpecies
                            typename CloudType::parcelType* parcelQ = sortedCellOccupancy[celli][maxSpecies][collisionSelectMax[l]];

                            ionizationModel.prepareIonization(*parcelP,*parcelQ);

                            this->coulombCollisionModel().collide(parcelQ,parcelP,debyeLength,nMaxSpecies,nMinSpecies,n12);

                            if(ionizationModel.ionize())
                                ionizations++;
                            k++;
                        }
                        coulombCollisions += nMax;
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
    Info << "    Coulomb Collisions              = " << coulombCollisions << nl
         << "    Ionizations                     = " << ionizations << nl
         << "    Average Debye length            = " << average_DebyeLength/this->nCellsGlobal() << nl
         << "    Average b0                      = " << this->coulombCollisionModel().average_impactParameter()/coulombCollisions << nl
         << "    Average Coulomb logarithm       = " << this->coulombCollisionModel().average_coulombLog()/coulombCollisions << nl << endl;

    //Reset
    this->coulombCollisionModel().average_coulombLog() = 0.0;
    this->coulombCollisionModel().average_impactParameter() = 0.0;
}


// ************************************************************************* //
