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

#include "SputterEmission.H"
#include "Random.H"
// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SputterEmission<CloudType>::SputterEmission
(
    const dictionary& dict,
    CloudType& cloud
)
:
    EmissionModel<CloudType>(dict,cloud,typeName),
    species_(),
    energies_(),
    probability_()
{
    label nSpecies = cloud.typeIdList().size();


    species_.setSize(nSpecies);
    energies_.setSize(nSpecies);
    probability_.setSize(nSpecies,scalarList(nSpecies,0.0));

    forAllConstIter(IDLList<entry>, this->coeffDict(), iter)
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            label id = findIndex(cloud.typeIdList(),iter().keyword());
            if(id == -1)
                FatalErrorInFunction << "Undefined typeId " << iter().keyword() << abort(FatalError);

            scalarList energies = subDict.lookup("sputterEnergies");
            energies_[id].transfer(energies);

            scalarList probability = subDict.lookup("sputterProbabilities");
            probability_[id].transfer(probability);

            wordList typeIdWords = subDict.lookup("sputterSpecies");

            if(typeIdWords.size() != energies_[id].size() || typeIdWords.size() != probability.size())
                FatalErrorInFunction << "Different number of entries" << nl << abort(FatalError);

            DynamicList<label> sputterSpecies;
            forAll(typeIdWords,i)
            {
                label sId = findIndex(cloud.typeIdList(), typeIdWords[i]);
                if(sId == -1)
                    FatalErrorInFunction << "Undefined typeId " << iter().keyword() << abort(FatalError);


                sputterSpecies.append(sId);
            }
            species_[id].transfer(sputterSpecies);
        }

    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SputterEmission<CloudType>::~SputterEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SputterEmission<CloudType>::emission()
{}

template<class CloudType>
void Foam::SputterEmission<CloudType>::collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    CloudType& cloud(this->cloud());
    Random& rndGen(cloud.rndGen());

    label typeId = p.typeId();
    if(p.patch() == this->patchId() && species_[typeId].size() > 0)
    {
        forAll(probability_[typeId],i)
        {
                scalar prop = probability_[typeId][i];
                if(rndGen.scalar01() < prop)
                {
                    scalar temperature;
                    if(energies_[typeId][i] < 0.0)
                    {
                        //This calculation includes the dirft velocity, that is an error...
                        temperature = (p.U()&p.U())*cloud.constProps(typeId).mass()/(3*constant::physicoChemical::k.value());
                    }
                    else
                        temperature = energies_[typeId][i]/*eV*/*constant::electromagnetic::e.value()/constant::physicoChemical::k.value();

                    this->emitParticleExplicit(temperature,typeId,p.stepFraction(),p.nParticle());
                }
        }
    }

}
// ************************************************************************* //
