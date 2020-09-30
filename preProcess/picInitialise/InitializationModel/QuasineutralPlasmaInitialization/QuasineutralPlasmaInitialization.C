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

#include "QuasineutralPlasmaInitialization.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::QuasineutralPlasmaInitialization<CloudType>::QuasineutralPlasmaInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName),
    ionSpecies_(this->coeffDict().lookup("ionSpecies")),
    ionTypeId_(-1),
    cellCount_(cloud.typeIdList().size(),1),
    useMaxwellian_(true)
{
    ionTypeId_ = findIndex(cloud.typeIdList(), ionSpecies_);
    if(ionTypeId_ == -1)
        FatalErrorInFunction << "Ion species " << ionSpecies_ << " is not defined in cloud" << abort(FatalError);


    if(cloud.electronTypeId() == -1)
        FatalErrorInFunction << "Electron species is not defined in cloud" << abort(FatalError);

    scalar eCharge = cloud.constProps()[cloud.electronTypeId()].charge();
    scalar iCharge = cloud.constProps()[ionTypeId_].charge();
    if((iCharge + eCharge) != 0.0)
        FatalErrorInFunction << "Model expects the ion and electron charge to be of same magnitude and opposite signs" << abort(FatalError);

    cellCount_ = readLabel(this->coeffDict().lookup("cellCount"));
    cellCount_[ionTypeId_] /= cloud.constProps()[ionTypeId_].nParticle();
    cellCount_[cloud.electronTypeId()] /= cloud.constProps()[cloud.electronTypeId()].nParticle();

    Info << "Parcels ion: " << cellCount_[ionTypeId_] << nl
         << "Parcels electron: " << cellCount_[cloud.electronTypeId()] << endl;

    useMaxwellian_ = readBool(this->coeffDict().lookup("useMaxwellian"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::QuasineutralPlasmaInitialization<CloudType>::~QuasineutralPlasmaInitialization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::QuasineutralPlasmaInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    Random& rndGen(cloud.rndGen());
    const fvMesh& mesh(cloud.mesh());

    scalar eMass = cloud.constProps()[cloud.electronTypeId()].mass();
    scalar iMass = cloud.constProps()[ionTypeId_].mass();

    //Ions
    forAll(mesh.cells(), celli)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
             mesh,
             celli
        );
        for(label i = 0; i < cellCount_[ionTypeId_]; i++)
        {
            label rndTet = static_cast<label>((static_cast<scalar>(cellTets.size())-1.0)*rndGen.scalar01());
            const tetIndices& cellTetIs = cellTets[rndTet];
            tetPointRef tet = cellTetIs.tet(mesh);
            point p = tet.randomPoint(rndGen);

            meshTools::constrainDirection(mesh, mesh.solutionD(), p);

            vector U;
            if(useMaxwellian_)
            {
                U = cloud.equipartitionLinearVelocity
                (
                       this->temperatures_[ionTypeId_],
                       iMass
                );
            }
            else
            {
                vector rndDir = rndGen.sample01<vector>();
                rndDir*= 1.0/mag(rndDir);


                U = rndDir*Foam::sqrt(this->temperatures_[ionTypeId_]*constant::physicoChemical::k.value()*3.0/iMass);
            }

            cloud.addNewParcel(p,celli,U,ionTypeId_);
        }


    }

    //Electrons
    forAll(mesh.cells(), celli)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
             mesh,
             celli
        );
        for(label i = 0; i < cellCount_[cloud.electronTypeId()]; i++)
        {
            label rndTet = static_cast<label>((static_cast<scalar>(cellTets.size())-1.0)*rndGen.scalar01());
            const tetIndices& cellTetIs = cellTets[rndTet];
            tetPointRef tet = cellTetIs.tet(mesh);
            point p = tet.randomPoint(rndGen);

            meshTools::constrainDirection(mesh, mesh.solutionD(), p);

            vector U;
            if(useMaxwellian_)
            {
                U = cloud.equipartitionLinearVelocity
                (
                       this->temperatures_[cloud.electronTypeId()],
                       eMass
                );
            }
            else
            {
                vector rndDir = rndGen.sample01<vector>();
                rndDir*= 1.0/mag(rndDir);

                U = rndDir*Foam::sqrt(this->temperatures_[cloud.electronTypeId()]*constant::physicoChemical::k.value()*3.0/eMass);
            }

            cloud.addNewParcel(p,celli,U,cloud.electronTypeId());
        }

    }



}

// ************************************************************************* //

