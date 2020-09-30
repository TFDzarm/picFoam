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

Application
    picInitialise

Description
    Initialise a case for picFoam by reading the initialisation dictionary
    system/picInitialise.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"
#include "Random.H"
#include "barycentric.H"
#include "picCloud.H"
#include "InitializationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




int main(int argc, char *argv[])
{
    argList::addOption("dict","file (def: picInitialiseDict)", "Dictionary to be read in particle initialization");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    word initDict;
    args.optionReadIfPresent<word>("dict", initDict, "picInitialiseDict");

    IOdictionary picInitialiseDict
    (
        IOobject
        (
            initDict,
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Initialising pic for Time = " << runTime.timeName() << nl << endl;

    picCloud pic("pic", mesh, picInitialiseDict);

    word model = picInitialiseDict.lookup("InitializationModel");
    autoPtr<InitializationModel<picCloud>> initializationModel(InitializationModel<picCloud>::New
                                                               (
                                                                   model,
                                                                   picInitialiseDict,
                                                                   pic
                                                               ));

    //init
    initializationModel->readDefaultSettings();

    Info<< nl << "Initialising particles..." << endl;
    initializationModel->initialiseParticles();

    label totalMolecules = pic.size();

    if (Pstream::parRun())
    {
        reduce(totalMolecules, sumOp<label>());
    }

    Info<< "    Total number of parcels added: " << totalMolecules
        << nl << endl;

    //Build this because we set N:<species> there
    pic.buildCellOccupancy();



    Info << "[Statistics]" << endl;
    const List<List<DynamicList<picCloud::parcelType*>>>& sortedCellOccupancy = pic.sortedCellOccupancy();
    scalarField pCaverage(pic.typeIdList().size(),0.0);
    List<label> pCmin(pic.typeIdList().size(),labelMax);
    List<label> pCmax(pic.typeIdList().size(),0);

    forAll(sortedCellOccupancy,celli)
    {
        forAll(pic.typeIdList(),id)
        {
            label size = sortedCellOccupancy[celli][id].size();
            pCaverage[id] += size;
            if(size < pCmin[id])
                pCmin[id] = size;
            if(size > pCmax[id])
                pCmax[id] = size;
        }
    }
    pCaverage /= sortedCellOccupancy.size();
    Info << "Cell occupancy:" << nl;
    forAll(pic.typeIdList(),id)
    {
        word species = pic.typeIdList()[id];
        Info << "    [" << species << "] avg: " << pCaverage[id] << ", max: " << pCmax[id] << ", min: " << pCmin[id] << nl;
    }
    Info << endl;


    //Calculate temperature
    {
        Field<vector> vDrift(pic.typeIdList().size(),Zero);
        Field<scalar> vSqr(pic.typeIdList().size(),0.0);
        Field<scalar> nParticleTypes(pic.typeIdList().size(),0.0);
        forAllIter(picCloud, pic, iter)
        {
            vDrift[iter().typeId()] += iter().U()*iter().nParticle();
            vSqr[iter().typeId()] += (iter().U()& iter().U())*iter().nParticle();
            nParticleTypes[iter().typeId()] += iter().nParticle();
        }
        if (Pstream::parRun())
        {
            Pstream::listCombineGather(vDrift, maxEqOp<vector>());
            Pstream::listCombineScatter(vDrift);
            Pstream::listCombineGather(vSqr, maxEqOp<scalar>());
            Pstream::listCombineScatter(vSqr);
            Pstream::listCombineGather(nParticleTypes, maxEqOp<scalar>());
            Pstream::listCombineScatter(nParticleTypes);
        }

        Info << "Drift velocity: " << nl;
        forAll(pic.typeIdList(), id){
            if(nParticleTypes[id] == 0.0)
                continue;

            vector vD = vDrift[id]/nParticleTypes[id];

            Info << "    [" << pic.typeIdList()[id] << "]: " << vD << "  | " << mag(vD) << " m/s"<< nl;
        }
        Info << endl;

        Info << "Calculated temperatures: " << nl;
        forAll(pic.typeIdList(), id){
            if(nParticleTypes[id] == 0.0)
                continue;

            vector vD = vDrift[id]/nParticleTypes[id];
            scalar T = pic.constProps(id).mass()/(3.0*constant::physicoChemical::k.value())*(vSqr[id]/nParticleTypes[id] - (vD&vD));

            Info << "    [" << pic.typeIdList()[id] << "]: " << T << " K | " << T*constant::physicoChemical::k.value()/constant::electromagnetic::e.value() << " eV"<< nl;
            if(initializationModel->calculateTemperatures())
                initializationModel->temperatures()[id] = T;
        }
        Info << endl;
    }

    //Calculate number densities
    {
        Field<scalar> numberDensities(pic.typeIdList().size(),0.0);
        forAllIter(picCloud, pic, iter)
        {
            numberDensities[iter().typeId()] += iter().nParticle();
        }

        scalar meshVolume = sum(mesh.cellVolumes());

        if (Pstream::parRun())
        {
            reduce(meshVolume, sumOp<scalar>());
            Pstream::listCombineGather(numberDensities, maxEqOp<scalar>());
            Pstream::listCombineScatter(numberDensities);
        }
        numberDensities /= meshVolume;

        if(initializationModel->calculateNumberDensities())
            initializationModel->numberDensities() = numberDensities;

        Info << "Calculated number densities: " << nl;
        forAll(pic.typeIdList(), id)
        {
            Info << "    [" << pic.typeIdList()[id] << "]: " << numberDensities[id] << " m^-3" << nl;
        }
        Info << endl;
    }

    pic.backgroundGas().initialize(initializationModel->temperatures(),initializationModel->numberDensities());
    if(readBool(picInitialiseDict.lookup("initalizeCollisionModels")))
    {
        Info << "Initializing collision models using" << nl
             << "    number density: " << (initializationModel->calculateNumberDensities() ? "calculated" : "default") << nl
             << "    temperature: " << (initializationModel->calculateTemperatures() ? "calculated" : "default") << endl;

        pic.binaryCollision().initialize(initializationModel->temperatures(),initializationModel->numberDensities());
        pic.electronNeutralCollision().initialize(initializationModel->temperatures(),initializationModel->numberDensities());
    }

    //Calculating fields
    pic.calculateFields();

    if(readBool(picInitialiseDict.lookup("solveMaxwellEquations"))) {
        Info << "Solving Maxwell equations" << endl;
        //pic.boundaryModels().preUpdate_Boundary();
        pic.maxwellSolver().solveFields();
        //pic.boundaryModels().postUpdate_Boundary();
    }

    if(readBool(picInitialiseDict.lookup("initalizeLeapFrog")))
    {
        Info << "Updating inital velocity" << endl;
        scalar dt = mesh.time().deltaTValue();
        forAllIter(picCloud, pic, iter)
        {
            picCloud::parcelType& p = iter();
            pic.particlePusher().updateVelocity(p,-0.5*dt);
        }
    }


//Info...
    Info << "Mesh volume:" << "\n    " << sum(mesh.cellVolumes()) << " m^3" << endl;



    Field<scalar> charge(pic.typeIdList().size(),0.0);
    Field<scalar> N(pic.typeIdList().size(),0.0);
    forAllIter(picCloud, pic, iter)
    {
        charge[iter().typeId()] += iter().charge()*iter().nParticle();
        N[iter().typeId()] += iter().nParticle();
    }
    forAll(pic.typeIdList(),i)
    {
        if(N[i] > 0.0)
            charge[i] /= N[i];
    }

    scalar debyeLength = 0.0;
    forAll(pic.chargedSpecies(),i)
    {
        label typeId = pic.chargedSpecies()[i];
        scalar T = initializationModel->temperatures()[typeId];
        scalar n = initializationModel->numberDensities()[typeId];
        scalar e = charge[typeId];
        if(n <= 0.0 || T <= 0.0)
            continue;

        debyeLength += n*e*e/constant::electromagnetic::epsilon0.value() * 1.0/(T*constant::physicoChemical::k.value());
    }
    if(debyeLength > 0.0)
        debyeLength = ::sqrt(1.0/debyeLength);

    Info << "Average cell size:" << "\n    " << pic.avgCellLengthScale() << " m" << nl
         << "Debye length:" << "\n    " << debyeLength << " m" << endl;

    if(pic.avgCellLengthScale() > debyeLength)
        Info << "    WARNING!!! Average cell size is larger than the debye length" << endl;

    if(pic.electronTypeId() > -1)
    {
        const picCloud::parcelType::constantProperties& cP = pic.constProps(pic.electronTypeId());

        scalar plasmaFreq = 0.0;
        scalar n = initializationModel->numberDensities()[pic.electronTypeId()];
        if(n > 0.0)
        {
            plasmaFreq = ::sqrt(n*cP.charge()*cP.charge()/constant::electromagnetic::epsilon0.value()/cP.mass());
        }
        Info << "Plasma frequency x deltaT:" << "\n    " << plasmaFreq*runTime.deltaTValue() << endl;
        if(plasmaFreq*runTime.deltaTValue() > 2.0)
            Info << "    WARNING!!! Plasma frequency times deltaT is larger than 2" << endl;
    }

    //Diagnostics
    pic.info();

///////////////////////////////////////////////////////////////////////////////////////////////////////

    IOstream::defaultPrecision(15);

    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing picCloud."
            << nl << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;


    return 0;
}


// ************************************************************************* //
