/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "picFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "picCloud.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(picFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        picFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::picFields::picFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    forceWrite_(false),
    typeIdList_()
{
    forceWrite_.readIfPresent("forceWrite", dict);
    dict.readIfPresent("typeIdList",typeIdList_);
    
    //Insert a ":" at the front
    forAll(typeIdList_, i)
    {
        word& typeId = typeIdList_[i];
        typeId.insert(typeId.begin(),':');
    }

    //Append empty string to list to calculate global fields
    typeIdList_.append("");

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::picFields::~picFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::picFields::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::picFields::execute()
{
    return true;
}


bool Foam::functionObjects::picFields::write()
{
    forAll(typeIdList_, j)
    {
        word ext = typeIdList_[j];

        const volScalarField& rhoNMean = obr_.lookupObject<volScalarField>
        (
            "rhoN"+ext+"Mean"
        );

        const volScalarField& rhoMMean = obr_.lookupObject<volScalarField>
        (
            "rhoM"+ext+"Mean"
        );

        const volVectorField& momentumMean = obr_.lookupObject<volVectorField>
        (
            "momentum"+ext+"Mean"
        );

        const volScalarField& linearKEMean = obr_.lookupObject<volScalarField>
        (
            "linearKE"+ext+"Mean"
        );

        const volVectorField& fDMean = obr_.lookupObject<volVectorField>
        (
            "fD"+ext+"Mean"
        );

        if (min(mag(rhoNMean)).value() > vSmall)
        {
            Info<< "Calculating picFields" << ext << '.' << endl;

            //Construct and calculate
            Info<< "    Calculating UMean field." << endl;
            volVectorField UMean
            (
                IOobject
                (
                    "UMean"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),
                momentumMean/rhoMMean
            );

            Info<< "    Calculating translationalT field." << endl;
            volScalarField translationalT
            (
                IOobject
                (
                    "translationalT"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),

                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
               *(linearKEMean - 0.5*rhoMMean*(UMean & UMean))
            );

            Info<< "    Calculating pressure field." << endl;
            volScalarField p
            (
                IOobject
                (
                    "p"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),
                physicoChemical::k.value()*rhoNMean*translationalT
            );

            //Calculate pressure at the boundaries
            volScalarField::Boundary& pBf = p.boundaryFieldRef();

            forAll(mesh_.boundaryMesh(), i)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[i];

                if (isA<wallPolyPatch>(patch))
                {
                    pBf[i] =
                        fDMean.boundaryField()[i]
                      & (patch.faceAreas()/patch.magFaceAreas());
                }
            }

            //Print info
            Info<< "    mag(UMean) max/min : "
                << max(mag(UMean)).value() << " "
                << min(mag(UMean)).value() << endl;

            Info<< "    translationalT max/min : "
                << max(translationalT).value() << " "
                << min(translationalT).value() << endl;

            Info<< "    p max/min : "
                << max(p).value() << " "
                << min(p).value() << endl;

            //Write fields
            UMean.write();

            translationalT.write();

            p.write();

            Info<< "picFields" << ext << " written." << nl << endl;
        }
        else if (forceWrite_)
        {
            Info<< "Force calculation of picFields" << ext << '.' << endl;

            //Construct zero fields
            Info<< "    Creating UMean field." << endl;
            volVectorField UMean
            (
                IOobject
                (
                    "UMean"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),
                mesh_,
                dimensionedVector("zero",  dimensionSet(0, 1, -1, 0, 0, 0, 0), Zero),
                calculatedFvPatchScalarField::typeName
            );

            Info<< "    Creating translationalT field." << endl;
            volScalarField translationalT
            (
                IOobject
                (
                    "translationalT"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),
                mesh_,
                dimensionedScalar("zero",  dimensionSet(0, 0, 0, 1, 0, 0, 0), Zero),
                calculatedFvPatchScalarField::typeName
            );

            Info<< "    Creating pressure field." << endl;
            volScalarField p
            (
                IOobject
                (
                    "p"+ext,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ
                ),
                mesh_,
                dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero),
                calculatedFvPatchScalarField::typeName
            );

            //Calculate
            Info<< "    Calculating picFields." << endl;
            forAll(mesh_.cells(),i)
            {
                if(rhoNMean[i] > vSmall) {
                    UMean[i] = momentumMean[i]/rhoMMean[i];
                    if(linearKEMean[i] > vSmall)
                    {
                        translationalT[i] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean[i])*(linearKEMean[i] - 0.5*rhoMMean[i]*(UMean[i] & UMean[i]));
                        p[i] = physicoChemical::k.value()*rhoNMean[i]*translationalT[i];
                    }
                }
            }

            //Calculate pressure at the boundaries
            volScalarField::Boundary& pBf = p.boundaryFieldRef();

            forAll(mesh_.boundaryMesh(), i)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[i];

                if (isA<wallPolyPatch>(patch))
                {
                    pBf[i] =
                        fDMean.boundaryField()[i]
                      & (patch.faceAreas()/patch.magFaceAreas());
                }
            }

            //Print info
            Info<< "    mag(UMean) max/min : "
                << max(mag(UMean)).value() << " "
                << min(mag(UMean)).value() << endl;

            Info<< "    translationalT max/min : "
                << max(translationalT).value() << " "
                << min(translationalT).value() << endl;

            Info<< "    p max/min : "
                << max(p).value() << " "
                << min(p).value() << endl;


            //Write fields
            UMean.write();

            translationalT.write();

            p.write();

            Info<< "picFields" << ext << " written." << nl << endl;
        }
        else
        {
            Info<< "Small value (" << min(mag(rhoNMean))
                << ") found in rhoNMean field. "
                << "Not calculating picFields to avoid division by zero."
                << endl;
        }
    }
    return true;
}


// ************************************************************************* //
