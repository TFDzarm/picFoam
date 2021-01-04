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
#include "DiagnosticsList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::DiagnosticsList<CloudType>::DiagnosticsList(CloudType& owner)
:
    PtrList<Foam::DiagnosticInfo<CloudType>>(),
    owner_(owner)
{}


template<class CloudType>
Foam::DiagnosticsList<CloudType>::DiagnosticsList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<Foam::DiagnosticInfo<CloudType>>(),
    owner_(owner)
{
    setupModels(dict,owner);
}

template<class CloudType>
void Foam::DiagnosticsList<CloudType>::setupModels(const dictionary& dict, CloudType& owner)
{
    wordList diagnostics = dict.lookup("modelList");
    if(diagnostics.empty())
        return;

    Info << "+ Selecting DiagnosticInfo:" << endl;
    this->setSize(diagnostics.size());
    forAllReverse(diagnostics, i)
    {
        const word& diagnosticInfo = diagnostics[i];

        this->set(i,DiagnosticInfo<CloudType>::New(
                            diagnosticInfo,
                            dict,
                            owner
                        ));
    }
}


template<class CloudType>
Foam::DiagnosticsList<CloudType>::DiagnosticsList
(
    const DiagnosticsList<CloudType>& iml
)
:
    PtrList<Foam::DiagnosticInfo<CloudType>>(iml),
    owner_(iml.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiagnosticsList<CloudType>::~DiagnosticsList()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
const CloudType& Foam::DiagnosticsList<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::DiagnosticsList<CloudType>::owner()
{
    return owner_;
}


//Called in PICCloud::info()
template<class CloudType>
void Foam::DiagnosticsList<CloudType>::info()
{
    if(this->empty())
        return;

    //Check if at least one model should execute.
    const CloudType& cloud(this->owner());
    bool runDiagnostic = false;
    forAll(*this, i)
    {
        if(this->operator[](i).shouldExecute()) {
            runDiagnostic = true;
            break;
        }
    }
    if(runDiagnostic)
    {
        forAllConstIter(typename CloudType, cloud, iter)//Go through all particles
        {
            const typename CloudType::parcelType& p = iter();
            forAll(*this, i)//Call gatherDiagnostic for all models in the list
            {
                if(this->operator[](i).shouldExecute()/* || this->operator[](i).shouldWrite()*/)
                    this->operator[](i).gatherDiagnostic(p);
            }
        }
        //Print the info for all models
        forAll(*this, i)
        {
            if(this->operator[](i).shouldExecute()) {
                this->operator[](i).info();
                Info << "----" << endl;//flush
            }
        }
    }

}

/*
template<class CloudType>
void Foam::DiagnosticsList<CloudType>::write()
{
    forAll(*this, i)
    {
        if(this->operator[](i).shouldWrite())
            this->operator[](i).write();
    }
}*/


// ************************************************************************* //

