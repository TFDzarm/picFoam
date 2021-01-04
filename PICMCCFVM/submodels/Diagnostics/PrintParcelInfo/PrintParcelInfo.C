/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "PrintParcelInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PrintParcelInfo<CloudType>::PrintParcelInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    printPosition_(false),//Default values...
    printVelocity_(false),
    printSums_(false),
    sumPosition_(),
    sumVelocity_()
{
    //Which information do we want to print?
    printPosition_.readIfPresent("position",this->coeffDict());
    printVelocity_.readIfPresent("velocity",this->coeffDict());
    printSums_.readIfPresent("sums",this->coeffDict());
}

template<class CloudType>
Foam::PrintParcelInfo<CloudType>::PrintParcelInfo(const PrintParcelInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    printPosition_(im.printPosition_),
    printVelocity_(im.printVelocity_),
    printSums_(im.printSums_),
    sumPosition_(im.sumPosition_),
    sumVelocity_(im.sumVelocity_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PrintParcelInfo<CloudType>::~PrintParcelInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PrintParcelInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const CloudType& cloud(this->owner());

    //Directly print the info for every particle
    if(printPosition_)
        Info << "(" << p.origId() << ") [position]: " << cloud.mesh().time().timeName() << "," << p.position().x() << "," << p.position().y() << "," << p.position().z() << "\n";
    if(printVelocity_) {
        vector corrU = cloud.particlePusher().correctedVelocity(p, 0.5*cloud.mesh().time().deltaTValue());
        Info << "(" << p.origId() << ") [velocity]: " << cloud.mesh().time().timeName() << "," << p.U().x() << "," << p.U().y() << "," << p.U().z() << nl
             << "                " << corrU.x() << "," << corrU.y() << "," << corrU.z()
             << "\n";

    }

    //Add up the sums
    if(printSums_) {
        sumPosition_ += p.position();
        sumVelocity_ += p.U();
    }

}

//- Print info
template<class CloudType>
void Foam::PrintParcelInfo<CloudType>::info()
{
    //Print the sums
    if(printSums_) {
        reduce(sumPosition_, sumOp<vector>());
        reduce(sumVelocity_, sumOp<vector>());
        Info << "Parcel position sum: " << sumPosition_.x() << "," << sumPosition_.y() << "," << sumPosition_.z() << "\n";
        Info << "Parcel velocity sum: " << sumVelocity_.x() << "," << sumVelocity_.y() << "," << sumVelocity_.z() << endl;
        sumPosition_ = Zero;
        sumVelocity_ = Zero;
    }

}

// ************************************************************************* //
