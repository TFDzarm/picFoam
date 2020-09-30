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

Class
    Foam::TestReaction

Description
    Test event stop simulation when the associated patch is hit.

SourceFiles
    TestReaction.C

\*---------------------------------------------------------------------------*/

#ifndef TestReaction_H
#define TestReaction_H

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TestReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TestReaction :
        public BoundaryEvent<CloudType>
{
public:

    //- Runtime type information
    TypeName("TestReaction");


    // Constructors

    TestReaction(const TestReaction<CloudType>& im);

    //- Construct from dictionary
        TestReaction
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~TestReaction();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new TestReaction<CloudType>(*this)
        );
    }

    // Member Functions

        //- Stop run if collision occurs
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
        {
            CloudType& cloud(this->owner());
            const polyMesh& mesh(cloud.mesh());
            const polyPatch& patch = mesh.boundaryMesh()[p.patch()];

            Info << "[TestReaction]\n    Patch: " << patch.name() << "\n    Vel: " << p.U() << " |" << mag(p.U()) << "|\n    preVel: " << td.boundaryVelocity() << " |" << mag(td.boundaryVelocity()) <<
                 "|\n" << endl;

            FatalErrorInFunction << "\n\n END SIMULATION \n\n" << exit(FatalError);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TestReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
