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

Class
    Foam::PICParcel

Description
    PIC parcel class

SourceFiles
    PICParcelI.H
    PICParcel.C
    PICParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef PICParcel_H
#define PICParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"
#include "PICCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class PICParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const PICParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class PICParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class PICParcel
:
    public ParcelType
{
    // Private member data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold PIC particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Viscosity index
            scalar omega_;

            //- Charge [A s] (inital charge)
            scalar charge_;

            //- Solve movement of this species?
            bool solveMovement_;

            //- Number of real atoms represented by a parcel (inital number)
            scalar nParticle_;

    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the particle mass [kg]
            inline scalar mass() const;

            //- Return const access to the hard sphere diameter [m]
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the viscosity index
            inline scalar omega() const;

            //- Return the charge
            inline scalar charge() const;

            //- Return the switch for solveMovement_
            inline bool solveMovement() const;

            //- Return the nParticle
            inline scalar nParticle() const;

    };


    //- Do not use base tracking data
    //typedef typename ParcelType::trackingData trackingData;
    //- Class used to pass kinematic tracking data to the trackToFace function
    class trackingData
    :
        public ParcelType::trackingData
    {
    public:
        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud
            ) :
                ParcelType::trackingData(cloud),
                boundaryVelocity_(zero()),
                requireResync_(false)
            {}

            inline vector boundaryVelocity() const;

            inline vector& boundaryVelocity();

            inline bool requireResync() const;

            inline bool& requireResync();

    private:
            vector boundaryVelocity_;//synced or unsynced
            bool requireResync_;
    };


protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Parcel type id
            label typeId_;

            //- Charge multiplier for constantProperties::charge (used for ionization, default = 1)
            label chargeModifier_;

            //- Charge of the parcel, is updated everytime the charge changeds (charge_ = constantProperties::charge * chargeModifier_)
            scalar charge_;

            //- Number of real atoms represented by a parcel
            scalar nParticle_;


public:

    //- Runtime type information
    TypeName("PICParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline PICParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& U,
            const label typeId
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        inline PICParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const vector& U,
            const label typeId
        );

        //- Construct from Istream
        PICParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new PICParcel<ParcelType>(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<PICParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<PICParcel<ParcelType>>
                (
                    new PICParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to charge number (Value makes only sense for ions ... for all other species always 1)
            inline label Zstar() const;

            //- Return const access to the charge multiplier
            inline label chargeModifier() const;

            //- Return const access to the charge of the parcel
            inline scalar charge() const;

            //- Return const access to the number of real the particle
            inline scalar nParticle() const;

        // Edit

            //- Return access to velocity
            inline vector& U();

            //-NEW- Return access to the charge multiplier
            inline label& chargeModifier();

            //-NEW- Return access to the charge of the parcel
            inline scalar& charge();

            //-NEW- Return access to the number of real the particle
            inline scalar& nParticle();


        // Main calculation loop

            // Tracking

                //- Prepare movement of the parcel
                template<class TrackCloudType>
                bool move
                (
                    TrackCloudType& cloud,
                    trackingData& td,
                    const scalar trackTime
                );

                //- Actually move the parcel
                template<class TrackCloudType>
                void moveForward
                (
                    TrackCloudType& cloud,
                    trackingData& td,
                    const scalar trackTime
                );

                template<class TrackCloudType>
                void syncVelocityAtBoundary
                (
                        TrackCloudType& cloud,
                        scalar fraction
                );


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackCloudType>
            bool hitPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackCloudType>
            void hitProcessorPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType&, trackingData&);

            //- Overridable function
            //template<class TrackCloudType>
            //void hitSymmetryPatch(TrackCloudType&, trackingData&);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const transformer&);

            //Do not update fields, simply reflect!
            void specularReflect();

            template<class TrackCloudType>
            void wallReflection(TrackCloudType&, trackingData&);

            template<class TrackCloudType>
            void wallAbsorption(TrackCloudType&, trackingData&);

        // I-O

            static void readFields(Cloud<PICParcel<ParcelType>>& c);

            static void writeFields(const Cloud<PICParcel<ParcelType>>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const PICParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PICParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PICParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
