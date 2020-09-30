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
    Foam::PerezIonizationModel

Description
    Ionization of charged particles in coulomb collision
    F. Pérez et al.: Improved modeling of relativistic collisions and collisional ionization in particle-in-cell
    codes (https://doi.org/10.1063/1.4742167)

    Code structure inspired by:
    http://www.maisondelasimulation.fr/projects/Smilei/html/index.html

    References for binding energies:
        A. Kramida, Yu. Ralchenko, J. Reader, and and NIST ASD Team.
        NIST Atomic Spectra Database (ver. 5.7.1),
        Available: https://physics.nist.gov/asd [2020, August 20].
        National Institute of Standards and Technology, Gaithersburg, MD., 2019.

        D. Thomas. Binding energies of electrons in atoms from h(z=1) to lw(z=103).
        Available: http://www.chembio.uoguelph.ca/educmat/atomdata/bindener/elecbind.htm [2020, August 20]
        , 1997.

SourceFiles
    PerezIonizationModel.C

\*---------------------------------------------------------------------------*/

#ifndef PerezIonizationModel_H
#define PerezIonizationModel_H

#include "IonizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PerezIonizationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PerezIonizationModel :
        public IonizationModel<CloudType>
{
    //- Electron typeId
    label electronTypeId_;

    //- Ionization data read from dictionary
    struct IonizationData
    {
        IonizationData();
        ~IonizationData() {}

        //- Does this species ionize?
        bool isActive(){ return ionSpecies_ != -1; }

        //- Calculate binding energy
        scalar bindingEnergy(label Zstar, label k);

        //- Calculate cross section
        void calculateIonizationCrossSection();

        //- Ion species
        label ionSpecies_;

        //- Neutral species... not used!
        label neutralSpecies_;

        //- Maximum ionization state we want
        label ionizationLimit_;

        //- Atomic number of the species
        label atomicNumber_;

        //- Binding energy of orbital
        scalarList bindingEneries_;

        //- Ionization energy needed for state Zstar
        scalarList ionizationEneries_;

        //- Mean ionization cross section for state Zstar and collision kinetic energy
        List<List<scalar>> ionizationCrossSection_;

        //- Mean energy transferd during ionization to the new electron
        List<List<scalar>> ionizationTransferredEnergy_;

        //- Mean lost during ionization of the colliding electron
        List<List<scalar>> ionizationLostEnergy_;

        //- Name of the dictonary data was read from
        word dictonaryName_;
    };

    //- List of processed data
    List<IonizationData> ionizationData_;
public:

    //- Runtime type information
    TypeName("Perez");


    // Constructors

    //- Construct from dictionary
        PerezIonizationModel
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PerezIonizationModel();


    // Member Functions

        //- Initialize
        virtual void initIonization(label pTypeId, label qTypeId, scalar nP, scalar nQ);

        //- Prepare density ratios before ionizations
        virtual void prepareNumberDensities(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ);

        //- Perpare ionization of the parcel
        virtual void prepareIonization(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ);

        //- Ionize the parcel
        virtual bool ionize();

private:
        //- Is a ionization possible?
        bool performIonization_;

        //- Stored electron
        typename CloudType::parcelType* electron_;

        //- Stored ion
        typename CloudType::parcelType* parcel_;

        //- Number densities of the species
        scalar np_;
        scalar ne_;

        //- Hybrid species
        scalar nei_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PerezIonizationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
