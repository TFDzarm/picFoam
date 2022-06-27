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

#include "PerezIonizationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PerezIonizationModel<CloudType>::PerezIonizationModel
(
    const dictionary& dict,
    CloudType& cloud
)
:
    IonizationModel<CloudType>(dict,cloud,"PerezIonization"),
    electronTypeId_(-1),
    ionizationData_(this->owner().typeIdList().size()),
    performIonization_(false),
    electron_(nullptr),
    parcel_(nullptr),
    np_(0.0),
    ne_(0.0),
    nei_(0.0)

{
    dictionary relations
    (
        cloud.particleProperties().subDict("moleculeProperties").subDict("SpeciesRelations")
    );

    //We need to know which typeId is the electron species
    word electronSpecies = relations.lookup("electronTypeId");
    electronTypeId_ = findIndex(this->owner().typeIdList(),electronSpecies);
    if(electronTypeId_ == -1)
        FatalErrorInFunction << "electronSpecies " << electronSpecies << " is not defined" << abort(FatalError);

    //Read ionization data
    forAllConstIter(IDLList<entry>, this->coeffDict(), iter)//for every subdict...
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();
            word ionSpecies = subDict.lookup("ionSpecies");
            label idIon = findIndex(this->owner().typeIdList(),ionSpecies);
            if(idIon == -1)
                FatalErrorInFunction << "ionSpecies " << ionSpecies << " is not defined" << abort(FatalError);

            label atomicNumber = readLabel(subDict.lookup("atomicNumber"));
            label ionizationLimit = subDict.lookupOrDefault<label>("ionizationLimit",atomicNumber);//Limit for the ionization, default is the atomic number
            if(ionizationLimit > atomicNumber)
                ionizationLimit = atomicNumber;

            if(ionizationLimit <= 1)
            {
                Info << "|->    IonizationModel: " << iter().keyword() << ".ionizationLimit is " << ionizationLimit << ", no ionization for ion species " << ionSpecies << endl;
                if(ionizationLimit < 1)//Go to the next subdict...
                    continue;
            }

            scalarList bindingEnergy = subDict.lookup("bindingEnergy");
            scalarList ionizationEnergy = subDict.lookup("ionizationEnergy");

            if(bindingEnergy.size() != atomicNumber)
                FatalErrorInFunction << "Number of entries in list \"bindingEnergy\" (" << bindingEnergy.size() << ") does not match the atomic number " << atomicNumber << "." << abort(FatalError);

            if(ionizationEnergy.size() != atomicNumber)
                FatalErrorInFunction << "Number of entries in list \"ionizationEnergy\" (" << ionizationEnergy.size() << ") does not match the atomic number " << atomicNumber << "." << abort(FatalError);

            ionizationData_[idIon].dictonaryName_ = iter().keyword();
            ionizationData_[idIon].ionSpecies_ = idIon;
            ionizationData_[idIon].atomicNumber_ = atomicNumber;
            ionizationData_[idIon].ionizationLimit_ = ionizationLimit;
            ionizationData_[idIon].bindingEneries_ = bindingEnergy;
            ionizationData_[idIon].ionizationEneries_ = ionizationEnergy;
            ionizationData_[idIon].ionizationCrossSection_.resize(atomicNumber);
            ionizationData_[idIon].ionizationTransferredEnergy_.resize(atomicNumber);
            ionizationData_[idIon].ionizationLostEnergy_.resize(atomicNumber);

            ionizationData_[idIon].calculateIonizationCrossSection();//Create the cross section table for this ion species

            word neutralSpecies = subDict.lookupOrDefault<word>("neutralSpecies","");
            label idNeutal = findIndex(this->owner().typeIdList(),neutralSpecies);
            if(idNeutal == -1) {
                if(ionizationLimit == 1)
                    Info << "IonizationModel: WARNING " << iter().keyword() << ".ionizationLimit is 1 and " << iter().keyword() << ".neutralSpecies is not provided" << endl;
            }
            else {
                ionizationData_[idIon].neutralSpecies_ = idNeutal;
                ionizationData_[idNeutal] = ionizationData_[idIon];
            }


        }
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PerezIonizationModel<CloudType>::~PerezIonizationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PerezIonizationModel<CloudType>::initIonization(label pTypeId, label qTypeId, scalar nP, scalar nQ)
{
    //Set up in a way that np_ is always the value of ion species
    ne_ = 0.0;
    nei_ = 0.0;
    if(pTypeId == electronTypeId_)
        np_ = nQ;
    else if(qTypeId == electronTypeId_)
        np_ = nP;
}

template<class CloudType>
void Foam::PerezIonizationModel<CloudType>::prepareNumberDensities(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ)
{
    //Paper: Only count pairs that can collide. Here: We do not need to count only the pairs that can collide.
    // ne and nei do not include newly created electrons in our case!
    label pType = pP.typeId();
    label qType = pQ.typeId();

    scalar E;
    scalar We, Wi;
    label Zstar, parcelTypeId;

    if(pType == electronTypeId_) {
        Zstar = pQ.Zstar();
        parcelTypeId = qType;

        We = pP.nParticle();
        Wi = pQ.nParticle();

        scalar mass = this->owner().constProps(pType).mass();
        scalar gamma = 1.0/sqrt(1.0-sqr(mag(pP.U())/constant::universal::c.value()));
        E = (gamma-1.0)*mass*constant::universal::c.value()*constant::universal::c.value()/constant::electromagnetic::e.value();//Energy of the electron

    }
    else if(qType == electronTypeId_) {
        Zstar = pP.Zstar();
        parcelTypeId = pType;

        We = pQ.nParticle();
        Wi = pP.nParticle();

        scalar mass = this->owner().constProps(qType).mass();
        scalar gamma = 1.0/sqrt(1.0-sqr(mag(pQ.U())/constant::universal::c.value()));
        E = (gamma-1.0)*mass*constant::universal::c.value()*constant::universal::c.value()/constant::electromagnetic::e.value();//Energy of the electron
    }
    else {
        return;
    }

    if( Zstar >= ionizationData_[parcelTypeId].ionizationLimit_)
        return;

    //Find the entry in the interpolation table
    scalar x = 6.142165*::log( E );

    if ( x < 0.0)
        x = 0.0;
    else if ( x > 99.0 )
        x = 99.0;

    scalar cs = ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][label(x)];//Get the cross section

    if( cs>0.0 ) {
        ne_ += We;
        nei_ += We<Wi ? We : Wi;
    }
}

template<class CloudType>
void Foam::PerezIonizationModel<CloudType>::prepareIonization(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ)
{
    //Set up internal pointers
    label pType = pP.typeId();
    label qType = pQ.typeId();
    if(pType == electronTypeId_ && ionizationData_[qType].isActive()) {
        electron_ = &pP;
        parcel_ = &pQ;
    }
    else if(qType == electronTypeId_ && ionizationData_[pType].isActive()) {
        electron_ = &pQ;
        parcel_ = &pP;
    }
    else {
        performIonization_ = false;
        return;
    }

    label parcelTypeId = parcel_->typeId();

    //Check if the ionization is allowed
    if(parcel_->Zstar() >= ionizationData_[parcelTypeId].ionizationLimit_) {
        performIonization_ = false;
        return;
    }
    performIonization_ = true;
}

template<class CloudType>
bool Foam::PerezIonizationModel<CloudType>::ionize()
{
    if(!performIonization_)//Can we try the ionize?
        return false;

    //Calculate particle properties
    label parcelTypeId = parcel_->typeId();
    label atomicNumber = ionizationData_[parcelTypeId].atomicNumber_;
    label Zstar = parcel_->Zstar();
    label initialZstar = Zstar;
    scalar gammaE = 1.0/::sqrt(1.0-::sqr(::mag(electron_->U())/constant::universal::c.value()));
    scalar gammaP = 1.0/::sqrt(1.0-::sqr(::mag(parcel_->U())/constant::universal::c.value()));
    scalar massE = this->owner().constProps(electron_->typeId()).mass();
    scalar massP = this->owner().constProps(parcel_->typeId()).mass();
    scalar deltaT = this->owner().mesh().time().deltaTValue();

    scalar wE = electron_->nParticle();
    scalar wP = parcel_->nParticle();
    scalar wEwP = wE/wP;
    scalar wPwE = 1./wEwP;
    // begin

    //Code inspired by Smilie (see header description): go from lab frame to ion frame
    scalar gamma_s = gammaP*gammaE*(1.0-mag(parcel_->U())*mag(electron_->U())/(constant::universal::c.value()*constant::universal::c.value()));
    if(gamma_s <= 1.0)//rare but can happen ...energy is low ... so ionization unlikely
    {//important reset vars
        parcel_ = nullptr;
        electron_ = nullptr;
        return false;
    }

    scalar ionCoeff_s = constant::universal::c.value()*::sqrt(gamma_s*gamma_s-1.0)/(gammaP)*(np_*ne_/nei_)*deltaT;
    scalar eVkinEnergy_s = (gamma_s-1.0)*massE*constant::universal::c.value()*constant::universal::c.value()/constant::electromagnetic::e.value();

    scalar rndU = this->owner().rndGen().scalar01();
    label ionizationMax = ionizationData_[parcelTypeId].ionizationLimit_-Zstar-1;
    scalar s, d, crossSection, transferredEnergy, lostEnergy;
    List<scalar> ionERate(atomicNumber,0.0);
    List<scalar> ionERate_inv(atomicNumber,0.0);
    List<scalar> ionProp(atomicNumber,0.0);
    scalar com_ionProp = 0.0;
    label l;

    //Get the cross section, energy transferred to the new electron and energy lost by the colliding electron
    for(int k = 0; k <= ionizationMax; k++)
    {
        s = 6.142165*::log(eVkinEnergy_s);
        if( s < 0.0 ) break;
        if( s < 99.0 )//interpolate
        {
           l = label(s);
           d = s-scalar(l);
           crossSection = (ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][l+1]-ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][l])*d + ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][l];
           transferredEnergy = (ionizationData_[parcelTypeId].ionizationTransferredEnergy_[Zstar][l+1]-ionizationData_[parcelTypeId].ionizationTransferredEnergy_[Zstar][l])*d + ionizationData_[parcelTypeId].ionizationTransferredEnergy_[Zstar][l];
           lostEnergy = (ionizationData_[parcelTypeId].ionizationLostEnergy_[Zstar][l+1]-ionizationData_[parcelTypeId].ionizationLostEnergy_[Zstar][l])*d + ionizationData_[parcelTypeId].ionizationLostEnergy_[Zstar][l];

        }
        else//extrapolate
        {
            d = s - 99.0;
            crossSection = (ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][99]-ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][98])*d + ionizationData_[parcelTypeId].ionizationCrossSection_[Zstar][99];
            transferredEnergy = ionizationData_[parcelTypeId].ionizationTransferredEnergy_[Zstar][99];
            lostEnergy = ionizationData_[parcelTypeId].ionizationLostEnergy_[Zstar][99];
        }
        if(crossSection <= 0.0)//Cannot ionize
           break;

        if(lostEnergy > gamma_s-1.0)
            break;

        ionERate[k] = ionCoeff_s*crossSection/gammaE;
        ionERate_inv[k] = 1.0/ionERate[k];
        ionProp[k] = ::exp(-ionERate[k]);

        scalar cp;
        //Nuter et al. 2011 https://doi.org/10.1063/1.3559494
        if(k==0)
            com_ionProp = ionProp[k];
        else if(k<ionizationMax) {
            int x,b;
            for(x = 0; x < k; x++)
            {
                cp = 1.0 - ionERate[k]*ionERate_inv[x];

                for( b=0  ; b<x; b++ )
                    cp *= 1.0 - ionERate[x]*ionERate_inv[b];
                for( b=x+1; b<k; b++ )
                    cp *= 1.0 - ionERate[x]*ionERate_inv[b];

                com_ionProp += (ionProp[k]-ionProp[x])/cp;
            }
        }
        else {
            int x,b;
            for(x = 0; x < k; x++)
            {
                cp = 1.0 - ionERate[k]*ionERate_inv[x];

                for( b=0  ; b<x; b++ )
                    cp *= 1.0 - ionERate[x]*ionERate_inv[b];
                for( b=x+1; b<k; b++ )
                    cp *= 1.0 - ionERate[x]*ionERate_inv[b];

                com_ionProp += (1.0-ionProp[k]+ionERate[k]*ionERate_inv[x]*(ionProp[x]-1.0))/cp;

            }
        }

        if( rndU < com_ionProp )
            break;

        vector pOld = electron_->U() * gammaE * massE;
        vector pI = parcel_->U() * gammaP * massP;
        scalar p2 = gamma_s*gamma_s-1.0;

        scalar rndU2 = this->owner().rndGen().scalar01();
        if(rndU2 < wEwP)
        {
            //Increase charge
            parcel_->chargeModifier()++;

            //New electron + momentum transfer
            typename CloudType::parcelType* newElectron = this->owner().addNewParcel(electron_->position(),electron_->cell(),electron_->U(),electronTypeId_);//copy electron
            newElectron->nParticle() = wP;//Ion weight

            scalar gammaW = transferredEnergy + 1.0;//Ekin = (gamma-1)*mc^2 -> gamma = Ekin/mc^2 +1 .... here transferredEnergy is already dimensionless
            scalar alpha_w = sqrt((gammaW*gammaW-1.0)/p2);

            vector pTransferNew = alpha_w*pOld+(gammaW-alpha_w*gamma_s)*pI*massE/massP;
            scalar gammaEnew = sqrt(1.0+::pow(mag(pTransferNew)/(massE*constant::universal::c.value()),2));
            newElectron->U() = pTransferNew/(massE*gammaEnew);
        }
        if(rndU2 < wPwE)
        {
            //Change of momentum of the colliding electron
            scalar alpha_e = sqrt((pow(gamma_s-lostEnergy,2)-1.0)/p2);
            vector pLostEnergyNew = alpha_e * pOld + ((1.0-alpha_e)*gamma_s-lostEnergy)*pI*massE/massP;
            scalar gammaLostEnew = sqrt(1.0+pow(mag(pLostEnergyNew)/(massE*constant::universal::c.value()),2));
            electron_->U() = pLostEnergyNew/(massE*gammaLostEnew);

            gammaE += ((1.0-alpha_e)*gamma_s-lostEnergy)*gammaP;

            gamma_s -= lostEnergy;//Decrease gamma_s
        }

        Zstar++;


    }

    if( initialZstar == Zstar )
    {//important reset vars
        parcel_ = nullptr;
        electron_ = nullptr;
        return false;
    }

    parcel_->charge() = this->owner().constProps(parcelTypeId).charge()*parcel_->chargeModifier();//update charge of the ion!!!


    // end reset vars
    parcel_ = nullptr;
    electron_ = nullptr;
    return true;
}

// ************************************************************************* //

template<class CloudType>
PerezIonizationModel<CloudType>::IonizationData::IonizationData() :
    ionSpecies_(-1),
    neutralSpecies_(-1),
    ionizationLimit_(-1),
    atomicNumber_(0),
    bindingEneries_(),
    ionizationEneries_(),
    ionizationCrossSection_(),
    ionizationTransferredEnergy_(),
    ionizationLostEnergy_(),
    dictonaryName_("")
{ }

template<class CloudType>
scalar PerezIonizationModel<CloudType>::IonizationData::bindingEnergy(label Zstar, label k)//by Carlson
{
    return ionizationEneries_[Zstar] - bindingEneries_[atomicNumber_-Zstar-1] + bindingEneries_[k];
}


template<class CloudType>
void PerezIonizationModel<CloudType>::IonizationData::calculateIonizationCrossSection()
{
    //Create the cross section table for interpolations
    namespace ca = constant::atomic;
    namespace cu = constant::universal;
    namespace cm = constant::mathematical;
    namespace ce = constant::electromagnetic;

    scalar r_e = 2.8179403227e-15;
    label tableSize = 100;
    //1eV-10MeV
    scalar tableStepSize = (tableSize-1)/::log(1e7/1);

    scalar E, Bk, eps, eps_dash, u_dash, b_dash,beta_eps2,beta_b2,beta_u2,sigma0,ediv2,A1,A2,A3,sigmak,wk,ek;
    label N;
    for(label Zstar = 0; Zstar < atomicNumber_; Zstar++)
    {
        ionizationCrossSection_[Zstar] = List<scalar>(tableSize,0.0);
        ionizationTransferredEnergy_[Zstar] = List<scalar>(tableSize,0.0);
        ionizationLostEnergy_[Zstar] = List<scalar>(tableSize,0.0);

        for(label i = 0; i < tableSize; i++)
        {
            E = ::exp( scalar(i)/ tableStepSize);
            N = 1;
            for(label k = 0; k < atomicNumber_-Zstar; k++)
            {
                Bk = bindingEnergy(Zstar,k);
                if(k < atomicNumber_-Zstar-1) {
                    if(Bk == bindingEnergy(Zstar,k+1)){
                        N++;
                        continue;
                    }
                }
                eps = E/Bk;
                //E*=ce::e.value();//in Joule -> crash rounding error
                //Bk*=ce::e.value();//in Joule
                if(eps > 1.0)
                {
                    b_dash = Bk/510997.37/*ce::e.value()/(ca::me.value()*cu::c.value()*cu::c.value())*/;
                    eps_dash = E/510997.37/*ce::e.value()/(ca::me.value()*cu::c.value()*cu::c.value())*/;
                    u_dash = b_dash;//mean kinetic energy of orbital k: unknown use Uk=Bk as in Smilie!

                    beta_eps2 = 1.0-1.0/((1.0+eps_dash)*(1.0+eps_dash)); // wrong -> beta_eps2 *= beta_eps2;
                    beta_b2 = 1.0-1.0/((1.0+b_dash)*(1.0+b_dash)); // wrong -> beta_b2 *= beta_b2;
                    beta_u2 = 1.0-1.0/((1.0+u_dash)*(1.0+u_dash)); // wrong -> beta_u2 *= beta_u2;

                    sigma0 = 2.0*cm::pi*r_e*r_e*N/(b_dash*(beta_eps2+beta_b2+beta_u2));

                    ediv2 = (1.0+eps_dash*0.5); ediv2 *= ediv2;

                    A1 = 1.0/(1.0+eps)*(1.0+2*eps_dash)/ediv2;
                    A2 = (eps-1.0)/2.0 * b_dash*b_dash/ediv2;
                    A3 = ::log(beta_eps2/(1.0-beta_eps2))-beta_eps2-::log(2.0*b_dash);

                    sigmak = sigma0*(0.5*A3*(1.0-1.0/(eps*eps))+1.0-1.0/eps+A2-A1*::log(eps));
                    wk = sigma0*(0.5*A3*(eps-1.0)*(eps-1.0)/(eps*(eps-1.0))+2.0*::log((eps+1.0)*0.5)-::log(eps)+A2*(eps-1.0)*0.25-A1*(::log(eps)-(eps+1.0)*::log(0.5*(eps+1.0))));
                    ek = sigmak+wk;

                    ionizationCrossSection_[Zstar][i] += sigmak;
                    ionizationTransferredEnergy_[Zstar][i] += wk*b_dash;//Perez uses Bk .... here use b_dash to make dimensionless
                    ionizationLostEnergy_[Zstar][i] += ek*b_dash; // see above^^

                }
                N=1;
            }
            // The transferred and lost energies are averages over the orbitals
            if( ionizationCrossSection_[Zstar][i] > 0.0 ) {
                ionizationTransferredEnergy_[Zstar][i] /= ionizationCrossSection_[Zstar][i];
                ionizationLostEnergy_[Zstar][i] /= ionizationCrossSection_[Zstar][i];
            }
        }

    }
}
