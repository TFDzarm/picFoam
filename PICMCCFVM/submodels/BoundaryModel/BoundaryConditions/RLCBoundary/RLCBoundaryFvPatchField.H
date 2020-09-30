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
    Foam::RLCBoundaryFvPatchField

Description
    Boundary condition for the RLC BoundaryModel, coefficents are set by the circuitRLC model.

SourceFiles
    RLCBoundaryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef RLCBoundaryFvPatchField_H
#define RLCBoundaryFvPatchField_H

#include "fixedGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RLCBoundaryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class RLCBoundaryFvPatchField
:
    public fvPatchField<scalar>
{

    // Private data

        //- Coefficients set by the BoundaryModel
        Field<scalar> coeff1_;
        Field<scalar> coeff2_;

        //- Charges d^(n+1)Q/dt^(n+1) calculated in the BoundaryModel
        scalar Q_;
        scalar Q1_;
        scalar Q2_;
        scalar Q3_;

        //- Surface charge
        scalar sigma_;
public:

    //- Runtime type information
    TypeName("RLCBoundary");


    // Constructors

        //- Construct from patch and internal field
        RLCBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RLCBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given RLCBoundaryFvPatchField onto a new patch
        RLCBoundaryFvPatchField
        (
            const RLCBoundaryFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        RLCBoundaryFvPatchField
        (
            const RLCBoundaryFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new RLCBoundaryFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RLCBoundaryFvPatchField
        (
            const RLCBoundaryFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new RLCBoundaryFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }

        // Return/Modify coefficient fields

            Field<scalar>& coeff1(){ return coeff1_; }
            const Field<scalar>& coeff1() const{ return coeff1_;}
            Field<scalar>& coeff2(){ return coeff2_; }
            const Field<scalar>& coeff2() const{ return coeff2_;}

        // Return/Modify charges
            scalar& Q(){ return Q_; }
            const scalar& Q() const{ return Q_; }

            scalar& Q1(){ return Q1_; }
            const scalar& Q1() const{ return Q1_; }

            scalar& Q2(){ return Q2_; }
            const scalar& Q2() const{ return Q2_; }

            scalar& Q3(){ return Q3_; }
            const scalar& Q3() const{ return Q3_; }

            scalar& sigma(){ return sigma_; }
            const scalar& sigma() const{ return sigma_; }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

        //- Update the coefficients associated with the patch field
        //virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
