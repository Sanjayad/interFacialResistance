/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    dynamicYokoiAlphaContactAngleFvPatchScalarField

Description
    Numerical Studies of the Influence of the Dynamic Contact Angle on a Droplet Impacting on a Dry Surface (2009)

SourceFiles
    dynamicYokoiAlphaContactAngleFvPatchScalarField.C

Author
    Sanjay Adhikari 2017

\*---------------------------------------------------------------------------*/

#ifndef dynamicYokoiAlphaContactAngleFvPatchScalarField_H
#define dynamicYokoiAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class dynamicYokoiAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleFvPatchScalarField
{
    // Private data

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receeding contact angle
        scalar thetaR_;
            
		//- Equilibrium contact angle (advancing or receding), as parameter
        scalar thetaE_;
		
		//- Surface constant for advancing contact angle
		scalar ka_;
		
		//- Surface constant for receding contact angle
		scalar kr_;

        //- Name of the dynamic viscosity field
        word muName_;

        //- Name of the surface tension
        word sigmaName_;		


public:

    //- Runtime type information
    TypeName("dynamicYokoiAlphaContactAngle");

    // Static data members

        //- Conversion factor for radians into degrees
        //static const scalar convertToDeg;

        //- Conversion factor for degrees into radians
        //static const scalar convertToRad;

    // Constructors

        //- Construct from patch and internal field
        dynamicYokoiAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicYokoiAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  dynamicYokoiAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        dynamicYokoiAlphaContactAngleFvPatchScalarField
        (
            const dynamicYokoiAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicYokoiAlphaContactAngleFvPatchScalarField
        (
            const dynamicYokoiAlphaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicYokoiAlphaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicYokoiAlphaContactAngleFvPatchScalarField
        (
            const dynamicYokoiAlphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicYokoiAlphaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
