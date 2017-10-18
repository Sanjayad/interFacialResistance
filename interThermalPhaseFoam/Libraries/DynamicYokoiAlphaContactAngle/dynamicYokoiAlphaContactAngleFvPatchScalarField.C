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

\*---------------------------------------------------------------------------*/

#include "dynamicYokoiAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFields.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

//const scalar dynamicYokoiAlphaContactAngleFvPatchScalarField::convertToDeg =
  //           180.0/constant::mathematical::pi;

//const scalar dynamicYokoiAlphaContactAngleFvPatchScalarField::convertToRad =
  //           constant::mathematical::pi/180.0;

//const scalar dynamicYokoiAlphaContactAngleFvPatchScalarField::theta0 = 70.0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicYokoiAlphaContactAngleFvPatchScalarField::
dynamicYokoiAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(0.0),
    thetaR_(0.0),
	thetaE_(0.0),
	ka_(1.0), //- To avoid division by zero
	kr_(1.0),
    muName_("undefined"),
    sigmaName_("undefined")
{}

dynamicYokoiAlphaContactAngleFvPatchScalarField::
dynamicYokoiAlphaContactAngleFvPatchScalarField
(
    const dynamicYokoiAlphaContactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(acpsf, p, iF, mapper),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
	thetaE_(acpsf.thetaE_),
	ka_(acpsf.ka_),
	kr_(acpsf.kr_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


dynamicYokoiAlphaContactAngleFvPatchScalarField::
dynamicYokoiAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
	thetaE_(readScalar(dict.lookup("thetaE"))),
	ka_(readScalar(dict.lookup("ka"))),
	kr_(readScalar(dict.lookup("kr"))),
    muName_(dict.lookup("muEffYokoi")),
    sigmaName_(dict.lookup("sigmaYokoi"))
{
    evaluate();
}


dynamicYokoiAlphaContactAngleFvPatchScalarField::
dynamicYokoiAlphaContactAngleFvPatchScalarField
(
    const dynamicYokoiAlphaContactAngleFvPatchScalarField& acpsf
)
:
    alphaContactAngleFvPatchScalarField(acpsf),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
	thetaE_(acpsf.thetaE_),
	ka_(acpsf.ka_),
	kr_(acpsf.kr_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


dynamicYokoiAlphaContactAngleFvPatchScalarField::
dynamicYokoiAlphaContactAngleFvPatchScalarField
(
    const dynamicYokoiAlphaContactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(acpsf, iF),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
	thetaE_(acpsf.thetaE_),
	ka_(acpsf.ka_),
	kr_(acpsf.kr_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> dynamicYokoiAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    
	//- Lookup and return the patchField of dynamic viscosity of mixture
    //     and surface tension
    if((muName_ != "muDynamic") || (sigmaName_ != "sigmaDynamic"))
    {
        FatalErrorIn
        (
            "dynamicYokoiAlphaContactAngleFvPatchScalarField"
        )   << " muEffYokoi or sigmaYokoi set inconsitently, muEffYokoi = " << muName_
            << ", sigmaYokoi = " << sigmaName_ << '.' << nl
            << "    Set both muEffYokoi and sigmaYokoi according to the "
           "definition of dynamicYokoiAlphaContactAngle"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }


    //RATTNER to change, we probably want to switch mu to muEff
    const fvPatchField<scalar>& mup =
    patch().lookupPatchField<volScalarField, scalar>(muName_);


    const fvPatchField<scalar>& sigmap =
    patch().lookupPatchField<volScalarField, scalar>(sigmaName_);

    vectorField nf = patch().nf();

    // Calculate the component of the velocity parallel to the wall
    vectorField Uwall = Up.patchInternalField() - Up; 
	//- Up.patchInternalField() - Up doesn't work for the total slip boundary condition as the velocity on the patch equals the internal field velocity
    
	Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall = nHat - (nf & nHat)*nf;

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the wall
    scalarField uwall = nWall & Uwall;

    //- Calculate local Capillary number
    scalarField Ca = mup*mag(uwall)/sigmap;
    //-    Calculate and return the value of contact angle on patch faces,
    //     a general approach: the product of Uwall and nWall is negative
    //     for advancing and positiv for receding motion.
    //     thetaDp is initialized to thetaE degrees corresponding to no wall adhesion
//	scalarField thetaDp(patch().size(), convertToRad*thetaE_);
	scalarField thetaDp(patch().size(), thetaE_);
	forAll(uwall, pfacei)
	{
	    if(uwall[pfacei] < 0.0)
	    {
//		        thetaDp[pfacei] = min ( thetaA_*convertToRad, thetaE_*convertToRad + pow( Ca[pfacei]/ka_, 1.0/3.0 ) );
//	        thetaDp[pfacei] = min ( thetaA_, thetaE_ + pow( Ca[pfacei]/ka_, 1.0/3.0 ) );
	        thetaDp[pfacei] = min ( thetaA_, thetaE_ + cbrt( Ca[pfacei]/ka_ ) );
	    }
	    else if (uwall[pfacei] > 0.0)
	    {
//		        thetaDp[pfacei] = max ( thetaR_*convertToRad, thetaE_*convertToRad - pow( Ca[pfacei]/ka_, 1.0/3.0 ) );
//	        thetaDp[pfacei] = max ( thetaR_, thetaE_ + pow( -Ca[pfacei]/ka_, 1.0/3.0 ) );
	        thetaDp[pfacei] = max ( thetaR_, thetaE_ + cbrt( -Ca[pfacei]/kr_ ) );
	    }
	}

	forAll(uwall, pfacei){		
		Info <<  uwall[pfacei] << tab << thetaDp[pfacei] << endl;
	}
/*
//////////////////// eb - Print out some data ////////////////////////
    Info << "pfacei: " << tab  << "nf: "<< tab << "Uwall: " << tab
         << "nHat: " << tab << "nWall: " << tab << "uwall: " << tab
         << "mup: " << tab<< "sigmap: " << tab << "Ca: " << tab << "thetaD: " << endl;

    forAll(uwall, pfacei)
    {
        Info << "[" << pfacei << "]" << tab << nf[pfacei] << tab << Uwall[pfacei] << tab
             << nHat[pfacei] << tab << nWall[pfacei] << tab << uwall[pfacei] << tab
             << mup[pfacei] << tab << sigmap[pfacei] << tab << Ca [pfacei]<< tab
             << convertToDeg*thetaDp[pfacei] << endl;
    }
//////////////////////////////////////////////////////////////////////
*/
//	return convertToDeg*thetaDp;
	return 1.0*thetaDp;
}

void dynamicYokoiAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaE") << thetaE_ << token::END_STATEMENT << nl;
    os.writeKeyword("ka") << ka_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr") << kr_ << token::END_STATEMENT << nl;
    os.writeKeyword("muEffYokoi") << muName_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaYokoi") << sigmaName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, dynamicYokoiAlphaContactAngleFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
