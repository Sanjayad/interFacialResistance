/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Alex Rattner
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

#include "interfacialResistance_SGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(interfacialResistance_SGS, 0);
    addToRunTimeSelectionTable(thermalPhaseChangeModel, interfacialResistance_SGS, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::interfacialResistance_SGS::interfacialResistance_SGS
(
		const word& name,
		const dictionary& thermalPhaseChangeProperties,
		const twoPhaseThermalMixture& twoPhaseProperties,
		const volScalarField& T,
		const volScalarField& alpha1
)
:
    thermalPhaseChangeModel(name, thermalPhaseChangeProperties, twoPhaseProperties, T, alpha1),
	mesh_(T.mesh()),
	Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    ),
   	Q_pc_sgs_
    (
        IOobject
        (
            "sgsPhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    ),

	qFlux_sgs_ // Declare the flux from integral as a volumeScalarField
    (
        IOobject
        (
            "sgsPhaseChangeHeatFlux",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,0,-3,0,0,0,0), 0 )
    ),
	wet
    (
        IOobject
        (
            "wet",
            T_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,0,0,0,0,0,0), 0 )
     ),
	faceTime
    (
        IOobject
        (
            "faceTime",
            T_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 15E-1 )
    ),
	T_sp_coeff_
    (
        IOobject
        (
            "T_sp_coeff",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-1,-3,-1,0,0,0), 0 )
    ),
	T_sc_coeff_
    (
        IOobject
        (
            "T_sc_coeff",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    ),
	InterfaceMeshGraph( mesh_, alpha1 ),
    InterfaceField_
    (
        IOobject
        (
            "InterfaceField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        scalar(0)
    ),
	WallField
    (
        IOobject
        (
            "WallField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0)
    ),
	interfaceArea //Is initialized to zero
	(
        IOobject
        (
            "interfaceArea",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(0,2,0,0,0,0,0), 0 )
	),
	threshold_ //Is initialized to zero
	(
        IOobject
        (
            "threshold",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(0,0,0,0,0,0,0), 0 )
	),
	R_g( thermalPhaseChangeProperties_.lookup("R_g") ),
	sigmaHat( thermalPhaseChangeProperties_.lookup("sigmaHat") ),
	v_lv( (32.0/twoPhaseProperties_.rho2().value()) - (1.0/twoPhaseProperties_.rho1().value()) ),
	hi( (2.0*sigmaHat.value()/(2.0-sigmaHat.value())) * (h_lv_.value()*h_lv_.value()/(T_sat_.value()*v_lv)) * pow(1.0/(2.0*3.1416*R_g.value()*T_sat_.value()),0.5) ),
	C_1 ( 8.379E5 ),
	C_2 ( -0.2356 ),
	C_3 ( 3E-4 )

{
	//Read in the cond/evap int. thresholds
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::interfacialResistance_SGS::calcQ_pc()
{
//Info << "hi" << tab << hi << endl;
		
	//Get the sets of interface cell face pairs for evaporation/condensation
	std::vector<MeshGraph::CellFacePair> CondIntCellFacePairs, EvapIntCellFacePairs;


	//Find internal interface cell pairs using graph traversal
	InterfaceMeshGraph.Reset();
	InterfaceMeshGraph.GetInterfaceCellFacePairs(CondIntCellFacePairs, CondThresh);
	InterfaceMeshGraph.Reset();
	InterfaceMeshGraph.GetInterfaceCellFacePairs(EvapIntCellFacePairs, EvapThresh);

	//Compute the interpolated T field to see which face pairs are actually evaporating/condensing:
	surfaceScalarField Tf = fvc::interpolate(T_);

	//Reset interface field, then interpolate
	InterfaceField_ = 0;

	//Loop through cond cells:
	for (std::vector<MeshGraph::CellFacePair>::iterator it = CondIntCellFacePairs.begin(); it != CondIntCellFacePairs.end(); it++)
	{
		//Check that temp is below T_sat for condensation
		if ( Tf[(*it).f] <= T_sat_.value() )
		{   InterfaceField_[(*it).c1] = 1;  InterfaceField_[(*it).c2] = 1;  }
	}

	//Loop through evap cells:
	for (std::vector<MeshGraph::CellFacePair>::iterator it = EvapIntCellFacePairs.begin(); it != EvapIntCellFacePairs.end(); it++)
	{
		//Check that temp is above T_sat for evaporation
		if ( Tf[(*it).f] >= T_sat_.value() )
		{   InterfaceField_[(*it).c1] = 1;  InterfaceField_[(*it).c2] = 1;  }
	}

	//Spit out internal interface cells count
	//Info<< "Internal interface cells: " << gSum(InterfaceField_) << endl;


	//Now add wall cells to the interfaceField:
	labelList WallCells;
	forAll( mesh_.boundary(), pI )
	{
		if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
		{  WallCells.append( mesh_.boundary()[pI].faceCells() );  }
	}
	WallField = 0;
	forAll( WallCells, cI )
	{   
		WallField[WallCells[cI]] = 1;
		//InterfaceField_[WallCells[cI]] = 1;
	}

	//List total int. cells
	//Info<< "Total interface cells: " << gSum(InterfaceField_) << endl;
	
	//Reset all Q_pc to 0
	Q_pc_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );
	Q_pc_sgs_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );

	//Reset all threshold to 0
	threshold_ = dimensionedScalar( "dummy", dimensionSet(0,0,0,0,0,0,0), 0 );

	//Compute some helpful props:
	//For some reason dT is dimensionless
	const dimensionedScalar& dT = alpha1_.time().deltaTValue() * dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 1.0 );
	const dimensionedScalar& rho1 = twoPhaseProperties_.rho1();
	const dimensionedScalar& rho2 = twoPhaseProperties_.rho2();

	interfaceArea.internalField() = mag(fvc::grad(alpha1_))*mesh_.V();
Info << "interfaceArea1 = " << gSum(interfaceArea.internalField()) << endl;
Info << "interfaceArea2 = " << gSum(InterfaceField_*interfaceArea.internalField()) << endl;
Info << "hi = " << hi << endl;
Info << "vlv = " << v_lv << endl;
Info << "dtUtilizedByTheThermalPhaseChangeModel = " << dT.value() << endl;

	//limited phase change heat
	//T_sp_coeff_.internalField() = hi*interfaceArea/(mesh_.V()*twoPhaseProperties_.cp());
	//T_sc_coeff_.internalField() = -hi*interfaceArea*T_sat_/(mesh_.V()*twoPhaseProperties_.cp());
	T_sp_coeff_.internalField() = hi*interfaceArea/mesh_.V();
	T_sc_coeff_.internalField() = -hi*interfaceArea*T_sat_/mesh_.V();
	Q_pc_.internalField() = hi*interfaceArea*(T_-T_sat_)/mesh_.V(); 
	forAll(mesh_.cells(),pI)
	{
		threshold_[pI] = ( (alpha1_[pI] >= 0.00) && (alpha1_[pI] <= 1.00) ) ? 1.0 : 0.0; 
	}
	

	//Now apply subgridscale model on the wall patches
	surfaceScalarField alpha1f = fvc::interpolate(alpha1_); //Alpha1 on faces

	//- Applying the subGrid Model
	forAll(mesh_.boundary(),pI)
	{
		if(isA<wallFvPatch>(mesh_.boundary()[pI]))
		{
			scalarField& faceTimePatch = faceTime.boundaryField()[pI];
			faceTimePatch += dT.value();
			scalarField& wetPatch = wet.boundaryField()[pI];  
			scalarField& alphaPatch = alpha1f.boundaryField()[pI];
			scalarField& qFlux_sgsPatch = qFlux_sgs_.boundaryField()[pI];
			
			forAll(wetPatch, fI) 
			{
				wetPatch[fI] = ( (alphaPatch[fI] > 0.9) || ((wetPatch[fI] == 1.0) && (alphaPatch[fI] > 0.1)) ) ? 1.0 : 0.0;
			}

			faceTimePatch = max(SMALL, (1.0-wetPatch)*faceTimePatch);
			
			forAll(faceTimePatch, fI)
			{					
				if(faceTimePatch[fI] <= 0.0)
				{

					qFlux_sgsPatch[fI] = (1.0-wetPatch[fI])*C_1*pow(SMALL,C_2);  // It's actually heat transfer coefficient
				}
				else if(faceTimePatch[fI] > 0.0 && faceTimePatch[fI] <= C_3)
				{
					qFlux_sgsPatch[fI] = (1.0-wetPatch[fI])*C_1*pow(faceTimePatch[fI],C_2);
				}
				else
				{
					qFlux_sgsPatch[fI] = 0.0;
				}			
			}
//Info << "patch time = " << faceTimePatch;
		}	
	}
	//Calculate volumetric SGS phase change rate in the cell	
	Q_pc_sgs_.internalField() = fvc::surfaceIntegrate( (1-alpha1f)*qFlux_sgs_*mesh_.magSf() )*(T_-T_sat_);

	//decaying Phase Change Heat per unit volume
	//Q_pc_.internalField() = threshold_*twoPhaseProperties_.rho()*twoPhaseProperties_.cp()*((1.0-exp(-hi*interfaceArea*dT.value()/(mesh_.V()*twoPhaseProperties_.rho()*twoPhaseProperties_.cp())))*(T_-T_sat_)/dT.value());
	//Q_pc_.internalField() = twoPhaseProperties_.rho()*twoPhaseProperties_.cp()*((1.0-exp(-hi*interfaceArea*dT.value()/(mesh_.V()*twoPhaseProperties_.rho()*twoPhaseProperties_.cp())))*(T_-T_sat_)/dT.value());

	//Unlimited phase change heat
	//Q_pc_ = InterfaceField_*twoPhaseProperties_.rho()*twoPhaseProperties_.cp()*((T_-T_sat_)/dT);
	

	//Fluid availability limits
	//Get cond/evap limits
	volScalarField LimCond = (1.0-alpha1_)*( rho2*h_lv_ / dT );
	//No evaporation on wall cells!
	volScalarField LimEvap = (1.0-WallField)*alpha1_*rho1*h_lv_ / dT;
	//volScalarField LimEvap = alpha1_*rho1*h_lv_ / dT;

	//Apply fluid limiting
	volScalarField Q_pc_fluid = neg(Q_pc_)*max(Q_pc_, -LimCond) + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

	//Volume-based limiting (i.e. relative phase change rate can't exceed |1| per time step
	volScalarField PCV_fac = dT*(Q_pc_ / h_lv_)*( (scalar(1.0)/twoPhaseProperties_.rho2()) - (scalar(1.0)/twoPhaseProperties_.rho1()) );

	//Again, don't allow evap on wall	
	volScalarField Q_pc_vol = Q_pc_ * mag( min( max(1.0/(PCV_fac+SMALL), -1.0), (1.0-WallField) ) );
	//volScalarField Q_pc_vol = Q_pc_ * mag( min( max(1.0/(PCV_fac+SMALL), -1.0), 1.0 ) );

	//Composite limit
	Q_pc_ = neg(Q_pc_)*max( max( Q_pc_, Q_pc_fluid ), Q_pc_vol) + pos(Q_pc_)*min( min( Q_pc_, Q_pc_fluid ), Q_pc_vol);
}


bool Foam::thermalPhaseChangeModels::interfacialResistance_SGS::read(const dictionary& thermalPhaseChangeProperties)
{
	thermalPhaseChangeModel::read(thermalPhaseChangeProperties);

	//Read in the cond/evap int. thresholds
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	return true;
}


// ************************************************************************* //
