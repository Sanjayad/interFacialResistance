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

Class
	Foam::thermalPhaseChangeModel::HiLoRelaxed

Description
    This model forces wall and interface cells to the equilibrium state,
    causing evaporation or condensation as appropriate. The phase change
    heat term is limited to not produce/sink more than the volume of a 
	cell in a single time step(like the courantAlpha limit).

    This model uses different interface thresholds for identifying
	evaporating and condensing cell pairs.

	This is similar to the previous HiLoFacePair code, but the phase change
	rate is numerically relaxed by some user-defined value (should be <= 1,
	where 0 is no phase change, and 1 is no relaxation). This seems to help
	with overshooting where the temperature in a cell can go from superheated
	to subcooled in one time step.


SourceFiles
    HiLoRelaxed.C

\*--------------------------------------------------------------------*/

#ifndef HiLoRelaxed_H
#define HiLoRelaxed_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                              Class HiLoRelaxed
\*--------------------------------------------------------------------*/

class HiLoRelaxed
:
    public thermalPhaseChangeModel
{

// Private Member Functions

	//- Calculate and return the phase change heat rate
	void calcQ_pc();

	//- This phase change model needs access to the mesh:
	const fvMesh& mesh_;

	//- Phase change heat
	volScalarField Q_pc_;

	//- Graph of mesh grid, used for finding interface cells
	MeshGraph InterfaceMeshGraph;

	//- Field for tracking Interface cells
    volScalarField InterfaceField_;

	//- Field for blocking evaporation on wall cells
	volScalarField WallField;

	//Scalar parameters for identifying condensation and evaporation int. cell pairs:
	scalar CondThresh;
	scalar EvapThresh;
	scalar RelaxFac; //The relaxation factor (0-1)

public:

    //- Runtime type information
    TypeName("HiLoRelaxed");


    // Constructors

	//- construct from components
	HiLoRelaxed
	(
		const word& name,
		const dictionary& thermalPhaseChangeProperties,
		const twoPhaseThermalMixture& twoPhaseProperties,
		const volScalarField& T,
		const volScalarField& alpha1
	);


    //- Destructor
    ~HiLoRelaxed()
    {}


    // Member Functions

        //- Return phase change heat
		tmp<volScalarField> Q_pc() const
		{
			return Q_pc_;
		}

        //- Correct the thermalPhaseChange model
        void correct()
		{
			calcQ_pc();
		}

		//- Return interface cells marker field:
		tmp<volScalarField> InterfaceField() const
		{
			return InterfaceField_;
		}
		

        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
