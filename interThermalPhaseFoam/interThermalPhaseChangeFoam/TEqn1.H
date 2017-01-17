//Use limited alpha1:    
const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );

//Update kEff
kEff = twoPhaseProperties.lambdaf() + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );

//Keep a reference to the global time step:
//dimensionedScalar totalDeltaT = runTime.deltaT();

//Latent heat accounted for in source term
T.correctBoundaryConditions();
phaseChangeModel->correct();			
//Form and solve the energy equation
fvScalarMatrix TEqn
(
	fvm::ddt(rhoCp, T)
	+ fvm::div(rhoPhiCp, T)
	- fvm::laplacian(kEff, T)
	==
	(
		- fvm::Sp(phaseChangeModel->T_sp_coeff(), T) 
		- phaseChangeModel->T_sc_coeff()
	)

);
TEqn.solve();
//Now reevaluate H for the updated temperature field
H = ( (T - T_0)*(limAlpha1*rho1*cp1 + (1.0-limAlpha1)*rho2*cp2) )/rho;
