	//Use limited alpha1:    
	const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );
	
	//Thermal conductivity
	surfaceScalarField Kf
    (
         "Kf",
		fvc::interpolate( limAlpha1*lambda1 + (1-limAlpha1)*lambda2 )
    );

	//Check how many Energy equation sub-cycles to use
	label nEnergySubCycles( pimple.dict().lookupOrDefault<int>("nEnergySubCycles", 1) );

	//Update phase change model at new temperature field values
	phaseChangeModel->correct();

	//Latent heat accounted for in source term
	T.correctBoundaryConditions();

	//Energy-temperature loop
	label nEnergyLoops(readLabel(pimple.dict().lookup("nEnergyLoops")));
	for (int TEqnCount=0; TEqnCount < nEnergyLoops; TEqnCount++)
	{
		phaseChangeModel->correct();		
		//Form and solve the energy equation
		fvScalarMatrix TEqn
		(
			  fvm::ddt(rhoCp, T) 
			+ fvm::div(rhoPhiCp, T)
			- fvm::laplacian(Kf, T)
			==
			-  fvm::Sp(phaseChangeModel->T_sp_coeff(), T)
            -         phaseChangeModel->T_sc_coeff()
		);
		TEqn.solve();
		//Now reevaluate H for the updated temperature field
		H = ( (T - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
	}

