//      Evaluate time step of Navier-Stokes equations

//      volScalarField muEff(turbulence->muEff())  ;
	
	rhoFil = fvc::surfaceSum
        (
        mesh.magSf()*fvc::interpolate(rho)
        )/fvc::surfaceSum(mesh.magSf());
	
	muEff = turbulence -> muEff() ;
        muArt = rho * 2.0 * mag((rho-rhoFil)/(rho+rhoFil)) * Chi * Foam::pow(cellVolu,2.0/3.0) *  Foam::sqrt(2.0)*mag(0.5*symm(fvc::grad(U))) ;

	volScalarField k = thermo.Cp()*(muEff+muArt)/Pr ;

        scalar CoNum = 0.0;
        

        
        if (mesh.nInternalFaces())
        {
         surfaceScalarField cave = fvc::interpolate(c) ;  
         surfaceScalarField un   = fvc::interpolate(mag(U));
         surfaceScalarField lambda = mag(un + cave);
         surfaceScalarField evmax    = lambda;
//       surfaceScalarField evmu     = fvc:: interpolate(muEff/rho);
         surfaceScalarField evmu     = fvc:: interpolate((muEff+muArt)/rho);
         surfaceScalarField eval     = fvc:: interpolate(k/rho/thermo.Cp());
         //dx is the distance between two centers projected along the face normal
         
         surfaceScalarField dx
         (
          1./mesh.surfaceInterpolation::deltaCoeffs()
         );

//       surfaceScalarField dx       = 1./mesh.surfaceInterpolation::nonOrthDeltaCoeffs();
         surfaceScalarField dx2      = sqr(dx);
         double dtdx_inv= max(evmax/dx).value() ;//inviscid limitation
         double dtdx_rmu= max(evmu/dx2).value() ;//viscous limitation
         double dtdx_ral= max(eval/dx2).value() ;//viscous limitation
 
         CoNum   = max(dtdx_inv,dtdx_rmu)     ;
         CoNum   = max(dtdx_ral,CoNum)        ;
         CoNum   = 1./CoNum;
        }
//
//
        runTime.setDeltaT
        (
         maxCo*CoNum
        );


