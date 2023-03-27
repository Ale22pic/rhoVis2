//  =========                 |
//  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//   \\    /   O peration     |
//    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
//     \\/     M anipulation  |
//-------------------------------------------------------------------------------
//    License
//
//    This file is part of OpenFOAM.
//
//    OpenFOAM is free software: you can redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

#include "fvCFD.H"
#include "psiThermo.H"
#include "hePsiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include <fstream>     

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main

int main(int argc, char *argv[])
{
    
    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "readThermophysicalProperties.H"
    #include "variables.H"



//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    turbulence->validate();
    
    
    Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s" << nl << endl;
    
    
    while (runTime.run()) 
    {
    
    
        Info<< "Time = " << runTime.timeName() << nl << endl;
    
        runTime++;

//     	Saving quantities at preavious time step

	rhoOld  = rho; 
       	rhoUOld = rhoU; 
       	rhoEOld = rhoE; 

//     	RK Time step
       	
	for (int cycle =0; cycle < gamCoeff.size(); cycle++)
       	{
    
//     	Speed of sound and Mach number

	c = Foam::sqrt(thermo.Cp()/thermo.Cv()/psi);

//      Interpolated quantities at cell faces     
        
	rhoave = fvc::interpolate(rho) ;
        Uave = fvc::interpolate(U) ;

//      Flux at the intercell
        
	phi = fvc::interpolate(U) & mesh.Sf() ;
        phit = fvc::interpolate(rhoU) & mesh.Sf() ; 
      
//      Enthalpy

	H    = (rhoE + p)/rho ;

//      Enthalpy at the intercell
        
	Have = fvc:: interpolate(H) ;

//      Pressure at the intercell

        pave = fvc::interpolate(p) ;    
 
//	Update sensonr only once every RK step
	

//      Evaluate viscous terms

	rhoFil = fvc::surfaceSum(mesh.magSf() * fvc::interpolate(rho)) / fvc::surfaceSum(mesh.magSf()) ;

	muEff = turbulence -> muEff() ;
	
        muArt = rho * 2.0 * mag((rho-rhoFil)/(rho+rhoFil)) * Chi * Foam::pow(cellVolu,2.0/3.0) *  Foam::sqrt(2.0)*mag(0.5*symm(fvc::grad(U))) ;
	
        tauMC = (muEff+muArt)*dev2(Foam::T(fvc::grad(U))) ;

	muave = fvc::interpolate(muEff+muArt) ; 

        k = thermo.Cp()*(muEff+muArt)/Pr ; 

        kave = fvc::interpolate(k) ; 
        
//	Momentum viscous flux
        
	momVisFlux = muave*(fvc::snGrad(U)*mesh.magSf()) ;
        
//	Energy viscous flux
        
	heatFlux = kave*fvc::snGrad(T)*mesh.magSf() ;
	
        visWork = (momVisFlux + fvc::dotInterpolate(mesh.Sf(),tauMC)) & Uave ;
        
        enVisFlux = heatFlux + visWork ;


// 	Total fluxes, Eulerian + viscous 

        
	rhoFlux = rhoave*phi ;
        
        momFlux = rhoave*Uave*phi + pave*mesh.Sf() - momVisFlux ;

	enFlux = rhoave*Have*phi - enVisFlux ;


	volScalarField rhoFln = -runTime.deltaT()*rhoCoeff[cycle]*rhoFl ;
	volVectorField momFln = -runTime.deltaT()*rhoCoeff[cycle]*momFl ;
	volScalarField enFln = -runTime.deltaT()*rhoCoeff[cycle]*enFl ;


//	This are quantities divided by the volume, since the divergence is
//	applied	

        rhoFl = fvc::div(rhoFlux) ;
        
        momFl = fvc::div(momFlux) - fvc::div(tauMC) ;
        
        enFl = fvc::div(enFlux) ;



// 	RK sub-step
        
	rho  = rho + rhoFln + gamCoeff[cycle]*runTime.deltaT()*(
                - rhoFl) ;         
        
        rhoU = rhoU + momFln + gamCoeff[cycle]*runTime.deltaT()*(
                - momFl) ;

        rhoE = rhoE + enFln + gamCoeff[cycle]*runTime.deltaT()*(
                -enFl ) ;

//	Update primitive variables and boundary conditions


	U.ref() = rhoU() / rho();
        
        U.correctBoundaryConditions();
	
        rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

        e = rhoE/rho - 0.5*magSqr(U);
        
        e.correctBoundaryConditions();

	
//	Thermodinamic library
        
	thermo.correct();


        rhoE.boundaryFieldRef() = rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()) );

	p.ref() = rho() / psi();
        
        p.correctBoundaryConditions();

	rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField(); //psi=1/(R*T)


       	}
	
//	End of RK time integration
       	
	runTime.write();
       	turbulence->correct(); //turbulence model


        #include "diagnostics.H" //print tke on diagnostics.dat
        #include "step.H"        //Evaluate Courant Number
        #include "setDeltaT.H"   //Andjust time step


    	}
        
	runTime.write();
    
	Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s" << nl << endl;

    	Info<< "End\n" << endl;

    	return 0;
}
// ************************************************************************* //
