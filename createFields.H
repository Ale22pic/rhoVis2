#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;


autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);

psiThermo& thermo = pThermo();

volScalarField& e = thermo.he();

//const volScalarField& cv = thermo.Cv();


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField cellVolu
(
    IOobject
    (
        "cellVolu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimVolume,0.0)
);

    
cellVolu.ref() = mesh.V();


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


#include "readInput.H"


volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);


surfaceScalarField surfArea 
(
    IOobject
    (
        "surfArea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.magSf() 
);


volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);


surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

volScalarField Chi
(
    IOobject
    (
        "Chi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Chi", dimensionSet(0,0,0,0,0,0,0), chi)
);


surfaceScalarField phi("phi", fvc::flux(U));

surfaceScalarField phit("phit", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;

autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phit,
        thermo
    )
);


// RK coefficient 

//scalarList rkCoeff(0) ;
//
//rkCoeff.append(1./4.) ;
//rkCoeff.append(1./3.) ;
//rkCoeff.append(1./2.) ;
//rkCoeff.append(1./1.) ;

//scalarList rkCoeff(0) ;
//
//rkCoeff.append(0.1100) ;
//rkCoeff.append(0.2766) ;
//rkCoeff.append(0.5000) ;
//rkCoeff.append(1.0000) ;

scalarList rhoCoeff(0) ;

rhoCoeff.append(0.) ;
rhoCoeff.append(-17./60.) ;
rhoCoeff.append(-5./12.) ;

scalarList gamCoeff(0) ;

gamCoeff.append(8./15.) ;
gamCoeff.append(5./12.) ;
gamCoeff.append(3./4.) ;

