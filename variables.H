volScalarField muEff
(
    IOobject
    (
        "muEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muEff",dimensionSet(1,-1,-1,0,0,0,0),scalar(0))
);

volScalarField muArt
(
    IOobject
    (
        "muArt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muArt",dimensionSet(1,-1,-1,0,0,0,0),scalar(0))
);

surfaceScalarField rhoave
(
    IOobject
    (
        "rhoave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoave",dimensionSet(1,-3,0,0,0,0,0),scalar(0))
);

surfaceVectorField Uave
(
    IOobject
    (
        "Uave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Uave",dimensionSet(0,1,-1,0,0,0,0),vector(0,0,0))
);

volVectorField rhoUOld
(
    IOobject
    (
        "rhoUOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("rhoUOld",dimensionSet(1,-2,-1,0,0,0,0),vector(0,0,0))
);     
      

volScalarField rhoOld
(
    IOobject
    (
        "rhoOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoOld",dimensionSet(1,-3,0,0,0,0,0),scalar(0))
);


volScalarField rhoEOld
(
    IOobject
    (
        "rhoEOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoEOld",dimensionSet(1,-1,-2,0,0,0,0),scalar(0))
);

volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("c",dimensionSet(0,1,-1,0,0,0,0),scalar(0))
);


volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("H",dimensionSet(0,2,-2,0,0,0,0),scalar(0))
);

surfaceScalarField Have
(
    IOobject
    (
        "Have",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Have",dimensionSet(0,2,-2,0,0,0,0),scalar(0))
);

surfaceScalarField pave
(
    IOobject
    (
        "pave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pave",dimensionSet(1,-1,-2,0,0,0,0),scalar(0))
);

volScalarField rhoFil
(
    IOobject
    (
        "rhoFil",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("rhoFil",dimensionSet(1,-3,0,0,0,0,0),scalar(0))
);

volTensorField tauMC
(
    IOobject
    (
        "tauMC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("tauMC",dimensionSet(1,-1,-2,0,0,0,0),tensor(0,0,0,0,0,0,0,0,0))
);

surfaceScalarField muave
(
    IOobject
    (
        "muave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muave",dimensionSet(1,-1,-1,0,0,0,0),scalar(0))
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("k",dimensionSet(1,1,-3,-1,0,0,0),scalar(0))
);

surfaceScalarField kave
(
    IOobject
    (
        "kave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("kave",dimensionSet(1,1,-3,-1,0,0,0),scalar(0))
);

surfaceVectorField momVisFlux
(
    IOobject
    (
        "momVisFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("momVisFlux",dimensionSet(1,1,-2,0,0,0,0),vector(0,0,0))
);

surfaceScalarField heatFlux
(
    IOobject
    (
        "heatFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("heatFlux",dimensionSet(1,2,-3,0,0,0,0),scalar(0))
);

surfaceScalarField visWork
(
    IOobject
    (
        "visWork",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("visWork",dimensionSet(1,2,-3,0,0,0,0),scalar(0))
);

surfaceScalarField rhoFlux
(
    IOobject
    (
        "rhoFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoFlux",dimensionSet(1,0,-1,0,0,0,0),scalar(0))
);

surfaceVectorField momFlux
(
    IOobject
    (
        "momFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("momFlux",dimensionSet(1,1,-2,0,0,0,0),vector(0,0,0))
);

surfaceScalarField enFlux
(
    IOobject
    (
        "enFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("enFlux",dimensionSet(1,2,-3,0,0,0,0),scalar(0))
);

surfaceScalarField enVisFlux
(
    IOobject
    (
        "enVisFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("enVisFlux",dimensionSet(1,2,-3,0,0,0,0),scalar(0))
);

volScalarField rhoFl
(
    IOobject
    (
        "rhoFl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoFl",dimensionSet(1,-3,-1,0,0,0,0),scalar(0))
);

volVectorField momFl
(
    IOobject
    (
        "momFl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("momFl",dimensionSet(1,-2,-2,0,0,0,0),vector(0,0,0))
);

volScalarField enFl
(
    IOobject
    (
        "enFl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("enFl",dimensionSet(1,-1,-3,0,0,0,0),scalar(0))
);


