/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    wallShearStress

Description
    Calculates and reports wall shear stress for all patches, for the
    specified times when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include <fstream> // 标准输出，到文件

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar scalarField_simpleStatistics(const scalarField& s)
{
    Info<< "Simple Statistics on field : " << endl;
    Info<< "size of field : " << s.size() << endl;
    Info<< "min : " << min(s) << endl;
    Info<< "max : " << max(s) << endl;
    scalar mean = sum(s)/s.size(); 
    Info<< "athemtic mean : " << mean << endl;
    Info<< nl;

    return mean;
}

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearStress
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volSymmTensorField Reff(model->devReff());

    forAll(wallShearStress.boundaryField(), patchI)
    {
        wallShearStress.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
    Foam::argList::addBoolOption
    (
        "noWrite",
        "no output data to file"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H" 
	#include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);
    bool noWriting = args.optionFound("noWrite");
    #include "createNamedMesh.H"

    if(!noWriting)
    {
        mkDir("postProcessing");
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                //IOobject::AUTO_WRITE
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "wallShearStress",
                sqr(dimLength)/sqr(dimTime),
                vector::zero
            )
        );

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus_"+fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "yPlus", // confirm IOobject overwrites this
                dimless,
                scalar(0.)
            )
        );

        Info << "reading nu_mean" << endl;
        volScalarField nu_mean
        (
            IOobject
            (
                "nu_mean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "nu_mean", // confirm IOobject overwrites this
                dimless,
                scalar(0.)
            )
        );
        Info << "read nu_mean complete" << endl;

        IOobject UHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);
            calcIncompressible(mesh, runTime, U, wallShearStress);

            const polyBoundaryMesh& pp = mesh.boundaryMesh();
            const label patchLabel = pp.findPatchID(patchName);
            volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y(); //#include "nearWallDist.H"

            if (patchLabel != -1)
            {
                vectorField& tauByRho = wallShearStress.boundaryField()[patchLabel];
                scalarField uTau = Foam::sqrt(mag(tauByRho));
                Info<< "On patch " << patchName << endl;
                Info<< "On patch " << patchName << endl;
                scalarField& d_ = d[patchLabel];
                Info<< "calculating yPlus using nu_mean at the boundaryField " << endl;
                forAll(nu_mean.boundaryField()[patchLabel], i)
                {
                    if (nu_mean.boundaryField()[patchLabel][i] == 0){
                        Info<< "i = " << endl;
                    }
                }
                yPlus.boundaryField()[patchLabel] = uTau * d_ / nu_mean.boundaryField()[patchLabel];
                Info<< "Finish calculating yPlus " << endl;
                Info << "d_ : " << endl;
                scalar dMean = scalarField_simpleStatistics(d_);
                Info << "uTau : " << endl;
                scalar uTauMean = scalarField_simpleStatistics(uTau);
                Info << "yPlus : " << endl;
                scalar yPlusMean = scalarField_simpleStatistics(yPlus.boundaryField()[patchLabel]);

                {
                    std::ofstream txtOutput
                    (
                        fileName(
                                string("postProcessing")/string(yPlus.name())
                                ).c_str(),
                        ios_base::app
                    );

                    forAll(wallShearStress.boundaryField()[patchLabel], i)
                    {
                        txtOutput
                            << runTime.timeName() << " "
                            << std::endl;
                    }
                }
            }
            else
            {
                Info<< "no patch named " << patchName << endl;
                Info<< "taking break ..." << endl;
                break;
            }

        }
        else
        {
            Info<< "    no U field" << endl;
        }

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
