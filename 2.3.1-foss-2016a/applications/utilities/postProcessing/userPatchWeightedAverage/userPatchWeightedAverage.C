#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

//#include "OFstream.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("scalarFieldName");
    argList::validArgs.append("patchName");

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    argList::validArgs.append("scalarFieldName");
    argList::validArgs.append("patchName");
    word scalarFieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);

    const polyBoundaryMesh& pp = mesh.boundaryMesh();
    const label patchLabel = pp.findPatchID(patchName);
    const surfaceScalarField& magSf = mesh.magSf();
    if (patchLabel != -1)
    {
        Info<< "Found patch named : " << patchName << endl;
	    IOList<scalar> dS 
	    (
	        IOobject
	        (
	            "dS",
	            "constant",
	            runTime,
	            IOobject::NO_READ,
	            IOobject::AUTO_WRITE
	        ),
            magSf.boundaryField()[patchLabel]
	    );

	    forAll(timeDirs, timeI)
	    {
	        runTime.setTime(timeDirs[timeI], timeI); 
	        Info<< "Time = " << runTime.timeName() << endl;
	
	        volScalarField sField
	        (
	            IOobject
	            (
	                scalarFieldName,
	                runTime.timeName(),
	                mesh,
	                IOobject::MUST_READ
	            ),
	            mesh
	        );
	
	        scalar average = 0.0;
	        forAll(sField.boundaryField()[patchLabel], facei)
            {
                average += sField.boundaryField()[patchLabel][facei] * dS[facei];
            }
            average /= sum(dS);

	        Info<< "Average " << scalarFieldName << " on patch " 
                << patchName << " : " << average << endl;
	    }
    }
    {
        Info<< "Cannot find patch named : " << patchName << endl;
    }


    return 0;
}

// ************************************************************************* //
