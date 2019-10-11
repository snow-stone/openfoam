#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

//#include "OFstream.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const word patchName = "Port2";
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
        Info<< "Proceeding dS output: " << patchName << endl;
        //OFstream path(dS.objectPath());
        dS.write();
    }
    else
    {
        Info<< "Cannot find patch named : " << patchName << endl;
    }

    
    return 0;
}

// ************************************************************************* //
