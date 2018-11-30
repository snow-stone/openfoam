#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"


	IOList<label> cellSet
	(
	    IOobject
	    (
		    "selectCells",
		    runTime.time().constant(),
			mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		)	
	);

	volVectorField U
	(
	    IOobject
		(
	        "U",
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		),
		mesh
	);

	/*
	forAll(cellSet, i)
	{
	    Info<< "U[" << cellSet[i] << "] =" << U.internalField()[cellSet[i]] << endl;
	}
	*/

	for (label i=0; i < 10; i++)
	{
	    Info<< "U[" << cellSet[i] << "] =" << U.internalField()[cellSet[i]] << endl;
	}

    return 0;
}

// ************************************************************************* //
