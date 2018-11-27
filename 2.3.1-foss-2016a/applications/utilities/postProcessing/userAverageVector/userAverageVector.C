#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::validArgs.append("vectorFieldName");

    #include "setRootCase.H"
    #include "createTime.H"

	word vectorFieldName(args.additionalArgs()[0]);

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

    Info<< "for the first timeName() : " << runTime.timeName() << endl;
    volVectorField dummyVField
    (
        IOobject
        (
            vectorFieldName,//read for its dimension 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    //Info<< "copy its dimension to create the mean field" << endl;
	volVectorField vFieldMean
	(
	    IOobject
	    (
	       vectorFieldName+"_mean",// name to be the mean field
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedVector
		(
		   vectorFieldName+"_mean",
		   dummyVField.dimensions(),
		   vector::zero
		)
	);

	label nFields = 0;

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

		volVectorField vField
		(
	       	IOobject
	       	(
	           	vectorFieldName,
	       	    runTime.timeName(),
	       		mesh,
	          	IOobject::MUST_READ
	      	),
			mesh
		);
		
		vFieldMean += vField; 
		nFields += 1;	
	}

	if (nFields >= 1)
	{
        Info<< "Iteration ended" << endl;
        Info<< "runTime.timeName() : " << runTime.timeName() << endl;
        Info<< "nFields : " << nFields << endl;
	    vFieldMean /= nFields;
	    vFieldMean.write(); // write only at the last time step
	}
	else
	{
	    Info<< "input field number : " << nFields << endl;
	}

    return 0;
}

// ************************************************************************* //
