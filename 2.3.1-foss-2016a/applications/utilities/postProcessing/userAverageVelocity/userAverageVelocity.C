#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::validArgs.append("velocityFieldName");

    #include "setRootCase.H"
    #include "createTime.H"

	word velocityFieldName(args.additionalArgs()[0]);

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

	volVectorField vFieldMean
	(
	    IOobject
	    (
	       velocityFieldName+"_mean",
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedVector
		(
		   velocityFieldName+"_mean",
		   dimVelocity,
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
	           	velocityFieldName,
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
