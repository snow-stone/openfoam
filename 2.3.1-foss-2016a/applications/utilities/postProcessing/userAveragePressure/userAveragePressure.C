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

	volScalarField sFieldMean
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
        dimensionedScalar
        (
            velocityFieldName,
            sqr(dimLength)/sqr(dimTime),
            0.0
        )
	);

	label nFields = 0;

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

		volScalarField sField
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

		sFieldMean += sField; // sField has no dimension while sFieldMean has "dimless"
		nFields += 1;	
	}

	if (nFields >= 1)
	{
        Info<< "Iteration ended" << endl;
        Info<< "runTime.timeName() : " << runTime.timeName() << endl;
        Info<< "nFields : " << nFields << endl;
	    sFieldMean /= nFields;
	    sFieldMean.write(); // write only at the last time step
	}
	else
	{
	    Info<< "input field number : " << nFields << endl;
	}

    return 0;
}

// ************************************************************************* //
