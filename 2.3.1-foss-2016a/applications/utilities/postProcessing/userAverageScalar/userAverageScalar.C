#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::validArgs.append("scalarFieldName");

    #include "setRootCase.H"
    #include "createTime.H"

	word scalarFieldName(args.additionalArgs()[0]);

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

    Info<< "for the first timeName() : " << runTime.timeName() << endl;
    volScalarField dummySField
    (
        IOobject
        (
            scalarFieldName, //read for its dimension 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    //Info<< "copy its dimension to create the mean field" << endl;
	volScalarField sFieldMean
	(
	    IOobject
	    (
	       scalarFieldName+"_mean", // name to be the mean field
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
        dimensionedScalar
        (
            scalarFieldName,
            dummySField.dimensions(),
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
	           	scalarFieldName,
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
