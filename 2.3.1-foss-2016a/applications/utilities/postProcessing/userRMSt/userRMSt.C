#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::validArgs.append("scalarFieldName");
    argList::validArgs.append("timeOfAverageField");

    #include "setRootCase.H"
    #include "createTime.H"

	word scalarFieldName(args.additionalArgs()[0]);
	word timeOfAverageField(args.additionalArgs()[1]);

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

    IOobject meanHeader
    (
	    string(scalarFieldName+"_mean"),
	    timeOfAverageField,
	    mesh,
	    IOobject::MUST_READ
    );

    volScalarField rmsField
    (
        IOobject
        (
            "tt_mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
		dimensionedScalar
	    (
            "tt_mean",
		    dimTemperature,
		    0.0	
		)
    );

	volScalarField varField = sqr(rmsField);

    if(meanHeader.headerOk())
	{
		volScalarField mean(meanHeader, mesh);
        label nFields = 0;

		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI); 
		    Info<< "Time = " << runTime.timeName() << endl;

			volScalarField scalarField
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

			volScalarField sPrime = scalarField - mean;
			varField += sqr(sPrime);

            nFields += 1;

	    }
        if (nFields >= 1)
        {
            Info<< "Iteration ended" << endl;
            Info<< "runTime.timeName() : " << runTime.timeName() << endl;
            Info<< "nFields : " << nFields << endl;
			varField = varField / nFields;
            rmsField = sqrt(varField);
            rmsField.write(); // write only at the last time step
        }
        else
        {
            FatalError
                << "nFields = 0"
                << nl << exit(FatalError);
        }
	}
	else
	{
		Info<< "No mean field" << endl;
	}

    return 0;
}

// ************************************************************************* //
