#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::validArgs.append("velocityFieldName");
    argList::validArgs.append("timeOfAverageField");

    #include "setRootCase.H"
    #include "createTime.H"

	word velocityFieldName(args.additionalArgs()[0]);
	word timeOfAverageField(args.additionalArgs()[1]);

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

    IOobject meanHeader
    (
	    string(velocityFieldName+"_mean"),
	    timeOfAverageField,
	    mesh,
	    IOobject::MUST_READ
    );

    if(meanHeader.headerOk())
	{
		volVectorField mean(meanHeader, mesh);

		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI); 
		    Info<< "Time = " << runTime.timeName() << endl;

			volVectorField velocityField
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

			volVectorField UPrime = velocityField - mean;
	
	        volTensorField reyTensor
	        (
	            IOobject
	            (
	                "reyTensor",
	                runTime.timeName(),
	                mesh,
	                IOobject::NO_READ,
	                IOobject::AUTO_WRITE
	            ),
				UPrime * UPrime
			);
	
			reyTensor.write();
	
	    }
	}
	else
	{
		Info<< "No mean field" << endl;
	}

    return 0;
}

// ************************************************************************* //
