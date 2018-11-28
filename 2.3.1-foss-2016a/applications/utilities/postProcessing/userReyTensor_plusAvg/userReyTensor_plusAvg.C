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

    volSymmTensorField reyTensorMean
    (
        IOobject
        (
           "reyTensor_mean",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
           "reyTensor_mean",
           sqr(dimVelocity),
           symmTensor::zero
        )
    );

    if(meanHeader.headerOk())
	{
		volVectorField mean(meanHeader, mesh);
        label nFields = 0;

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
	
	        volTensorField reyTensor0
	        (
	            IOobject
	            (
	                "reyTensor0",
	                runTime.timeName(),
	                mesh,
	                IOobject::NO_READ,
	                IOobject::NO_WRITE
	            ),
				UPrime * UPrime
			);

            //convert to symmTensor for less storage. And indeed reynTensor0 is mathematically symmetric
            //volSymmTensorField reyTensor(symm(reyTensor)); short but output name cannot be overwritten
            //taking the longer version

	        volSymmTensorField reyTensor
	        (
	            IOobject
	            (
	                "reyTensor",
	                runTime.timeName(),
	                mesh,
	                IOobject::NO_READ,
	                IOobject::AUTO_WRITE
	            ),
                symm(reyTensor0) // = UPrime * UPrime
			);

            //check for difference
			//Info<< "diff : " << mag(reyTensor0.internalField() - reyTensor.internalField()) << endl;

            reyTensorMean += reyTensor;
            nFields += 1;

			reyTensor.write();
	    }
        if (nFields >= 1)
        {
            Info<< "Iteration ended" << endl;
            Info<< "runTime.timeName() : " << runTime.timeName() << endl;
            Info<< "nFields : " << nFields << endl;
            reyTensorMean /= nFields;
            reyTensorMean.write(); // write only at the last time step
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
