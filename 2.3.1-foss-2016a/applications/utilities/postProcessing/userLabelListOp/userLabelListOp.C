#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "cellSet.H"

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

	IOList<label> slice
	(
	    IOobject
	    (
		    "selectCells",
		    "constant",
			mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		)	
	);

    volVectorField UMean
    (
        IOobject
        (
           velocityFieldName+"_mean",
		   timeOfAverageField,
           mesh,
           IOobject::MUST_READ,
		   IOobject::NO_WRITE
        ),
        mesh
    );

    vectorField slice_U(slice.size());
    vectorField slice_U_mean(slice.size());

    vectorField sum(slice.size(), vector::zero);
    Info<<" sum = " << sum << endl;
	label nField = 0;

    forAll(slice, i)
    {
        slice_U_mean[i] = UMean.internalField()[slice[i]];
    	if (i < 3)
    	{
    	    Info<< "UMean[" << slice[i] <<"] : "<< slice_U_mean[i] << endl;
    	}
   	}

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI);
		Info<< "Time = " << runTime.timeName() << endl;

	    volVectorField U
	    (
	        IOobject
		    (
	            velocityFieldName,
	            runTime.timeName(),
	            mesh,
	            IOobject::MUST_READ,
		        IOobject::NO_WRITE	
		    ),
		    mesh
	    );

    	//Info<< "slice : " << slice << endl;
    	Info<< "size(slice) : " << slice.size() << endl;
    
    	forAll(slice, i)
    	{
    	    slice_U[i] = U.internalField()[slice[i]];
    		if (i < 3)
    		{
    		    Info<< "U[" << slice[i] <<"] : "<< slice_U[i] << endl;
    		}
    	}

		sum += slice_U;
		nField++;
	}

	if(nField >= 1)
	{
	    Info<<"labelList op mean : " << endl;
    	forAll(slice, i)
		{
    		if (i < 3)
    		{
    	    	Info<< "U[" << slice[i] <<"] : "<< sum[i]/nField << endl;
    		}
    	}
	}

    return 0;
}

// ************************************************************************* //
