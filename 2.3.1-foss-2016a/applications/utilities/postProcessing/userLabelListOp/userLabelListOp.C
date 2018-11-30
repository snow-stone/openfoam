#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "cellSet.H"

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

    volScalarField TMean
    (
        IOobject
        (
           scalarFieldName+"_mean",
		   timeOfAverageField,
           mesh,
           IOobject::MUST_READ,
		   IOobject::NO_WRITE
        ),
        mesh
    );

	label cellN  = slice.size();
    scalarField slice_T(cellN);
    scalarField slice_T_mean(cellN);

    scalarField TMeanTemp(cellN, 0.0);
    Info<<" TMeanTemp = " << TMeanTemp << endl;
	label nField = 0;
	scalar TMeanSpatial;
	scalar TVARSpatial;
	scalar TRMSSpatial;

    forAll(slice, i)
    {
        slice_T_mean[i] = TMean.internalField()[slice[i]];
    	if (i < 3)
    	{
    	    Info<< "TMean[" << slice[i] <<"] : "<< slice_T_mean[i] << endl;
    	}
   	}

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI);
		Info<< "Time = " << runTime.timeName() << endl;

	    volScalarField T
	    (
	        IOobject
		    (
	            scalarFieldName,
	            runTime.timeName(),
	            mesh,
	            IOobject::MUST_READ,
		        IOobject::NO_WRITE	
		    ),
		    mesh
	    );

		TMeanSpatial = 0.0;
		TVARSpatial = 0.0;
		TRMSSpatial = 0.0;
    
    	forAll(slice, i)
    	{
    	    slice_T[i] = T.internalField()[slice[i]];
    		if (i < 3)
    		{
    		    Info<< "T[" << slice[i] <<"] : "<< slice_T[i] << endl;
    		}
    	}

		forAll(slice, i)
		{
			TMeanSpatial += slice_T[i];
		}
		TMeanSpatial /= cellN;
		Info<< "TMeanSpatial : " << TMeanSpatial << endl;
		forAll(slice, i)
		{
			TVARSpatial += sqr(slice_T[i] - TMeanSpatial);
		}
		TVARSpatial /= cellN;
		//TRMSSpatial = sqrt(TVARSpatial); // error: call of overloaded ‘sqrt(Foam::scalar&)’ is ambiguous
		TRMSSpatial = Foam::sqrt(TVARSpatial);
		Info<< "TRMSSpatial : " << TRMSSpatial << endl;

		TMeanTemp += slice_T;
		nField++;
	}

	if(nField >= 1)
	{
	    Info<<"labelList op mean : " << endl;
    	forAll(slice, i)
		{
    		if (i < 3)
    		{
    	    	Info<< "T[" << slice[i] <<"] : "<< TMeanTemp[i]/nField << endl;
    		}
    	}
	}

    return 0;
}

// ************************************************************************* //
