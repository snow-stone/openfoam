#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "cellSet.H"

#include <fstream>

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("scalarFieldName");
	//argList::validArgs.append("timeOfAverageField");
	argList::validArgs.append("sliceStore");
	argList::validArgs.append("sliceName");
	Foam::argList::addBoolOption
	(
	    "noWriteLog",
	    "suppress writting to userDefinedLog/TKE"
    );

    #include "setRootCase.H"
    #include "createTime.H"

	word scalarFieldName(args.additionalArgs()[0]);
	//word timeOfAverageField(args.additionalArgs()[1]);
	word sliceStore(args.additionalArgs()[1]);
	word sliceName(args.additionalArgs()[2]);
	bool writeLog = !args.optionFound("noWriteLog");

	scalar D = 0.008;
    scalar n = 130;
    scalar surface = D * D;
	scalar dS = surface / (n * n);

    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

	IOList<label> slice
	(
	    IOobject
	    (
		    sliceName,
			sliceStore,
			mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		)	
	);

/*
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
*/

	label cellN  = slice.size();
    scalarField slice_T(cellN);
    scalarField slice_T_mean(cellN);

    scalarField TMeanTemp(cellN, 0.0);
    Info<<" TMeanTemp = " << TMeanTemp << endl;
	label nField = 0;
	scalar TMeanSpatial;
	scalar TVARSpatial;
	//scalar TRMSSpatial;

    forAll(slice, i)
    {
        //slice_T_mean[i] = TMean.internalField()[slice[i]];
        /*
    	if (i < 3)
    	{
    	    Info<< "TMean[" << slice[i] <<"] : "<< slice_T_mean[i] << endl;
    	}
		*/
   	}

	mkDir("userDefinedLog");
	std::ofstream sliceLog
	(
		fileName(string("userDefinedLog")/string(sliceName+"_mean_rms_Dai")).c_str(),
		ios_base::app
	);

	forAll(timeDirs, timeI) { runTime.setTime(timeDirs[timeI], timeI);
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
		//TRMSSpatial = 0.0;
    
    	forAll(slice, i)
    	{
    	    slice_T[i] = T.internalField()[slice[i]];
			/*
    		if (i < 3)
    		{
    		    Info<< "T[" << slice[i] <<"] : "<< slice_T[i] << endl;
    		}
			*/
    	}

		TMeanSpatial = 0.5;
		Info<< "TMeanSpatial : " << TMeanSpatial << endl;
		Info<< "Version Dai : " << TMeanSpatial << endl;
		forAll(slice, i)
		{
			TVARSpatial += sqr(slice_T[i] - TMeanSpatial) * dS ;
		}
		//TVARSpatial /= cellN;
		TVARSpatial *= 1.0/(0.5 * 0.5 * surface);
		//TRMSSpatial = sqrt(TVARSpatial); // error: call of overloaded ‘sqrt(Foam::scalar&)’ is ambiguous
		//TRMSSpatial = Foam::sqrt(TVARSpatial);
		//Info<< "TRMSSpatial : " << TRMSSpatial << endl;

		TMeanTemp += slice_T;
		nField++;

		if (writeLog)
		{
		    sliceLog << runTime.timeName() << " "
				     << TMeanSpatial << " "
				     << TVARSpatial << " "
				     << std::endl;	 
		}
	}

	if(nField >= 1)
	{
		/*
	    Info<<"labelList op mean : " << endl;
    	forAll(slice, i)
		{
    		if (i < 3)
    		{
    	    	Info<< "T[" << slice[i] <<"] : "<< TMeanTemp[i]/nField << endl;
    		}
    	}
		*/
	}

    return 0;
}

// ************************************************************************* //
