#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream>

int main(int argc, char *argv[])
{
    Info << "Hello world" << endl;

    timeSelector::addOptions();
    argList::validArgs.append("timeOfAverageField");
	Foam::argList::addBoolOption
	(
	    "noWriteField",
	    "suppress writting for scalarField TKE"
    );
	Foam::argList::addBoolOption
	(
	    "noWriteLog",
	    "suppress writting to userDefinedLog/TKE"
    );
    //argList::validArgs.append("noWriteField");
    //argList::validArgs.append("noWriteLog");

    #include "setRootCase.H"
    #include "createTime.H"

	word timeOfAverageField(args.additionalArgs()[0]);
    bool writeField = !args.optionFound("noWriteField");
    bool writeLog = !args.optionFound("noWriteLog");

	instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"

    IOobject U_meanHeader
    (
        "U_mean",
       	timeOfAverageField, 
        mesh,
        IOobject::MUST_READ
    );

	if(U_meanHeader.headerOk())
	{
		Info<< "Reading U_mean" << " @ time step : " << timeOfAverageField << endl;
		volVectorField U_mean(U_meanHeader, mesh);	

		mkDir("userDefinedLog");
		std::ofstream TKELog
		(
	   		fileName(string("userDefinedLog")/string("TKE")).c_str(),
			ios_base::app
		);

		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI); 
	    	Info<< "Time = " << runTime.timeName() << endl;

        	IOobject Uheader
        	(
            	"U",
            	runTime.timeName(),
            	mesh,
            	IOobject::MUST_READ
        	);

			if (Uheader.headerOk())
			{
				volVectorField U(Uheader, mesh);

				volScalarField TKE
				(
					IOobject
					(
						"TKE",
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					0.5 * (U - U_mean) & (U - U_mean)
				);

				Info<< "gSum(TKE) = " << gSum(TKE) << endl;
				Info<< "weightedAverage(TKE) = " << TKE.weightedAverage(mesh.V()).value() << endl;

				if (writeLog)
                {
				    TKELog
						<< runTime.timeName() << " "
                        << gSum(TKE) << " "
                        << TKE.weightedAverage(mesh.V()).value() << std::endl;
				}
				if (writeField)
				{
				    TKE.write();
                                }
			}
			else
			{
				Info<< "No field U @ Time = " << runTime.timeName() << endl;
			}
    	}
	}
	else
	{
		Info<< "No U_mean input" << endl;
	}


    return 0;
}
