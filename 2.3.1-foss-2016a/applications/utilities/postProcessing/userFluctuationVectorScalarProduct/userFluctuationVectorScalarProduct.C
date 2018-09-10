#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"


int main(int argc, char *argv[])
{
    Info << "Hello world" << endl;

    timeSelector::addOptions();
    argList::validArgs.append("fieldNameV");
    argList::validArgs.append("fieldNameS");

	#include "setRootCase.H"
    #include "createTime.H"

	word fieldNameV(args.additionalArgs()[0]);
	word fieldNameS(args.additionalArgs()[1]);
	//word outputFileName(string(fieldNameV+fieldNameS));
	word outputFieldName("ut");

	instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"


		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI); 
	    	Info<< "Time = " << runTime.timeName() << endl;

    		IOobject vHeader
    		(
        		string(fieldNameV+"_fluctuation"),
       			runTime.timeName(), 
        		mesh,
        		IOobject::MUST_READ
    		);

        	IOobject sHeader
        	(
        		string(fieldNameS+"_fluctuation"),
            	runTime.timeName(),
            	mesh,
            	IOobject::MUST_READ
        	);

			if (vHeader.headerOk() && sHeader.headerOk())
			{
				volVectorField v(vHeader, mesh);
				volScalarField s(sHeader, mesh);

				volVectorField vs
				(
					IOobject
					(
						outputFieldName,
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					//s * v                     doesn't work
					//cmptMultiply(v, v)        vector by vector 
					mesh,
					Foam::vector::zero
				);

				forAll(mesh.C(), celli)
				{
					vs[celli] = v[celli] * s[celli];
				}

				Info<< "gSum(fluctuation) = " << gSum(vs) << endl;
				Info<< "weightedAverage(fluctuation) = " << vs.weightedAverage(mesh.V()).value() << endl;

				vs.write();
			}
			else
			{
				Info<< "No field @ Time = " << runTime.timeName() << endl;
			}
    	}


    return 0;
}
