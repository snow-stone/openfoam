#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream> // 标准输出，到文件

int main(int argc, char *argv[])
{
    Info << "Hello world" << endl;

    timeSelector::addOptions();
    argList::validArgs.append("fieldName");

	#include "setRootCase.H"
    #include "createTime.H"

	word fieldName(args.additionalArgs()[0]);

	// create list of time steps based on the -time argument
	instantList timeDirs = timeSelector::select0(runTime, args);
    //set the runTime at the first time step of the given range
    //runTime.setTime(timeDirs[0], 0);

    // now: create mesh based on the first time step of the given range 
	#include "createMesh.H"

	mkDir("userDefinedLog");
	std::ofstream rangeLog
	(
	   	fileName(string("userDefinedLog")/string("rangeForField_"+fieldName)).c_str(),
		ios_base::app
	);

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

        IOobject header
        (
            fieldName,
            runTime.timeName(), 
            mesh,
            IOobject::MUST_READ
        );

	    if (header.headerOk())
	    {
		    Info<< "headerOk for field : " << fieldName << endl;

		    volVectorField v(header, mesh);

			volVectorField vorticity
			(
                IOobject
                (
                    string("vorticity_"+fieldName),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
				fvc::curl(v)
			);

		    Info<< "mag("
				<< fieldName 
				<< ")"
				<< " max/min : "
			    << max(mag(vorticity)).value() << " "
			    << min(mag(vorticity)).value() << endl;

			rangeLog
				<< runTime.timeName() << " "
				<< max(mag(vorticity)).value() << " "
				<< min(mag(vorticity)).value() << std::endl;

			vorticity.write();
	    }
	    else
	    {
		    Info<< "No " << fieldName << " input" << endl;
	    }
    }

    return 0;
}
