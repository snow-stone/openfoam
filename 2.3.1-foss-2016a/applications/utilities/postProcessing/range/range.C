#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream> // 标准输出，到文件
//using std::ofstream; // 不加这个class ofstream找不到

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
		runTime.setTime(timeDirs[timeI], timeI); // 重新启用runTime来做输出，不选timeDirs[timeI]作为输出
	    Info<< "Time = " << runTime.timeName() << endl;

        IOobject header
        (
            fieldName,
            runTime.timeName(), // 从这里看出上面runTime要定义的必要性，因为要读每个timeName()里面的数据，否则缺省值就是runTime.setTime(timeDirs[0], 0);
            mesh,
            IOobject::MUST_READ
        );

	    if (header.headerOk())
	    {
		    Info<< "headerOk for field : " << fieldName << endl;

		    volScalarField s(header, mesh);

		    Info<< fieldName << " max/min : "
			    << max(s).value() << " "
			    << min(s).value() << endl;

			rangeLog
				<< runTime.timeName() << " "
				<< max(s).value() << " "
				<< min(s).value() << std::endl;
	    }
	    else
	    {
		    Info<< "No " << fieldName << " input" << endl;
	    }
    }

    return 0;
}
