#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream> // 标准输出，到文件

int main(int argc, char *argv[])
{
    Info << "Hello world" << endl;

    timeSelector::addOptions();
    argList::validArgs.append("velocityFieldName");

	#include "setRootCase.H"
    #include "createTime.H"

	word velocityFieldName(args.additionalArgs()[0]);

	// const
	word nuFieldName = "nu";
	word outputFieldName = string("LocalRe_")
			             + string("mag")
						 + velocityFieldName;
	scalar D = 0.008;

	instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"

	mkDir("userDefinedLog");
	std::ofstream rangeLog
	(
	   	fileName(string("userDefinedLog")/string("rangeForField_"+outputFieldName)).c_str(),
		ios_base::app
	);

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

        IOobject velocityHeader
        (
            velocityFieldName,
            runTime.timeName(), 
            mesh,
            IOobject::MUST_READ
        );

        IOobject nuHeader
        (
            nuFieldName,
            runTime.timeName(), 
            mesh,
            IOobject::MUST_READ
		);

	    if (velocityHeader.headerOk() && nuHeader.headerOk())
	    {
		    Info<< "headerOk for field : "
				<< velocityFieldName
				<< " and "
				<< nuFieldName
				<< endl;

		    volVectorField U(velocityHeader, mesh);
		    volScalarField nu(nuHeader, mesh);

			volScalarField LocalRe
			(
			    IOobject
			    (
					outputFieldName,
				    runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE	
				),
				mag(U)*D/nu	
			);

		    Info<< outputFieldName << " max/min : "
			    << max(LocalRe).value() << " "
			    << min(LocalRe).value() << endl;

			rangeLog
				<< runTime.timeName() << " "
				<< max(LocalRe).value() << " "
				<< min(LocalRe).value() << std::endl;

			LocalRe.write();
	    }
	    else
	    {
		    Info<< "No " 
				<< velocityFieldName 
				<< " or "
				<< nuFieldName
				<< " input" << endl;
	    }
    }

    return 0;
}
