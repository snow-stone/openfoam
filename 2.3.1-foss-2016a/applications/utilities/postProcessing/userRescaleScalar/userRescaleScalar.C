#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("scalingFactor");

	#include "setRootCase.H"
    #include "createTime.H"

	word fieldName(args.additionalArgs()[0]);
	scalar scalingFactor = atof(args.additionalArgs()[1].c_str());

	instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"

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

		    volScalarField s(header, mesh);

            Info<< "before scaling : " << nl
                << "max = " << max(s) << nl
                << "min = " << min(s) << nl
                << endl;

			volScalarField s_nonD
			(
                IOobject
                (
                    string(fieldName+"_nonD"),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
				s/scalingFactor
			);

            Info<< "after scaling : " << nl
                << "max = " << max(s_nonD) << nl
                << "min = " << min(s_nonD) << nl
                << endl;

            s_nonD.write();

	    }
	    else
	    {
		    Info<< "No " << fieldName << " input" << endl;
	    }
    }

    return 0;
}
