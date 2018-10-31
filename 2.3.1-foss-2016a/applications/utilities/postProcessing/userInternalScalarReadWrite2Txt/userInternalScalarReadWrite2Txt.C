#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream> // 标准输出，到文件
//using std::ofstream; // 不加这个class ofstream找不到

void scalarField_simpleStatistics(const scalarField& s)
{
    Info<< "Simple Statistics on field : " << endl;
	Info<< "size of field : " << s.size() << endl;
	Info<< "min : " << min(s) << endl;
	Info<< "max : " << max(s) << endl;
	Info<< "athemtic mean : " << sum(s)/s.size() << endl;
	Info<< nl;
}

int main(int argc, char *argv[])
{
    Info << "Hello world" << endl;

    timeSelector::addOptions();
	Foam::argList::addBoolOption
	(
	    "noWrite",
	    "no output data to file"
    );
    argList::validArgs.append("fieldName");

	#include "setRootCase.H"
    #include "createTime.H"

	word fieldName(args.additionalArgs()[0]);
	bool noWriting = args.optionFound("noWrite");

	instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"

	if(!noWriting)
	{
	    mkDir("postProcessing");
	}

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

			Info<< "internalField " << fieldName << endl;
			scalarField_simpleStatistics(s.internalField());

			Info<< fieldName 
				<< " volumic average :"
				//<< sum(mesh.V()*s)/sum(mesh.V()) << endl;
				<< s.weightedAverage(mesh.V()) << endl;

			if (!noWriting)
			{
    	        std::ofstream txtOutput
    	        (
    	   	        fileName(
    				    string("postProcessing")/string("txtInternal_"+fieldName+"_"+runTime.timeName())
    				).c_str(),
    		        ios_base::app
    	        );
    
    			forAll(s.internalField(), i)
    			{
    			    txtOutput
    				    << s.internalField()[i]
    				    << std::endl;
    			}
			}
	    }
	    else
	    {
		    Info<< "No " << fieldName << " input" << endl;
	    }
    }

    return 0;
}
