#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#include <fstream> // 标准输出，到文件
//using std::ofstream; // 不加这个class ofstream找不到
//

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
    argList::validArgs.append("patchName");

	#include "setRootCase.H"
    #include "createTime.H"

	word fieldName(args.additionalArgs()[0]);
	word patchName(args.additionalArgs()[1]);
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
		    const polyBoundaryMesh& pp = mesh.boundaryMesh();
			const label patchLabel = pp.findPatchID(patchName);

			if (patchLabel != -1)
			{
				// Downgraded to scalarField. But... boundaryField() has no dimension so
				// it is not really a downgrade. (Meaning it's different when taking max
				// of an internalField(). The output will give its dimension etc...
				scalarField& s_FieldOnPatch = s.boundaryField()[patchLabel];
				Info<< "On patch " << patchName << endl;

				Info<< "field " << fieldName << endl;
				scalarField_simpleStatistics(s_FieldOnPatch);

				const surfaceScalarField& magSf = mesh.magSf();
				//Info<< magSf.boundaryField()[patchLabel] << endl;
				Info<< "field magSf " << endl;
				scalarField_simpleStatistics(magSf.boundaryField()[patchLabel]);

				scalar weightedAverage = sum(s_FieldOnPatch * magSf.boundaryField()[patchLabel])
											/
										 sum(magSf.boundaryField()[patchLabel]);
    			Info<< fieldName 
    				<< " surface weighted average :" << weightedAverage << endl;

				//const surfaceVectorField& Sf = mesh.Sf();
				//Info<< Sf.boundaryField()[patchLabel] << endl;
				
				//const surfaceScalarField& phi = runTime.db().lookupObject<surfaceScalarField>("phi"); Grammar works runTime fail
				//const volVectorField& U = runTime.db().lookupObject<volVectorField>("U"); same...
				//I suppose .... need to read and register...
				//Info << s.boundaryField()[patchLabel] << endl;
				//Info<< "field phi " << endl;
				scalarField_simpleStatistics(s.boundaryField()[patchLabel]);
				//Info<< "sum(phi) " << sum(s.boundaryField()[patchLabel]) << endl;
			}

			if (!noWriting)
			{

				if (patchLabel != -1)
				{
        	        std::ofstream txtOutput
        	        (
        	   	        fileName(
        				    string("postProcessing")/string("txtBoundary_"+patchName+"_"+fieldName+"_"+runTime.timeName())
        				).c_str(),
        		        ios_base::app
        	        );
        
        			forAll(s.boundaryField()[patchLabel], i)
        			{
        			    txtOutput
        				    << s.boundaryField()[patchLabel][i]
        				    << std::endl;
        			}
				}
				else
				{
				    Info<< "No patch named " << patchName << endl;
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
