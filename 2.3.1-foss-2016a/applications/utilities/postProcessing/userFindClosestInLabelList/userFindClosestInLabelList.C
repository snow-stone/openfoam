#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
	argList::validArgs.append("sliceStore");
	argList::validArgs.append("slicePrefix");

    #include "setRootCase.H"
    #include "createTime.H"

	word sliceStore(args.additionalArgs()[0]);
	word slicePrefix(args.additionalArgs()[1]);

    instantList timeDirs = timeSelector::select0(runTime, args);

	WarningIn("User defined here : only one reference vector !") << endl;
	IOList<vector> refVectors
	(
        IOobject
		(
			"refVectors",
			sliceStore,
			runTime,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	vector ref = refVectors[0];
	Info << "ref : " << ref << endl;

    IOList<word> sliceNumber
	(
	    IOobject
		(
		    "sliceNumber",
		    sliceStore,
		    runTime,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		)
	);
	Info << sliceNumber << endl;

    #include "createMesh.H"

	forAll(sliceNumber, i)
	{
	
		IOList<label> slice
		(
		    IOobject
		    (
			    slicePrefix+sliceNumber[i],
				sliceStore,
				mesh,
			    IOobject::MUST_READ,
			    IOobject::NO_WRITE	
			)	
		);
	
		label cellN  = slice.size();
		scalarField distance(cellN, 0.0);
	
		//WarningIn("only first 3 elements here") << endl;
		forAll(slice, i)
		{
			distance[i] = Foam::sqrt
					(
					 (mesh.C()[slice[i]] - ref) & (mesh.C()[slice[i]] - ref)
					);
		}
	
		label l_min = findMin(distance);
	    Info << "index for min " << slice[l_min] << nl
			 << "coord : "
			 << mesh.C()[slice[l_min]] << nl
			 << "distance : "
			 << distance[l_min]
			 << endl;
	}

    return 0;
}

// ************************************************************************* //
