#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "cellSet.H"

#include <fstream>

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
	argList::validArgs.append("sliceStore");
	argList::validArgs.append("sliceName");

    #include "setRootCase.H"
    #include "createTime.H"

	word sliceStore(args.additionalArgs()[0]);
	word sliceName(args.additionalArgs()[1]);

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

	label cellN  = slice.size();
	vector ref(0.0, 0.0, 0.0);
	scalarField distance(cellN, 0.0);

	WarningIn("only first 3 elements here") << endl;
	forAll(slice, i)
	{
		distance[i] = Foam::sqrt
				(
				 (mesh.C()[slice[i]] - ref) & (mesh.C()[slice[i]] - ref)
				);
	}

	label l_min = findMin(distance);
    Info << "index for min " << slice[l_min] << "\n"
		 << "coord : "
		 << mesh.C()[slice[l_min]] << "\n"
		 << "distance : "
		 << distance[l_min]
		 << endl;

    return 0;
}

// ************************************************************************* //
