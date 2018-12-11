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

    IOList<word> sliceNumberList
	(
	    IOobject
		(
		    "sliceNumberList",
		    sliceStore,
		    runTime,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE	
		)
	);
	Info << sliceNumberList << endl;

    List<label> labelGroup;

    #include "createMesh.H"
	forAll(sliceNumberList, slicei)
	{
	
		IOList<label> slice
		(
		    IOobject
		    (
			    slicePrefix+sliceNumberList[slicei],
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
	
		label label4min = findMin(distance);

		Info<< "For slicei = " << slicei << endl;
		Info<< "    @slice" << sliceNumberList[slicei] << endl;
		Info<< "    label to memorise : " << slice[label4min] << endl;
		Info<< "    coords : " << mesh.C()[slice[label4min]] << endl;
		Info<< "    offset relative to the first element of slice : " 
			<< (slice[label4min] - slice[0]) << endl;
        Info<< nl;

		labelGroup.append(slice[label4min]);
	}

	Info<< "labelGroup : "
        << labelGroup << endl;

	IOList<label> IOlabelGroup
	(
	    IOobject
		(
		    "labelGroup",
			"constant",   // cannot write to sliceStore or "." ...
			runTime,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
        labelGroup
	);

	//write as a labelList. To be read by another program
	//Info<< "IOlabelGroup : " << IOlabelGroup << endl;
	IOlabelGroup.write();

    return 0;
}

// ************************************************************************* //
