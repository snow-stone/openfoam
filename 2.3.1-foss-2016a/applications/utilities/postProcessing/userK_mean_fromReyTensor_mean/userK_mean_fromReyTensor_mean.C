#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

	volScalarField kMean
	(
	    IOobject
	    (
	       "k_mean",
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedScalar
		(
		    "k_mean",
		    sqr(dimVelocity),
		    0.0
		)
    );

    forAll(timeDirs, timeI)
	{
	    runTime.setTime(timeDirs[timeI], timeI); 
		Info<< "Time = " << runTime.timeName() << endl;

        volSymmTensorField reyTensorMean
	    (
       	    IOobject
	   	    (
	       	    "reyTensor_mean",
	   	        runTime.timeName(),
	   		    mesh,
	       	    IOobject::MUST_READ
	   	    ),
		    mesh
	    );

		kMean = 0.5 * tr(reyTensorMean);
		Info<<"Writing kMean"<<endl;

        kMean.write();
    }

    return 0;
}

// ************************************************************************* //
