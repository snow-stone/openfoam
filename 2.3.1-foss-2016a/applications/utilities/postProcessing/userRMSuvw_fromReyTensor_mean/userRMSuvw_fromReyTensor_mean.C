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

	volScalarField uuMean
	(
	    IOobject
	    (
	       "uu_mean",
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedScalar
		(
		    "uu_mean",
		    sqr(dimVelocity),
		    0.0
		)
    );

	volScalarField vvMean
	(
	    IOobject
	    (
	       "vv_mean",
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedScalar
		(
		    "vv_mean",
		    sqr(dimVelocity),
			0.0
		)
	);

	volScalarField wwMean
	(
	    IOobject
	    (
	       "ww_mean",
	       runTime.timeName(),
	       mesh,
	       IOobject::NO_READ,
	       IOobject::AUTO_WRITE
	    ),
		mesh,
		dimensionedScalar
		(
		    "ww_mean",
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

		uuMean = reyTensorMean.component(symmTensor::XX);
		vvMean = reyTensorMean.component(symmTensor::YY);
		wwMean = reyTensorMean.component(symmTensor::ZZ);

        uuMean.write();
        vvMean.write();
        wwMean.write();
    }

    return 0;
}

// ************************************************************************* //
