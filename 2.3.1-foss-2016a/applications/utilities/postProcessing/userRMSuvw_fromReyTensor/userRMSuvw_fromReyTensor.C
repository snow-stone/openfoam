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

	label nFields = 0;

    forAll(timeDirs, timeI)
	{
	    runTime.setTime(timeDirs[timeI], timeI); 
		Info<< "Time = " << runTime.timeName() << endl;

        volSymmTensorField reyTensor
	    (
       	    IOobject
	   	    (
	       	    "reyTensor",
	   	        runTime.timeName(),
	   		    mesh,
	       	    IOobject::MUST_READ
	   	    ),
		    mesh
	     );

	     volScalarField uu
	     (
	        IOobject
	        (
	            "uu",
	            runTime.timeName(),
	            mesh,
	            IOobject::NO_READ,
	            IOobject::AUTO_WRITE
	        ),
		    reyTensor.component(tensor::XX)
		 );

	     volScalarField vv
	     (
	        IOobject
	        (
	            "vv",
	            runTime.timeName(),
	            mesh,
	            IOobject::NO_READ,
	            IOobject::AUTO_WRITE
	        ),
		    reyTensor.component(tensor::YY)
		);

	    volScalarField ww
	    (
	        IOobject
	        (
	            "ww",
	            runTime.timeName(),
	            mesh,
	            IOobject::NO_READ,
	            IOobject::AUTO_WRITE
	        ),
		    reyTensor.component(tensor::ZZ)
		);

		uuMean += uu;
		vvMean += vv;
		wwMean += ww;

        nFields += 1;

    }

    if (nFields >= 1)
	{
	     Info<< "input field number : " << nFields << endl;
	     uuMean /= nFields;
	     vvMean /= nFields;
	     wwMean /= nFields;
		 uuMean.write();
		 vvMean.write();
		 wwMean.write();
	}
	else
	{
	     Info<< "input field number : 0 " << endl;
	}

    return 0;
}

// ************************************************************************* //
