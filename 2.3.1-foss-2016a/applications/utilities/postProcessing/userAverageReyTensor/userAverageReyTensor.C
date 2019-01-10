#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

	instantList timeDirs = timeSelector::select0(runTime, args);
    Info<< " before creating mesh " << endl;
	#include "createMesh.H"

    volSymmTensorField reyTensorMean
    (
        IOobject
        (
           "reyTensor_mean",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
           "reyTensor_mean",
           sqr(dimVelocity),
           symmTensor::zero
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

            reyTensorMean += reyTensor;
            nFields += 1;
	    }
        
        if (nFields >= 1)
        {
            Info<< "Iteration ended" << endl;
            Info<< "runTime.timeName() : " << runTime.timeName() << endl;
            Info<< "nFields : " << nFields << endl;
            reyTensorMean /= nFields;
            reyTensorMean.write(); // write only at the last time step
        }

    return 0;
}

// ************************************************************************* //
