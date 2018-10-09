#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"


//volScalarField userCalcNu(volScalarField strainRate) 
//--> FOAM FATAL ERROR: 
//LHS and RHS of + have different dimensions
//     dimensions : [0 0 0 0 0 0 0] + [0 0 -2 0 0 0 0]
scalarField userCalcNu(scalarField strainRate)
{
	scalar nuInf_ = 2e-6;
	scalar nu0_   = 3e-4;
	scalar k_     = 1;
	scalar n_     = 0.326;
	scalar a_     = 2.0;  // defaut BirdCarreau

    return
		nuInf_
	  + (nu0_ - nuInf_)
	    *pow(scalar(1) + pow(k_*strainRate, a_), (n_ - 1.0)/a_);
}

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
	argList::validArgs.append("velocityFieldName");

	#include "setRootCase.H"
    #include "createTime.H"

	word velocityFieldName(args.additionalArgs()[0]);
	word outputFieldName(string("nu_"+velocityFieldName));

	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

		volVectorField velocityField
		(
        	IOobject
        	(
            	velocityFieldName,
            	runTime.timeName(),
           		mesh,
            	IOobject::MUST_READ
        	),
			mesh
		);

        volScalarField strainRate
        (
           IOobject
            (
                "strainRate_"+velocityFieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
			Foam::sqrt(2.0)*mag(symm(fvc::grad(velocityField)))
		);

		volScalarField nu
		(
		    IOobject
		    (
                outputFieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
			),
            //userCalcNu(strainRate) Not a volScalarField !
		    mesh,
			dimensionedScalar
			(
				outputFieldName,
			    dimless,//sqr(dimLength)/dimTime,
			    scalar(0.)
			)
		);

		nu.internalField()=userCalcNu(strainRate);

		nu.write();

    }

    return 0;
}


// ************************************************************************* //
