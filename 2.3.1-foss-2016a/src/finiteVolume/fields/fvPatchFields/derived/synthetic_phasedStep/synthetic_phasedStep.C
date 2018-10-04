/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "synthetic_phasedStep.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

synthetic_phasedStep::synthetic_phasedStep
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    turbulenceSwitch_(false),
    U_tau_(0),
	t0_(0),
	a_(0),
	f_(0),
	phi0_(0),
	s_(1.0),
    n_(1, 0, 0),
    y_(0, 1, 0),
	z_(0, 0, 1)
{}


synthetic_phasedStep::synthetic_phasedStep
(
    const synthetic_phasedStep& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    turbulenceSwitch_(ptf.turbulenceSwitch_),
    U_tau_(ptf.U_tau_),
	t0_(ptf.t0_),
	a_(ptf.a_),
	f_(ptf.f_),
	phi0_(ptf.phi0_),
	s_(ptf.s_),
    n_(ptf.n_),
    y_(ptf.y_),
	z_(ptf.z_)
{}


synthetic_phasedStep::synthetic_phasedStep
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
	turbulenceSwitch_(readBool(dict.lookup("turbulenceSwitch"))),
    U_tau_(readScalar(dict.lookup("U_tau"))),
    t0_(readScalar(dict.lookup("t0"))),
    a_(readScalar(dict.lookup("a"))),
    f_(readScalar(dict.lookup("f"))),
    phi0_(readScalar(dict.lookup("phi0"))),
    s_(readScalar(dict.lookup("s"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
	z_(dict.lookup("z"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn("synthetic_phasedStep(dict)")
            << "n or y or z given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);
	z_ /= mag(z_);
	
	Info<< "Version synthetic_phasedStep." << endl;
	Info<< "TurbulenceSwitch : " << turbulenceSwitch_ << endl;
	Info<< "U_tau = " << U_tau_ << endl;
	Info<< "t0 = " << t0_ << endl;
	Info<< "a = " << a_ << endl;
	Info<< "f = " << f_ << endl;
	Info<< "phi0 = " << phi0_ << endl;
	Info<< "scalable factor = " << s_ << endl;
	Info<<" patchName : "<< patch().name() << endl;

	//Info<<" timeName() : "<< runTime.timeName() << endl; runTime not in this scope
	Info<<" time : "<< this->db().time().timeOutputValue()<< endl;

    evaluate();
}


synthetic_phasedStep::synthetic_phasedStep
(
    const synthetic_phasedStep& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
	turbulenceSwitch_(fcvpvf.turbulenceSwitch_),
    U_tau_(fcvpvf.U_tau_),
	t0_(fcvpvf.t0_),
	a_(fcvpvf.a_),
	f_(fcvpvf.f_),
	phi0_(fcvpvf.phi0_),
	s_(fcvpvf.s_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
	z_(fcvpvf.z_)
{}


// * * * * * * * * * * * * * * * non-Member Functions  * * * * * * * * * * * * * //
	scalar viscousLayer(scalar yPlus){
		return yPlus;
	}

	scalar bufferLayer(scalar yPlus){
		return 5.0 * log(yPlus) - 3.05;
	}

	scalar logLayer(scalar yPlus){
		return 2.5 * log(yPlus) + 5.5;
	}

	scalarField envelopeUrRMS(scalarField yPlus){
		return 		
				(
					-6.38e-09*pow(yPlus,4) 
					+3.14e-06*pow(yPlus,3)
					-5.48e-04*pow(yPlus,2)
					+3.81e-02*pow(yPlus,1)
					-6.20e-02
				 );
	} 

	scalar envelopeUrRMS(scalar yPlus){  // this is slightly negative when yPlus=1
		return 		
				(
					-6.38e-09*pow(yPlus,4) 
					+3.14e-06*pow(yPlus,3)
					-5.48e-04*pow(yPlus,2)
					+3.81e-02*pow(yPlus,1)
					-6.20e-02
				 );
	} 

	scalarField envelopeUzRMS(scalarField yPlus){
		return 		
				(
				 /* requires much more precision when you are a polynome degree 6
				  *
					-2.71e-11*pow(yPlus,6) 
					+1.49e-08*pow(yPlus,5)
					-3.16e-06*pow(yPlus,4)
					+3.23e-04*pow(yPlus,3)
					-1.63e-02*pow(yPlus,2)
					+3.41e-01*pow(yPlus,1)
					+3.37e-01
				*/
					-2.71469e-11*pow(yPlus,6) 
					+1.48853e-08*pow(yPlus,5)
					-3.15553e-06*pow(yPlus,4)
					+3.23651e-04*pow(yPlus,3)
					-1.62823e-02*pow(yPlus,2)
					+3.40719e-01*pow(yPlus,1)
					+3.37271e-01
				 );
	} 

	scalar envelopeUzRMS(scalar yPlus){
		return 	
				(
					-2.71469e-11*pow(yPlus,6) 
					+1.48853e-08*pow(yPlus,5)
					-3.15553e-06*pow(yPlus,4)
					+3.23651e-04*pow(yPlus,3)
					-1.62823e-02*pow(yPlus,2)
					+3.40719e-01*pow(yPlus,1)
					+3.37271e-01
				 );
	} 

    scalar sinusPhased_uTau(scalar U_base, scalar a, scalar omega, scalar t, scalar phi0){
		return U_base * (1 + a * sin(omega * t + phi0));
	}

    scalar stepPhased_uTau(scalar U_base, scalar a, scalar period, scalar t, scalar phi0){
		// imagine when 2pi = phi0, delta_t = period
		// phi0 always positive
		scalar delta_t = phi0 / constant::mathematical::twoPi * period;
		t = t + delta_t;
		while ( t >= period ) {
			t = t - period;
		}
		if (t >= 0.0 && t < period) { // t >= 0 is always true
			if (t < period / 2.0){
				return U_base * (1 + a);
			}
			else
			{
				return U_base * (1 - a);
			}
		}
		else
		{
			Info << " Warning : something has gone wrong " << endl;
			Info << " make it zero " << endl;
			return 0;
		}
	}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void synthetic_phasedStep::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // const
	scalar R(0.004);
    scalar  yesOrNo;
	if (turbulenceSwitch_)
	{
	    yesOrNo = 1.0;
	}
	else
	{
	    yesOrNo = 0.0;
	}
	Info<< " yesOrNo : " << yesOrNo << endl;

	//scalar nu(1.0e-6);
	//scalar U_tau_(0.0218); // to make the bulk velocity to 0.3
	//scalar U_tau_(0.0363); // to make the bulk velocity to 0.5
	//scalar U_tau_(0.0945); // to make the bulk velocity to 1.3
/*
	Info << "envelope Rr(1)" << envelopeUrRMS(1) << endl;
	Info << "envelope Rr(20)" << envelopeUrRMS(20) << endl;
	Info << "envelope Rr(100)" << envelopeUrRMS(100) << endl;
	Info << "envelope Rz(1)" << envelopeUzRMS(1) << endl;
	Info << "envelope Rz(20)" << envelopeUzRMS(20) << endl;
	Info << "envelope Rz(100)" << envelopeUzRMS(100) << endl;
*/

    // oscillating U_tau
	// time relative to t0_
	const scalar t = this->db().time().timeOutputValue() - t0_;
	const scalar a = a_;
	const scalar period = 1.0 / f_;
	const scalar phi0 = constant::mathematical::twoPi * phi0_;
	//U_tau_ = U_tau_ * (1 + a * sin(omega * t + phi0));       this is wrong. U_tau_ is modified and saved every time step...
	//const scalar uTau = U_tau_ * (1 + a * sin(omega * t + phi0));
	//const scalar uTau = sinusPhased_uTau( U_tau_, a, omega, t, phi0);
	const scalar uTau = stepPhased_uTau( U_tau_, a, period, t, phi0);
    Info<< "patch : " << patch().name() << endl;
	Info<< "uTau = " << uTau << endl;

    // Get range and orientation
    // Maybe better to have a fixed const value for ctr
    boundBox bb(patch().patch().localPoints(), true);
    vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

	scalarField coord_pow2 = 4*( sqr((c - ctr) & y_)+sqr((c - ctr) & z_) )/(sqr((bb.max() - bb.min()) & y_));

	scalarField y1 = R * (1 - sqrt(coord_pow2));
	scalarField yPlus = y1 * uTau / 1.0e-6;  // consider replace with nu
	//Info << "here is y : " << y1 << endl;
	//Info << "here is yPlus : " << yPlus << endl;
	Info << "Time : "<< this->db().time().timeOutputValue()<< endl;
	Info << "max of yPlus = " << max(yPlus) << endl; // many lines of output when reconstructing
	Info << "min of yPlus = " << min(yPlus) << endl;

	scalarField envelope_Rz = uTau * envelopeUzRMS(yPlus);
	scalarField envelope_Rr = uTau * envelopeUrRMS(yPlus);

    // Random part
	// Random ranGen_(label(0))
	Random ranGen_(label(this->db().time().timeIndex()));

    Field<vector>& patchField = *this; // change <Type> to vector 

    Field<vector> randomField(this->size()); // change <Type> to vector 
    Field<vector> envelopeField(this->size()); // change <Type> to vector 


    forAll(patchField, facei)
    {
        ranGen_.randomise(randomField[facei]);
	    envelopeField[facei]
			= envelope_Rz[facei] * n_ * s_ ;

		if (yPlus[facei] < 5)
		{
			//Info << "yPlus = " << yPlus[facei] << endl;
			patchField[facei] = n_ * uTau * viscousLayer(yPlus[facei]);
		}
		else if ( 5 <= yPlus[facei] && yPlus[facei] < 30)
		{
			patchField[facei] = n_ * uTau * bufferLayer(yPlus[facei]);
		}
		else if (yPlus[facei] >= 30) // no upper limit here, there should be one
		{
			patchField[facei] = n_ * uTau * logLayer(yPlus[facei]);
		}
		else
		{
			Info << "Exception ! Wrong value of yPlus !" << endl;;
		}
    }

    //ToDo: consider relaxation and large structure (in contrary to point to point perturbation)
    // Correction-factor to compensate for the loss of RMS fluctuation
    // due to the temporal correlation introduced by the alpha parameter.

	randomField = 2*(randomField - 0.5*pTraits<vector>::one);

    patchField = patchField
	    + yesOrNo * cmptMultiply
		  (
		    envelopeField,			
		    randomField
		  );

    fixedValueFvPatchField<vector>::updateCoeffs(); // change <Type> to vector 
}


// Write : this is important ! if new data member is add, this must be be adapted. 
// DecomposePar will write to processor*/timeDir/U using this method.
// You will not want in the decomposed case member list are not complete.
void synthetic_phasedStep::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("turbulenceSwitch")
        << turbulenceSwitch_ << token::END_STATEMENT << nl;
    os.writeKeyword("U_tau")
        << U_tau_ << token::END_STATEMENT << nl;
    os.writeKeyword("t0")
        << t0_ << token::END_STATEMENT << nl;
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("f")
        << f_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi0")
        << phi0_ << token::END_STATEMENT << nl;
    os.writeKeyword("s")
        << s_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, synthetic_phasedStep); // this is important too. maybe creating the real class member

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
