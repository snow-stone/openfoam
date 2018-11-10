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

#include "synthetic_phased_squareDuct.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

synthetic_phased_squareDuct::synthetic_phased_squareDuct
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    turbulenceSwitch_(false),
	t0_(0.),
	a_(0.),
	f_(0.),
	phi0_(0.),
	G_(1.),
    n_(1, 0, 0),
    y_(0, 1, 0),
	z_(0, 0, 1),
	N_(label(1)),
	//w_(1.),
	//h_(1.),
	ctr_(0., 0., 0.)
{}


synthetic_phased_squareDuct::synthetic_phased_squareDuct
(
    const synthetic_phased_squareDuct& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    turbulenceSwitch_(ptf.turbulenceSwitch_),
	t0_(ptf.t0_),
	a_(ptf.a_),
	f_(ptf.f_),
	phi0_(ptf.phi0_),
	G_(ptf.G_),
    n_(ptf.n_),
    y_(ptf.y_),
	z_(ptf.z_),
	N_(ptf.N_),
	//w_(ptf.w_),
	//h_(ptf.h_),
	ctr_(ptf.ctr_)
{}


synthetic_phased_squareDuct::synthetic_phased_squareDuct
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
	turbulenceSwitch_(readBool(dict.lookup("turbulenceSwitch"))),
    t0_(readScalar(dict.lookup("t0"))),
    a_(readScalar(dict.lookup("a"))),
    f_(readScalar(dict.lookup("f"))),
    phi0_(readScalar(dict.lookup("phi0"))),
    G_(readScalar(dict.lookup("G"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
	z_(dict.lookup("z")),
	//w_(dict.lookup("w")),
	//h_(dict.lookup("h")),
	N_(readLabel(dict.lookup("N"))),
	ctr_(dict.lookup("ctr"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn("synthetic_phased_squareDuct(dict)")
            << "n or y or z given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);
	z_ /= mag(z_);
	
	Info<< "Version synthetic_phased_squareDuct." << endl;
	Info<< "TurbulenceSwitch : " << turbulenceSwitch_ << endl;
	Info<< "t0 = " << t0_ << endl;
	Info<< "a = " << a_ << endl;
	Info<< "f = " << f_ << endl;
	Info<< "phi0 = " << phi0_ << endl;
	Info<< "scalable factor G = " << G_ << endl;
	Info<< "geometry of the patch : " << endl;
	//Info<< "  w = " << w_ << endl;
	//Info<< "  h = " << h_ << endl;
	Info<< "patch center : " << ctr_ << endl;
	Info<<" patchName : "<< patch().name() << endl;

	Info<<" time : "<< this->db().time().timeOutputValue()<< endl;

    evaluate();
}


synthetic_phased_squareDuct::synthetic_phased_squareDuct
(
    const synthetic_phased_squareDuct& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
	turbulenceSwitch_(fcvpvf.turbulenceSwitch_),
	t0_(fcvpvf.t0_),
	a_(fcvpvf.a_),
	f_(fcvpvf.f_),
	phi0_(fcvpvf.phi0_),
	G_(fcvpvf.G_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
	z_(fcvpvf.z_),
	N_(fcvpvf.N_),
	//w_(fcvpvf.w_),
	//h_(fcvpvf.h_),
	ctr_(fcvpvf.ctr_)
{}

// * * * * * * * * * * * * * * * non-Member Functions  * * * * * * * * * * * * * //

// i-th term of the series of the non-dimensioned part for series for Hagen-poiseuille square duct
scalar squareDuctUz_i(label i, scalar w, scalar h, scalar x, scalar y)
{
	label n = 2*i -1;
	scalar Uz_i =  1.0/pow(n,3)
			   * ( 1.0
	             - cosh(n * constant::mathematical::pi * x / h)
				   /
				   cosh(n * constant::mathematical::pi * w / (2.0 * h))
				 )
			   * sin(n * constant::mathematical::pi * (y+0.5) / w);
	
	return Uz_i;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void synthetic_phased_squareDuct::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

	// const
	//label N = 20;
	scalar w = 0.008;
	scalar h = 0.008;
	//vector ctr(-0.08, 0.0, 0.0);
	//scalar G = 7.372 * 1.0 / sqr(w);

	// time relative to t0_
	const scalar t = this->db().time().timeOutputValue() - t0_;
	const scalar a = a_;
	const scalar omega = constant::mathematical::twoPi * f_;
	const scalar phi0 = constant::mathematical::twoPi * phi0_;
	const scalar sinusPart =  1 + a * sin(omega * t + phi0);

	const vectorField& c = patch().Cf();
	vectorField coord_local = c - ctr_;
    vectorField& patchField = *this;

    forAll(patchField, facei)
    {
		scalar Uz = 0.0;
		scalar x  = coord_local[facei] & y_; // direction vector for x 
		scalar y  = coord_local[facei] & z_; // direction vector for y

		// sum of the series until order N_
		for(label i = 0; i <= N_; i++)
		{
			//Info << "i = " << i << endl;
		    Uz = Uz + squareDuctUz_i(i, w, h, x, y);
		}
	    patchField[facei] = n_ 
				* G_ * 4 * sqr(h) / pow(constant::mathematical::pi,3) * Uz 
				* sinusPart ;
    }

    fixedValueFvPatchField<vector>::updateCoeffs(); // change <Type> to vector 
}


// Write : this is important ! if new data member is add, this must be be adapted. 
// DecomposePar will write to processor*/timeDir/U using this method.
// You will not want in the decomposed case member list are not complete.
void synthetic_phased_squareDuct::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("turbulenceSwitch")
        << turbulenceSwitch_ << token::END_STATEMENT << nl;
    os.writeKeyword("t0")
        << t0_ << token::END_STATEMENT << nl;
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("f")
        << f_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi0")
        << phi0_ << token::END_STATEMENT << nl;
    os.writeKeyword("G")
        << G_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("N")
        << N_ << token::END_STATEMENT << nl;
    //os.writeKeyword("w")
    //    << w_ << token::END_STATEMENT << nl;
    //os.writeKeyword("h")
    //    << h_ << token::END_STATEMENT << nl;
    os.writeKeyword("ctr")
        << ctr_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, synthetic_phased_squareDuct); // this is important too. maybe creating the real class member

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
