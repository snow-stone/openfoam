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

Class
    synthetic_phasedStep

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    synthetic_phasedStep.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef synthetic_phasedStep_H
#define synthetic_phasedStep_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "Random.H" //userDefined

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class synthetic_phasedStep Declaration
\*---------------------------------------------------------------------------*/

class synthetic_phasedStep
:
    public fixedValueFvPatchVectorField
{
    // Private data
	
		//- turbulent switch
		bool turbulenceSwitch_;

        //- U_tau
        scalar U_tau_;

		//- startTime
		scalar t0_;

		//- amplitude
		scalar a_;  // fraction of U_bulk

		//- frequency
		scalar f_;

		//- initial phase
		scalar phi0_; // phi0=phi0_*2Pi

		//- scale for turbulence level
		scalar s_;

        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;

		//- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("synthetic_phasedStep");


    // Constructors

        //- Construct from patch and internal field
        synthetic_phasedStep
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        synthetic_phasedStep
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given synthetic_phasedStep
        //  onto a new patch
        synthetic_phasedStep
        (
            const synthetic_phasedStep&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new synthetic_phasedStep(*this)
            );
        }

        //- Construct as copy setting internal field reference
        synthetic_phasedStep
        (
            const synthetic_phasedStep&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new synthetic_phasedStep(*this, iF)
            );
        }


    // Member functions

        //- Return U_tau value
        scalar& U_tau()
        {
            return U_tau_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }

        //- Return z direction
        vector& z()
        {
            return z_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
