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
    parabolicVelocityFvPatchVectorField2Dpf

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    parabolicVelocityFvPatchVectorField2Dpf.C
	a 2D version

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef parabolicVelocityFvPatchVectorField2Dpf_H
#define parabolicVelocityFvPatchVectorField2Dpf_H

#include "Random.H" // 2Dp
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class parabolicVelocityFvPatchVectorField2Dpf
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- U_tau
        scalar U_tau_;

        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;

		//- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("parabolicVelocity2Dpf");


    // Constructors

        //- Construct from patch and internal field
        parabolicVelocityFvPatchVectorField2Dpf
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicVelocityFvPatchVectorField2Dpf
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicVelocityFvPatchVectorField2Dpf
        //  onto a new patch
        parabolicVelocityFvPatchVectorField2Dpf
        (
            const parabolicVelocityFvPatchVectorField2Dpf&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicVelocityFvPatchVectorField2Dpf(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicVelocityFvPatchVectorField2Dpf
        (
            const parabolicVelocityFvPatchVectorField2Dpf&,
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
                new parabolicVelocityFvPatchVectorField2Dpf(*this, iF)
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
