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
    turbulentFvPatchVectorField_switch_ScalableCmptStream

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    turbulentFvPatchVectorField_switch_ScalableCmptStream.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef turbulentFvPatchVectorField_switch_ScalableCmptStream_H
#define turbulentFvPatchVectorField_switch_ScalableCmptStream_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "Random.H" //userDefined

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class turbulentFvPatchVectorField_switch_ScalableCmptStream Declaration
\*---------------------------------------------------------------------------*/

class turbulentFvPatchVectorField_switch_ScalableCmptStream
:
    public fixedValueFvPatchVectorField
{
    // Private data
	
		//- turbulent switch
		bool turbulenceSwitch_;

        //- U_tau
        scalar U_tau_;

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
    TypeName("turbulentFvPatchVectorField_switch_ScalableCmptStream");


    // Constructors

        //- Construct from patch and internal field
        turbulentFvPatchVectorField_switch_ScalableCmptStream
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentFvPatchVectorField_switch_ScalableCmptStream
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentFvPatchVectorField_switch_ScalableCmptStream
        //  onto a new patch
        turbulentFvPatchVectorField_switch_ScalableCmptStream
        (
            const turbulentFvPatchVectorField_switch_ScalableCmptStream&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentFvPatchVectorField_switch_ScalableCmptStream(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentFvPatchVectorField_switch_ScalableCmptStream
        (
            const turbulentFvPatchVectorField_switch_ScalableCmptStream&,
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
                new turbulentFvPatchVectorField_switch_ScalableCmptStream(*this, iF)
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
