/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef convertVolFields_H
#define convertVolFields_H

#include "convertVolFieldItems.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
class fvMesh;
    
void convertVolField
(
    const volVectorField& ccCyl,
    const fvMesh& mesh,
    const HashSet<word>& selectedFields
)
{
    // Search for list of objects for this time
    IOobjectList objects(mesh, mesh.time().timeName());

        // Convert volFields
        // ~~~~~~~~~~~~~
    
        convertScalarFieldItem//<scalar>
        (
            ccCyl,
            mesh,
            objects,
            selectedFields
        );

        convertVectorFieldItem//<vector>
        (
            ccCyl,
            mesh,
            objects,
            selectedFields
        );
        
        convertSphericalTensorFieldItem//<sphericalTensor>
        (
            ccCyl,
            mesh,
            objects,
            selectedFields
        );
        
        convertSymmTensorFieldItem//<symmTensor>
        (
            ccCyl,
            mesh,
            objects,
            selectedFields
        );
        
        convertTensorFieldItem//<tensor>
        (
            ccCyl,
            mesh,
            objects,
            selectedFields
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
