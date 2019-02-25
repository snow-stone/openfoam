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

#ifndef convertVolFieldItems_H
#define convertVolFieldItems_H

#include "GeometricField.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void convertScalarFieldItem
(
    const volVectorField& cCyl,
    const fvMesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
        {
//            const word& fieldName = fieldIter()->name();

            if (selectedFields.found(fieldIter()->name()))
            {
                Info << " This functionality is still not active."
                             " Cannot convert scalar field " << fieldIter()->name()
                       << " into cylindrical coordinates " << endl;
            }
        }
}

void convertVectorFieldItem
(
    const volVectorField& cCyl,
    const fvMesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    typedef GeometricField<vector, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.found(fieldName))
        {
            Info<< " Reading field " << fieldName << endl;

            volVectorField fieldToConvert
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
             );

            const word& fieldNameCC = fieldName + "cyl";

                volVectorField fieldCC
                (
                    IOobject
                    (
                        fieldNameCC,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    1*fieldToConvert
                );

                Info<< "    converting " << fieldName
                      << " to cylindrical coordinates" << endl;


                // transformation of velocity field U from
                // cartesian -> cylindrical
                //   [Ux Uy]*[cos(thetha)    -sin(theta)]
                //                [sin(thetha)     cos(theta)]


                // the radial component
                fieldCC.replace
                (
                    vector::X,
                    (
                        fieldToConvert.component(vector::X)
                       *cos(cCyl.component(vector::Y))
                     )
                    +
                    (
                        fieldToConvert.component(vector::Y)
                       *sin(cCyl.component(vector::Y))
                    )
                );
                // the azimuthal component
                fieldCC.replace
                (
                    vector::Y,
                    (-1*fieldToConvert.component(vector::X)
                        *sin(cCyl.component(vector::Y)))
                    +(fieldToConvert.component(vector::Y)
                        *cos(cCyl.component(vector::Y)))
                );
                // the axial component
                fieldCC.replace
                (
                    vector::Z,
                    fieldToConvert.component(vector::Z)
                );

                // write the field in cyclindrical coordinates

                Info<< "    writing " << fieldNameCC
                      << " to the time directory" << endl;

                fieldCC.write();
        }
    }
}

void convertSymmTensorFieldItem
(
    const volVectorField& cCyl,
    const fvMesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    // TODO
}

void convertSphericalTensorFieldItem
(
    const volVectorField& cCyl,
    const fvMesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    typedef GeometricField<sphericalTensor, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
        {
//            const word& fieldName = fieldIter()->name();
            if (selectedFields.found(fieldIter()->name()))
            {
                Info << " This functionality is still not active."
                             " Cannot convert spherical tensor field " << fieldIter()->name()
                       << " into cylindrical coordinates " << endl;
            }
        }
}

void convertTensorFieldItem
(
    const volVectorField& cCyl,
    const fvMesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    typedef GeometricField<tensor, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
        {
//            const word& fieldName = fieldIter()->name();

            if (selectedFields.found(fieldIter()->name()))
            {
                Info << " This functionality is still not active."
                             " Cannot convert tensor field " << fieldIter()->name()
                       << " into cylindrical coordinates " << endl;
            }
        }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
