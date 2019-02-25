#include "createMesh.H"

// defining cylindrical coordinate system

    Info<<" creating cylindrical system (r, teta, z) from cartesian system\n"
               " with r = x and z = z\n" << endl;

    cylindricalCS cyl
    (
        "cylindricalCS",
        point(0, 0, 0),   // center point of ccs
        vector(0, 0, 1), // axis of ccs
        vector(1, 0, 0), // base axis for cylindrical angle or start of angle
        false               // keyword false:
    );                         // Use radians since cos and sin work with radians

//Read the vectors pointing from the origin of the mesh to the cell centers
volVectorField cc
    (
        IOobject
        (
            "cc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C()
    );

//The vector field cc will now be transformed into cylindrical coordinates
volVectorField cCyl
    (
           IOobject
           (
              "cCyl",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
           ),
           mesh,
           vector (0,0,0)
    );

    cCyl.internalField() = cyl.localVector(cc.internalField());
    forAll (cCyl.boundaryField(), patchI)
    {
           cCyl.boundaryField()[patchI] =
               cyl.localVector(cc.boundaryField()[patchI]);
    }

// write the vector field cc of cylindrical coordinates

    Info<< " writing the cylindrincal coordinate system"
                " to the time directory" << endl;

    cCyl.write();

