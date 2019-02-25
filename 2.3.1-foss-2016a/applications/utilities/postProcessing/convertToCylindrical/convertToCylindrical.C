/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    convertToCylindrical

Description
    Converts a velocity field from Cartesian coordinates to cylindrical
    coordinates

Author
     Satish MALIK (Satish.Malik@ideel-factory.fr, svmalik)
     Built on the lines of Bryan Lewis, Penn State University
     Assembled from forum threads by Hrvoje Jasak and Hakkan Nilsson

Usage
     After the simulation has completed, run this application to convert the
     velocity field to cylindrical coordinates (r,theta,z)

     This utility is built to convert "U" and "UMean" into cyclindrical
     coordinates

     The model must be oriented with the x-y plan at the r-theta plane
     and the z-axis must be the center axis of rotation

     For different orientation look at the class cylindricalCS for details
 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cylindricalCS.H"
#include "convertVolFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //- Active flag
            Switch active=true;

    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be converted to cylindrical coordinates  Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );


    timeSelector::addOptions();

#include "setRootCase.H"

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    if (selectedFields.empty())
     {
         Info<< "Please input the fields to be converted into "
               << "cylindrical coordinates.\n"
               << "For help with usage type:  "
               << "convertToCylindrical -help \n\n"
               << "Disabling conversion of fields.\n" << endl;

         active=false;
     }

    if (active)
    {
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createCoordinateSystem.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();


         convertVolField
         (
             cCyl,
             mesh,
             selectedFields
         );
    }
    return 0;
    }
    Info<< "\nEnd\n" << endl;
}



// ************************************************************************* //
