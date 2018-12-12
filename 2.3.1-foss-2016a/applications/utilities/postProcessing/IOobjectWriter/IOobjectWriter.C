/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    vorticity

Description
    Calculates and writes the vorticity of velocity field U.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"

#include "IOList.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main( int argc, char *argv[])
{
    argList args(argc, argv);

	Time runTime
	(
	    Time::controlDictName,
		args.rootPath(),
		args.caseName()
	);

    IOList<scalar> IOScalarList 
    (
        IOobject
        (
            "IOSList",
            "constant",
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

	List<scalar> sList(8);
    sList[0] = 7.0;
    sList[1] = 9.0;
    sList[2] = 1.0;
    sList[3] = 2.1;
    sList[4] = 4.0;
    sList[5] = 7.0;
    sList[6] = 4.0;
    sList[7] = 0.0;
	IOList<scalar> IOScalarList1
	(
		IOobject
	   	(
			"IOSList1",
		    "constant",
			runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
		),
		sList
	);

    Info<< "Writing " << IOScalarList.name() << " to " << IOScalarList.objectPath() << endl;

    OFstream os(IOScalarList.objectPath());
    OFstream os1(IOScalarList1.objectPath());

	//writeHeader only
	IOScalarList.writeHeader(os);
	IOScalarList1.writeHeader(os1);

	//output as list
	//os << IOScalarList;
	//os1 << IOScalarList1;

	//output and header
	//IOScalarList.write();
	//IOScalarList1.write();
	
    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
