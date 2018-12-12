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
#include "fvMesh.H"
#include "timeSelector.H"

#include "fvc.H"

#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main( int argc, char *argv[])
{
    timeSelector::addOptions();
	argList::validArgs.append("vectorFieldName");
    argList::validArgs.append("timeOfAverageField");
    argList::validArgs.append("position");

    #include "setRootCase.H"
    #include "createTime.H"

	word vectorFieldName(args.additionalArgs()[0]);
	word timeOfAverageField(args.additionalArgs()[1]);
	word position(args.additionalArgs()[2]);

	instantList timeDirs = timeSelector::select0(runTime, args);

    IOList<label> labelGroup
    (
        IOobject
        (
            "labelGroup",
            position,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	//output and header
	//labelGroup.write();
	Info<< labelGroup << endl;

	mkDir("userDefinedLog");
	std::ofstream fluctuationLog_x
	(
	   	fileName(string("userDefinedLog")/string("fluctuation_labelGroup_"+vectorFieldName+"_x")).c_str(),
		ios_base::app
	);
	std::ofstream fluctuationLog_y
	(
	   	fileName(string("userDefinedLog")/string("fluctuation_labelGroup_"+vectorFieldName+"_y")).c_str(),
		ios_base::app
	);
	std::ofstream fluctuationLog_z
	(
	   	fileName(string("userDefinedLog")/string("fluctuation_labelGroup_"+vectorFieldName+"_z")).c_str(),
		ios_base::app
	);

    #include "createMesh.H"

	volVectorField mean
	(
	    IOobject
	    (
		    vectorFieldName+"_mean",
		    timeOfAverageField,
		    mesh,
		    IOobject::MUST_READ
	    ),
		mesh
	);

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI); 
	    Info<< "Time = " << runTime.timeName() << endl;

		volVectorField vectorField
		(
        	IOobject
        	(
            	vectorFieldName,
           	    runTime.timeName(),
           		mesh,
            	IOobject::MUST_READ
        	),
			mesh
		);

		volVectorField sPrime = vectorField - mean;

        fluctuationLog_x << runTime.timeName() << " ";
        fluctuationLog_y << runTime.timeName() << " ";
        fluctuationLog_z << runTime.timeName() << " ";
        forAll(labelGroup, i)
        {
            //fluctuationLog << labelGroup[i] << " " ;
            fluctuationLog_x << sPrime.internalField()[labelGroup[i]].component(vector::X) << " " ;
            fluctuationLog_y << sPrime.internalField()[labelGroup[i]].component(vector::Y) << " " ;
            fluctuationLog_z << sPrime.internalField()[labelGroup[i]].component(vector::Z) << " " ;
        }
        fluctuationLog_x << std::endl;
        fluctuationLog_y << std::endl;
        fluctuationLog_z << std::endl;
	}
	
    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
