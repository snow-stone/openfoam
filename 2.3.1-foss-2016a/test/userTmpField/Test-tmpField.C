/*---------------------------------------------------------------------------*\
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
    tmpFieldTest

Description
    Tests for possible memory leaks in the tmp<Field> algebra.

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    scalarField f1(5, 1.0), f2(5, 2.0), f3(3, 3.0);
	//scalarField x(1.0, 2.0, 3.0);
	scalarField x(3);
	forAll(x , i)
	{
		x[i] = i;
	}

	Info<< "f1 " << f1 << endl;
	Info<< "f2 " << f2 << endl;
	Info<< "f3 " << f3 << endl;
	Info<< "x "  << x  << endl;

	scalarField f12 = f1 * f2;
	Info<< "f12 " << f12 << endl;

	Info<< "f3*x " << f3 * x << endl;
	Info<< "cmptMultiply(f3,x) " << cmptMultiply(f3, x) << endl;

    Info<< "end" << endl;
}


// ************************************************************************* //
