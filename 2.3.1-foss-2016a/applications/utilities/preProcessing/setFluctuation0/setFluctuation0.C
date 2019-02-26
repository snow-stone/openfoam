/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
------------------------------------------------------------------------------- License This file is part of OpenFOAM.

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
    averageAlongAxis

Description
    Average field along a chosen axis.

    Assuming that the mesh is periodic in the direction of the axis.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ReadFields.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"
#include "OFstream.H"

#include "OSspecific.H"
#include "Random.H"

// * * * * * * * * * * * * * * * non-Member Functions  * * * * * * * * * * * * * //
namespace Foam
{

	scalar viscousLayer(scalar yPlus){
		return yPlus;
	}

	scalar bufferLayer(scalar yPlus){
		return 5.0 * log(yPlus) - 3.05;
	}

	scalar logLayer(scalar yPlus){
		return 2.5 * log(yPlus) + 5.5;
	}

	scalarField envelopeUrRMS(scalarField yPlus){
		return 		
				(
					-6.38e-09*pow(yPlus,4) 
					+3.14e-06*pow(yPlus,3)
					-5.48e-04*pow(yPlus,2)
					+3.81e-02*pow(yPlus,1)
					-6.20e-02
				 );
	} 

	scalar envelopeUrRMS(scalar yPlus){  // this is slightly negative when yPlus=1
		return 		
				(
					-6.38e-09*pow(yPlus,4) 
					+3.14e-06*pow(yPlus,3)
					-5.48e-04*pow(yPlus,2)
					+3.81e-02*pow(yPlus,1)
					-6.20e-02
				 );
	} 

	scalarField envelopeUzRMS(scalarField yPlus){
		return 		
				(
				 /* requires much more precision when you are a polynome degree 6
				  *
					-2.71e-11*pow(yPlus,6) 
					+1.49e-08*pow(yPlus,5)
					-3.16e-06*pow(yPlus,4)
					+3.23e-04*pow(yPlus,3)
					-1.63e-02*pow(yPlus,2)
					+3.41e-01*pow(yPlus,1)
					+3.37e-01
				*/
					-2.71469e-11*pow(yPlus,6) 
					+1.48853e-08*pow(yPlus,5)
					-3.15553e-06*pow(yPlus,4)
					+3.23651e-04*pow(yPlus,3)
					-1.62823e-02*pow(yPlus,2)
					+3.40719e-01*pow(yPlus,1)
					+3.37271e-01
				 );
	} 

	scalar envelopeUzRMS(scalar yPlus){
		return 	
				(
					-2.71469e-11*pow(yPlus,6) 
					+1.48853e-08*pow(yPlus,5)
					-3.15553e-06*pow(yPlus,4)
					+3.23651e-04*pow(yPlus,3)
					-1.62823e-02*pow(yPlus,2)
					+3.40719e-01*pow(yPlus,1)
					+3.37271e-01
				 );
	} 

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("wallModel");
    argList::addBoolOption("noMeanFields");
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"
#   include "readTransportProperties.H"


    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "Computing wall quantities for time " << runTime.timeName() << endl;
        
        volScalarField r
        (
           IOobject
           (
               "r",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
               ),
               mesh,
			   dimensionedScalar("r", dimLength, scalar(0.0))
         );

        volVectorField U1
        (
           IOobject
           (
               "U1",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
               ),
               mesh,
			   dimensionedVector("U1", dimLength/dimTime, vector(0, 0, 0))
         );
/*	
	//This block reads a labelListList
	//face is a list of faces and face it self is a list so polyMesh/faces is an example of
	//labelListList

	//labelListIOList proximityCells  // this gives error labelListIOList is not in this scope
	IOList<labelList> proximityCells  // while this class is well recognized
	(
		IOobject
		(
			"proximityCells",
			runTime.time().constant(), 
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	Info << "size : " << proximityCells.size() << endl;
	Info << "size : " << proximityCells[0] << endl;
	Info << "size : " << proximityCells[1] << endl;
*/

//	IOList<label> myCells(
//    	IOobject(
//                 "myCells",
//                 runTime.time().constant(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//            )
//     );
        // Get index of patch

        const word w0("movingWall");
		const word w1("outlet0");
		const word w2("wall");

		label nb_fvPatch = mesh.boundaryMesh().size();

		const label w0PatchID = mesh.boundaryMesh().findPatchID(w0);
		const label w1PatchID = mesh.boundaryMesh().findPatchID(w1);
		const label w2PatchID = mesh.boundaryMesh().findPatchID(w2);

		Info<< "Size of fvBoundryMesh = "
			<< nb_fvPatch << "\n"<< endl;
		Info<< "PatchID for "
			<< w0
			<< " is "
			<< w0PatchID << "\n" << endl;
		Info<< "PatchID for "
			<< w1
			<< " is "
			<< w1PatchID << "\n" << endl;
		Info<< "PatchID for "
			<< w2
			<< " is "
			<< w2PatchID << "\n" << endl;


        // Distance to the wall for the near-wall cells
        volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();

        const fvPatchList& patches = mesh.boundary();

		//scalarField r_in = r.internalField();
		//Info << "r_in = " << r_in << endl;

		//forAll(r_in, celli)
		//Info << "Nb of cells : " << mesh.nCells() << endl;
		Info << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//Info << "cell center : " << mesh.C() << endl;
		//Info << "cell center : " << mesh.boundaryMesh().C() << endl;
		Info << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        vector ctr(0, 0, 0);
		Info << "Calculating distance from the axis z" << endl;
		Info << "Calculating distance from the axis z" << endl;
		Info << "Calculating distance from the axis z" << endl;
		forAll(r.mesh().C(),celli)
		{
			scalar x,y;
			//Info << "cell" << celli << " : x = " << r.mesh().C()[celli].x() << endl;;
			x = r.mesh().C()[celli].x();
			y = r.mesh().C()[celli].y();
			//z = r.mesh().C()[celli].z();
            //r_in = mag(ctr - vector(x, y, z));
            r.internalField()[celli] = mag(ctr - vector(x, y, 0));
		}
		//Info << "r_in = " << r_in << endl;
		//Info << "r = " << r << endl;
		scalar R(0.004);
		scalar uTau(0.045);
		//scalarField yPlus = R - r.internalField();
		scalarField yPlus = (R - r.internalField()) * uTau / 1.0e-6;

        Info << "max of yPlus = " << max(yPlus) << endl;
		Info << "min of yPlus = " << min(yPlus) << endl;

		//scalarField envelope_Rz = Foam::sqrt(3.0) * uTau * Foam::envelopeUzRMS(yPlus); // Foam::envelope... ????
		//scalarField envelope_Rr = Foam::sqrt(3.0) * uTau * Foam::envelopeUrRMS(yPlus);
		scalarField envelope_Rz = Foam::sqrt(3.0) * uTau * envelopeUzRMS(yPlus); // Foam::envelope... ????
		scalarField envelope_Rr = Foam::sqrt(3.0) * uTau * envelopeUrRMS(yPlus);

		vector n_ (0, 0, 0);
		Info << "this is the non-mean flow version !! n_ = " << n_ << endl;
		vectorField& internalU = U1.internalField();
		//Random ranGen_(label(this->db().time().timeIndex()));
		Random ranGen_(label(runTime.time().timeIndex()));
		vectorField randomField(U1.internalField().size());
		vectorField envelopeField(U1.internalField().size());
		//Info << "runTime.time() : " << runTime.time() << endl;
		Info << "runTime.time().timeIndex() : " << runTime.time().timeIndex() << endl;
		forAll(U1.mesh().C(),celli)
		{
			//internalU[celli] = vector(0, 0, 20); // assignment failed. no write on U1
			///U1.internalField()[celli] = vector(0, 0, 20); // assignment succeed. write uniform internalField on U1
			//U1.internalField()[celli] = vector(0, 0, 20); // assignment succeed. write uniform internalField on U1

			if (yPlus[celli] < 5)
			{
				//Info << "yPlus = " << yPlus[facei] << endl;
				U1.internalField()[celli] = n_ * uTau * viscousLayer(yPlus[celli]);
			}
			else if ( 5 <= yPlus[celli] && yPlus[celli] < 30)
			{
				U1.internalField()[celli] = n_ * uTau * bufferLayer(yPlus[celli]);
			}
			else if (yPlus[celli] >= 30) // no upper limit here, there should be one
			{
				U1.internalField()[celli] = n_ * uTau * logLayer(yPlus[celli]);
			}
			else
			{
				Info << "Exception ! Wrong value of yPlus !" << endl;;
			}

			ranGen_.randomise(randomField[celli]);
			envelopeField[celli] = vector(	envelope_Rr[celli],
							envelope_Rr[celli],
							envelope_Rz[celli]       );
			/*
			envelopeField[celli] = vector(	envelope_Rr[celli]/2.0,
							envelope_Rr[celli]/2.0,
							envelope_Rz[celli]       );
			*/

		}

		randomField = 2*(randomField - 0.5*pTraits<vector>::one);

		U1.internalField() =
				U1.internalField()
			+	cmptMultiply
				(
				 	envelopeField,
					randomField
				 );

        fileName writePathRoot("./V3_csv/computeWallQuantities"/runTime.timeName());
        mkDir(writePathRoot);

        // Go through patches
        forAll(patches, patchI)
        {
            const fvPatch& currPatch = patches[patchI];

            if (currPatch.name() == w0)
            {
            	Info << "entering patch " << w0
            			<< " size = " << currPatch.size() << endl;
            	Info << "entering patch.Cf() " << w0
            			<< " size = " << currPatch.Cf().size() << endl;

        		labelList myList;
/*
				// This block works, well... the initialization for a list is nasty
        		labelList lst;
				//wrong answers
				//lst = {1, 2, 3, 4};
				//lst = {label 1, label 4};
				//lst = (label 1, label 2, label 3, label 4); // wrong
				//lst = (1, 2, 3, 4);  // gives warning only, check size of lst => 0. Wrong answser
				label a = 279;
				label b = 289;
				label c = 299;
				lst.append(a);
				lst.append(b);
				lst.append(c);
				Info << "lst size : " << lst.size() << endl;
				forAll(lst, i)
				{
      				Info<< "lst : " << lst[i] << " " <<  currPatch.Cf()[lst[i]].component(vector::Y) << endl;
				}
*/
        		scalar x,y,z;
        		vector ctr(0, -0.2, 0);

            	forAll(currPatch, faceI)
            	{
            		// this is for internal mesh and uses the global label list
            		// Info << "coordinates : " << mesh.C()[faceI] << endl;

            		// this is for boundary mesh
					// get the labelList of cells satisfying the criteria below
            		if (currPatch.Cf()[faceI].component(vector::X) < 0.0001
            				&&
						currPatch.Cf()[faceI].component(vector::X) > 0	)
            		{
                		//Info << "coordinates : " << faceI << " " << currPatch.Cf()[faceI] << endl;
                		myList.append(faceI);
            		}

            		x = currPatch.Cf()[faceI].component(vector::X);
            		y = currPatch.Cf()[faceI].component(vector::Y);
            		z = currPatch.Cf()[faceI].component(vector::Z);

            		r.boundaryField()[patchI][faceI] = mag(ctr - vector(x, y, z));
            	}
            	//Info<< "myList = " << myList << endl;

            	mkDir(writePathRoot/currPatch.name());
                //OFstream yCoorFile(
                //        fileName(writePathRoot/currPatch.name()/"yCoor"));
                OFstream xCoorFile(
                        fileName(writePathRoot/currPatch.name()/"yCoor"));


				/*
                forAll(myList, i)
                {
                	Info << "label and coordinates :" << myList[i] << " " << currPatch.Cf()[myList[i]] << endl;
                	yCoorFile << currPatch.Cf()[myList[i]].component(vector::X) <<
                        " " << currPatch.Cf()[myList[i]].component(vector::Y) <<
                        " " << currPatch.Cf()[myList[i]].component(vector::Z) << endl;
                }
				*/

//				labelList order;
//				scalarList a(myCells.size());
//				Info << "myCells : " << myCells << endl;
//				Info << "myCells xcoor : " << endl;
//
//				forAll(myCells, i)
//				{
//      				Info<< currPatch.Cf()[myCells[i]].component(vector::X) << endl;
//					//a.append(currPatch.Cf()[myCells[i]].component(vector::X));  // eventually size of a is two times the size of myCells which is not right at all !
//					a[i] = currPatch.Cf()[myCells[i]].component(vector::X);
//				}
//
//				sortedOrder(a, order);
//				sort(a);
//
//				forAll(order, i)
//				{
//					Info<< "order[" << i << "]" << order[i] << "  ordered myCells : " << myCells[order[i]] << " xcoord : " << currPatch.Cf()[myCells[order[i]]].component(vector::X)<< endl;
//                	xCoorFile << currPatch.Cf()[myCells[order[i]]].component(vector::X) <<
//                        " " << currPatch.Cf()[myCells[order[i]]].component(vector::Y) <<
//                        " " << currPatch.Cf()[myCells[order[i]]].component(vector::Z) << endl;
//				}

            	/*
            	// if vect is one random point, findCell will return global label
        		vector vect(0.0004 , -2.96525e-21 , 0.004);
        		Info << "point label : " << mesh.findCell(vect) << endl;
        		Info << "coordinates : " << mesh.C()[mesh.findCell(vect)] << endl;
        		*/
        		//Info << "point label : " << currPatch.findCell(vect) << endl;  // doesn't work
            }

		}
        r.write();
        U1.write();

    }
    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
