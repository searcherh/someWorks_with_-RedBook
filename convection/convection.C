/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    convection

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    //初始化场
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

        volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    const volVectorField& C = mesh.C();
    const fvBoundaryMesh& bMesh = mesh.boundary();
    
    forAll(C, cellI)
    {
        scalar x = C[cellI].component(0);
        scalar y = C[cellI].component(1);

        if(y > 0.5*x)
        {
            T[cellI] = 1;
        }        
    }

    forAll(bMesh, patchI)
    {
        const fvPatch& patch = bMesh[patchI];
        fvPatchField<scalar>& Tpatch = T.boundaryFieldRef()[patchI];

        forAll(patch, faceI)
        {
            const vectorField& Cf = patch.Cf();
            scalar x = Cf[faceI].component(0);
            scalar y = Cf[faceI].component(1);

            if (y > 0.5*x)
            {
                Tpatch[faceI] = 1;
            }
        }
    }

    T.write();
    //初始化场完成

    //建立方程并求解，div的插值格式在fvSchemes中更改
    while(runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"

        fvScalarMatrix TEqn
        (          
          fvm::div(phi, T)
        );

        TEqn.solve();        
        T.correctBoundaryConditions();
        runTime.write();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    
    Info<< "End\n" << endl;
    }
    return 0;
    
}


// ************************************************************************* //
