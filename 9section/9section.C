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
    9section

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    //初始化给定分布函数的标量场T

    const volVectorField& C = mesh.C();

    volScalarField T 
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE   
        ),
        mesh,
        dimensionedScalar(dimensionSet(0, 0, 0, 1, 0, 0, 0), 0)
    );

    forAll(C,cellI)
    {
        scalar x = C[cellI].component(0);
        scalar y = C[cellI].component(1);

        scalar phi_xyz = 10*(x*x)*(y*y) + 3*y; 
        T[cellI] = phi_xyz;
    }

    const fvBoundaryMesh& bMesh = mesh.boundary();
    forAll(bMesh, patchI)
    {
        const fvPatch& patch = bMesh[patchI];
        fvPatchField<scalar>& Tpatch = T.boundaryFieldRef()[patchI];

        forAll(patch,faceI)
        {
            vector FC = patch.Cf()[faceI];
            scalar x = FC.component(0);
            scalar y = FC.component(1);

            scalar phi_Bxyz = 10*(x*x)*(y*y) + 3*y;
            Tpatch[faceI] = phi_Bxyz;          
        }
    }
    //


    //计算梯度，在实际算例的fvSchemes中更改格式
    volVectorField Tgrad = fvc::grad(T);
    //

    //建立精确梯度场,并使用interpolate函数插值得到界面上的数值梯度值
    volVectorField realVolGradT 
    (
        IOobject
        (
            "realGradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE   
        ),
        mesh,
        dimensionedVector(dimensionSet(0, -1, 0, 1, 0, 0, 0), vector::zero)
    );

    //const volVectorField& C = mesh.C();
    forAll(C,cellI)
    {
        scalar x = C[cellI].component(0);
        scalar y = C[cellI].component(1);

        scalar gradx = 20*x*(y*y);
        scalar grady = 20*(x*x)*y + 3; 
        
        realVolGradT[cellI] = vector(gradx, grady, 0);
    }

    //const fvBoundaryMesh& bMesh = mesh.boundary();
    forAll(bMesh,patchI)
    {
        const fvPatch& patch = bMesh[patchI];
        fvPatchField<vector>& realVolGradT_BGrad = realVolGradT.boundaryFieldRef()[patchI];

        forAll(patch,faceI)
        {
            vector FC = patch.Cf()[faceI];
            scalar x = FC.component(0);
            scalar y = FC.component(1);

            scalar gradx = 20*x*(y*y);
            scalar grady = 20*(x*x)*y + 3;
            
            realVolGradT_BGrad[faceI] = vector(gradx, grady, 0);          
        }
    }
    //

    surfaceVectorField calcSurGradT
    (
        IOobject
        (
            "calcSurGradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Foam::fvc::interpolate(realVolGradT)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
