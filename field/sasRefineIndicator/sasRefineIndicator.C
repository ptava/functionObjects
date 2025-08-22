/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
#include "sasRefineIndicator.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "IOobjectList.H"
#include "fvMeshSubset.H"

namespace Foam
{
namespace functionObjects
{

// Enum mapping
const Foam::Enum<sasRefineIndicator::focusRegion> sasRefineIndicator::focusRegionNames_
({
    { sasRefineIndicator::focusRegion::periphery, "periphery" },
    { sasRefineIndicator::focusRegion::core,      "core" },
    { sasRefineIndicator::focusRegion::combined,  "combined" }
});

const Foam::Enum<sasRefineIndicator::functionType> sasRefineIndicator::functionTypeNames_
({
    { sasRefineIndicator::functionType::constant, "constant" },
    { sasRefineIndicator::functionType::gaussian, "gaussian" }
});

defineTypeNameAndDebug(sasRefineIndicator, 0);
addToRunTimeSelectionTable(functionObject, sasRefineIndicator, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sasRefineIndicator::sasRefineIndicator
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    // Read configuration first (sets zoneSubSetPtr_, resultName_, etc.)
    read(dict);

    auto* fldPtr = new volScalarField
    (
        IOobject
        (
            resultName_,
            mesh_.time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, -GREAT)
    );
    mesh_.objectRegistry::store(fldPtr);
}

// * * * * * * * * * * * * * * * *  Helpers  * * * * * * * * * * * * * * * * //

tmp<volScalarField::Internal> sasRefineIndicator::markCoreConstant
(
    const labelList& cellLabels,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar coreWeight
) const
{
    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );

    auto& G = tG.ref();
    for (const label i : cellLabels)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply constant value to all positive values.
        scalar d = c2[i] - c1[i];
        G[i] = d > 0.0 ? coreWeight : -GREAT;
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreOddScaler
(
    const labelList& cellLabels,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar coreWeight,
    const scalar sigma
) const
{
    const scalar invTwoSigma = 0.5/(sigma*sigma);

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );

    auto& G = tG.ref();
    for (const label i : cellLabels)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply an odd, monotonic, sign-preserving function
        scalar d = c2[i] - c1[i];
        G[i] = d * (1 + coreWeight * (1 - exp(invTwoSigma * d * d)));
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markPeripheryGaussSink
(
    const labelList& cellLabels,
    const volScalarField::Internal& Lvk,
    const volScalarField::Internal& c2,
    const scalar peripheryWeight1,
    const scalar peripheryWeight2,
    const scalar sigma
) const
{
    const scalar invTwoSigma = 0.5/(sigma*sigma);

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );

    auto& G = tG.ref();
    for (const label i : cellLabels)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        scalar nLvk = Lvk[i] / c2[i];
        G[i] = peripheryWeight1 * exp(-invTwoSigma * pow(nLvk - 1, 2))
             - peripheryWeight2 * pow(nLvk - 1, 2);
    }

    return tG;
}

// * * * * * * * * * * * * * *  Calculation  * * * * * * * * * * * * * * * * //

void sasRefineIndicator::calcIndicator()
{
    // Mandatory fields on main mesh
    if
    (
        !obr_.foundObject<volScalarField>("Lvk")
     || !obr_.foundObject<volScalarField>("C1")
     || !obr_.foundObject<volScalarField>("C2")
    )
    {
        FatalErrorInFunction
            << "Required fields 'Lvk', 'C1', 'C2' must exist in objectRegistry '"
            << obr_.name() << "'. None are optional." << nl
            << "Available volScalarField objects: "
            << obr_.lookupClass<volScalarField>().toc()
            << exit(FatalError);
    }

    const volScalarField& Lvk = obr_.lookupObject<volScalarField>("Lvk");
    const volScalarField& C1  = obr_.lookupObject<volScalarField>("C1");
    const volScalarField& C2  = obr_.lookupObject<volScalarField>("C2");

    const labelList& cellLabels = mesh_.cellZones()[regionName_];
    auto& fld = mesh_.lookupObjectRef<volScalarField>(resultName_);
    auto& fldI = fld.internalFieldRef();

    const auto& LvkI = Lvk.internalField();
    const auto& C1I  = C1.internalField();
    const auto& C2I  = C2.internalField();

    switch (focusRegion_)
    {
        case focusRegion::core:
            // compute fldI based on the function name selected
            fldI = functionType_ == functionType::constant
                ? markCoreConstant(cellLabels, C1I, C2I, coreWeight_)
                : markCoreOddScaler(cellLabels, C1I, C2I, coreWeight_, sigma_);
            break;
        case focusRegion::periphery:
            fldI = markPeripheryGaussSink
            (
                cellLabels,
                LvkI,
                C2I,
                peripheryWeight1_,
                peripheryWeight2_,
                sigma_
            );
            break;
        case focusRegion::combined:
            break;
    }

    if (debug_)
    {
        Info<< type() << " '" << name() << "': indicator range = ["
            << gMin(fld.internalField()) << ", "
            << gMax(fld.internalField()) << "]" << nl;
    }
}

// * * * * * * * * * * * * * *  Read/Execute/Write  * * * * * * * * * * * * * //

bool sasRefineIndicator::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    resultName_ = dict.getOrDefault<word>("result", "sasRefineIndicator");
    sigma_ = dict.getOrDefault<scalar>("sigma", 0.05);
    coreWeight_ = dict.getOrDefault<scalar>("coreWeight", 10.0);
    peripheryWeight1_ = dict.getOrDefault<scalar>("peripheryWeight1", 10.0);
    peripheryWeight2_ = dict.getOrDefault<scalar>("peripheryWeight2", 1000.0);
    debug_ = dict.getOrDefault<Switch>("debug", false);
    focusRegion_ = focusRegionNames_.get("focusRegion", dict);

    // Testing
    if (focusRegion_ == focusRegion::core)
    {
        functionType_ = functionTypeNames_.getOrDefault(
            "functionType", dict, functionType::constant
        );
    }

    // Optional sub-mesh selection (like momentumError)
    if (dict.found("regionName"))
    {
        regionName_ = dict.get<word>("regionName");
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Required 'cellZone' not found in dictionary." << nl
            << "Do not execute on entire mesh" << nl
            << exit(FatalIOError);
    }

    if (sigma_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "'sigma' must be > 0 (got " << sigma_ << ")." << exit(FatalIOError);
    }

    if (this->log)
    {
        Info<< type() << ' ' << name() << ':' << nl
            << "  focusRegion      : " << focusRegionNames_[focusRegion_] << nl
            << "  sigma            : " << sigma_ << nl
            << "  coreWeight       : " << coreWeight_ << nl
            << "  peripheryWeight1 : " << peripheryWeight1_ << nl
            << "  peripheryWeight2 : " << peripheryWeight2_ << nl
            << "  result           : " << resultName_ << nl
            << "  regionName       : " << regionName_ << nl
            << endl;
    }

    return true;
}


bool sasRefineIndicator::execute()
{
    calcIndicator();
    return true;
}


bool sasRefineIndicator::write()
{
    Log << " functionObjects::" << type() << ' ' << name();

    Log << " writing field: " << resultName_ << nl
        << " for region: " << regionName_ << endl;
    const auto& fld = mesh_.lookupObject<volScalarField>(resultName_);
    fld.write();

    return true;
}

} // End namespace functionObjects
} // End namespace Foam
