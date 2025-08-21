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

Class
    Foam::functionObjects::sasRefineIndicator

Description
    This function object calculates a field for mesh refinement based for
    kOmegaSSTSAS turbulence model, looking at its formulation of the von Karman
    length scale.

    Lvk = max(c1, c2)
    where c2 stands for the `high wave number damper`, i.e. grid filter

    Procecess Lvk in order to obtain the refinement field that leads to refine
    where the grid filter is active

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    sasRefineIndicator.C

\*---------------------------------------------------------------------------*/
#ifndef functionObjects_sasRefineIndicator_H
#define functionObjects_sasRefineIndicator_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "Switch.H"
#include "Enum.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class sasRefineIndicator
\*---------------------------------------------------------------------------*/

class sasRefineIndicator
:
    public fvMeshFunctionObject
{
public:

    // Focus region selection
    enum focusRegion
    {
        core,
        periphery,
        combined
    };
    static const Enum<focusRegion> focusRegionNames_;

    enum functionType
    {
        constant,
        gaussian
    };
    static const Enum<functionType> functionTypeNames_;

protected:

    // ---- User options
    focusRegion  focusRegion_;
    scalar       sigma_;
    scalar       coreWeight_;
    scalar       peripheryWeight_;
    functionType functionType_;
    word         resultName_;
    Switch       debug_;

    // ---- Region definition 
    // TODO: address 'momentumError' procedure beacause isn't working
    word regionName_;

    // ---- Created field names (scoped)
    word meshFieldName_;   // field held on the main mesh
    word subFieldName_;    // field on the sub-mesh (if any)

    // ---- Helpers

    //- Transform field via constant value (applyed only on some cells)
    tmp<volScalarField::Internal> markCoreConstant
    (
        const labelList& cellLabels,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar coreWeight
    ) const;

    //- Transform field via gaussian function (applyed on whole region)
    tmp<volScalarField::Internal> markCoreGaussian
    (
        const labelList& cellLabels,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar coreWeight,
        const scalar sigma
    ) const;

    //- Compute the indicator (writes into created fields)
    void calcIndicator();

public:

    TypeName("sasRefineIndicator");

    // Constructors

    sasRefineIndicator
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Read dictionary, with deprecations handled
    virtual bool read(const dictionary&);

    //- Execute: compute
    virtual bool execute();

    //- Write
    virtual bool write();

    //- No-op end
    virtual bool end() { return true; }
};

} // End namespace functionObjects
} // End namespace Foam

#endif
