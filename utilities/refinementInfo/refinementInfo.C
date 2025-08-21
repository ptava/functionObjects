/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "refinementInfo.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(refinementInfo, 0);
    addToRunTimeSelectionTable(functionObject, refinementInfo, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::refinementInfo::writeFileHeader(Ostream& os)
{
    if (headerDone_)
    {
        writeBreak(file());
    }
    else
    {
        writeHeader(os, "Refinement Info");
    }

    writeCommented(os, "Time");

    for (const word& name : parameters_)
    {
        writeTabbed(os, name);
    }

    os << endl;

    headerDone_ = true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::refinementInfo::refinementInfo(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
  : fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    headerDone_(false),
    parameters_(dict.lookup("refinementParameters"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::refinementInfo::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        writeFile::read(dict);
    }

    return true;
}


bool Foam::functionObjects::refinementInfo::execute()
{
    if (!headerDone_)
    {
        writeFileHeader(file());
    }

    writeCurrentTime(file());

    for (const word& name : parameters_)
    {
        if (obr().foundObject<uniformDimensionedScalarField>(name))
        {
            const auto& parameter =
                obr().lookupObject<uniformDimensionedScalarField>(name);
            const string& value = Foam::name(parameter.value());
            writeTabbed(file(), value);
        }
        else
        {
            writeTabbed(file(), "N/A");
        }
    }

    file() << endl;

    return true;
}

bool Foam::functionObjects::refinementInfo::write()
{
    return true;
}

// ************************************************************************* //
