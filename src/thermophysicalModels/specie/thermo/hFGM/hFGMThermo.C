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

\*---------------------------------------------------------------------------*/

#include "hFGMThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hFGMThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("hFGMThermo<EquationOfState>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorIn("hFGMThermo<EquationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorIn("hFGMThermo<EquationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hFGMThermo<EquationOfState>::hFGMThermo(Istream& is)
:
    EquationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
    Tcommon_(readScalar(is))
{
    checkInputData();

    forAll(highCpCoeffs_, i)
    {
        is >> highCpCoeffs_[i];
    }

    forAll(lowCpCoeffs_, i)
    {
        is >> lowCpCoeffs_[i];
    }

    // Check state of Istream
    is.check("hFGMThermo::hFGMThermo(Istream& is)");
}


template<class EquationOfState>
Foam::hFGMThermo<EquationOfState>::hFGMThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tlow_(readScalar(dict.subDict("thermodynamics").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("thermodynamics").lookup("Thigh"))),
    Tcommon_(readScalar(dict.subDict("thermodynamics").lookup("Tcommon"))),
    highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
    lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hFGMThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("Tcommon", Tcommon_);
    dict.add("highCpCoeffs", highCpCoeffs_);
    dict.add("lowCpCoeffs", lowCpCoeffs_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hFGMThermo<EquationOfState>& hft
)
{
    os  << static_cast<const EquationOfState&>(hft) << nl
        << "    " << hft.Tlow_
        << tab << hft.Thigh_
        << tab << hft.Tcommon_;

    os << nl << "    ";

    forAll(hft.highCpCoeffs_, i)
    {
        os << hft.highCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(hft.lowCpCoeffs_, i)
    {
        os << hft.lowCpCoeffs_[i] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const hFGMThermo<EquationOfState>& hft)"
    );

    return os;
}


// ************************************************************************* //
