/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "rhoFGM.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoFGM<Specie>::rhoFGM(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoFGM<Specie>::rhoFGM
(
    const word& name,
    const rhoFGM<Specie>& rf
)
:
    Specie(name, rf)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoFGM<Specie> >
Foam::rhoFGM<Specie>::clone() const
{
    return autoPtr<rhoFGM<Specie> >(new rhoFGM<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::rhoFGM<Specie> >
Foam::rhoFGM<Specie>::New(Istream& is)
{
    return autoPtr<rhoFGM<Specie> >(new rhoFGM<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::rhoFGM<Specie> > Foam::rhoFGM<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<rhoFGM<Specie> >(new rhoFGM<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoFGM<Specie>::rho(scalar p, scalar T) const
{
    return p/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::rhoFGM<Specie>::psi(scalar, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::rhoFGM<Specie>::Z(scalar, scalar) const
{
    return 1.0;
}


template<class Specie>
inline Foam::scalar Foam::rhoFGM<Specie>::cpMcv(scalar, scalar) const
{
    return this->RR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoFGM<Specie>::operator+=(const rhoFGM<Specie>& rf)
{
    Specie::operator+=(rf);
}


template<class Specie>
inline void Foam::rhoFGM<Specie>::operator-=(const rhoFGM<Specie>& rf)
{
    Specie::operator-=(rf);
}


template<class Specie>
inline void Foam::rhoFGM<Specie>::operator*=(const scalar f)
{
    Specie::operator*=(f);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoFGM<Specie> Foam::operator+
(
    const rhoFGM<Specie>& rf1,
    const rhoFGM<Specie>& rf2
)
{
    return rhoFGM<Specie>
    (
        static_cast<const Specie&>(rf1)
      + static_cast<const Specie&>(rf2)
    );
}


template<class Specie>
inline Foam::rhoFGM<Specie> Foam::operator-
(
    const rhoFGM<Specie>& rf1,
    const rhoFGM<Specie>& rf2
)
{
    return rhoFGM<Specie>
    (
        static_cast<const Specie&>(rf1)
      - static_cast<const Specie&>(rf2)
    );
}


template<class Specie>
inline Foam::rhoFGM<Specie> Foam::operator*
(
    const scalar f,
    const rhoFGM<Specie>& rf
)
{
    return rhoFGM<Specie>(f*static_cast<const Specie&>(rf));
}


template<class Specie>
inline Foam::rhoFGM<Specie> Foam::operator==
(
    const rhoFGM<Specie>& rf1,
    const rhoFGM<Specie>& rf2
)
{
    return rf2 - rf1;
}


// ************************************************************************* //
