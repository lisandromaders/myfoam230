/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo>::FGMTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::FGMTransport<Thermo>::FGMTransport
(
    const word& name,
    const FGMTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::FGMTransport<Thermo> >
Foam::FGMTransport<Thermo>::clone() const
{
    return autoPtr<FGMTransport<Thermo> >
    (
        new FGMTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::FGMTransport<Thermo> >
Foam::FGMTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<FGMTransport<Thermo> >
    (
        new FGMTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::FGMTransport<Thermo> >
Foam::FGMTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<FGMTransport<Thermo> >
    (
        new FGMTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cpv(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo>& Foam::FGMTransport<Thermo>::operator=
(
    const FGMTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;

    return *this;
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator+=
(
    const FGMTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ + molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator-=
(
    const FGMTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ - molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator+
(
    const FGMTransport<Thermo>& ct1,
    const FGMTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return FGMTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ + molr2*ct2.mu_,
        1.0/(molr1/ct1.rPr_ + molr2/ct2.rPr_)
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator-
(
    const FGMTransport<Thermo>& ct1,
    const FGMTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return FGMTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ - molr2*ct2.mu_,
        1.0/(molr1/ct1.rPr_ - molr2/ct2.rPr_)
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator*
(
    const scalar s,
    const FGMTransport<Thermo>& ct
)
{
    return FGMTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator==
(
    const FGMTransport<Thermo>& ct1,
    const FGMTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
