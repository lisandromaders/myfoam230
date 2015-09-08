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

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::FGMTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo>::FGMTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts)
{}


template<class Thermo>
inline Foam::FGMTransport<Thermo>::FGMTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::FGMTransport<Thermo>::FGMTransport
(
    const word& name,
    const FGMTransport& ft
)
:
    Thermo(name, ft),
    As_(ft.As_),
    Ts_(ft.Ts_)
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

//Pay attention to the fact I am working with the mu viscosity instead of nu viscosity.. Compressible (mu) x Incompressible(nu)

template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar mu_tcc;
    mu_tcc = 4.00e-05;
    return mu_tcc;
    //return As_*::sqrt(T)/(1.0 + Ts_/T);
}

// kappa is the thermal conductivity. Units: W/(mK)

template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}

// alphah is the thermal diffusivity. Units: m^2/s

template<class Thermo>
inline Foam::scalar Foam::FGMTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo>&
Foam::FGMTransport<Thermo>::operator=
(
    const FGMTransport<Thermo>& ft
)
{
    Thermo::operator=(ft);

    As_ = ft.As_;
    Ts_ = ft.Ts_;

    return *this;
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator+=
(
    const FGMTransport<Thermo>& ft
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(ft);

    molr1 /= this->nMoles();
    scalar molr2 = ft.nMoles()/this->nMoles();

    As_ = molr1*As_ + molr2*ft.As_;
    Ts_ = molr1*Ts_ + molr2*ft.Ts_;
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator-=
(
    const FGMTransport<Thermo>& ft
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(ft);

    molr1 /= this->nMoles();
    scalar molr2 = ft.nMoles()/this->nMoles();

    As_ = molr1*As_ - molr2*ft.As_;
    Ts_ = molr1*Ts_ - molr2*ft.Ts_;
}


template<class Thermo>
inline void Foam::FGMTransport<Thermo>::operator*=
(
    const scalar f
)
{
    Thermo::operator*=(f);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator+
(
    const FGMTransport<Thermo>& ft1,
    const FGMTransport<Thermo>& ft2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ft1) + static_cast<const Thermo&>(ft2)
    );

    scalar molr1 = ft1.nMoles()/t.nMoles();
    scalar molr2 = ft2.nMoles()/t.nMoles();

    return FGMTransport<Thermo>
    (
        t,
        molr1*ft1.As_ + molr2*ft2.As_,
        molr1*ft1.Ts_ + molr2*ft2.Ts_
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator-
(
    const FGMTransport<Thermo>& ft1,
    const FGMTransport<Thermo>& ft2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ft1) - static_cast<const Thermo&>(ft2)
    );

    scalar molr1 = ft1.nMoles()/t.nMoles();
    scalar molr2 = ft2.nMoles()/t.nMoles();

    return FGMTransport<Thermo>
    (
        t,
        molr1*ft1.As_ - molr2*ft2.As_,
        molr1*ft1.Ts_ - molr2*ft2.Ts_
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator*
(
    const scalar f,
    const FGMTransport<Thermo>& ft
)
{
    return FGMTransport<Thermo>
    (
        f*static_cast<const Thermo&>(ft),
        ft.As_,
        ft.Ts_
    );
}


template<class Thermo>
inline Foam::FGMTransport<Thermo> Foam::operator==
(
    const FGMTransport<Thermo>& ft1,
    const FGMTransport<Thermo>& ft2
)
{
    return ft2 - ft1;
}


// ************************************************************************* //
