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

#include "hFGMThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState>::hFGMThermo
(
    const EquationOfState& ft,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename hFGMThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename hFGMThermo<EquationOfState>::coeffArray& lowCpCoeffs
)
:
    EquationOfState(ft),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
        lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
    }
}


template<class EquationOfState>
inline const typename Foam::hFGMThermo<EquationOfState>::coeffArray&
Foam::hFGMThermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState>::hFGMThermo
(
    const word& name,
    const hFGMThermo& hft
)
:
    EquationOfState(name, hft),
    Tlow_(hft.Tlow_),
    Thigh_(hft.Thigh_),
    Tcommon_(hft.Tcommon_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = hft.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = hft.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningIn
        (
            "hFGMThermo<EquationOfState>::limit(const scalar T) const"
        )   << "attempt to use hFGMThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class EquationOfState>
inline const typename Foam::hFGMThermo<EquationOfState>::coeffArray&
Foam::hFGMThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::hFGMThermo<EquationOfState>::coeffArray&
Foam::hFGMThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::cp
(
    const scalar p,
    const scalar T
) const
{
    return 2200;
//   const coeffArray& a = coeffs(T);
//   return this->RR*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return this->RR*
    (
        ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
      + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    const scalar Tstd = specie::Tstd;
    return this->RR*
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::hFGMThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    this->RR*
    (
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*::log(T)
      + a[6]
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hFGMThermo<EquationOfState>::operator+=
(
    const hFGMThermo<EquationOfState>& hft
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(hft);

    molr1 /= this->nMoles();
    scalar molr2 = hft.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, hft.Tlow_);
    Thigh_ = min(Thigh_, hft.Thigh_);

    if (hFGMThermo<EquationOfState>::debug && notEqual(Tcommon_, hft.Tcommon_))
    {
        FatalErrorIn
        (
            "hFGMThermo<EquationOfState>::operator+="
            "(const hFGMThermo<EquationOfState>& hft) const"
        )   << "Tcommon " << Tcommon_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << hft.Tcommon_ << " for "
            << (hft.name().size() ? hft.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<hFGMThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs_[coefLabel] =
            molr1*highCpCoeffs_[coefLabel]
          + molr2*hft.highCpCoeffs_[coefLabel];

        lowCpCoeffs_[coefLabel] =
            molr1*lowCpCoeffs_[coefLabel]
          + molr2*hft.lowCpCoeffs_[coefLabel];
    }
}


template<class EquationOfState>
inline void Foam::hFGMThermo<EquationOfState>::operator-=
(
    const hFGMThermo<EquationOfState>& hft
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(hft);

    molr1 /= this->nMoles();
    scalar molr2 = hft.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, hft.Tlow_);
    Thigh_ = min(Thigh_, hft.Thigh_);

    if (hFGMThermo<EquationOfState>::debug && notEqual(Tcommon_, hft.Tcommon_))
    {
        FatalErrorIn
        (
            "hFGMThermo<EquationOfState>::operator-="
            "(const hFGMThermo<EquationOfState>& hft) const"
        )   << "Tcommon " << Tcommon_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << hft.Tcommon_ << " for "
            << (hft.name().size() ? hft.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<hFGMThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs_[coefLabel] =
            molr1*highCpCoeffs_[coefLabel]
          - molr2*hft.highCpCoeffs_[coefLabel];

        lowCpCoeffs_[coefLabel] =
            molr1*lowCpCoeffs_[coefLabel]
          - molr2*hft.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState> Foam::operator+
(
    const hFGMThermo<EquationOfState>& hft1,
    const hFGMThermo<EquationOfState>& hft2
)
{
    EquationOfState eofs = hft1;
    eofs += hft2;

    scalar molr1 = hft1.nMoles()/eofs.nMoles();
    scalar molr2 = hft2.nMoles()/eofs.nMoles();

    typename hFGMThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename hFGMThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        register label coefLabel=0;
        coefLabel<hFGMThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            molr1*hft1.highCpCoeffs_[coefLabel]
          + molr2*hft2.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            molr1*hft1.lowCpCoeffs_[coefLabel]
          + molr2*hft2.lowCpCoeffs_[coefLabel];
    }

    if
    (
        hFGMThermo<EquationOfState>::debug
     && notEqual(hft1.Tcommon_, hft2.Tcommon_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const hFGMThermo<EquationOfState>& hft1,"
            " const hFGMThermo<EquationOfState>& hft2)"
        )   << "Tcommon " << hft1.Tcommon_ << " for "
            << (hft1.name().size() ? hft1.name() : "others")
            << " != " << hft2.Tcommon_ << " for "
            << (hft2.name().size() ? hft2.name() : "others")
            << exit(FatalError);
    }

    return hFGMThermo<EquationOfState>
    (
        eofs,
        max(hft1.Tlow_, hft2.Tlow_),
        min(hft1.Thigh_, hft2.Thigh_),
        hft1.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState> Foam::operator-
(
    const hFGMThermo<EquationOfState>& hft1,
    const hFGMThermo<EquationOfState>& hft2
)
{
    EquationOfState eofs = hft1;
    eofs -= hft2;

    scalar molr1 = hft1.nMoles()/eofs.nMoles();
    scalar molr2 = hft2.nMoles()/eofs.nMoles();

    typename hFGMThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename hFGMThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        register label coefLabel=0;
        coefLabel<hFGMThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            molr1*hft1.highCpCoeffs_[coefLabel]
          - molr2*hft2.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            molr1*hft1.lowCpCoeffs_[coefLabel]
          - molr2*hft2.lowCpCoeffs_[coefLabel];
    }

    if
    (
        hFGMThermo<EquationOfState>::debug
     && notEqual(hft1.Tcommon_, hft2.Tcommon_)
    )
    {
        FatalErrorIn
        (
            "operator-"
            "(const hFGMThermo<EquationOfState>& hft1,"
            " const hFGMThermo<EquationOfState>& hft2)"
        )   << "Tcommon " << hft1.Tcommon_ << " for "
            << (hft1.name().size() ? hft1.name() : "others")
            << " != " << hft2.Tcommon_ << " for "
            << (hft2.name().size() ? hft2.name() : "others")
            << exit(FatalError);
    }

    return hFGMThermo<EquationOfState>
    (
        eofs,
        max(hft1.Tlow_, hft2.Tlow_),
        min(hft1.Thigh_, hft2.Thigh_),
        hft1.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState> Foam::operator*
(
    const scalar f,
    const hFGMThermo<EquationOfState>& hft
)
{
    return hFGMThermo<EquationOfState>
    (
        f*static_cast<const EquationOfState&>(hft),
        hft.Tlow_,
        hft.Thigh_,
        hft.Tcommon_,
        hft.highCpCoeffs_,
        hft.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::hFGMThermo<EquationOfState> Foam::operator==
(
    const hFGMThermo<EquationOfState>& hft1,
    const hFGMThermo<EquationOfState>& hft2
)
{
    return hft2 - hft1;
}


// ************************************************************************* //
