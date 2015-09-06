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

Class
    Foam::FGMTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    FGMTransportI.H
    FGMTransport.C

\*---------------------------------------------------------------------------*/

#ifndef FGMTransport_H
#define FGMTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class FGMTransport;

template<class Thermo>
inline FGMTransport<Thermo> operator+
(
    const FGMTransport<Thermo>&,
    const FGMTransport<Thermo>&
);

template<class Thermo>
inline FGMTransport<Thermo> operator-
(
    const FGMTransport<Thermo>&,
    const FGMTransport<Thermo>&
);

template<class Thermo>
inline FGMTransport<Thermo> operator*
(
    const scalar,
    const FGMTransport<Thermo>&
);

template<class Thermo>
inline FGMTransport<Thermo> operator==
(
    const FGMTransport<Thermo>&,
    const FGMTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const FGMTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class FGMTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class FGMTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline FGMTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline FGMTransport(const word&, const FGMTransport&);

        //- Construct from Istream
        FGMTransport(Istream&);

        //- Construct from dictionary
        FGMTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<FGMTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<FGMTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<FGMTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "FGM<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline FGMTransport& operator=(const FGMTransport&);

        inline void operator+=(const FGMTransport&);

        inline void operator-=(const FGMTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend FGMTransport operator+ <Thermo>
        (
            const FGMTransport&,
            const FGMTransport&
        );

        friend FGMTransport operator- <Thermo>
        (
            const FGMTransport&,
            const FGMTransport&
        );

        friend FGMTransport operator* <Thermo>
        (
            const scalar,
            const FGMTransport&
        );

        friend FGMTransport operator== <Thermo>
        (
            const FGMTransport&,
            const FGMTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const FGMTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FGMTransportI.H"

#ifdef NoRepository
#   include "FGMTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
