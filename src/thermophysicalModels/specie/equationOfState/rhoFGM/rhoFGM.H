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

Class
    Foam::rhoFGM

Description
    Perfect gas equation of state.

SourceFiles
    rhoFGMI.H
    rhoFGM.C

\*---------------------------------------------------------------------------*/

#ifndef rhoFGM_H
#define rhoFGM_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class rhoFGM;

template<class Specie>
inline rhoFGM<Specie> operator+
(
    const rhoFGM<Specie>&,
    const rhoFGM<Specie>&
);

template<class Specie>
inline rhoFGM<Specie> operator-
(
    const rhoFGM<Specie>&,
    const rhoFGM<Specie>&
);

template<class Specie>
inline rhoFGM<Specie> operator*
(
    const scalar,
    const rhoFGM<Specie>&
);

template<class Specie>
inline rhoFGM<Specie> operator==
(
    const rhoFGM<Specie>&,
    const rhoFGM<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rhoFGM<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class rhoFGM Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rhoFGM
:
    public Specie
{

public:

    // Constructors

        //- Construct from components
        inline rhoFGM(const Specie& sp);

        //- Construct from Istream
        rhoFGM(Istream&);

        //- Construct from dictionary
        rhoFGM(const dictionary& dict);

        //- Construct as named copy
        inline rhoFGM(const word& name, const rhoFGM&);

        //- Construct and return a clone
        inline autoPtr<rhoFGM> clone() const;

        // Selector from Istream
        inline static autoPtr<rhoFGM> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<rhoFGM> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rhoFGM<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const rhoFGM&);
        inline void operator-=(const rhoFGM&);

        inline void operator*=(const scalar);


    // Friend operators

        friend rhoFGM operator+ <Specie>
        (
            const rhoFGM&,
            const rhoFGM&
        );

        friend rhoFGM operator- <Specie>
        (
            const rhoFGM&,
            const rhoFGM&
        );

        friend rhoFGM operator* <Specie>
        (
            const scalar f,
            const rhoFGM&
        );

        friend rhoFGM operator== <Specie>
        (
            const rhoFGM&,
            const rhoFGM&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rhoFGM&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoFGMI.H"

#ifdef NoRepository
#   include "rhoFGM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
