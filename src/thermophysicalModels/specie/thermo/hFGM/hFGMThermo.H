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
    Foam::hFGMThermo

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    hFGMThermoI.H
    hFGMThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hFGMThermo_H
#define hFGMThermo_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hFGMThermo;

template<class EquationOfState>
inline hFGMThermo<EquationOfState> operator+
(
    const hFGMThermo<EquationOfState>&,
    const hFGMThermo<EquationOfState>&
);

template<class EquationOfState>
inline hFGMThermo<EquationOfState> operator-
(
    const hFGMThermo<EquationOfState>&,
    const hFGMThermo<EquationOfState>&
);

template<class EquationOfState>
inline hFGMThermo<EquationOfState> operator*
(
    const scalar,
    const hFGMThermo<EquationOfState>&
);

template<class EquationOfState>
inline hFGMThermo<EquationOfState> operator==
(
    const hFGMThermo<EquationOfState>&,
    const hFGMThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hFGMThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hFGMThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hFGMThermo
:
    public EquationOfState
{

public:

    // Public data

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline hFGMThermo
        (
            const EquationOfState& ft,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highCpCoeffs,
            const coeffArray& lowCpCoeffs
        );

        //- Construct from Istream
        hFGMThermo(Istream&);

        //- Construct from dictionary
        hFGMThermo(const dictionary& dict);

        //- Construct as a named copy
        inline hFGMThermo(const word&, const hFGMThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hFGM<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the common temperature
            inline scalar Tcommon() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowCpCoeffs() const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hFGMThermo&);
        inline void operator-=(const hFGMThermo&);


    // Friend operators

        friend hFGMThermo operator+ <EquationOfState>
        (
            const hFGMThermo&,
            const hFGMThermo&
        );

        friend hFGMThermo operator- <EquationOfState>
        (
            const hFGMThermo&,
            const hFGMThermo&
        );

        friend hFGMThermo operator* <EquationOfState>
        (
            const scalar,
            const hFGMThermo&
        );

        friend hFGMThermo operator== <EquationOfState>
        (
            const hFGMThermo&,
            const hFGMThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hFGMThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hFGMThermoI.H"

#ifdef NoRepository
#   include "hFGMThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
