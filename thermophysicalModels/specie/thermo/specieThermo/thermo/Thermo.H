/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Thermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    ThermoI.H
    Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef Thermo_H
#define Thermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class Thermo;

template<class thermo>
inline Thermo<thermo> operator+
(
    const Thermo<thermo>&,
    const Thermo<thermo>&
);

template<class thermo>
inline Thermo<thermo> operator-
(
    const Thermo<thermo>&,
    const Thermo<thermo>&
);

template<class thermo>
inline Thermo<thermo> operator*
(
    const scalar,
    const Thermo<thermo>&
);

template<class thermo>
inline Thermo<thermo> operator==
(
    const Thermo<thermo>&,
    const Thermo<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const Thermo<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class Thermo
:
    public thermo
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max temperature jump of energy -> temperature inversion functions
        static const scalar TJump_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private member functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(T) and dF(T)/dT
        inline scalar T
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (Thermo::*F)(const scalar, const scalar) const,
            scalar (Thermo::*dFdT)(const scalar, const scalar) const,
            scalar (Thermo::*limit)(const scalar) const
        ) const;


public:

    // Constructors

        //- construct from components
        inline Thermo(const thermo& sp);

        //- Construct from Istream
        Thermo(Istream&);

        //- Construct as named copy
        inline Thermo(const word& name, const Thermo&);


    // Member Functions

        // Fundamental properties
        // (These functions must be provided in derived types)

            // Heat capacity at constant pressure [J/(kmol K)]
            //scalar cp(const scalar) const;

            // Enthalpy [J/kmol]
            //scalar h(const scalar) const;

            // Sensible enthalpy [J/kmol]
            //scalar hs(const scalar) const;

            // Chemical enthalpy [J/kmol]
            //scalar hc(const scalar) const;

            // Entropy [J/(kmol K)]
            //scalar s(const scalar) const;


        // Calculate and return derived properties
        // (These functions need not provided in derived types)

            // Mole specific properties

   

                //- Heat capacity at constant volume [J/(kmol K)]
                inline scalar cv(const scalar, const scalar) const;

                //- gamma = cp/cv []
                inline scalar gamma(const scalar, const scalar) const;

                //- Internal energy [J/kmol]
                inline scalar e(const scalar, const scalar) const;

                //- Sensible internal energy [J/kmol]
                inline scalar es(const scalar, const scalar) const;

                //- Gibbs free energy [J/kmol]
                inline scalar g(const scalar, const scalar) const;

                //- Helmholtz free energy [J/kmol]
                inline scalar a(const scalar, const scalar) const;


            // Mass specific properties

                //- Heat capacity at constant pressure [J/(kg K)]
                inline scalar Cp(const scalar, const scalar) const;

                //- Heat capacity at constant volume [J/(kg K)]
                inline scalar Cv(const scalar, const scalar) const;

                //- Enthalpy [J/kg]
                inline scalar H(const scalar, const scalar) const;

                //- Sensible enthalpy [J/kg]
                inline scalar Hs(const scalar, const scalar) const;

                //- Chemical enthalpy [J/kg]
                inline scalar Hc() const;

                //- Entropy [J/(kg K)]
                inline scalar S(const scalar, const scalar) const;

                //- Internal energy [J/kg]
                inline scalar E(const scalar, const scalar) const;

                //- Gibbs free energy [J/kg]
                inline scalar G(const scalar, const scalar) const;

                //- Helmholtz free energy [J/kg]
                inline scalar A(const scalar, const scalar) const;


        // Equilibrium reaction thermodynamics

            //- Equilibrium constant [] i.t.o fugacities
            //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar, const scalar) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar, const scalar) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar, const scalar) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx(const scalar T, const scalar p) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar T,
                const scalar n
            ) const;


        // Energy->temperature  inversion functions

            //- Temperature from Enthalpy given an initial temperature T0
            inline scalar TH(const scalar, const scalar, const scalar) const;

            //- Temperature from sensible Enthalpy given an initial T0
            inline scalar THs(const scalar, const scalar, const scalar) const;

            //- Temperature from internal energy given an initial temperature T0
            inline scalar TE(const scalar, const scalar, const scalar) const;


    // Member operators

        inline void operator+=(const Thermo&);
        inline void operator-=(const Thermo&);

        inline void operator*=(const scalar);


    // Friend operators

        friend Thermo operator+ <thermo>
        (
            const Thermo&,
            const Thermo&
        );

        friend Thermo operator- <thermo>
        (
            const Thermo&,
            const Thermo&
        );

        friend Thermo operator* <thermo>
        (
            const scalar s,
            const Thermo&
        );

        friend Thermo operator== <thermo>
        (
            const Thermo&,
            const Thermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const Thermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoI.H"

#ifdef NoRepository
#   include "Thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
