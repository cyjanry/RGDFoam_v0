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

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

\*---------------------------------------------------------------------------*/

#include "Thermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

template<class thermo>
const Foam::scalar Foam::Thermo<thermo>::tol_
(
    debug::tolerances("speciesThermoTol", 1e-4)
);


template<class thermo>
const Foam::scalar Foam::Thermo<thermo>::TJump_
(
    debug::tolerances("speciesThermoTJump", 20)
);


template<class thermo>
const int Foam::Thermo<thermo>::maxIter_
(
    debug::optimisationSwitch("speciesThermoMaxIter", 100)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
Foam::Thermo<thermo>::Thermo(Istream& is)
:
    thermo(is)
{
    is.check("Thermo::Thermo(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const Thermo<thermo>& st)
{
    os  << static_cast<const thermo&>(st);

    os.check("Ostream& operator<<(Ostream& os, const Thermo& st)");
    return os;
}


// ************************************************************************* //
