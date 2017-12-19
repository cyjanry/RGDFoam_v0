/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "tabularEOS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::tabularEOS<Specie>::tabularEOS
(
    const Specie& sp
)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::tabularEOS<Specie>::tabularEOS
(
    const word& name,
    const tabularEOS<Specie>& ico
)
:
    Specie(name, ico)
{

    densityTable = extrapolation2DTable<scalar>("constant/densityTable");
    // extrapolation2DTable<scalar>::ERROR is used to protect the interpolation range
    densityTable.outOfBounds(extrapolation2DTable<scalar>::ERROR);
}


template<class Specie>
inline Foam::autoPtr<Foam::tabularEOS<Specie> >
Foam::tabularEOS<Specie>::clone() const
{
    return autoPtr<tabularEOS<Specie> >(new tabularEOS<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::tabularEOS<Specie> >
Foam::tabularEOS<Specie>::New(Istream& is)
{
    return autoPtr<tabularEOS<Specie> >(new tabularEOS<Specie>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::rho(scalar p, scalar T) const
{

     //Info << "rho at T=" << T << " and p=" << p << " : " << densityTable(T,p) << endl;
     return densityTable(T, p);
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::psi(scalar p, scalar T) const
{
    return rho(p, T) / (p + SMALL);
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p,T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::cpMcv(scalar p, scalar T) const
{
    return cpMcvTable(T,p)*this->W();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::tabularEOS<Specie>::operator+=(const tabularEOS<Specie>& ico)
{
    //scalar molr1 = this->nMoles();

    //Specie::operator+=(ico);
/*
    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ico.rho_;
*/
}


template<class Specie>
inline void Foam::tabularEOS<Specie>::operator-=(const tabularEOS<Specie>& ico)
{
/*
    scalar molr1 = this->nMoles();
*/
    //Specie::operator-=(ico);
/*
    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ico.rho_;
*/
}


template<class Specie>
inline void Foam::tabularEOS<Specie>::operator*=(const scalar s)
{
    //Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::tabularEOS<Specie> Foam::operator+
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{
    return tabularEOS<Specie>
    (
        static_cast<const Specie&>(ico1)
    );
}


template<class Specie>
inline Foam::tabularEOS<Specie> Foam::operator-
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{

    return tabularEOS<Specie>
    (
        static_cast<const Specie&>(ico1)
    );
}


template<class Specie>
inline Foam::tabularEOS<Specie> Foam::operator*
(
    const scalar s,
    const tabularEOS<Specie>& ico
)
{
    return tabularEOS<Specie>(s*static_cast<const Specie&>(ico));
}


template<class Specie>
inline Foam::tabularEOS<Specie> Foam::operator==
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{
    return ico2 - ico1;
}


// ************************************************************************* //
