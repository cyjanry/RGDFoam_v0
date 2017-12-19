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
inline Foam::tabularTransport<Thermo>::tabularTransport
(
    const Thermo& t,
    const extrapolation2DTable<scalar> mu,
    const extrapolation2DTable<scalar> kappa
)
:
    Thermo(t),
    mu_(mu),
    kappa_(kappa)
{
    mu_.outOfBounds(extrapolation2DTable<scalar>::EXTRAPOLATE);
    kappa_.outOfBounds(extrapolation2DTable<scalar>::EXTRAPOLATE);
}


template<class Thermo>
inline Foam::tabularTransport<Thermo>::tabularTransport
(
    const word& name,
    const tabularTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::clone() const
{
    return autoPtr<tabularTransport<Thermo> >
    (
        new tabularTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<tabularTransport<Thermo> >
    (
        new tabularTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabularTransport<Thermo> >
    (
        new tabularTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabularTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_(T, p);
}


template<class Thermo>
inline Foam::scalar Foam::tabularTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_(T, p); // this->Cpv(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
//inline Foam::scalar Foam::tabularTransport<Thermo>::alphah
inline Foam::scalar Foam::tabularTransport<Thermo>::alpha
(
    const scalar p,
    const scalar T
) const
{
    //return kappa(p, T) / this->Cpv(p, T); // mu(p, T)*rPr_;
     return kappa(p, T) / this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabularTransport<Thermo>& Foam::tabularTransport<Thermo>::operator=
(
    const tabularTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    kappa_ = ct.kappa_;

    return *this;
}


template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator+=
(
    const tabularTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);
    /*
    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ + molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    */
}


template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator-=
(
    const tabularTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);
    /*
    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ - molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    */
}


template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator+
(
    const tabularTransport<Thermo>& ct1,
    const tabularTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return tabularTransport<Thermo>
    (
        t,
        ct1.mu,
	ct1.kappa
    );
}


template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator-
(
    const tabularTransport<Thermo>& ct1,
    const tabularTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return tabularTransport<Thermo>
    (
        t,
        ct1.mu,
	ct1.kappa
    );
}


template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator*
(
    const scalar s,
    const tabularTransport<Thermo>& ct
)
{
    return tabularTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.kappa_
    );
}


template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator==
(
    const tabularTransport<Thermo>& ct1,
    const tabularTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
