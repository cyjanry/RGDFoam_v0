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

\*---------------------------------------------------------------------------*/

#include "Thermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class thermo>
inline Foam::Thermo<thermo>::Thermo
(
    const thermo& sp
)
:
    thermo(sp)
{}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (Thermo<thermo>::*F)(const scalar, const scalar) const,
    scalar (Thermo<thermo>::*dFdT)(const scalar, const scalar) const,
    scalar (Thermo<thermo>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Tmin = 0;
    scalar Tmax = 0;
    scalar Ttol = T0*tol_;
    int    iter = 0;
    bool newton = true;
    //Info << "T(), T0:" << T0 << ", p:" << p << ", f: " << f << endl;

    //return (this->*F)(p,Test)/(this->*dFdT)(p,Test);

    



    do
    {
        Test = Tnew;
	if(newton == true)
	{
	  scalar A = (this->*F)(p, Test);
	  scalar B = (this->*dFdT)(p, Test);
	  
	  Tnew =
	      (this->*limit)
	      (Test - (A - f)/B);
	  
	  
	  //Info << "Test: " << Test << ", F:" << A << ", dFdT: " << B << endl;
	  //Info << "h(Test) - h0:" << A-f << endl;
	  //Info << "(h(Test) - h0)/cp(Test)" << 0.5 * (A-f)/B << endl;
	  //Info << Test << " --> " << Tnew << endl;
	  
	}
	else
	{
	  Tnew = (Tmin + Tmax) / 2; // new midpoint
	  scalar A = (this->*F)(p, Tnew);
	  //new interval
	  if(sign(A - f) == -1)
	    Tmin = Tnew;
	  else
	    Tmax = Tnew;
	}
        if (iter++ > maxIter_)
        {
	    //Info << "max iter, Tmin: " << Tmin << ", Tmax: " << Tmax << endl;
	    //When two values of the function were found with different signs, switch over to the slower, but more robust secant method
	    if(newton == true)
	    {
	      //Info << "Switching from newton method to bisection, now searching for proper boundaries" << endl;
	      newton = false;
	      //find min / max values
	      Tmin = T0;
	      int i = 0;
	      for(i = 0; i < 1000; i++)
	      {
		Tmin *= 0.9;
		if(sign((this->*F)(p, Tmin) - f) == -1)
		  break;
	      }
	      if(i == 999)
		FatalErrorIn
		(
		    "thermo<Thermo, Type>::T(scalar f, scalar T0, "
		    "scalar (thermo<Thermo, Type>::*F)"
		    "(const scalar) const, "
		    "scalar (thermo<Thermo, Type>::*dFdT)"
		    "(const scalar) const, "
		    "scalar (thermo<Thermo, Type>::*limit)"
		    "(const scalar) const"
		    ") const"
		)   << "Failed to find minimum temperature for bisection method."
		    << abort(FatalError);
		      
	      Tmax = T0;
	      for(i = 0; i < 1000; i++)
	      {
		Tmax *= 1.1;
		if(sign((this->*F)(p, Tmax) - f) == 1)
		  break;
	      }
	      if(i == 999)
		FatalErrorIn
		(
		    "thermo<Thermo, Type>::T(scalar f, scalar T0, "
		    "scalar (thermo<Thermo, Type>::*F)"
		    "(const scalar) const, "
		    "scalar (thermo<Thermo, Type>::*dFdT)"
		    "(const scalar) const, "
		    "scalar (thermo<Thermo, Type>::*limit)"
		    "(const scalar) const"
		    ") const"
		)   << "Failed to find maximum temperature for bisection method."
		    << abort(FatalError);
	    }
	    else // In this case there are no known values from both sides of the root, so it can't be found with the bisection method
	      FatalErrorIn
	      (
		  "thermo<Thermo, Type>::T(scalar f, scalar T0, "
		  "scalar (thermo<Thermo, Type>::*F)"
		  "(const scalar) const, "
		  "scalar (thermo<Thermo, Type>::*dFdT)"
		  "(const scalar) const, "
		  "scalar (thermo<Thermo, Type>::*limit)"
		  "(const scalar) const"
		  ") const"
	      )   << "Maximum number of iterations (" << maxIter_ << ") exceeded"
		  << abort(FatalError);
	     iter = 0; // Start again with bisection method
        }

    } while (mag(Tnew - Test) > Ttol);
    //Info << "New temperature: " << Tnew << endl;
    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::Thermo<thermo>::Thermo
(
    const word& name,
    const Thermo& st
)
:
    thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::cv(const scalar p, const scalar T) const
{
  return this->cp(p, T) - this->cpMcv(p, T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::gamma(const scalar p, const scalar T) const
{
  scalar CP = this->cp(p,T);
  return CP/(CP - this->cpMcv(p, T));
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::e(const scalar p, const scalar T) const
{
  return this->h(p, T) - p*this->W()/this->rho(p,T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::es(const scalar p, const scalar T) const
{
  return this->hs(p, T) - p*this->W()/this->rho(p,T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::g(const scalar p, const scalar T) const
{
  return this->h(p, T) - T*this->s(p, T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::a(const scalar p, const scalar T) const
{
  return this->e(p, T) - T*this->s(p, T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Cp(const scalar p, const scalar T) const
{
  return this->cp(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Cv(const scalar p, const scalar T) const
{
  return this->cv(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::H(const scalar p, const scalar T) const
{
  return this->h(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Hs(const scalar p, const scalar T) const
{
  return this->hs(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Hc() const
{
    return this->hc()/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::S(const scalar p, const scalar T) const
{
  return this->s(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::E(const scalar p, const scalar T) const
{
  return this->e(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::G(const scalar p, const scalar T) const
{
  return this->g(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::A(const scalar p, const scalar T) const
{
  return this->a(p, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::K(const scalar p, const scalar T) const
{
  scalar arg = -this->nMoles()*this->g(p, T)/(this->RR*T);

    if (arg < 600.0)
    {
        return ::exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Kp(const scalar p, const scalar T) const
{
  return K(p, T);
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Kc(const scalar p, const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
      return Kp(p, T);
    }
    else
    {
      return Kp(p, T)*pow(this->Pstd/(this->RR*T), this->nMoles());
    }
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Kx
(
    const scalar p,
    const scalar T
) const
{
    if (equal(this->nMoles(), SMALL))
    {
      return Kp(p, T);
    }
    else
    {
      return Kp(p, T)*pow(this->Pstd/p, this->nMoles());
    }
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
      return Kp(p, T);
    }
    else
    {
      return Kp(p, T)*pow(n*this->Pstd/p, this->nMoles());
    }
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::TH
(
    const scalar h,
    const scalar p,
    const scalar T0
) const
{
    return T
      (
       h,
       p,
       T0,
       &Thermo<thermo>::H,
       &Thermo<thermo>::Cp,
       &Thermo<thermo>::limit
       );
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
  return T
    (
     hs,
     p,
     T0,
     &Thermo<thermo>::Hs,
     &Thermo<thermo>::Cp,
     &Thermo<thermo>::limit
     );
}


template<class thermo>
inline Foam::scalar Foam::Thermo<thermo>::TE
(
    const scalar e,
    const scalar p,
    const scalar T0
) const
{
    return T
      (
       e,
       T0,
       &Thermo<thermo>::E,
       &Thermo<thermo>::Cv,
       &Thermo<thermo>::limit
       );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::Thermo<thermo>::operator+=
(
    const Thermo<thermo>& st
)
{
    thermo::operator+=(st);
}

template<class thermo>
inline void Foam::Thermo<thermo>::operator-=
(
    const Thermo<thermo>& st
)
{
    thermo::operator-=(st);
}

template<class thermo>
inline void Foam::Thermo<thermo>::operator*=(const scalar s)
{
    thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::Thermo<thermo> Foam::operator+
(
    const Thermo<thermo>& st1,
    const Thermo<thermo>& st2
)
{
    return Thermo<thermo>
    (
        static_cast<const thermo&>(st1) + static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::Thermo<thermo> Foam::operator-
(
    const Thermo<thermo>& st1,
    const Thermo<thermo>& st2
)
{
    return Thermo<thermo>
    (
        static_cast<const thermo&>(st1) - static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::Thermo<thermo> Foam::operator*
(
    const scalar s,
    const Thermo<thermo>& st
)
{
    return Thermo<thermo>
    (
        s*static_cast<const thermo&>(st)
    );
}


template<class thermo>
inline Foam::Thermo<thermo> Foam::operator==
(
    const Thermo<thermo>& st1,
    const Thermo<thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
