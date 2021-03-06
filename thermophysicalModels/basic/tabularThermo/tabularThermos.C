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
#include "makeBasicPsiThermo.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "hPsiThermo.H"
#include "pureMixture.H"






//#include "hpsiThermo.H"
//#include "makeThermo.H"

//#include "specie.H"
//#include "perfectGas.H"
#include "tabularEOS.H"
//#include "hConstThermo.H"
//#include "hPolynomialThermo.H"
//#include "hTabularThermo.H"
//#include "eConstThermo.H"
//#include "janafThermo.H"
//#include "sensibleEnthalpy.H"
//#include "sensibleInternalEnergy.H"
//#include "thermo.H"

//#include "constTransport.H"
//#include "sutherlandTransport.H"
//#include "polynomialTransport.H"
#include "tabularTransport.H"

//#include "hePsiThermo.H"
#include "heTabularThermo.H"
#include "hTabularThermo.H"
//#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */
//Some simple combinations for testing purposes
/*
makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);
*/
makeBasicPsiThermo
(
    //hPsiThermo,
    heTabularThermo,
    pureMixture,
    constTransport,
    //sensibleInternalEnergy,
    hConstThermo,//hTabularThermo,
    perfectGas
    //specie
);
//makeBasicPsiThermo
//(
//    hPsiThermo,
//    pureMixture,
//    sutherlandTransport,
//    hConstThermo,
//    perfectGas
//);
/*
makeThermo
(
    //hpsiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    perfectGas,
    specie
);

makeThermo
(
    //hpsiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    tabularEOS,
    specie
);

makeThermo
(
    //hpsiThermo,
    hePsiThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hConstThermo,
    tabularEOS,
    specie
);


//One of the last two combinations is usually used for plasmas, using tabulated data
//This one uses two tables for h(p, T) and T(p, h)
makeThermo
(
    //hpsiThermo,
    heTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    tabularEOS,
    specie
);

//This one only uses one table for h(p,T) and calculates T(p, h) manually by root finding
makeThermo
(
    //hpsiThermo,
    hePsiThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    tabularEOS,
    specie
);

*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
