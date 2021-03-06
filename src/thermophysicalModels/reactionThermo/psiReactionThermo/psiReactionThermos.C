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

#include "makeReactionThermo.H"

#include "psiReactionThermo.H"
#include "hePsiThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "rhoFGM.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "hFGMThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"
#include "FGMTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// constTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport(and FGMTransport), hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport(or FGMTransport), janafThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hFGMThermo,
    rhoFGM,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hFGMThermo,
    rhoFGM,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    FGMTransport,
    sensibleEnthalpy,
    hFGMThermo,
    rhoFGM,
    specie
);


// Multi-component thermo for sensible enthalpy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    FGMgasHThermoPhysics
);


// Multi-component thermo for internal energy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    FGMgasEThermoPhysics
);


// Multi-component reaction thermo for sensible enthalpy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    FGMgasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    FGMgasHThermoPhysics
);


// Multi-component reaction thermo for internal energy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    FGMgasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    FGMgasEThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
