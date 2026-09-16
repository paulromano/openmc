.. _methods_recoil:

==================
Recoil Production
==================

.. currentmodule:: openmc

When :attr:`Settings.recoil_production` is enabled, supported
continuous-energy neutron reactions create the residual nucleus left by the
reaction as a secondary particle. OpenMC can also create emitted protons,
deuterons, tritons, helium-3 nuclei, and alpha particles. A
:class:`ParticleProductionFilter` can score the energy and direction of
these secondary particles to obtain primary knock-on atom (PKA) and gas-
production spectra.

Recoil production has two complementary parts:

1. OpenMC uses exact two-body kinematics when the processed nuclear data
   describe the emitted neutron or identify a discrete two-body final state.
2. When the processed data omit a required product distribution, OpenMC
   reconstructs a physically admissible event. In particular, a fitted
   **surrogate model** supplies the energy and direction of emitted light
   charged particles. Momentum conservation then determines the heavy
   residual.

The surrogate is necessary because the ACE neutron format used by OpenMC
retains neutron and photon distributions but not the evaluated charged-particle
distributions. It is intended to recover useful PKA and light-ion spectra
without adding all ENDF File 6 charged-particle data to OpenMC's HDF5 format.
Its functional form is motivated by statistical nuclear-reaction theory, while
its parameters are fitted to charged-particle distributions read directly from
ENDF evaluations. The model is described in
:ref:`methods_recoil_light_ions`.

These secondary particles are scored but not transported. OpenMC has no
stopping-power model for heavy ions, so it removes them from the secondary bank
after scoring. Bound thermal scattering represented by
:math:`S(\alpha,\beta)` data or NCrystal does not identify a unique free
recoiling nucleus and is outside the present model. Fission fragments are also
excluded.

--------------------
Kinematic Foundation
--------------------

Conservation of momentum
------------------------

Consider a neutron incident on a stationary target nucleus. Target thermal
motion, which requires a different interpretation, is discussed separately in
:ref:`methods_recoil_elastic`. Let
:math:`\mathbf{p}_{n,\mathrm{in}}` be the incident neutron momentum,
:math:`\mathbf{p}_i` the momentum of emitted particle :math:`i`,
:math:`\mathbf{p}_R` the momentum of the residual nucleus, and
:math:`M_R` the ground-state mass of the residual. Momentum conservation
gives

.. math::
    :label: recoil-momentum

    \mathbf{p}_R =
    \mathbf{p}_{n,\mathrm{in}} - \sum_i \mathbf{p}_i .

For the neutron energies considered here, the massive products are
nonrelativistic. In Newtonian mechanics, momentum is :math:`p=mv` and
kinetic energy is :math:`E=mv^2/2`; eliminating the speed :math:`v` gives
:math:`E=p^2/(2m)`. Applying that relation to the residual gives

.. math::
    :label: recoil-energy

    E_R = \frac{|\mathbf{p}_R|^2}{2M_R}.

These equations do not require a particular system of units. Energies, masses,
momenta, and the speed of light :math:`c` must simply be expressed
consistently. OpenMC uses the relativistic relation
:math:`p_\gamma=E_\gamma/c` for photons.

Conservation of charge and nucleon number determines the residual nuclide. Let
:math:`Z_T` and :math:`A_T` denote the atomic number and mass
number of the target, and let :math:`Z_i` and :math:`A_i` denote
those quantities for emitted particle :math:`i`. The incident neutron
adds one nucleon and no charge, so

.. math::
    :label: recoil-identity

    Z_R = Z_T - \sum_i Z_i ,

.. math::

    A_R = A_T + 1 - \sum_i A_i .

Here :math:`Z_R` and :math:`A_R` identify the residual nucleus.
OpenMC represents that nucleus in its ground state.

Energy release and excitation
-----------------------------

Reconstructing a recoil requires knowing how much energy the reaction products
may share. The incident neutron supplies kinetic energy
:math:`E_{\mathrm{in}}`, while conversion between the initial and final rest
masses can release additional energy or consume some of it. The Q value
measures that rest-mass contribution. A positive Q value adds energy to the
products; a negative Q value uses part of the incident energy and gives the
reaction an energy threshold. OpenMC therefore needs Q both to limit sampled
product energies and to decide whether a reconstructed final state is
physically possible. It enters the event energy balance in
:eq:`recoil-budget`.

The ground-state Q value is the rest-mass energy released when all products
are left in their ground states. If :math:`M_T` is the target nuclear mass,
:math:`m_n` is the neutron mass, :math:`M_R` is the ground-state mass of
the final residual, and :math:`m_j` are the masses of all other final
particles, mass-energy balance gives

.. math::
    :label: recoil-qm

    Q_M =
    \left(M_T + m_n - M_R - \sum_j m_j\right)c^2 .

The tabulated target and residual masses are based on AME2020 [AME2020]_.
OpenMC uses bare nuclear masses consistently in the Q value and the subsequent
kinematic calculations.

ENDF File 3 distinguishes the mass-difference value ``QM`` from the reaction
value ``QI`` [ENDF102]_. For a lumped channel (reactions that represent the sum
over all final states, e.g., MT=103), ``QI`` and ``QM`` generally agree. For a
continuum member of a split representation (e.g., MT=649), ``QI`` is an
effective threshold value and can be several MeV below the ground-state energy
release. OpenMC therefore uses :math:`Q_M` for continuum, lumped, and
multiparticle channels. It uses the evaluated ``QI`` for an MT number that names
one residual level, because that value includes the excitation energy of the
named level.

The total ground-state rest mass of the final particles is

.. math::
    :label: recoil-final-mass

    M_{\mathrm{final}} = M_R + \sum_i m_i ,

where the sum includes every emitted particle with mass. For a stationary
target, the incident neutron and target have total laboratory momentum
:math:`\mathbf{p}_{n,\mathrm{in}}`. The kinetic energy associated with motion of
the final system's center of mass is therefore
:math:`|\mathbf{p}_{n,\mathrm{in}}|^2/(2M_{\mathrm{final}})`. Subtracting that
translational energy from the laboratory energy balance gives the energy
available for relative motion and excitation in the center-of-mass frame:

.. math::
    :label: recoil-u0

    U_0 = E_{\mathrm{in}} + Q
        - \frac{|\mathbf{p}_{n,\mathrm{in}}|^2}
               {2M_{\mathrm{final}}}.

Thus, :math:`U_0` initializes the energy available for recoil reconstruction.
As products are assigned kinetic energies and momenta, OpenMC tracks the
corresponding energy :math:`U` of the system that remains to be divided into
products. If :math:`\mathcal{S}` is the set of products already assigned,
:math:`\mathbf{P}_{\mathrm{rem}}` is the momentum of the remaining system, and
:math:`M_{\mathrm{rem}}` is its total rest mass, the laboratory energy balance
gives

.. math::
    :label: recoil-running-budget

    U = E_{\mathrm{in}} + Q
        - \sum_{i\in\mathcal{S}} E_i
        - \frac{|\mathbf{P}_{\mathrm{rem}}|^2}{2M_{\mathrm{rem}}}.

Before any product is assigned, :math:`\mathcal{S}` is empty,
:math:`\mathbf{P}_{\mathrm{rem}}=\mathbf{p}_{n,\mathrm{in}}`, and
:math:`M_{\mathrm{rem}}=M_{\mathrm{final}}`, so
:eq:`recoil-running-budget` reduces to the definition of :math:`U_0`. Each
subsequent emission reduces :math:`U`, which then determines the kinematic
endpoint for the next emission.

After every emitted particle has been assigned, the remaining system is simply
the residual nucleus: the reaction product left after the neutron, photons, or
light ions have departed. OpenMC does not sample its momentum or kinetic energy
independently. Momentum conservation, :eq:`recoil-momentum`, fixes its momentum
:math:`\mathbf{p}_R`; its kinetic energy then follows from :eq:`recoil-energy`
as :math:`E_R=|\mathbf{p}_R|^2/(2M_R)`. OpenMC creates the residual nucleus as a
secondary particle with that energy and direction.

At this final stage,
:math:`\mathbf{P}_{\mathrm{rem}}=\mathbf{p}_R` and
:math:`M_{\mathrm{rem}}=M_R`. Substituting these values into
:eq:`recoil-running-budget` gives the residual excitation energy

.. math::
    :label: recoil-budget

    E_x = E_{\mathrm{in}} + Q - \sum_i E_i - E_R,

where the sum now includes every emitted particle. Thus, :math:`E_x` is simply
the part of the reaction energy not carried as kinetic energy by an emitted
particle or by the residual nucleus. OpenMC uses :math:`E_x` for energy
bookkeeping but does not create an excited nuclear state or simulate its later
decay. Requiring :math:`E_x\geq0` prevents a reconstructed event from assigning
more kinetic energy than the reaction provides. In an exact two-body channel,
the two product kinetic energies exhaust the available energy and
:math:`E_x=0`, apart from numerical roundoff.

Emission endpoint
-----------------

When OpenMC models an emitted light particle, it must first determine the
largest kinetic energy that the current event can give that particle. This
upper limit is the **kinematic endpoint**. It prevents the surrogate model from
sampling an energy that would leave too little energy for the remaining
products and their recoil.

Suppose a light particle :math:`b` of mass :math:`m_b` is emitted
from a system with current available energy :math:`U` from
:eq:`recoil-running-budget`. Let
:math:`M_D` be the combined mass of the daughter and any products not yet
emitted. In the two-body center-of-mass frame, the two sides have equal and
opposite momentum :math:`p`. Their kinetic energies are

.. math::

    E_b^{\mathrm{cm}} = \frac{p^2}{2m_b},
    \qquad
    E_D^{\mathrm{cm}} = \frac{p^2}{2M_D}
      = \frac{m_b}{M_D}E_b^{\mathrm{cm}} .

Requiring
:math:`E_b^{\mathrm{cm}}+E_D^{\mathrm{cm}}\leq U` gives

.. math::
    :label: recoil-endpoint

    E_{b,\max} = U\frac{M_D}{m_b+M_D}.

For an exact two-body channel, the emitted particle and residual must use all
of :math:`U`, so OpenMC assigns :math:`E_b^{\mathrm{cm}}=E_{b,\max}`. If the
MT number identifies an excited residual level, the energy needed to populate
that level is already included in the level-specific Q value; no additional
unaccounted excitation remains.

For a continuum or multiparticle channel, the surrogate may sample an energy
below the endpoint. After the emission, the internal energy left in the
remaining system is

.. math::

    U_D =
    U - E_b^{\mathrm{cm}}\left(1+\frac{m_b}{M_D}\right) \geq 0 .

This energy is not discarded. If another product still needs to be emitted,
:math:`U_D` becomes the value of :math:`U` used to calculate that product's
endpoint. If no products remain, the remaining system is the residual nucleus
and :math:`U_D` becomes its final bookkeeping excitation :math:`E_x`. Thus,
the endpoint equation both bounds each sample and updates the running energy
budget for a sequence of emissions.

-----------------
Reaction Families
-----------------

.. _methods_recoil_elastic:

Elastic scattering
------------------

For elastic scattering, define the momentum transferred from the neutron to the
target as

.. math::
    :label: elastic-momentum-transfer

    \mathbf{q} =
    \mathbf{p}_{n,\mathrm{in}}-\mathbf{p}_{n,\mathrm{out}},

where :math:`\mathbf{p}_{n,\mathrm{out}}` is the sampled outgoing neutron
momentum. OpenMC assigns the recoil direction
:math:`\widehat{\mathbf{q}}` and the stationary-target-equivalent recoil
energy

.. math::
    :label: elastic-target-rest-energy

    E_R = \frac{|\mathbf{q}|^2}{2M_T}.

For a target initially at rest, this is exactly its final kinetic energy. It is
also the momentum-transfer energy conventionally used in PKA and displacement
models.

OpenMC samples target motion in its free-gas and resonance-scattering models.
If the target initially has laboratory momentum :math:`\mathbf{p}_T`,
momentum conservation gives

.. math::

    \mathbf{p}_{T,\mathrm{out}} = \mathbf{p}_T+\mathbf{q}.

Its final laboratory kinetic energy is consequently

.. math::

    E_T^\mathrm{out} =
    \frac{|\mathbf{p}_T+\mathbf{q}|^2}{2M_T},

and its change in laboratory kinetic energy is

.. math::

    \Delta E_T =
    \frac{|\mathbf{q}|^2}{2M_T}
    + \frac{\mathbf{p}_T\mathbin{\cdot}\mathbf{q}}{M_T}.

Neither quantity generally equals :eq:`elastic-target-rest-energy`.
The energy change can even be negative when a neutron gains energy from target
motion. It therefore cannot be used as the nonnegative energy of a produced
PKA. OpenMC reports :eq:`elastic-target-rest-energy`; target motion still
affects it through the sampled outgoing neutron.

For a stationary target of neutron-mass ratio
:math:`A=M_T/m_n`, the greatest transfer occurs in a head-on collision in
which the neutron reverses direction in the center-of-mass frame. Solving
momentum and kinetic-energy conservation for that limiting case gives

.. math::

    E_{R,\max} = \frac{4A}{(A+1)^2}E_{\mathrm{in}} .

One-neutron inelastic scattering
--------------------------------

For discrete inelastic levels (MT=51-90), the evaluated outgoing-neutron energy
and angle specify a two-body final state. OpenMC obtains the residual momentum
by inserting the sampled neutron in :eq:`recoil-momentum`. This is exact
for the nuclear level and neutron distribution represented by the evaluation.
Momentum from later de-excitation photons is neglected.

Continuum inelastic scattering (MT=91) uses the same momentum balance against
the sampled outgoing neutron. The evaluated distribution is inclusive over
unresolved residual states, so the inferred residual excitation varies from
event to event. OpenMC does not sample an explicit recoil subsection that may
also be stored in the evaluation; it reports the residual implied by momentum
balance against the sampled neutron.

Reactions emitting several neutrons
-----------------------------------

For a reaction such as (n,2n), processed data provide a one-neutron marginal
distribution and an average multiplicity, not a joint distribution for all
neutrons in one event. OpenMC transports one sampled neutron. Its usual
secondary-production treatment duplicates that neutron as needed to reproduce
the expected neutron yield. This preserves the evaluated single-neutron
distribution and expected multiplicity, but the duplicated particles do not
constitute a physically correlated final state.

Recoil reconstruction leaves the transported neutron unchanged and samples
each additional neutron independently from the same evaluated marginal.
Independent samples can request more energy than the reaction provides.
OpenMC therefore accepts an additional neutron only if enough energy remains
for all unmodeled products and for the minimum translational kinetic energy of
the residual system.

Reconstruction fails when the transported neutron is already incompatible with
the energy budget or when no acceptable auxiliary sample is found. In that
case, neutron transport is unchanged and no residual secondary particle is
created. OpenMC never creates a residual whose nuclide identity assumes an
emitted neutron that was omitted from its momentum balance.

The approximate mean recoil energy can be derived directly from
:eq:`recoil-momentum`. Let

.. math::

    \mathbf{V}_{\mathrm{cm}} =
    \frac{\mathbf{p}_{n,\mathrm{in}}}{m_n+M_T}

be the velocity of the entrance-channel center of mass, and let
:math:`\mathbf{k}_i` be emitted-neutron momentum :math:`i` in that
frame. The residual laboratory momentum is

.. math::

    \mathbf{p}_R =
    M_R\mathbf{V}_{\mathrm{cm}}-\sum_i\mathbf{k}_i .

If the emission directions are independent and have zero mean, then
:math:`\langle\mathbf{k}_i\rangle=0` and
:math:`\langle\mathbf{k}_i\mathbin{\cdot}\mathbf{k}_j\rangle=0` for
:math:`i\ne j`. Expanding :math:`|\mathbf{p}_R|^2` in
:eq:`recoil-energy` then gives

.. math::

    \langle E_R\rangle \simeq
    \frac{M_Rm_nE_{\mathrm{in}}}{(m_n+M_T)^2}
    + \frac{m_n}{M_R}
      \sum_i\langle E_i^{\mathrm{cm}}\rangle .

The first term is translation of the center of mass. The second is the average
recoil caused by the independently directed emitted neutrons, where
:math:`E_i^{\mathrm{cm}}` is the center-of-mass kinetic energy of neutron
:math:`i`. Correlations missing from the evaluated marginal would contribute
additional cross terms.

Radiative capture
-----------------

For radiative capture, momentum conservation gives

.. math::

    \mathbf{p}_R =
    \mathbf{p}_{n,\mathrm{in}}-\sum_k\mathbf{p}_{\gamma,k},

where :math:`\mathbf{p}_{\gamma,k}` is the momentum of cascade photon
:math:`k`. Processed reaction data provide inclusive photon distributions
and average yields, not the event-by-event joint distribution of a cascade.
OpenMC converts each yield to a stochastic integer multiplicity and samples
photon energies and directions independently. A common scale factor is then
applied to the sampled photon energies so that

.. math::

    \sum_k E_{\gamma,k}+E_R=E_{\mathrm{in}}+Q .

This procedure preserves the sampled multiplicity and enforces energy and
momentum conservation. It is exact for a fully specified single-photon final
state. A multiphoton recoil spectrum remains model dependent because the
processed data do not supply correlations among cascade photons.

With survival biasing, OpenMC creates the implicit nonfission absorption
secondary with the corresponding absorbed weight and then continues the
neutron to a scattering event. A :class:`ReactionFilter` consequently
reports the scattering MT rather than the implicit absorption MT. Use analog
absorption (``settings.survival_biasing = False``) for reaction-resolved
absorption PKA tallies.

.. _methods_recoil_light_ions:

Light charged particles: fitted surrogate model
------------------------------------------------

Purpose and scope
~~~~~~~~~~~~~~~~~

ACE-derived neutron libraries do not contain the evaluated energy-angle
distributions of emitted protons, deuterons, tritons, helium-3 nuclei, or alpha
particles. OpenMC therefore cannot sample these products directly even though
the original ENDF evaluation may describe them in File 6. The light-ion model
is explicitly a **surrogate** for those missing distributions.

The surrogate has two goals:

- reproduce the broad energy and angular trends in evaluated charged-particle
  distributions with a small, evaluation-independent set of coefficients; and
- produce a complete final state that satisfies charge, nucleon-number, energy,
  and momentum conservation event by event.

The model is not a replacement for an optical-model, Hauser-Feshbach, exciton,
or direct-reaction calculation. It uses compact, physically motivated
functional forms whose fixed coefficients were fitted to evaluated ENDF
charged-particle distributions.

When several products are missing, OpenMC emits them sequentially in the rest
frame of the undecayed system. The transported neutron, if present, is used
first because its evaluated sample must remain unchanged. The missing products
are then placed in random order so that the algorithm does not systematically
give the first product more of the available energy. Each accepted emission
reduces the internal energy according to :eq:`recoil-endpoint`.

Energy distribution
~~~~~~~~~~~~~~~~~~~

The starting point is the Weisskopf-Ewing statistical-emission spectrum
[Weisskopf1940]_,

.. math::
    :label: recoil-weisskopf

    P(E)\mathop{}\!\mathrm{d}E
    \propto
    E\,\sigma_{\mathrm{inv}}(E)\,\rho_D(E_x)
    \mathop{}\!\mathrm{d}E .

Here :math:`E` is the light ion's center-of-mass kinetic energy,
:math:`\sigma_{\mathrm{inv}}(E)` is the cross section for the inverse
reaction in which the ion is absorbed by the daughter, and
:math:`\rho_D(E_x)` is the density of daughter states at remaining
excitation :math:`E_x`. The factor :math:`E` follows from phase
space and detailed balance in the Weisskopf-Ewing derivation.

A transport collision kernel cannot evaluate the inverse cross section and
detailed daughter level density for every event. OpenMC replaces them with a
Coulomb transmission factor and an empirical endpoint factor:

.. math::
    :label: recoil-light-ion

    P(E) \propto
    E\,T_C(E)
    \left(1-\frac{E}{E_{\max}^{\mathrm{shape}}}\right)^\nu ,

with

.. math::

    0\leq E\leq E_{\max}^{\mathrm{event}} .

The two endpoints have different purposes.
:math:`E_{\max}^{\mathrm{shape}}` is the two-body endpoint for the channel
that emits this ion alone and determines the shape of the distribution.
:math:`E_{\max}^{\mathrm{event}}` is calculated from the energy remaining
in the current reconstructed event and limits the sampled energy. The endpoints
are equal for a one-ion channel and differ after another product has consumed
part of the available energy.

For a charged ion, the inverse reaction is strongly suppressed below the
Coulomb barrier. Quantum tunneling through a Coulomb potential gives the Gamow
factor :math:`\exp[-2\pi\eta(E)]` [Gamow1928]_, where the Sommerfeld
parameter is

.. math::
    :label: recoil-sommerfeld

    \eta(E) =
    \alpha Z_bZ_D
    \sqrt{\frac{\mu c^2}{2E}} .

In this expression, :math:`\alpha` is the fine-structure constant,
:math:`Z_b` and :math:`Z_D` are the atomic numbers of the ion and
daughter, and :math:`\mu` is their reduced mass,

.. math::

    \mu = \frac{m_bM_D}{m_b+M_D}.

The approximate height of the Coulomb barrier at the touching radius is

.. math::
    :label: recoil-coulomb-barrier

    V_C =
    \frac{\alpha\hbar c\,Z_bZ_D}
         {r_0\left(A_b^{1/3}+A_D^{1/3}\right)} ,

where :math:`A_b` and :math:`A_D` are the ion and daughter mass
numbers, :math:`\hbar` is the reduced Planck constant, and :math:`r_0`
is an effective nuclear-radius coefficient. This expression is the
electrostatic potential energy of charges :math:`Z_be` and :math:`Z_De`
separated by the assumed touching distance
:math:`r_0(A_b^{1/3}+A_D^{1/3})`.

OpenMC uses a smooth surrogate that has the WKB energy dependence and equals
one half at the nominal barrier:

.. math::
    :label: recoil-barrier

    T_C(E) =
    \left[
      1+\exp\left(
        2\pi g\,[\eta(E)-\eta(V_C)]
      \right)
    \right]^{-1}.

The dimensionless fitted coefficient :math:`g` adjusts the strength of
the idealized Coulomb exponent to account empirically for effects omitted by a
one-dimensional barrier. Equation :eq:`recoil-barrier` is therefore
inspired by WKB tunneling, not asserted to be an exact optical-model
transmission coefficient.

At very low energy, the exponential argument can exceed the floating-point
range even though ratios of probabilities remain meaningful. The
implementation therefore calculates :math:`\ln T_C` directly using a
numerically stable evaluation of :math:`\ln(1+\exp x)`. This avoids
artificially flattening the sub-barrier spectrum. The term sometimes called
the *softplus* function is only a numerical device; it does not add another
physical assumption.

The endpoint power in :eq:`recoil-light-ion` approximates the effect of
the daughter level density. It is empirical rather than a literal nuclear
level-density formula. Evaluated charged-particle spectra combine compound,
pre-equilibrium, and direct emission [Koning2012]_, so applying a pure
compound-nucleus level density to the entire spectrum makes it too soft. The
fitted exponent :math:`\nu` provides a compact compromise among those
components.

Energy-model parameters
~~~~~~~~~~~~~~~~~~~~~~~

The parameter values used by transport are

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   * - Parameter
     - Value
     - Role
   * - :math:`r_0`
     - 1.36093 fm
     - Effective radius in :eq:`recoil-coulomb-barrier`
   * - :math:`g`
     - 0.48566
     - Scale on the WKB-inspired exponent
   * - :math:`\nu`
     - 1.18221
     - Exponent of the empirical endpoint factor

Angular distribution
~~~~~~~~~~~~~~~~~~~~

Continuum light-ion directions use the Kalbach-Mann distribution
[Kalbach1988]_,

.. math::
    :label: recoil-kalbach

    f(\mu) =
    \frac{a}{2\sinh a}
    \left[
      \cosh(a\mu)+r\sinh(a\mu)
    \right].

Here :math:`\mu` is the cosine of the center-of-mass emission angle
relative to the incident neutron, :math:`a` controls the overall angular
slope, and :math:`r` is the pre-equilibrium fraction. The even
:math:`\cosh` term describes forward-backward-symmetric emission, while
the odd :math:`\sinh` term introduces forward bias. ENDF File 6 uses this
same functional form for LANG=2 distributions [ENDF102]_.

OpenMC calculates :math:`a` from Kalbach's published 1988 systematics and
multiplies it by one fitted scale factor. The processed ACE data do not retain
the evaluated :math:`r`, so OpenMC represents it with a bounded logistic
surrogate,

.. math::
    :label: recoil-angular-r

    r = \frac{1}{1+\exp(-u)} .

The logistic transformation ensures :math:`0<r<1`. Its input is

.. math::

    u = c_0+c_1x
        +c_2\ln\left(1+\frac{E_{\mathrm{in}}}{E_0}\right)

.. math::

    \phantom{u={}}
        +c_3A_D^{-1/3}
        +c_4\frac{N_D-Z_D}{A_D},

where

.. math::

    x = \frac{E_b^{\mathrm{cm}}}
             {E_{\max}^{\mathrm{shape}}},
    \qquad E_0=10\ \mathrm{MeV}.

The predictors describe outgoing-energy fraction :math:`x`, incident
energy :math:`E_{\mathrm{in}}`, daughter size :math:`A_D`, and
daughter neutron excess :math:`(N_D-Z_D)/A_D`. Here :math:`N_D`
and :math:`Z_D` are the daughter's neutron and proton numbers.
:math:`E_0` is a fixed scale that makes the logarithm dimensionless. This
logistic equation is empirical; it is not derived from Kalbach's theory.

The parameter values used by transport are

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Parameter
     - Value
     - Predictor
   * - :math:`c_0`
     - -8.42909
     - Intercept
   * - :math:`c_1`
     - 4.89255
     - Outgoing-energy fraction
   * - :math:`c_2`
     - 9.85086
     - Incident-energy dependence
   * - :math:`c_3`
     - -7.85305
     - Daughter-size dependence
   * - :math:`c_4`
     - 5.69894
     - Daughter neutron excess
   * - slope scale
     - 0.99102
     - Multiplier on Kalbach's 1988 slope

The angular surrogate is used only for continuum charged-particle channels.
An MT number in the discrete charged-particle bands (MT=600-849) identifies a
two-body residual level, but its charged-particle angular distribution is also
absent from ACE. OpenMC samples those directions isotropically in the
center-of-mass frame rather than applying a continuum pre-equilibrium model.

Sampling and event completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenMC samples the energy distribution by constructing a fixed-grid cumulative
distribution from :eq:`recoil-light-ion` in logarithmic probability
space and inverting it. The number of function evaluations per emitted ion is
bounded and does not depend on how narrow the spectrum is.

After sampling the energy and angle, OpenMC transforms the ion momentum from
the center-of-mass frame to the laboratory frame. The emission is accepted only
if the remaining system can satisfy :eq:`recoil-budget`. Once all emitted
particles have been constructed, :eq:`recoil-momentum` gives the heavy
residual. If any required product cannot be placed within the energy budget,
OpenMC creates no partial residual secondary particle.

Setting ``light_ion_model`` to ``'none'`` disables the light-ion
surrogate. For a charged-particle reaction, the residual then balances only
the products available from the processed library; emitted light ions are not
created.

Fission
-------

Fission-fragment production is outside the present capability, so a fission
event creates no recoil secondary particles.

-----------
Limitations
-----------

- :math:`S(\alpha,\beta)` and NCrystal thermal-scattering events create no
  recoil secondary particle; only free-gas elastic scattering is supported.
- The residual nucleus is treated as free. Lattice binding and subsequent
  slowing down are not modeled.
- De-excitation photon momentum is neglected except for radiative capture.
- Fission fragments are not produced.
- MT=5 and any other reaction whose exclusive exit channel cannot be inferred
  from the MT number create no recoil secondary particle. These reactions may
  contain appreciable charged-particle production in some evaluations.
- Multi-neutron and multiparticle final states are reconstructed from
  independently sampled marginal distributions rather than a correlated joint
  distribution.
- An event whose required exit channel cannot be completed creates no recoil
  secondary particle. OpenMC does not provide a runtime coverage counter;
  compare reaction-event and residual-production tally weights when measuring
  coverage.
- The light-ion surrogate omits optical-model transmission coefficients,
  explicit level densities, channel competition, direct-reaction amplitudes,
  and evaluated sequential decay.
- Residual excitation, isomeric state, and subsequent gamma recoil are not
  represented as transported states.
- Neglecting atomic electron binding introduces a reaction-dependent difference
  between the approximate and exact nuclear mass balance. This is small relative
  to the MeV-scale energies considered here but is retained as a model
  limitation.
- Massive-particle kinematics are nonrelativistic. The leading neutron
  correction is approximately :math:`E_{\mathrm{in}}/(2m_nc^2)`, below
  one percent at 14 MeV.
- Only continuous-energy transport produces recoil secondary particles.

----------
References
----------

.. [AME2020] M. Wang, W. J. Huang, F. G. Kondev, G. Audi, and S. Naimi,
   "The AME 2020 Atomic Mass Evaluation (II). Tables, Graphs and References,"
   *Chinese Physics C* **45**, 030003 (2021).
   `<https://doi.org/10.1088/1674-1137/abddaf>`_

.. [ENDF102] A. Trkov, M. Herman, and D. A. Brown, eds., *ENDF-6 Formats
   Manual: Data Formats and Procedures for the Evaluated Nuclear Data Files
   ENDF/B-VI, ENDF/B-VII and ENDF/B-VIII*, CSEWG Document ENDF-102,
   BNL-203218-2018-INRE (2018).
   `<https://www.nndc.bnl.gov/endfdocs/ENDF-102-2018.pdf>`_

.. [Gamow1928] G. Gamow, "Zur Quantentheorie des Atomkernes,"
   *Zeitschrift für Physik* **51**, 204-212 (1928).
   `<https://doi.org/10.1007/BF01343196>`_

.. [Weisskopf1940] V. F. Weisskopf and D. H. Ewing, "On the Yield of Nuclear
   Reactions with Heavy Elements," *Physical Review* **57**, 472-485 (1940).
   `<https://doi.org/10.1103/PhysRev.57.472>`_

.. [Kalbach1988] C. Kalbach, "Systematics of Continuum Angular Distributions:
   Extensions to Higher Energies," *Physical Review C* **37**, 2350-2370
   (1988). `<https://doi.org/10.1103/PhysRevC.37.2350>`_

.. [Koning2012] A. J. Koning and D. Rochman, "Modern Nuclear Data Evaluation
   with the TALYS Code System," *Nuclear Data Sheets* **113**, 2841-2934
   (2012). `<https://doi.org/10.1016/j.nds.2012.11.002>`_
