.. _methods_recoil:

==================
Recoil Production
==================

.. currentmodule:: openmc

When :attr:`Settings.recoil_production` is enabled, every continuous-energy
neutron collision creates additional entries in the secondary bank describing
the *recoil nucleus* left by the reaction and, optionally, the light
ions the reaction emits. Scoring those entries with a
:class:`ParticleProductionFilter` gives the primary knock-on atom (PKA)
spectrum, which is the starting point for displacement-damage estimates and for
hydrogen and helium production spectra.

These entries are production records, not transported particles. OpenMC has no
stopping-power model for heavy ions, so the records are removed from the
secondary bank as soon as the collision tallies have been scored.

--------------------
Kinematic Foundation
--------------------

Every reaction family uses the same nonrelativistic momentum balance. Writing
:math:`\mathbf{p}_{n,\text{in}}` for the incident neutron momentum and
:math:`\mathbf{p}_i` for the momentum of each emitted particle,

.. math::
    :label: recoil-momentum

    \mathbf{p}_R = \mathbf{p}_{n,\text{in}} - \sum_i \mathbf{p}_i , \qquad
    E_R = \frac{|\mathbf{p}_R|^2}{2 M_R} .

Masses and energies are in eV and momenta in eV with :math:`c = 1`, so
:math:`p = \sqrt{2mE}` for a massive particle and :math:`p = E` for a photon.

The identity of the recoil follows from the MT number,

.. math::

    Z_R = Z_T - \sum_i Z_i , \qquad A_R = A_T + 1 - \sum_i A_i ,

and the recoil is always represented in its ground state. Recoil excitation
energy is implicit: it is whatever remains of

.. math::
    :label: recoil-budget

    E_x = E_\text{in} + Q - \sum_i E_i - E_R \ge 0

after all kinetic energies are accounted for.

Products that OpenMC samples from evaluated data are used exactly as sampled.
Products for which the library carries no distribution are modelled, and the
model is constrained by :eq:`recoil-budget` so that no event can violate energy
conservation.

Masses
------

Masses are **nuclear**, obtained from the tabulated atomic mass as
:math:`M_\text{nuc}(Z,A) = M_\text{atom}(Z,A) - Z m_e`, with the proton,
deuteron, triton, helion and alpha taken from their CODATA nuclear values.
Electron *binding* energy is neglected: it does not cancel exactly in a
charged-particle Q value, but the residue is 13.6 eV for a proton channel,
79 eV for an alpha, and a few keV on a heavy target.

Mixing the two conventions is the error worth guarding against. A Q value built
from atomic targets and nuclear light ions is wrong by :math:`Z_b m_e c^2`,
which is 1.02 MeV for an alpha channel. A nuclide with no tabulated mass yields
no budget and therefore no recoil record, rather than a fabricated mass and a Q
value wrong by tens of MeV.

Energy budget
-------------

The Q value in :eq:`recoil-budget` is the **rest-mass energy release** of the
exit channel,

.. math::
    :label: recoil-qm

    Q_M = \left[M_T + m_n - M_D - \sum_j m_j\right] c^2 ,

computed from mass excesses so that the mass numbers cancel identically rather
than numerically. It is *not*, in general, the Q value stored in the nuclear
data library. ENDF MF=3 carries two: ``QM``, the mass-difference Q, and ``QI``,
the Q of the lowest state the MT represents — or, when the MT names no unique
state, an effective value chosen to place the threshold correctly. ENDF-102
warns that such a value cannot be relied on for energy-release calculations,
and it is ``QI`` that both the ACE and the direct-ENDF readers store in
:attr:`Reaction.q_value`.

For a lumped channel the two agree. For the continuum member of a split
representation — MT = 649, 699, 749, 799, 849 — ``QI`` lies below the release by
a median 3.5 keV and by as much as 7.0 MeV, measured over 190 such channels in
ENDF/B-VIII.1, JEFF-4.0, JENDL-5 and TENDL-2025. Using it as the budget
truncates the modelled light-ion spectrum well inside the evaluated one.

OpenMC therefore uses :eq:`recoil-qm` for continuum, lumped and multiparticle
channels, and the evaluated ``QI`` only for MT numbers that name one residual
level (MT = 51-90 and the discrete charged-particle bands), where it is the
level-specific Q and is exactly what a two-body channel needs. A level Q above
the ground-state release implies a negative excitation; an excess within
0.25 MeV is absorbed as mass-table disagreement, and a larger one means the
evaluation cannot be reconciled and the event produces no recoil.

Entrance and exit energetics
----------------------------

For a stationary target the energy available in the compound system's rest
frame is

.. math::
    :label: recoil-u0

    U_0 = E_\text{in} + Q - \frac{|\mathbf{p}_n|^2}{2 (M_T + m_n)}
        = E_\text{in} \frac{M_T}{M_T + m_n} + Q ,

the subtracted term being the kinetic energy of the centre of mass, which no
exit channel can spend. Emitting ion :math:`b` from a parent with remaining
internal energy :math:`U` leaves the daughter recoiling against it, so

.. math::
    :label: recoil-endpoint

    E_{b,\max} = U \frac{M_D}{m_b + M_D} , \qquad
    E_x = U - E_b^\text{cm}\left(1 + \frac{m_b}{M_D}\right) \ge 0 ,

with :math:`M_D` the daughter plus any product not yet emitted. These
expressions are shared with the offline calibration through a checked-in
fixture of 374 cases, so the distribution the transport kernel samples is the
one the surrogate was fitted to.

-----------------
Reaction Families
-----------------

Elastic scattering
------------------

Subtracting the sampled outgoing neutron momentum from the incident one gives
the recoil exactly.

The target's own thermal momentum is deliberately left out of
:eq:`recoil-momentum`. A PKA energy is the energy a collision *imparts* to an
atom, which is what displacement models and NJOY recoil matrices are built on;
adding the target's pre-collision momentum would instead report the atom's
total kinetic energy, which below a few eV is dominated by thermal motion and
has nothing to do with damage. The free-gas treatment still shapes the result
through the sampled outgoing neutron, so up-scattering broadens the recoil
distribution slightly past the stationary-target endpoint
:math:`4 A E / (A+1)^2` while leaving the mean at
:math:`2 A E (1 - \bar\mu_\text{cm}) / (A+1)^2`.

Discrete inelastic scattering
-----------------------------

For MT = 51-90 the outgoing neutron determines the recoil completely, so
subtracting the sampled neutron momentum is exact for the two-body level
transition described by the evaluated angular distribution. Momentum carried by
the de-excitation photons is neglected; for a recoil of mass :math:`M_R` and
uncorrelated photon directions they add
:math:`\overline{\sum_k E_{\gamma,k}^2} / 2 M_R c^2` to the mean recoil energy,
which for the Fe-56 continuum inelastic channel at 14 MeV — 3.7 photons
averaging 1.3 MeV — is of order 100 eV.

Continuum inelastic scattering
------------------------------

MT = 91 is treated the same way: the recoil recoils against the sampled
outgoing neutron. This is exact given the evaluated neutron energy-angle
distribution. Note that it can differ substantially from the explicit recoil
array some evaluations store in MF=6; see :ref:`methods_recoil_validation`.

Reactions emitting several neutrons
-----------------------------------

OpenMC's transport samples one outgoing neutron and duplicates it for integral
yields, which is correct for transport but not for an event-by-event momentum
balance: duplicating one momentum vector makes the emitted momenta perfectly
correlated, which broadens the recoil spectrum and shifts its mean.

Instead, each additional neutron of a multiplicity :math:`\nu > 1` channel is
sampled independently from the same evaluated distribution. Because the ENDF
distribution is inclusive and carries no joint final state, independent samples
can overrun the event's energy budget; a sample that would do so is rejected and
redrawn. That keeps every event kinematically possible while leaving the
marginal spectrum close to the evaluated one.

Rejection can fail. When it does the event produces **no record at all**,
rather than omitting the neutron and banking a recoil whose label claims it
left: the residual's mass number would then be one higher than the momentum
balance it carries. The same rule applies to a light ion that cannot be
emitted. Every such event is counted, so a channel that fails often is visible
rather than merely underrepresented.

With uncorrelated emission directions the cross terms in :eq:`recoil-momentum`
average out and the mean recoil approaches

.. math::

    \langle E_R \rangle \simeq
      \frac{M_R m_n E_\text{in}}{(m_n + M_T)^2}
      + \frac{m_n}{M_R} \sum_i \langle E_i^\text{cm} \rangle ,

the first term being the motion of the centre of mass and the second the recoil
against the emitted neutrons.

Radiative capture
-----------------

The recoil recoils against the emitted photons,
:math:`\mathbf{p}_R = \mathbf{p}_{n,\text{in}} - \sum_k \mathbf{p}_{\gamma,k}`.
Photon multiplicity and energies are resampled from the capture reaction's own
photon distribution, and the yield is converted to a stochastic integer
multiplicity so that the recoil belongs to the reaction a
:class:`ReactionFilter` reports.

This is an event-by-event model, so it produces a spectrum rather than the
average kick

.. math::

    \overline{E_R} = \frac{E}{A+1}
      + \frac{\overline{\sum_k E_{\gamma,k}^2}}{2 (A+1) m_n c^2}

that NJOY's HEATR module [MacFarlane2016]_ uses when explicit recoil data are
absent. The two agree in the mean but not in shape.

.. _methods_recoil_light_ions:

Light charged particles
-----------------------

ACE-derived libraries carry only neutron and photon products, so the proton,
deuteron, triton, helium-3, and alpha ions of channels such as (n,p),
(n,\ :math:`\alpha`), and (n,np) have no evaluated distribution to sample. When
``light_ion_model`` is ``'statistical'``, OpenMC emits them sequentially in the
rest frame of the system that has not yet decayed.

This is a **conditional independent-emission construction**, not a model of the
event's decay chain. Evaluated files supply marginal distributions and no joint
final state, so no ordering of the products is derivable from the data; the
transported neutron is taken first because it is the one product the library
does describe, and the modelled products follow in a randomized order so that
none is systematically favoured. Randomizing does not restore the correlations
the marginals omit.

Each emission follows :eq:`recoil-endpoint`: the ion can carry at most
:math:`E_{b,\max}`, and it removes :math:`E_b^\text{cm}(1 + m_b/M_D)` from the
parent's internal energy, leaving the daughter with
:math:`E_x = U\left(1 - E_b^\text{cm}/E_{b,\max}\right)`. Emission order is
randomized so no ion is systematically favoured.

Energy
~~~~~~

The centre-of-mass energy follows the Weisskopf-Ewing form for statistical
particle emission [Weisskopf1940]_,

.. math::
    :label: recoil-weisskopf

    P(E) \propto E\, \sigma_\text{inv}(E)\, \rho_D(E_x) ,

where the leading :math:`E` comes from detailed balance and phase space,
:math:`\sigma_\text{inv}` is the cross section of the inverse reaction
:math:`b + D \rightarrow` parent, and :math:`\rho_D` is the level density of the
daughter. OpenMC evaluates it as

.. math::
    :label: recoil-light-ion

    P(E) \propto E\, T_C(E)
        \left(1 - \frac{E}{E_{b,\max}^\text{shape}}\right)^{\nu} .

*Inverse cross section.* For a charged light ion :math:`\sigma_\text{inv}` is
governed almost entirely by the Coulomb barrier, so it is replaced by a WKB
Coulomb penetrability,

.. math::
    :label: recoil-barrier

    T_C(E) = \left[1 + e^{2\pi g\,[\eta(E) - \eta(V_C)]}\right]^{-1} ,
    \qquad
    \eta(E) = \frac{Z_b Z_D}{137.036}\sqrt{\frac{\mu c^2}{2E}} ,
    \qquad
    V_C = \frac{1.44\ \text{MeV fm}\ Z_b Z_D}{r_0 (A_b^{1/3} + A_D^{1/3})} ,

with :math:`\mu` the reduced mass of the ion and the daughter. It is normalized
so that :math:`T_C = 1/2` at the barrier top.

The decisive feature is the :math:`E^{-1/2}` inside the exponent. Because
:math:`\eta` grows as :math:`E` falls, so does the decay constant of the
transmission — which is the widening of the barrier at lower energy. A
Hill-Wheeler transmission [HillWheeler1953]_ with a fixed diffuseness falls at
one rate everywhere and cannot reproduce it; the Gamow form replaced it for
that reason.

The expression is evaluated in log space, as
:math:`\ln T_C = -\operatorname{softplus}(2\pi g[\eta - \eta_B])`, with no bound
on the exponent. Written directly, :math:`[1+e^x]^{-1}` underflows to exactly
zero once :math:`x` passes 745 — only tens of keV below the barrier for an
alpha against a heavy target — and a spectrum that is identically zero cannot
be normalized.

*Endpoint factor.* The :math:`(1-E/E_{b,\max}^\text{shape})^{\nu}` factor stands
in for the level density of the daughter. It is **not** a nuclear level density,
which would rise roughly exponentially with :math:`E_x`. Charged-particle
spectra in TALYS-based evaluations [Koning2012]_ mix compound-nucleus
evaporation with much harder pre-equilibrium and direct emission, and applying a
compound level density to the whole spectrum makes it far too soft. The
exponent is an empirical compromise between the two components, and
:math:`\nu + 2` should not be read as an exciton number unless a model with
that derivation is actually being fitted.

*Constants.* Three numbers are fitted: the effective barrier radius
:math:`r_0`, the WKB strength :math:`g`, and the exponent :math:`\nu`. Their
current values are the defaults of ``LightIonParams`` in
``include/openmc/recoil.h``, which is the single authoritative record; they are
deliberately not restated here, because a documented copy of a fitted constant
is a copy that goes stale. They were obtained by matching the evaluated ENDF
MF=6 centre-of-mass spectra of the charged-particle channels across several
general-purpose libraries.

Equation :eq:`recoil-light-ion` should be read as a calibrated surrogate whose
*form* is borrowed from statistical-model theory and whose *parameters* come
from evaluated data — not as an implementation of any published
nuclear-reaction model. The parameters are also strongly correlated: :math:`r_0`
and :math:`\nu` trade against each other along a shallow valley, so neither
should be interpreted on its own as a measured physical quantity.

*Sampling.* The spectrum is sampled by inverting a tabulated cumulative built in
log space on a fixed grid. The work per emitted ion is therefore bounded and
independent of how peaked the spectrum is, and no draw is ever replaced by a
fallback value.

Endpoint and direction
~~~~~~~~~~~~~~~~~~~~~~

Two endpoints appear in :eq:`recoil-light-ion` and they play different roles.
The *shape* endpoint is the one belonging to the channel that emits this ion
**alone**, built from :eq:`recoil-qm` and :eq:`recoil-u0` for that channel and
not from any evaluated Q; it sets the spectrum's shape, because the endpoint
factor stands in for the level density of the daughter that ion would leave on
its own. The *event* endpoint is :eq:`recoil-endpoint` evaluated on whatever
internal energy this particular event has left, and it truncates the sample.
The two coincide unless an earlier product has already spent part of the
budget. Discrete charged-particle levels (MT = 600-849) are exactly two-body,
so their light-ion energy is fixed at the event endpoint rather than sampled.

The direction uses the Kalbach-Mann form of evaluated MF=6 LANG=2 data,

.. math::

    f(\mu) = \frac{a}{2 \sinh a}
             \left[\cosh(a\mu) + r \sinh(a\mu)\right] ,

with the slope :math:`a` from Kalbach's systematics [Kalbach1988]_ — the same
expression the evaluations themselves use, taken as published — and the
pre-equilibrium fraction from a logistic in the outgoing energy fraction, the
incident energy and the size and neutron excess of the recoil nucleus,

.. math::

    r = \left[1 + e^{-u}\right]^{-1} , \qquad
    u = c_0 + c_1 \frac{E_b^\text{cm}}{E_b^\text{max}}
        + c_2 \ln\!\left(1 + \frac{E}{10\ \text{MeV}}\right)
        + c_3 A_D^{-1/3} + c_4 \frac{N_D - Z_D}{A_D} ,

calibrated against evaluated MF=6 LANG=2 distributions. Only :math:`r` is a
substitution.

This applies to continuum channels only. Kalbach's systematics describe a
channel fed by pre-equilibrium emission, and a named level is not one:
measured against 25,000 evaluated discrete distributions, carrying the
systematics onto MT = 600-849 doubles the recoil Wasserstein error relative to
assuming nothing and biases the first Legendre moment by :math:`+0.12`.
**Discrete charged-particle levels are therefore sampled isotropically in the
centre of mass**, which leaves that bias at :math:`-0.04`. No fitted angular
correction on top of isotropy survived being held out across libraries, so
none is applied; the residual error there is smaller than the disagreement
between libraries evaluating the same channel.

Setting ``light_ion_model`` to ``'none'`` skips this model entirely, in which
case the recoil of a charged-particle channel recoils against the incident
neutron alone.

Fission
-------

Fission is excluded. Fission-fragment recoil is not modelled.

.. _methods_recoil_validation:

--------------------------------
Comparison With NJOY Group Data
--------------------------------

The usual reference for PKA spectra is a group-wise recoil matrix produced by
NJOY, as consumed by codes such as SPECTRA-PKA [Gilbert2015]_. OpenMC agrees
closely with those matrices where both sides derive the recoil from the same
two-body kinematics, and differs in three understood places.

**Two-body channels agree.** For elastic scattering and for the discrete
inelastic levels MT = 51-90, NJOY computes the recoil from the evaluated MF=4
angular distribution and OpenMC subtracts the sampled neutron momentum. Mean
recoil energies agree to within a few percent over the whole range of
structural nuclides and incident energies, and the spectra agree in shape over
six decades, including the diffraction structure near the endpoint.

**Nuclear-data representation causes a small elastic offset.** OpenMC's
ACE-derived library and NJOY's ENDF MF=4 Legendre coefficients do not describe
the same forward-peaked angular distribution to better than a few percent. For
Fe-56 at 13.9 MeV the mean cosine is 0.8535 from the ACE-derived tabulated
distribution and 0.8462 as implied by the NJOY recoil matrix built from the same
evaluation, a 4% difference in :math:`1 - \bar\mu` and therefore in the mean
elastic recoil energy. This is a data-processing difference upstream of the
transport code, not a difference between the recoil models.

**Some evaluations store recoil arrays that violate momentum conservation.**
Evaluations generated with TALYS write an explicit recoil subsection in
MF=6, and NJOY uses it in preference to computing the recoil. For MT = 91 in
TENDL those arrays are tabulated on a coarse grid of about twenty uniform bins
starting at zero recoil energy, and they are systematically softer than the
evaluation's own neutron distribution requires. For Fe-56 at 14 MeV the stored
array has a mean of 0.175 MeV where momentum balance against the same file's
Kalbach-Mann neutron spectrum gives 0.286 MeV, and it places 43% of the recoils
below 47 keV even though the kinematic minimum is 8 keV. The same pattern holds
for every TENDL nuclide examined, with stored-to-balanced ratios of 0.28 to
0.71. The ENDF/B-VIII.1 recoil arrays for the same reaction are consistent with
momentum balance to within a few percent. OpenMC reports the momentum-conserving
result in both cases.

Because the light-ion model of :eq:`recoil-light-ion` is a calibrated surrogate
rather than evaluated data, agreement for the charged-particle channels should
be treated as approximate. Across the calibration set the mean light-ion energy
reproduces the evaluated value with a root-mean-square scatter of about 12% and
no significant bias, but individual nuclide-energy-channel combinations can
differ by 20% or more, and near threshold by considerably more. Note also that
the calibration and the comparison draw on the same family of evaluations, so
this is a measure of consistency with TENDL rather than of accuracy against
measured spectra.

-----------
Limitations
-----------

- S(:math:`\alpha,\beta`) and NCrystal thermal scattering produce no recoil
  record; only free-gas elastic scattering does.
- The recoil of a bound atom is treated as if it were free; no lattice binding
  energy is subtracted.
- De-excitation photon momentum is neglected for every reaction except capture.
- Fission fragments are not produced.
- Reactions whose exit channel cannot be determined from the MT number produce
  no record at all. In practice this means MT = 5, the ENDF catch-all, which
  carries inclusive product yields rather than a single recoil. Some
  evaluations put a substantial part of the charged-particle production there:
  ENDF/B-VIII.1 Fe-56 has 0.073 b in MT = 5 at 14 MeV against 0.114 b in (n,p),
  and its (n,alpha) cross section is a tenth of the value other libraries give
  because the rest of that channel sits in MT = 5.
- Multi-neutron final states are sampled independently rather than from a
  correlated joint distribution, and no ordering of the emitted products is
  implied by the order in which they are constructed.
- An event whose exit channel cannot be completed produces no record, so a
  channel with a low sampling success rate is underrepresented in the tally by
  the amount the counters report.
- The light-ion model omits optical-model transmission coefficients, explicit
  level densities, channel competition, direct reactions, and evaluated
  sequential decay.
- Recoil excitation, isomeric state, and subsequent gamma recoil are not
  represented.
- Atomic electron binding energy is neglected in the mass budget, which leaves
  a few keV unaccounted for in a charged-particle channel on a heavy target.
- All massive-particle kinematics are nonrelativistic, which understates the
  recoil energy by roughly :math:`E / 2 m_n c^2` — under 1% at 14 MeV.
- Only the continuous-energy transport mode produces recoils.

----------
References
----------

.. [Weisskopf1940] V. F. Weisskopf and D. H. Ewing, "On the Yield of Nuclear
   Reactions with Heavy Elements," *Physical Review* **57**, 472-485 (1940).
   `<https://doi.org/10.1103/PhysRev.57.472>`_

.. [HillWheeler1953] D. L. Hill and J. A. Wheeler, "Nuclear Constitution and
   the Interpretation of Fission Phenomena," *Physical Review* **89**,
   1102-1145 (1953). `<https://doi.org/10.1103/PhysRev.89.1102>`_

.. [Kalbach1988] C. Kalbach, "Systematics of continuum angular distributions:
   Extensions to higher energies," *Physical Review C* **37**, 2350-2370 (1988).
   `<https://doi.org/10.1103/PhysRevC.37.2350>`_

.. [Koning2012] A. J. Koning and D. Rochman, "Modern Nuclear Data Evaluation
   with the TALYS Code System," *Nuclear Data Sheets* **113**, 2841-2934 (2012).
   `<https://doi.org/10.1016/j.nds.2012.11.002>`_

.. [Gilbert2015] M. R. Gilbert, J. Marian, and J.-Ch. Sublet, "Energy spectra of
   primary knock-on atoms under neutron irradiation," *Journal of Nuclear
   Materials* **467**, 121-134 (2015).
   `<https://doi.org/10.1016/j.jnucmat.2015.09.023>`_

.. [MacFarlane2016] R. E. MacFarlane, D. W. Muir, R. M. Boicourt, A. C. Kahler,
   and J. L. Conlin, "The NJOY Nuclear Data Processing System, Version 2016,"
   Los Alamos National Laboratory report LA-UR-17-20093 (2016).
   `<https://doi.org/10.2172/1338791>`_

With survival biasing, the implicit absorption recoil is created with the
absorbed weight but the neutron then continues to a scattering event, so
:class:`ReactionFilter` reports the scattering MT rather than the absorption
one. Reaction-resolved absorption PKA tallies should use analog absorption
(``settings.survival_biasing = False``).
