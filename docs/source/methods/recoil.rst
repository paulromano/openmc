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
Nuclide masses come from the tabulated atomic masses; a nuclide missing from
that table falls back on :math:`A` atomic mass units.

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

For a parent of mass :math:`M` with internal energy :math:`U` emitting ion
:math:`b` and leaving daughter :math:`D`, the two-body kinematic maximum is

.. math::

    E_b^\text{max} = U \frac{M_D}{m_b + M_D} ,

and the emission removes :math:`E_b^\text{cm}(1 + m_b/M_D)` from :math:`U`, so
the daughter is left with excitation
:math:`E_x = U\left(1 - E_b^\text{cm}/E_b^\text{max}\right)`. Emission order is
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

    P(E) \propto E\, T_C(E) \sqrt{1 - E/E_b^\text{max}} .

*Inverse cross section.* For a charged light ion :math:`\sigma_\text{inv}` is
governed almost entirely by the Coulomb barrier, so it is replaced by a barrier
transmission coefficient,

.. math::
    :label: recoil-barrier

    T_C(E) = \left[1 + e^{(V_C - E)/\Delta}\right]^{-1} , \qquad
    V_C = \frac{1.44\ \text{MeV fm}\ Z_b Z_D}{r_0 (A_b^{1/3} + A_D^{1/3})} .

This is the Hill-Wheeler transmission through a parabolic barrier
[HillWheeler1953]_, for which the diffuseness is set by the barrier curvature,
:math:`\Delta = \hbar\omega / 2\pi`. The value used, :math:`\Delta = 0.8` MeV,
corresponds to :math:`\hbar\omega = 5.0` MeV, at the permeable end of the usual
few-MeV range. The reduced radius :math:`r_0 = 1.8` fm is larger than a
geometric touching radius of about 1.4 fm and therefore lowers :math:`V_C` by a
uniform factor of 0.78 — from 5.35 to 4.16 MeV for p + Mn-55, and from 9.24 to
7.19 MeV for :math:`\alpha` + Cr-53. That is the same kind of empirical barrier
reduction as the :math:`k_j` factors of the Dostrovsky inverse-cross-section
parameterization used in evaporation codes [Dostrovsky1959]_, which are likewise
below unity and typically of order 0.7-0.9 for medium-mass nuclei. Here it came
out of the fit described below rather than being imposed.

*Level density.* The square root is **not** a nuclear level density, which would
rise roughly exponentially with :math:`E_x`. Charged-particle spectra in
TALYS-based evaluations [Koning2012]_ mix compound-nucleus evaporation with much
harder pre-equilibrium and direct emission, and applying a compound level
density to the whole spectrum makes it far too soft. The weakly rising
:math:`\sqrt{E_x}` factor is an empirical compromise between the two components.

Three numbers were fitted: the barrier radius :math:`r_0`, the diffuseness
:math:`\Delta`, and the exponent of the :math:`E_x` factor. They were obtained
by matching the mean centre-of-mass light-ion energy of :eq:`recoil-light-ion`
against the evaluated ENDF MF=6 spectra of MT = 103-107 in TENDL, for thirteen
nuclides between beryllium and tantalum from 5 to 20 MeV. Equation
:eq:`recoil-light-ion` should therefore be read as a calibrated surrogate whose
*form* is borrowed from statistical-model theory and whose *parameters* come
from evaluated data — not as an implementation of any published
nuclear-reaction model.

Endpoint and direction
~~~~~~~~~~~~~~~~~~~~~~

The endpoint :math:`E_b^\text{max}` used in :eq:`recoil-light-ion` is the one
belonging to the channel that emits this ion *alone* — for a proton, the Q value
of (n,p) — because in a channel such as (n,np) the charged particle is
physically emitted first, from the hot compound nucleus, and only then does the
neutron follow. The sample is then truncated to what this particular event can
still afford. Discrete charged-particle levels (MT = 600-849) are exactly
two-body, so their light-ion energy is fixed at :math:`E_b^\text{max}` rather
than sampled.

The direction uses the Kalbach-Mann form of evaluated MF=6 LANG=2 data,

.. math::

    f(\mu) = \frac{a}{2 \sinh a}
             \left[\cosh(a\mu) + r \sinh(a\mu)\right] ,

with the slope :math:`a` from Kalbach's systematics [Kalbach1988]_ — the same
expression the evaluations themselves use — and the pre-equilibrium fraction
taken as :math:`r = E_b^\text{cm} / E_b^\text{max}`. The systematics for
:math:`a` are used as published; only :math:`r` is a substitution, and it
reproduces the qualitative behaviour of evaluated :math:`r` values, which rise
from nearly zero at low outgoing energy to 0.5-0.9 near the kinematic maximum.

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
  correlated joint distribution.
- The light-ion model omits optical-model transmission coefficients, explicit
  level densities, channel competition, direct reactions, and evaluated
  sequential decay.
- Recoil excitation, isomeric state, and subsequent gamma recoil are not
  represented.
- All massive-particle kinematics are nonrelativistic, which understates the
  recoil energy by roughly :math:`E / 2 m_n c^2` — under 1% at 14 MeV.
- Atomic masses are used as a proxy for nuclear masses.
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

.. [Dostrovsky1959] I. Dostrovsky, Z. Fraenkel, and G. Friedlander, "Monte Carlo
   Calculations of Nuclear Evaporation Processes. III. Applications to
   Low-Energy Reactions," *Physical Review* **116**, 683-702 (1959).
   `<https://doi.org/10.1103/PhysRev.116.683>`_

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
