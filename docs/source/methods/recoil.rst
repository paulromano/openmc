.. _methods_recoil:

==================
Recoil Production
==================

.. currentmodule:: openmc

When :attr:`Settings.recoil_production` is enabled, every continuous-energy
neutron collision creates additional entries in the secondary bank describing
the *heavy residual nucleus* left by the reaction and, optionally, the light
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

The identity of the residual follows from the MT number,

.. math::

    Z_R = Z_T - \sum_i Z_i , \qquad A_R = A_T + 1 - \sum_i A_i ,

and the residual is always represented in its ground state. Residual excitation
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

For MT = 51-90 the outgoing neutron determines the residual completely, so
subtracting the sampled neutron momentum is exact for the two-body level
transition described by the evaluated angular distribution. Momentum carried by
the de-excitation photons is neglected; for a residual of mass :math:`M_R` the
photons change the mean recoil energy by
:math:`\overline{|\sum_k \mathbf{p}_{\gamma,k}|^2} / 2M_R`, which is well under
1 keV for structural nuclides.

Continuum inelastic scattering
------------------------------

MT = 91 is treated the same way: the residual recoils against the sampled
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

The residual recoils against the emitted photons,
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

that NJOY's HEATR module uses when explicit recoil data are absent. The two
agree in the mean but not in shape.

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

and the emission removes :math:`E_b^\text{cm}(1 + m_b/M_D)` from :math:`U`.
Emission order is randomized so no ion is systematically favoured.

The centre-of-mass energy is sampled from

.. math::
    :label: recoil-light-ion

    P(E) \propto E\, T_C(E) \sqrt{1 - E/E_b^\text{max}} , \qquad
    T_C(E) = \left[1 + e^{(V_C - E)/\Delta}\right]^{-1} , \qquad
    V_C = \frac{1.44\ \text{MeV fm}\ Z_b Z_D}{r_0 (A_b^{1/3} + A_D^{1/3})} .

The factor :math:`E\,T_C(E)` stands in for the inverse-reaction cross section of
a Weisskopf-Ewing evaporation spectrum, and the square root for the level
density of the residual. The effective barrier radius :math:`r_0 = 1.8` fm and
diffuseness :math:`\Delta = 0.8` MeV are *not* optical-model quantities: the
radius is larger and the barrier softer than a geometric one so that the single
smooth transmission factor also absorbs sub-barrier tunnelling. They were
calibrated by matching the mean centre-of-mass ejectile energy of
:eq:`recoil-light-ion` against the evaluated ENDF MF=6 spectra of MT = 103-107
for thirteen nuclides between beryllium and tantalum from 5 to 20 MeV.

The endpoint :math:`E_b^\text{max}` used in :eq:`recoil-light-ion` is the one
belonging to the channel that emits this ion *alone* — for a proton, the Q value
of (n,p) — because in a channel such as (n,np) the charged particle is
physically emitted first, from the hot compound nucleus, and only then does the
neutron follow. The sample is then truncated to what this particular event can
still afford. Discrete charged-particle levels (MT = 600-849) are exactly
two-body, so their ejectile energy is fixed at :math:`E_b^\text{max}` rather
than sampled.

The direction uses the Kalbach-Mann form of evaluated MF=6 LANG=2 data,

.. math::

    f(\mu) = \frac{a}{2 \sinh a}
             \left[\cosh(a\mu) + r \sinh(a\mu)\right] ,

with the slope :math:`a` from Kalbach's systematics and the pre-equilibrium
fraction taken as :math:`r = E_b^\text{cm} / E_b^\text{max}`. That reproduces
the qualitative behaviour of evaluated :math:`r` values, which rise from nearly
zero at low outgoing energy to 0.5-0.9 near the kinematic maximum.

Setting ``light_ion_model`` to ``'none'`` skips this model entirely, in which
case the residual of a charged-particle channel recoils against the incident
neutron alone.

Fission
-------

Fission is excluded. Fission-fragment recoil is not modelled.

.. _methods_recoil_validation:

--------------------------------
Comparison With NJOY Group Data
--------------------------------

The usual reference for PKA spectra is a group-wise recoil matrix produced by
NJOY, as consumed by codes such as SPECTRA-PKA. OpenMC agrees closely with those
matrices where both sides derive the recoil from the same two-body kinematics,
and differs in three understood places.

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
distribution and 0.8473 as implied by the NJOY recoil matrix built from the same
evaluation, a 4% difference in :math:`1 - \bar\mu` and therefore in the mean
elastic recoil energy. This is a data-processing difference upstream of the
transport code, not a difference between the recoil models.

**Some evaluations store recoil arrays that violate momentum conservation.**
Evaluations generated with TALYS write an explicit heavy-residual subsection in
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
be treated as approximate. Across the calibration set the mean ejectile energy
reproduces the evaluated value with a root-mean-square scatter of about 12% and
no significant bias, but individual nuclide-energy-channel combinations can
differ by 20% or more, and near threshold by considerably more.

-----------
Limitations
-----------

- S(:math:`\alpha,\beta`) and NCrystal thermal scattering produce no recoil
  record; only free-gas elastic scattering does.
- The recoil of a bound atom is treated as if it were free; no lattice binding
  energy is subtracted.
- De-excitation photon momentum is neglected for every reaction except capture.
- Fission fragments are not produced.
- Multi-neutron final states are sampled independently rather than from a
  correlated joint distribution.
- The light-ion model omits optical-model transmission coefficients, explicit
  level densities, channel competition, direct reactions, and evaluated
  sequential decay.
- Residual excitation, isomeric state, and subsequent gamma recoil are not
  represented.
- All massive-particle kinematics are nonrelativistic, which understates the
  recoil energy by roughly :math:`E / 2 m_n c^2` — under 1% at 14 MeV.
- Atomic masses are used as a proxy for nuclear masses.
- Only the continuous-energy transport mode produces recoils.

With survival biasing, the implicit absorption recoil is created with the
absorbed weight but the neutron then continues to a scattering event, so
:class:`ReactionFilter` reports the scattering MT rather than the absorption
one. Reaction-resolved absorption PKA tallies should use analog absorption
(``settings.survival_biasing = False``).
