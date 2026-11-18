# Reference implementation of the anaesthetised-adult digital twin

> **If you use this code, or any part of it, please cite:**
>
> Birs, I., Mihai, M., Hegedus, E. & Muresan, C. Pharmacokinetic–pharmacodynamic
> blueprints for the anaesthetized adult digital twin: a coupled five-layer
> specification with an open reference implementation. *Sci. Rev.* **1**, 28
> (2026). https://doi.org/10.1038/s44573-026-00030-0

```bibtex
@article{Birs2026blueprints,
  author  = {Birs, Isabela and Mihai, Marcian and Hegedus, Erwin and Muresan, Cristina},
  title   = {Pharmacokinetic--pharmacodynamic blueprints for the anaesthetized adult digital twin: a coupled five-layer specification with an open reference implementation},
  journal = {Sci. Rev.},
  volume  = {1},
  pages   = {28},
  year    = {2026},
  doi     = {10.1038/s44573-026-00030-0}
}
```

Supplementary material for *Pharmacokinetic–Pharmacodynamic Blueprints for the
Anaesthetized Adult Digital Twin: A Coupled Five-Layer Specification with an
Open Reference Implementation*.

This package is the executable form of the specification in Section 5 of the
article. Any ambiguity in the text is resolved by the corresponding function
here. It is plain MATLAB.

The package is self-contained: distributed on its own it bundles a copy of the
`Patient` class, and inside the parent repository it uses the host definition.

## Running it

```matlab
cd pkpd_blueprint
run_all
```

`run_all` first checks every layer against a number reported in its source
publication, then runs the three demonstrations and writes the figures into
`figures/`. Total runtime is about ten minutes.

To use the twin on its own:

```matlab
setupPaths();
p    = Patient(1, 45, 172, 75, 'M');     % the host simulator's Patient class
twin = twinInit(p);
sim  = twinSimulate(twin, 60, 5/60, @(t,o) regimenCase(t,o,75,[]));
```

## Layout

| folder | contents |
|---|---|
| `models/` | one file per published model, each citing its source |
| `core/` | state map, coupled dynamics, output map, integrator, estimation |
| `experiments/` | the three demonstrations and the dosing regimens |
| `figures/` | output |

The state vector has 25 entries and is enumerated in `core/twinIndex.m` with
units and initial conditions. Drug quantities are carried as amounts rather
than concentrations; rocuronium and sugammadex are carried in micromoles
because their interaction is a 1:1 chemical reaction; norepinephrine is
carried in nanomoles, the unit of both its kinetic and its dynamic model.

## What each layer comes from

| layer | model | source |
|---|---|---|
| Pharmacokinetics | propofol, remifentanil | Eleveld 2018 *Br J Anaesth* 120:942–59; Eleveld 2017 *Anesthesiology* 126:1005–18 |
| L1 hypnosis | propofol–remifentanil surface, latent state, BIS observation model | Bouillon 2004 *Anesthesiology* 100:1353–72; confounders from Dahaba 2005 *Anesth Analg* 101:765–73 |
| L2 nociception | hierarchical tolerance surface | Bouillon 2004 |
| L3 respiration | CO2 kinetics with separate CO2 effect sites; joint depression | Bouillon 2003 *Anesthesiology* 99:779–87; Bouillon 2004 *Anesthesiology* 100:240–50; Nieuwenhuijs 2003 *Anesthesiology* 98:312–22; Koo 2015 *Yonsei Med J* 56:1408–14 |
| L4 haemodynamics | turnover model with baroreflex feedback; vasopressor | Su 2022 *Br J Anaesth* 128:806–16; Su 2023 *Br J Anaesth* 131:222–33; Li 2024 *Clin Pharmacokinet* 63:1597–608; Joachim 2024 *Br J Clin Pharmacol* 90:2861–9 |
| L5 neuromuscular | rocuronium–sugammadex with 1:1 encapsulation; arousal gate | Kleijn 2011 *Br J Clin Pharmacol* 72:415–33; Kim 2023 *Sci Rep* 13:20776 |
| Stimulus | multiplicative perturbation of the production rates | Garraza-Obaldia 2024 *Pharmaceutics* 16:1615 |

## Verification

`verifyAgainstLiterature` checks fourteen quantities against their sources,
including the Eleveld reference individual, the single-drug ventilatory
depressions of Nieuwenhuijs, the fixed point and baseline pressure of the
haemodynamic turnover model, the train-of-four at the reported potency, and
the norepinephrine clearance reduction reported by Li. All fourteen pass.
Each check is labelled. The seven `[source]` checks require the covariate
equations to reproduce a value reported by the publication (the Eleveld
reference individual and the Li clearance reduction). The seven
`[structure]` checks are true by construction; they guard against coding
errors but are not independent evidence.

One informational line is printed rather than tested. The Discussion of
Su 2023 quotes a fall in pressure of 11.7% and in heart rate of 2.7% at
1 µg mL⁻¹ propofol with 1 ng mL⁻¹ remifentanil, obtained from a simulation of
the full model; the infinite-time steady state of the published equations
gives −17% and +6%. The difference is in how the figure is computed, not in
the equations.

## Correction relative to the article

The endogenous norepinephrine production of Li et al. (Clin Pharmacokinet
2024;63:1597-1608, Table 3) is 497.7 ng/min. The version used for the
article read this value as nmol/h, which put the endogenous baseline at
4.1 nmol/L instead of the 1.4 nmol/L reported by Li et al. The age
coefficient on clearance has also been corrected to the published -0.344.
Both corrections are in `models/pkLiNorepinephrine.m`.

Propofol reduces norepinephrine clearance, so endogenous norepinephrine,
and through it MAP, responds to propofol even when no vasopressor is given.
The correction therefore changes the numbers that use MAP. Every other
number in the article, and every conclusion, is unchanged.

| quantity | article | this package |
|---|---|---|
| Fisher RSE with BIS and MAP: C50, ke0, CL (%) | 8.0, 5.5, 6.8 | 9.6, 6.2, 8.2 |
| max absolute correlation, BIS and MAP | 0.983 | 0.988 |
| condition number, BIS and MAP | 287 | 411 |
| Demonstration 2, MAP MDAPE (%), population / individualised | 20.0 / 4.5 | 20.3 / 3.6 |
| Demonstration 2, MAP RMSE (mmHg), population / individualised | 9.0 / 2.4 | 8.9 / 1.9 |
| Demonstration 2, BIS MDAPE interquartile range, individualised | 2.3-9.6 | 2.5-10.2 |
| Demonstration 2, BIS RMSE interquartile range, individualised | 1.4-5.4 | 1.3-5.9 |
| Demonstration 2, propofol Ce MDAPE (%), individualised | 57.4 | 57.9 |
| Demonstration 2, C50 recovery error (%), fitted | 70.8 | 70.5 |

## Known discrepancies and limitations

The package reproduces the article as published, apart from the
norepinephrine correction above. The following points are known and are
listed so that anyone reusing the code can take them into account.

### Differences between the code and the article's text

1. **Respiratory depression (layer 3).** The article specifies the two
   single-agent depression terms of Bouillon et al. (remifentanil C50
   0.92 ng/mL, gamma 1.25; propofol C50 1.33 ug/mL, gamma 1.68) combined
   multiplicatively. The code instead uses, by default, the synergistic
   response surface of Nieuwenhuijs et al. (`respSynergy = true` in
   `core/twinInit.m`). The arterial CO2 of Demonstration 1 printed in the
   article (peak 52 mmHg) was produced by this default. Setting
   `respSynergy = false` does not give the article's equation: it combines
   single-agent terms derived from Nieuwenhuijs, not those of Bouillon.
   With the Bouillon terms as written in the article, the peak is about
   66-67 mmHg. The combined depression of the Nieuwenhuijs surface has not been
   checked against the original publication beyond the two single-agent
   values.
2. **Molar masses.** Rocuronium and sugammadex are converted with 529.8 and
   2178 g/mol, as stated in the article. These are the masses of the
   rocuronium cation and of sugammadex sodium. Clinical doses are labelled
   as rocuronium bromide (609.7 g/mol) and as sugammadex (2002.1 g/mol), so
   the code overstates the molar dose of rocuronium by about 15% and
   understates that of sugammadex by about 8%. With the label masses, the
   train-of-four ratio before reversal in Demonstration 1 is about 0.8
   instead of 0.61.
3. **Rocuronium effect-site rate constant.** Kleijn et al. scale it with
   body weight, 0.134 (BW/70)^-0.25 1/min; the code uses 0.134 at every
   weight. The difference is small: the train-of-four ratio before reversal
   is 0.60 instead of 0.61.
4. **Arterial CO2 under controlled ventilation.** The article states 40 mmHg;
   the code starts at, and holds, the awake baseline of 41.6 mmHg. With
   40 mmHg the post-extubation peak is 51 instead of 52 mmHg.
5. **Haemodynamic baselines.** The code uses rounded values (TPR 0.016,
   HR 56.0, drift terms 0.090 and 0.12) where the article's state table
   lists 0.0163, 56.1, 0.0899 and 0.121. Baseline MAP is 89.9 instead of
   91.8 mmHg.

### Realism of the coupled simulation

6. **Cardiac output and remifentanil.** Cardiac output scales every
   clearance of propofol and remifentanil. The coupling was established for
   propofol and fentanyl (Bienert et al.); remifentanil is cleared by blood
   and tissue esterases and its clearance is unlikely to follow cardiac
   output. In Demonstration 1 the coupling raises the remifentanil
   effect-site concentration from 3.3 to 4.4 ng/mL at 30 min and from 5.7
   to 8.7 ng/mL at 55 min. The cardiac-output ratio also reaches the lower
   clip of 0.6 (`coBounds`) during haemorrhage and in Demonstration 2, so
   the scaling sits at that bound for long periods.
7. **Depth of hypotension.** The haemodynamic layer predicts larger falls in
   pressure than are typical in practice. In Demonstration 1 MAP falls from
   90 to about 65 mmHg within five minutes and norepinephrine is needed
   throughout maintenance; in Demonstration 2, which gives no vasopressor,
   the population twin settles at about 50 mmHg. Contributing factors are
   the published equations themselves (see the informational line under
   Verification), the remifentanil exposure of item 6, and the loop by which
   a lower cardiac output raises drug concentrations and lowers cardiac
   output further.
8. **Response to surgical stimuli.** Laryngoscopy and incision barely change
   MAP or heart rate, because the turnover of the haemodynamic states is
   slow (half-life of about 10 min) compared with stimuli of a few minutes.
9. **Ceiling of the arousal probability.** The two-point logistic of
   `core/arousalProbability.m` saturates at 0.728, so even a patient at
   BIS 90 has an arousal probability of 73%.

## Assumptions, stated plainly

1. **Stimulus intensity.** No published surface parameterises stimulus
   intensity as a continuous variable. The scale factor relating surgical
   stimuli to laryngoscopy (`kappa` in `models/pdNociception.m`) is assumed.
   It is the least constrained quantity in the specification.
2. **Arousal at reversal.** `core/arousalProbability.m` is a two-point
   logistic. Two points do not identify a curve; it carries the right
   qualitative behaviour, not an estimated relationship.
3. **Confounder weights.** The bias terms in `models/obsBis.m` default to
   zero because they are not established quantitatively. The structure is
   explicit so that they can be filled in when they are. The exception is
   the electromyographic term, which defaults to 6 BIS units at full
   stimulus and full neuromuscular recovery (`emgBias` in
   `core/twinInit.m`). This value is assumed.
4. **Stimulus on the haemodynamic layer.** Garraza-Obaldia et al. estimated
   the magnitude 0.284 for a decaying pulse on MAP production during
   hysteroscopy and found no effect on heart rate. Here a sustained
   stimulus profile multiplies the production of resistance and heart
   rate, as specified in Eq. (l4) of the article. Only the magnitude is
   taken from the publication.

## Note

These are in-silico demonstrations. They establish that the specification is
complete, self-consistent and executable, and they quantify properties of the
model architecture. They are not clinical validation, they provide no evidence
that any constituent model is correct, and no patient data were used. Section
6.5 of the article sets out what validation would require.
