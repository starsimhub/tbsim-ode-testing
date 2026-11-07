# 19: Incidence in v4 (first entry to asymptomatic TB per episode), and a fast ODE for calibration

**Status:** draft design, for the user to confirm before anything runs.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`).

## Background

- **The definition:** incidence is the first entry to asymptomatic TB in each episode, as tb_LAI_TPT counts it (CONVENTIONS.md 5; user, 2026-11-05).
  - **An episode** starts when a latent infection progresses (to non-infectious or asymptomatic TB). It ends with recovery or death.
  - **Not counted:** someone who goes asymptomatic → non-infectious → asymptomatic within one episode, which is a re-entry.
- **What v4 counts now:** its only incidence output (`v4/common.py`, `outputs()`) counts every entry to asymptomatic TB, re-entries included. v4's ODE has no state that knows whether a non-infectious person has already been asymptomatic in this episode.
- **Why now:** exp 20's calibration has two targets built on incidence (mortality ÷ incidence and prevalence ÷ incidence). It needs incidence by this definition in the ODE, and it needs a fast copy of the ODE (an equilibrium in milliseconds, not a second), for IMIS.

## Questions

1. Does the ODE's incidence, counted as the first entry to asymptomatic TB per episode (tb_LAI_TPT's definition), match the agents' count by the same definition?
2. How much does that differ from v4's current output, which counts every entry to asymptomatic TB, re-entries included?
3. Does a standalone copy of the ODE (our own code, outside Starsim's solver), with this incidence, reproduce v4's ODE exactly, at equilibrium and over time?

## Design

**The ODE (standalone copy, in `tbode4`):** v4's ODE is Starsim's compartmental solver applied to tbsim's declared model. It steps a sim at fixed `dt`, about a second per equilibrium. The copy writes the same equations as one function of the state and the parameters, and finds the equilibrium directly by root-finding, in milliseconds, as v0.10's `TBODE` did alongside tbsim's ODE. That speed is what exp 20's IMIS needs. It is a second implementation of the model, so it is used only if it reproduces v4's ODE (question 3). v4's own ODE stays the reference for ODE↔agents comparisons.

- **The equations:** the same as `tbsim.TB` as v4 declares them, with the two latent states and homogeneous mixing, and births replacing deaths.
- **One extra state, for counting only:** non-infectious TB is split into NON0 ("not yet asymptomatic this episode") and NON1 ("has been"), with identical rates, as `TBODE` did in v0.10.
  - **Where people enter:** progression from latent infection (recent or remote) enters NON0; asymptomatic → non-infectious enters NON1.
  - **How incidence is counted:** latent → asymptomatic, plus NON0 → asymptomatic.
  - **Effect on dynamics:** none, since NON0 + NON1 behave exactly as v4's single non-infectious state.
- **Equilibrium:** found by root-finding on the right-hand side, from the ODE's own run as a start.
- **Check against v4's ODE:** every state's share and every flow, at equilibrium and over 100 years from a seed, for set A and set B (two-state, φ = 0.01). Within about 1e-6 at equilibrium; within the step-size differences over time (v4's ODE also runs at fixed steps).

**Agents (v4 `tbsim.TB`, homogeneous mixing):** an analyzer counts incidence from transitions, with a per-agent flag for "has been asymptomatic this episode":

- **Always counted:** latent → asymptomatic (recent or remote).
- **Counted only if the flag is off:** non-infectious → asymptomatic (`worsening`).
- **Never counted:** symptomatic → asymptomatic (`remission`).
- **The flag:** set on any entry to asymptomatic TB; cleared when an episode ends (recovery, death) or a new one starts.
- **Re-entries:** v4's every-entry count is also recorded, to measure the difference.

**Runs:** as exp 17's question 1:

- **Model:** sets A and B (two-state, φ = 0.01, latent reinfection on), `beta` calibrated to 500/100k.
- **Agents:** 50,000, started at the ODE's equilibrium, 1900–2000, at monthly and weekly `dt`, 40 reps each (new seeds).
- **Measured:** incidence per person-year, agents ÷ ODE, mean over 1950–2000 (mean ± SE, 10–90% range), for both definitions. Also the share of entries to asymptomatic TB that are re-entries, in both.

## What would answer the questions

- **Question 1:** agents ÷ ODE within about 2 SE of 1. Exp 13 found every-entry incidence within 0.3% per prevalent case.
- **Question 2:** the size of the re-entry share. It matters for exp 20's targets, and for comparing with v4 outputs that use every entry.
- **Question 3:** the standalone ODE matches v4's ODE to the stated tolerance; otherwise it isn't used.

**Found while building it:** v4's existing incidence output (`v4/common.py`, `outputs()`) counts latent → asymptomatic and non-infectious → asymptomatic, but not remote latent → asymptomatic. With the two-state model it misses progression from the remote state. Here "every entry" includes it; v4's count as written is also recorded (`every_v4`).

## Files

- `v4/tbode4/fastode.py`: the standalone ODE, with NON0/NON1 and both incidence counts (`flows`, `summary`), its equilibrium by root-finding, and trajectories.
- `v4/tests/test_fastode.py`: question 3, a standing test. The standalone equilibrium is a fixed point of v4's ODE (within 1e-9 over 50 years), and trajectories from a seed agree within 1e-8 at weekly steps, for sets A and B. Rerun after bumping the v4 pins.
- `run.py`: `beta` calibrated on the standalone ODE; the agents' analyzer (`Incidence`); the runs → `results.obj` (calibrations, ODE flows, per-rep yearly means), `summary.csv`, `run.log`, `versions.json`.
- `SUMMARY.md`: findings, written after running.
