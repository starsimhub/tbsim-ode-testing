# Proposal 07: agents' starting states given by age group

**For:** Starsim v4 (`Agents.init_post()` in `starsim/solvers/agents.py`).
**Status:** draft, not sent.
**Track:** starsim `origin/v4` @ `d83f4f87` (4.0.0rc4).
**Evidence:** [exp 30](https://github.com/starsimhub/tbsim-ode-testing/blob/main/experiments/30_tb_corrected_demography/SUMMARY.md) (Finding 3), and [`repro.py`](https://github.com/starsimhub/tbsim-ode-testing/blob/main/proposals/07_starting_states/repro.py) here, which uses Starsim alone.

## Summary

We started a TB model's agents and its compartmental (ODE) twin from the same table: each state's share by single year of age and sex (`init=` with `ss.by` tables). The compartmental solver started exactly at the table. Agents didn't, in three ways:

| What happens | When | Effect in our runs | Proposed change |
|---|---|---|---|
| Each state's count in each group is rounded to whole agents | Always, with shares by group | Rare states vanish in small groups: 65+ started with 0.49–0.73 of the intended prevalence (100,000 agents, single years by sex) | 1: draw the counts, or round stochastically |
| States are drawn one after another from the agents left in the first state, each with its share unadjusted | A table read per agent: bins finer than the agents' groups, or a property such as sex without a stratification by it | Later states get too few agents: most adults started susceptible | 2: draw each agent's state once, from all the shares together |
| A share for the first state can't be read per agent | As above, with every state given | Refused: "the initial values sum to 0" | 2 (the same change) |

Agents and the ODE matched closely once we worked around all three (exp 30).

## 1. Counts by group are rounded to whole agents

**What happens now.** For shares given by group, `init_post()` gives each state an exact count in each group, keeping the group's total:

```python
plan[g] = (inds, share*len(inds), Compartmental.round_preserving(share*len(inds))) # The number in each state is exact, as it is in a compartmental run
```

`round_preserving()` floors every count and gives the leftover agents to the largest remainders. In a small group, a rare state's expected count is below 1. Its remainder loses to the larger states', so it gets nobody, in every seed. A later step puts one agent in a state that would otherwise have none at all, but only once, in one group.

**What it does** (`repro.py`, part 1): 0.3% infected at every age, 20,000 agents stratified by single years of age.

| Infected agents at the start | 0–14 | 15–49 | 50–64 | 65–79 | 80+ |
|---|---|---|---|---|---|
| Expected | 23.2 | 29.9 | 5.0 | 1.7 | 0.2 |
| Seeds 1, 2, 3 | 21, 21, 20 | 33, 34, 35 | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 |

**In our runs:** with 100,000 agents stratified by single years and sex, the groups at older ages hold about 50–150 agents. Active TB (0.2–0.35% for each of its two states), recent infection and susceptibility mostly rounded to nobody there. Agents' 65+ prevalence started at 0.49–0.73 of the compartmental solver's, and took decades to recover.

**Proposal.** Keep the counts right on average, rather than exact in each group:

- **The simplest:** draw each group's counts from a multinomial with its shares. In distribution, that's the same as drawing each agent's state from its group's shares.
- **Or round stochastically:** floor each count, then give the leftover agents to states with probability proportional to their remainders. That keeps each group's total and the variance small.
- **Either way, as an option or the default.** Exact counts suit large groups and comparisons with a compartmental run at the start. But they bias rare states downward in small groups, and the bias is the same in every seed, so reps don't average it out.

**Our workaround:** while agents' states are set, `Compartmental.round_preserving` is swapped for a multinomial draw seeded by the run's seed (`_drawn_start` in [exp 30's `run.py`](https://github.com/starsimhub/tbsim-ode-testing/blob/main/experiments/30_tb_corrected_demography/run.py)). With it, every state matched the table at every age within sampling noise, at 100,000 and 400,000 agents.

## 2. Shares read per agent are drawn one state at a time

**What happens now.** A share that agents read from each person's own properties (`init.on_people(sim)`) is drawn after the group step, one state at a time:

```python
first = self.members(self.graph.default_state)
share = ... init.lookup(self.sim, self.convert)[first]
chosen = first[self.draw(f'init_{state.name}', first) < share]
```

That happens when a table's bins of age are finer than the agents' age groups, or when a table is by a property that isn't a dimension of the sim, such as sex.

- **States drawn in turn:** each state is drawn from the agents still in the first state, with its share of everyone. So the second state gets share × (1 − the first state's share), the third less again, and the first state keeps the rest.
- **The first state can't be given this way:** its share isn't in the group shares, which then sum to 0, so the start is refused.

**What it does** (`repro.py`, parts 2 and 3): I = 0.1 and R = 0.5 at every age, by single year, with agents stratified by 5-year groups.

| Shares at the start | S | I | R |
|---|---|---|---|
| Expected | 0.40 | 0.10 | 0.50 |
| Agents stratified by 5-year groups (table read per agent) | 0.43–0.46 | 0.09–0.10 | 0.44–0.48 |
| Agents stratified by single years (table read per group) | 0.40 | 0.10 | 0.50 |

R comes out at about 0.5 × (1 − 0.1) = 0.45. With S = 0.4 given too, the 5-year version is refused ("The initial values of "sir" sum to 0, not 1, and it is not clear where the rest start"), while the single-year one runs.

**In our runs:** the TB model has seven starting states, and in adults the first (susceptible) holds about 0.2%. With an eighth-year table and agents stratified by WPP's age groups, each state after the first was drawn from a shrinking pool. So most adults stayed susceptible, and active TB was a fraction of its share.

**Proposal.** For shares read per agent, draw each agent's state once from all the states' shares together: the given shares, with the first state (or the one left without a value) taking the rest. That also lets the first state's share be given, as it can be by group.

**Our workaround:** a table by single year of age for both solvers, agents stratified by single years so the table is read by group, every state given but susceptible, and the sex split taken from a stratification by sex.

## Not proposed

- **Documentation:** the docstring of `ss.by` says a table by single year of age "works as agents whatever the model's age groups". That's true for one state, but not for several (part 2). A line saying how agents read starting shares (by group, rounded; or per agent) would help until either change is made.
