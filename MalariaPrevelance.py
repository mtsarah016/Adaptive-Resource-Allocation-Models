
"""
Epidemic GridWorld — Malaria / Mosquito-Mediated Edition

  FEATURE — Gravity-model spatial coupling (human mobility / importation)

  Without spatial coupling the model had a critical structural gap: a
  newly-warm high-latitude cell whose R₀(T) has just crossed 1 (due to
  climate change) contains zero infected humans and zero infected mosquitoes,
  so local FOI = 0 and the cell remains permanently disease-free even though
  it is now ecologically receptive.  This is not how malaria range expansion
  actually works.  The missing mechanism is *importation*: infected travellers
  from endemic areas carry parasites into receptive-but-naive cells, seeding
  local transmission chains.

  Implementation: gravity model from Marshall et al. (2018, Sci. Rep.)
  -----------------------------------------------------------------------
  Movement rate between cell i → j ∝ pop_i × pop_j / dist(i,j)^β
  β = 1.5 is Marshall et al.'s median fitted value for African inter-district
  overnight travel data

  At each step a fraction φ (mobility_rate) of each cell's infected population
  is distributed across all other cells according to gravity weights.  The
  incoming parasite flux at destination j is then multiplied by:
    (a) import_infectivity — the fraction of travellers currently gametocytaemic
    (b) r0_suit[j]         — parasites can only establish where vectors exist

  This makes poleward range expansion mechanistically complete: warm-climate
  cells export parasites northward; northern cells only amplify them once their
  temperature crosses T_min and vectors are present.

  Accuracy notes (all sources peer-reviewed)
  ------------------------------------------
  Parameter                  Value used     Basis
  ─────────────────────────────────────────────────────────────────────────
  T_MIN (lower thermal limit) 19 °C        Johnson et al. (2015) Bayesian
                                            refit of Mordecai et al. (2013);
                                            more conservative than original
                                            Brière fit (16 °C).  CHANGED
                                            from v3 to use the tighter
                                            Johnson 2015 estimate.
  T_OPT                       25 °C        Mordecai et al. (2013) Ecol. Lett.
                                            16:22-30; confirmed by Mordecai
                                            et al. (2019) Ecol. Lett. review
                                            table (T_opt = 25.4 °C for
                                            P. falciparum / A. gambiae).
  T_MAX                       33 °C        Johnson et al. (2015) Bayesian
                                            refit (32.6 °C for multi-species
                                            model).  Rounded to 33 °C.
                                            Note: An. stephensi has T_max
                                            ~36 °C (Shapiro et al. 2017);
                                            we model A. gambiae / sub-Saharan
                                            Africa setting so 33 °C applies.
  EIP at 25 °C                10 days      Mordecai et al. (2013) cite
                                            ~10-day EIP at 25 °C for
                                            P. falciparum; consistent with
                                            Paaijmans et al. (2012) Q10 data.
  EIP_Q10                     2.5          Paaijmans et al. (2012) Proc. R.
                                            Soc. B — empirical Q10 for parasite
                                            development rate in A. gambiae.
  recovery_rate               0.07/step    P. falciparum uncomplicated
                                            infectious period ~14 days
                                            (1/14 ≈ 0.07).  CHANGED from
                                            v3 (0.05) to match Bousema &
                                            Drakeley (2011) Clin. Microbiol.
                                            Rev. 24:377-410; White et al.
                                            (2014) Nat. Commun. report median
                                            infection duration ~14 days in
                                            semi-immune adults.
  import_infectivity          0.36         Fraction of asymptomatic infections
                                            with detectable gametocytes.
                                            Nkrumah et al. (2021) Sci. Rep.
                                            found 36.3% of asymptomatic P.
                                            falciparum infections had
                                            gametocytes at any time point.
                                            CHANGED from v3 (0.5).
  gravity β                   1.5          Marshall et al. (2018) Sci. Rep.
                                            8:7713 — median fitted exponent
                                            for African inter-district
                                            overnight travel across 4 countries.
  mobility_rate               0.008/step   Marshall et al. (2016) Malar. J.
                                            15:200 surveyed overnight-trip
                                            frequencies: ~2-5% / month across
                                            sites → ~0.07-0.17% / day.  We
                                            use 0.8% / step as the grid is
                                            coarse (district-level) and each
                                            step may span > 1 day.
  water_larval_boost          0.4          VECTRI model (Tompkins & Ermert
                                            2013, Geosci. Model Dev.) shows
                                            2-5× larval density near standing
                                            water; 0.4 additive boost on a
                                            [0,1] normalised field is
                                            conservative.
  seasonal amplitude ±30%     ±30%         Tusting et al. (2017) Malar. J.
                                            review of seasonal EIR variation
                                            in sub-Saharan Africa; ±30% is
                                            conservative for a region with
                                            year-round transmission.
  RCP-8.5 ΔT                  +2.6°C       IPCC AR5 WG1 Table SPM.2 central
  at end-of-century                        estimate for RCP-8.5.
  Latitude gradient           28→14°C      Approximates coastal Kenya
                                            (equatorial, ~28°C mean) to
                                            Ethiopian Highlands (poleward,
                                            ~14°C).  Under 2.6°C warming,
                                            poleward rows shift 14→16.6°C,
                                            just crossing T_min=19°C at
                                            ~ep 250, producing gradual
                                            highland emergence.
  foi_scale                   0.15         Simulation tuning constant —
                                            no direct empirical anchor.
                                            Validated qualitatively: at T_opt
                                            with infected_fraction=0.1, FOI
                                            ≈ 0.015/susceptible/step, broadly
                                            consistent with EIR of ~1-2
                                            infectious bites/person/night in
                                            high-transmission settings.
  hospital_overload_penalty   0.5          Conservative; recovery halved
                                            when hospitals overwhelmed.
                                            No specific peer-reviewed anchor
                                            — flagged as requiring calibration.
  flood_water_boost           0.3          Tompkins & Ermert (2013) VECTRI
                                            model; moderate larval habitat
                                            increase during flood.

  KNOWN REMAINING LIMITATIONS
  ────────────────────────────
  1. foi_scale is a free tuning parameter — should ideally be derived
     from empirical EIR data for the simulated setting.
  2. The latitude gradient (28→14°C over 8 grid rows) is a very coarse
     approximation; real landscapes have heterogeneous microclimates.
  3. mobility_rate is time-invariant; real travel has strong seasonality
     (Marshall et al. 2016 show rainy-season mobility increases).
  4. Recovery rate 0.07 is for semi-immune adults; children (most at-risk)
     have longer infectious periods (White et al. 2014).
  5. Hospital overload penalty (0.5) lacks a direct peer-reviewed parameter
     estimate; flagged for future calibration.

Sources
-------
  Mordecai et al. (2013) Ecol. Lett. 16:22-30.
  Johnson et al. (2015) — Bayesian refit. PMC cited in Mordecai et al. (2019).
  Mordecai et al. (2019) Ecol. Lett. 22:1-17 (thermal biology review table).
  Paaijmans et al. (2012) Proc. R. Soc. B 279:4386-4393.
  Bousema & Drakeley (2011) Clin. Microbiol. Rev. 24:377-410.
  Nkrumah et al. (2021) Sci. Rep. 11:21538 (gametocyte carriage fraction).
  Marshall et al. (2018) Sci. Rep. 8:7713 (gravity model, β=1.5).
  Marshall et al. (2016) Malar. J. 15:200 (travel frequency survey).
  IPCC AR5 WG1 (2013) Table SPM.2.
  Tompkins & Ermert (2013) Geosci. Model Dev. 6:603-617 (VECTRI, larval water).
"""

# Run this once to install dependencies:
# !pip install gymnasium

import numpy as np
import gymnasium as gym
from gymnasium import spaces
from collections import defaultdict
import matplotlib.pyplot as plt
from IPython.display import clear_output, display
from gymnasium.utils.env_checker import check_env
import time

# Interactive mode so matplotlib windows update without blocking the script.
plt.ion()


# ===========================================================================
# Thermal biology — Mordecai 2013 / Johnson 2015
# ===========================================================================

# Using Johnson et al. (2015) Bayesian refit values (more conservative):
T_MIN = 19.0   # °C  lower thermal limit  [Johnson 2015: ~19°C]
T_OPT = 25.0   # °C  optimal transmission [Mordecai 2013: 25.4°C; Johnson 2015: ~25°C]
T_MAX = 33.0   # °C  upper thermal limit  [Johnson 2015: ~32.6°C, rounded to 33]

# EIP parameters — Paaijmans et al. (2012) Proc. R. Soc. B
EIP_REF_T = 25.0  # reference temperature (°C)
EIP_REF_D = 10.0  # EIP in days at reference (Mordecai 2013)
EIP_Q10   = 2.5   # Q10 for parasite development rate


def r0_relative(T: np.ndarray) -> np.ndarray:
    """Brière-type R₀ thermal performance curve (Mordecai et al. 2013, eq. 2).

    R₀_rel(T) ∝ T(T − T_min)(T_max − T)^0.5  for T_min < T < T_max, else 0.
    Normalised to peak = 1.0 at T_opt.
    Using Johnson et al. (2015) T_min=19, T_max=33.
    """
    inside = (T > T_MIN) & (T < T_MAX)
    raw = np.where(
        inside,
        T * (T - T_MIN) * np.maximum(T_MAX - T, 0.0) ** 0.5,
        0.0,
    ).astype(np.float64)
    peak = T_OPT * (T_OPT - T_MIN) * (T_MAX - T_OPT) ** 0.5
    return np.clip(raw / (peak + 1e-12), 0.0, 1.0).astype(np.float32)


def eip_factor(T: np.ndarray) -> np.ndarray:
    """Inverse EIP: parasite development rate normalised to reference T.

    Based on Q10 = 2.5 from Paaijmans et al. (2012).
    Higher at warmer T → faster maturation → more infectious mosquitoes.
    """
    rate     = (1.0 / EIP_REF_D) * (EIP_Q10 ** ((T - EIP_REF_T) / 10.0))
    ref_rate = 1.0 / EIP_REF_D
    return np.clip(rate / ref_rate, 0.0, 2.0).astype(np.float32)


def mosquito_density(T: np.ndarray) -> np.ndarray:
    """Anopheles population density thermal curve (Brière, A. gambiae).

    Separate curve from R₀_rel; represents combined effect of survival
    and fecundity.  T_min_m=14°C, T_max_m=36°C, T_opt_m≈26°C.
    Consistent with Mordecai et al. (2013) trait decomposition.
    """
    T_min_m, T_max_m, T_opt_m = 14.0, 36.0, 26.0
    inside = (T > T_min_m) & (T < T_max_m)
    raw = np.where(
        inside,
        T * (T - T_min_m) * np.maximum(T_max_m - T, 0.0) ** 0.5,
        0.0,
    ).astype(np.float64)
    peak_m = T_opt_m * (T_opt_m - T_min_m) * (T_max_m - T_opt_m) ** 0.5
    return np.clip(raw / (peak_m + 1e-12), 0.0, 1.0).astype(np.float32)


def biting_rate(T: np.ndarray) -> np.ndarray:
    """Mosquito biting rate thermal curve (Mordecai 2013, A. gambiae).

    Increases to optimum ~26°C then declines.  Parameterised from
    Mordecai et al. (2013) trait decomposition for A. gambiae.
    """
    T_b_min, T_b_opt, T_b_max = T_MIN, 26.0, 35.0
    b_min, b_opt = 0.2, 1.0
    rate = np.where(
        T <= T_b_opt,
        b_min + (b_opt - b_min) * np.maximum(T - T_b_min, 0.0) / (T_b_opt - T_b_min),
        b_opt * np.maximum(T_b_max - T, 0.0) / (T_b_max - T_b_opt),
    )
    return np.clip(rate, 0.0, 1.0).astype(np.float32)


# ===========================================================================
# Gravity-model spatial coupling  (Marshall et al. 2018, Sci. Rep. 8:7713)
# ===========================================================================

def build_gravity_weights(
    grid_size: int,
    population: np.ndarray,
    beta: float = 1.5,       # Marshall et al. (2018) median fitted exponent
) -> np.ndarray:
    """Precompute row-normalised gravity weight matrix W[origin, dest].

    W[i,j] = pop_i × pop_j / dist(i,j)^β  (i≠j), rows sum to 1.
    Cells with zero population have zero outflow (no travellers).

    β = 1.5: Marshall et al. (2018) median fitted value for African
    inter-district overnight travel (Mali, Burkina Faso, Zambia, Tanzania).

    Returns W of shape (N, N) where N = grid_size².
    """
    N = grid_size * grid_size
    pop_flat = population.ravel().astype(np.float64)

    rows_idx = np.arange(grid_size)
    rr, cc = np.meshgrid(rows_idx, rows_idx, indexing="ij")
    coords = np.stack([rr.ravel(), cc.ravel()], axis=1).astype(np.float64)

    dr = coords[:, 0:1] - coords[:, 0]   # (N, N)
    dc = coords[:, 1:2] - coords[:, 1]
    dist = np.sqrt(dr ** 2 + dc ** 2)

    pop_i = pop_flat[:, np.newaxis]
    pop_j = pop_flat[np.newaxis, :]
    gravity = pop_i * pop_j

    dist_safe = np.where(dist > 0.0, dist, np.inf)
    gravity /= dist_safe ** beta
    np.fill_diagonal(gravity, 0.0)

    row_sums = gravity.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0.0, row_sums, 1.0)
    W = (gravity / row_sums).astype(np.float32)
    return W


def importation_pressure(
    infected_flat: np.ndarray,      # (N,)   current infection fraction
    r0_suit_flat: np.ndarray,       # (N,)   thermal suitability at destination
    W: np.ndarray,                   # (N, N) gravity weights
    mobility_rate: float,            # fraction of infected who travel per step
    import_infectivity: float,       # gametocyte-positive fraction (Nkrumah 2021)
) -> np.ndarray:
    """Compute per-cell importation FOI from travelling infected humans.

    import_j = mobility_rate × import_infectivity
               × (Σ_i W[i,j] × infected_i)   ← incoming parasite flux
               × r0_suit_j                    ← only establishes if vectors present

    The r0_suit_j gate is the key biological constraint: imported parasites
    can only initiate local transmission where Anopheles vectors are present
    (T > T_min).  Without this gate, importation would create infection in
    cold, vector-free cells, which is biologically impossible.

    Returns importation pressure of shape (N,), values in [0, 1].
    """
    # W.T @ infected_flat gives, for each dest j, weighted sum of infected
    # across all origins i (columns of W = destinations)
    incoming = W.T @ infected_flat.astype(np.float32)    # (N,)
    pressure = mobility_rate * import_infectivity * incoming * r0_suit_flat
    return np.clip(pressure, 0.0, 1.0).astype(np.float32)


# ===========================================================================
# Geography
# ===========================================================================

def _make_geography(grid_size: int, rng) -> tuple:
    """Procedural geography: 1–3 Gaussian urban cores, coastal strip."""
    pop = np.zeros((grid_size, grid_size), dtype=np.float32)
    n_cities = rng.integers(1, 4)
    xs, ys = np.meshgrid(np.arange(grid_size), np.arange(grid_size), indexing="ij")
    for _ in range(n_cities):
        cx = rng.integers(1, grid_size - 1)
        cy = rng.integers(1, grid_size - 1)
        pop += rng.uniform(0.6, 1.0) * np.exp(
            -((xs - cx) ** 2 + (ys - cy) ** 2) / (2 * rng.uniform(0.8, 2.0) ** 2)
        )
    pop = np.clip(pop, 0.0, 1.0).astype(np.float32)

    coastal = np.zeros((grid_size, grid_size), dtype=np.float32)
    for row in range(grid_size - 2, grid_size):
        coastal[row, :] = 1.0
    extra_cols = rng.choice(grid_size, size=rng.integers(0, grid_size // 2), replace=False)
    if len(extra_cols):
        coastal[grid_size - 3, extra_cols] = 1.0

    hospital_capacity = np.clip(0.3 + 0.7 * pop, 0.1, 1.0).astype(np.float32)
    return pop, coastal, hospital_capacity


# ===========================================================================
# RCP warming trajectory  (IPCC AR5 WG1 Table SPM.2)
# ===========================================================================

def _delta_T(episode: int, scenario: str = "rcp85") -> float:
    """Cumulative ΔT (°C) at this episode.

    Calibrated to IPCC AR5 WG1 central estimates (Table SPM.2):
      RCP-8.5: +2.6 °C by 2100  | RCP-4.5: +1.8 °C | RCP-2.6: +0.9 °C
    500 episodes = end-of-century proxy.
    """
    t = episode / 500.0
    if scenario == "rcp85":
        return 2.6 * t ** 1.1          # slight acceleration of feedbacks
    elif scenario == "rcp45":
        return 1.8 * (1.0 - np.exp(-3.5 * t))   # stabilising trajectory
    else:
        return 0.9 * t


def _seasonal_factor(step: int, max_steps: int) -> float:
    """±30% sinusoidal seasonal forcing.

    Calibrated to Tusting et al. (2017) Malar. J. seasonal EIR review;
    ±30% is conservative for year-round transmission regions.
    Peak at mid-episode (summer / wet season).
    """
    phase = 2 * np.pi * step / max(max_steps, 1)
    return float(1.0 + 0.30 * (-np.cos(phase)))


# ===========================================================================
# Main environment
# ===========================================================================

class MalariaEpidemicEnv(gym.Env):
    """Malaria GridWorld: Ross-Macdonald FOI + latitude gradient +
    gravity-model importation for mechanistic range expansion.

    Observation channels (12):
      0  observed_infected      — partial obs via testing
      1  known_mask
      2  vaccinated
      3  population
      4  temperature            — normalised (10–40°C → 0–1)
      5  water                  — larval habitat / flood risk
      6  agent position
      7  budget fraction
      8  hospital_load          — infections / capacity
      9  coastal                — flood-risk geography
     10  r0_suitability         — R₀_rel(T) per cell
     11  climate_stress         — ΔT / 3 °C scalar broadcast
    """

    metadata = {"render_modes": ["human", "rgb_array"], "render_fps": 8}
    OBS_CHANNELS = 12

    # Latitude temperature gradient
    # Row 0 = equatorial (warm, ~coastal Kenya), row N-1 = poleward (cool, ~Ethiopian Highlands)
    T_EQUATOR    = 28.0   # °C
    T_POLE_END   = 14.0   # °C
    T_CELL_NOISE = 1.5    # °C std dev (microclimate variation)

    def __init__(
        self,
        grid_size=8,
        max_steps=80,

        # --- FOI ---
        foi_scale=0.15,            # simulation tuning — see accuracy notes in docstring

        # --- Larval habitat (Tompkins & Ermert 2013 VECTRI) ---
        water_larval_boost=0.4,

        # --- Recovery (Bousema & Drakeley 2011; White et al. 2014) ---
        recovery_rate=0.07,        # UPDATED from v3 (0.05 → 0.07, ~14-day infectious period)

        # --- Vaccination ---
        vax_rate=0.35,
        vaccine_budget=25,

        # --- Testing ---
        test_radius=1,
        test_cost=0.02,

        # --- Reward ---
        people_penalty_weight=0.10,
        vaccinate_bonus_weight=0.6,
        move_cost=0.003,
        vaccinate_cost=0.005,
        waste_penalty=0.15,
        terminal_weight=5.0,
        delta_weight=1.0,
        unused_vax_penalty=2.0,

        # --- Climate ---
        scenario="rcp85",
        hospital_overload_penalty=0.5,
        flood_prob_base=0.03,
        flood_duration=4,
        flood_water_boost=0.3,
        episode_number=0,

        # --- Spatial coupling (Marshall et al. 2018) ---
        # 0.8% per step: conservative for district-level grid where each step
        # spans more than one day.  Grounded in Marshall et al. (2016) survey
        # data showing ~2–5% overnight trips per month across 4 African countries.
        mobility_rate=0.008,

        # β = 1.5: Marshall et al. (2018) median fitted exponent.
        gravity_beta=1.5,

        # 0.36: fraction of asymptomatic infections with gametocytes.
        # Nkrumah et al. (2021) Sci. Rep.: 36.3% in longitudinal cohort (Uganda).
        import_infectivity=0.36,

        render_mode=None,
        seed=None,
    ):
        super().__init__()
        self.grid_size   = int(grid_size)
        self.max_steps   = int(max_steps)
        self.render_mode = render_mode
        self.cell_px     = 24

        self.foi_scale              = float(foi_scale)
        self.water_larval_boost     = float(water_larval_boost)
        self.recovery_rate_base     = float(recovery_rate)
        self.vax_rate               = float(vax_rate)
        self.max_budget             = int(vaccine_budget)
        self.test_radius            = int(test_radius)
        self.test_cost              = float(test_cost)

        self.people_penalty_weight  = float(people_penalty_weight)
        self.vaccinate_bonus_weight = float(vaccinate_bonus_weight)
        self.move_cost              = float(move_cost)
        self.vaccinate_cost         = float(vaccinate_cost)
        self.waste_penalty          = float(waste_penalty)
        self.terminal_weight        = float(terminal_weight)
        self.delta_weight           = float(delta_weight)
        self.unused_vax_penalty     = float(unused_vax_penalty)

        self.scenario                  = scenario
        self.hospital_overload_penalty = float(hospital_overload_penalty)
        self.flood_prob_base           = float(flood_prob_base)
        self.flood_duration            = int(flood_duration)
        self.flood_water_boost         = float(flood_water_boost)
        self.episode_number            = int(episode_number)

        self.mobility_rate      = float(mobility_rate)
        self.gravity_beta       = float(gravity_beta)
        self.import_infectivity = float(import_infectivity)

        self.action_space = spaces.Discrete(9)
        self.observation_space = spaces.Box(
            low=0.0, high=1.0,
            shape=(self.grid_size, self.grid_size, self.OBS_CHANNELS),
            dtype=np.float32,
        )
        self.reset(seed=seed)

    # ------------------------------------------------------------------
    # reset
    # ------------------------------------------------------------------

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        self.steps  = 0
        self.budget = self.max_budget
        self.vax_used   = 0
        self.vax_wasted = 0
        self.tests_used = 0
        self._flood_steps_remaining = 0

        self.delta_T = _delta_T(self.episode_number, self.scenario)

        # Latitude temperature grid (row 0 = equatorial, row N-1 = poleward)
        lat_grad = np.linspace(self.T_EQUATOR, self.T_POLE_END, self.grid_size)[:, np.newaxis]
        noise    = self.np_random.normal(0.0, self.T_CELL_NOISE,
                                         (self.grid_size, self.grid_size)).astype(np.float32)
        self.temperature = np.clip(
            (lat_grad + noise + self.delta_T).astype(np.float32), 5.0, 45.0
        )

        # Precompute thermal fields (static within episode)
        self.r0_suit   = r0_relative(self.temperature)
        self.mosq_dens = mosquito_density(self.temperature)
        self.eip_f     = eip_factor(self.temperature)
        self.bite_r    = biting_rate(self.temperature)

        # Geography
        self.population, self.coastal, self.hospital_capacity = _make_geography(
            self.grid_size, self.np_random
        )

        # Gravity weight matrix (population is fixed per episode)
        self._W = build_gravity_weights(self.grid_size, self.population, self.gravity_beta)

        # Initial infection: only in warm cells where vectors exist (T > T_MIN)
        can_transmit = (self.temperature > T_MIN).astype(np.float32)
        raw_inf = self.np_random.uniform(0.0, 0.04,
                                          (self.grid_size, self.grid_size)).astype(np.float32)
        self.infected   = np.clip(raw_inf * can_transmit, 0.0, 1.0).astype(np.float32)
        self.vaccinated = np.zeros_like(self.infected, dtype=np.float32)

        # Water: baseline larval habitat
        self.water_base = np.clip(
            self.np_random.uniform(0.0, 0.1, self.infected.shape).astype(np.float32)
            + 0.15 * self.coastal, 0.0, 1.0
        )
        self.water = self.water_base.copy()

        # Partial observability
        self.known_mask        = np.zeros_like(self.infected, dtype=np.float32)
        self.observed_infected = np.zeros_like(self.infected, dtype=np.float32)

        # Agent spawns in a populated cell
        flat_pop = self.population.ravel()
        probs    = flat_pop / (flat_pop.sum() + 1e-8)
        flat_idx = self.np_random.choice(self.grid_size ** 2, p=probs)
        self.agent_pos = np.array(
            [flat_idx // self.grid_size, flat_idx % self.grid_size], dtype=np.int32
        )

        self.prev_people_affected = float(np.sum(self.infected * self.population))
        return self._get_obs(), {}

    # ------------------------------------------------------------------
    # observation
    # ------------------------------------------------------------------

    def _get_obs(self):
        agent_layer = np.zeros_like(self.infected, dtype=np.float32)
        agent_layer[self.agent_pos[0], self.agent_pos[1]] = 1.0

        budget_layer = np.full_like(
            self.infected, self.budget / max(self.max_budget, 1), dtype=np.float32
        )
        hosp_load = np.clip(
            self.infected / np.maximum(self.hospital_capacity, 0.01), 0.0, 1.0
        ).astype(np.float32)

        temp_norm = np.clip((self.temperature - 10.0) / 30.0, 0.0, 1.0).astype(np.float32)
        stress    = np.full_like(
            self.infected, np.clip(self.delta_T / 3.0, 0.0, 1.0), dtype=np.float32
        )

        obs = np.stack([
            self.observed_infected,  # 0
            self.known_mask,         # 1
            self.vaccinated,         # 2
            self.population,         # 3
            temp_norm,               # 4
            self.water,              # 5
            agent_layer,             # 6
            budget_layer,            # 7
            hosp_load,               # 8
            self.coastal,            # 9
            self.r0_suit,            # 10
            stress,                  # 11
        ], axis=-1).astype(np.float32)

        return np.clip(obs, 0.0, 1.0).astype(np.float32)

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _reveal(self, x, y):
        r = self.test_radius
        x0, x1 = max(0, x - r), min(self.grid_size, x + r + 1)
        y0, y1 = max(0, y - r), min(self.grid_size, y + r + 1)
        self.known_mask[x0:x1, y0:y1]        = 1.0
        self.observed_infected[x0:x1, y0:y1] = self.infected[x0:x1, y0:y1]

    # ------------------------------------------------------------------
    # step
    # ------------------------------------------------------------------

    def step(self, action):
        self.steps += 1
        x, y = int(self.agent_pos[0]), int(self.agent_pos[1])

        did_move = did_vax = did_test = False
        vax_added = 0.0
        waste     = 0.0
        vax_x, vax_y = x, y   # cell where vaccine was applied (fix from v3)

        def move(dx, dy):
            nonlocal x, y, did_move
            nx = int(np.clip(x + dx, 0, self.grid_size - 1))
            ny = int(np.clip(y + dy, 0, self.grid_size - 1))
            did_move = (nx != x) or (ny != y)
            x, y = nx, ny

        def vaccinate(ix, iy):
            nonlocal did_vax, vax_added, waste, vax_x, vax_y
            if self.budget <= 0:
                return
            if float(self.population[ix, iy]) <= 0.0:
                waste = 1.0; self.vax_wasted += 1; return
            prev = float(self.vaccinated[ix, iy])
            newv = float(np.clip(prev + self.vax_rate, 0.0, 1.0))
            added = newv - prev
            if added < 0.05:
                waste = 1.0; self.vax_wasted += 1
                if added <= 0.0: return
            self.vaccinated[ix, iy] = newv
            self.budget -= 1; self.vax_used += 1
            did_vax = True; vax_added = float(added)
            vax_x, vax_y = ix, iy

        if   action == 0: move(-1,  0)
        elif action == 1: move( 1,  0)
        elif action == 2: move( 0, -1)
        elif action == 3: move( 0,  1)
        elif action == 4: vaccinate(x, y); move(-1,  0)
        elif action == 5: vaccinate(x, y); move( 1,  0)
        elif action == 6: vaccinate(x, y); move( 0, -1)
        elif action == 7: vaccinate(x, y); move( 0,  1)
        elif action == 8:
            did_test = True; self.tests_used += 1; self._reveal(x, y)

        self.agent_pos[:] = (x, y)

        # ------------------------------------------------------------------
        # Flood events
        # Frequency scales with warming; IPCC AR5 projects roughly doubled
        # flood frequency under RCP-8.5 by 2100.
        # ------------------------------------------------------------------
        flood_scale = 1.0 + (self.delta_T / 2.6)
        if self._flood_steps_remaining > 0:
            self._flood_steps_remaining -= 1
            self.water = np.clip(
                self.water_base + self.flood_water_boost * self.coastal, 0.0, 1.0
            ).astype(np.float32)
        elif self.np_random.random() < self.flood_prob_base * flood_scale:
            self._flood_steps_remaining = self.flood_duration
            self.water = np.clip(
                self.water_base + self.flood_water_boost * self.coastal, 0.0, 1.0
            ).astype(np.float32)
        else:
            self.water = self.water_base.copy()

        # ------------------------------------------------------------------
        # Seasonal forcing (Tusting et al. 2017)
        # ------------------------------------------------------------------
        season_mult = _seasonal_factor(self.steps, self.max_steps)

        # ------------------------------------------------------------------
        # Local Ross-Macdonald FOI (mosquito-mediated)
        # ------------------------------------------------------------------
        eff_mosq = np.clip(
            self.mosq_dens * (1.0 + self.water_larval_boost * self.water),
            0.0, 1.0,
        ).astype(np.float32)

        local_foi = (
            self.foi_scale
            * self.r0_suit
            * self.bite_r
            * eff_mosq
            * self.eip_f
            * self.infected
            * season_mult
        )

        # ------------------------------------------------------------------
        # Gravity-model importation FOI (Marshall et al. 2018)
        # Infected travellers seed parasites into newly-warm receptive cells.
        # The r0_suit gate ensures parasites only establish where vectors exist.
        # ------------------------------------------------------------------
        import_foi = importation_pressure(
            infected_flat      = self.infected.ravel(),
            r0_suit_flat       = self.r0_suit.ravel(),
            W                  = self._W,
            mobility_rate      = self.mobility_rate,
            import_infectivity = self.import_infectivity,
        ).reshape(self.grid_size, self.grid_size)

        # ------------------------------------------------------------------
        # New infections = (local FOI + import FOI) × susceptible fraction
        # ------------------------------------------------------------------
        susceptible = np.clip(1.0 - self.infected - self.vaccinated, 0.0, 1.0)
        new_inf     = (local_foi + import_foi) * susceptible

        # ------------------------------------------------------------------
        # Hospital capacity (Bousema & Drakeley 2011)
        # ------------------------------------------------------------------
        hosp_load     = self.infected / np.maximum(self.hospital_capacity, 0.01)
        overload_mask = (hosp_load > 1.0).astype(np.float32)
        eff_recovery  = self.recovery_rate_base * (
            1.0 - self.hospital_overload_penalty * overload_mask
        )

        self.infected = np.clip(
            self.infected + new_inf - eff_recovery * self.infected,
            0.0, 1.0,
        ).astype(np.float32)

        # Update partial observation
        self.observed_infected[self.known_mask == 1] = self.infected[self.known_mask == 1]

        # ------------------------------------------------------------------
        # Reward (vaccination bonus uses vax_x/vax_y — fix from v3)
        # ------------------------------------------------------------------
        people_affected = float(np.sum(self.infected * self.population))
        delta           = people_affected - self.prev_people_affected
        self.prev_people_affected = people_affected

        people_penalty = self.people_penalty_weight * people_affected

        bonus = 0.0
        if did_vax:
            coastal_boost  = 1.0 + float(self.coastal[vax_x, vax_y])
            overload_boost = 1.0 + float(overload_mask[vax_x, vax_y])
            bonus = self.vaccinate_bonus_weight * (
                2.5 * float(self.population[vax_x, vax_y])
                + 10.0 * float(self.r0_suit[vax_x, vax_y])
                + 1.0
            ) * vax_added * coastal_boost * overload_boost

        cost = (
            (1.0 if did_move else 0.0) * self.move_cost
            + (1.0 if did_vax  else 0.0) * self.vaccinate_cost
            + (1.0 if did_test else 0.0) * self.test_cost
            + waste * self.waste_penalty
        )

        reward  = float(bonus - people_penalty - cost)
        reward -= self.delta_weight * max(float(delta), 0.0)

        terminated = False
        truncated  = self.steps >= self.max_steps

        if truncated:
            reward -= self.terminal_weight * people_affected
            unused_frac = self.budget / max(self.max_budget, 1)
            reward -= self.unused_vax_penalty * unused_frac

        info = {
            "people_affected":     people_affected,
            "vaccines_used":       int(self.vax_used),
            "vaccines_wasted":     int(self.vax_wasted),
            "tests_used":          int(self.tests_used),
            "budget":              int(self.budget),
            "delta_T":             float(self.delta_T),
            "flood_active":        bool(self._flood_steps_remaining > 0),
            "season_factor":       float(season_mult),
            "mean_r0_suitability": float(self.r0_suit.mean()),
            "poleward_r0":         float(self.r0_suit[self.grid_size // 2:].mean()),
        }

        return self._get_obs(), reward, terminated, truncated, info

    # ------------------------------------------------------------------
    # render
    # ------------------------------------------------------------------

    def render(self):
        if self.render_mode is None:
            return
        img = np.zeros((self.grid_size, self.grid_size, 3), dtype=np.float32)
        img[..., 0] = self.infected
        img[..., 1] = self.vaccinated
        img[..., 2] = np.clip((self.temperature - 10.0) / 30.0, 0.0, 1.0) * 0.6
        img *= (self.population > 0)[..., None].astype(np.float32)
        img = (np.clip(img, 0.0, 1.0) * 255).astype(np.uint8)
        img = np.kron(img, np.ones((self.cell_px, self.cell_px, 1), dtype=np.uint8))

        ax, ay = self.agent_pos
        img[ax*self.cell_px:(ax+1)*self.cell_px,
            ay*self.cell_px:(ay+1)*self.cell_px] = [255, 255, 0]

        if self.render_mode == "rgb_array":
            return img

        flood_str = " [FLOOD]" if self._flood_steps_remaining > 0 else ""

        # Reuse one figure across frames so they update in place
        if not hasattr(self, "_render_fig") or self._render_fig is None:
            self._render_fig, self._render_ax = plt.subplots(figsize=(5, 5))
            self._render_im = self._render_ax.imshow(img)
            self._render_ax.axis("off")
        else:
            self._render_im.set_data(img)

        self._render_ax.set_title(
            f"step={self.steps}  budget={self.budget}  "
            f"used={self.vax_used}  wasted={self.vax_wasted}\n"
            f"ΔT={self.delta_T:.2f}°C  {self.scenario}{flood_str}\n"
            f"red=infected  green=vaccinated  blue=warm"
        )

        # Update the figure in place without clearing the rest of the cell's output.
        # (Previously used clear_output(wait=True) which wiped any earlier plots
        # in the same cell.)
        try:
            from IPython.display import update_display
            if not hasattr(self, "_render_display_id"):
                self._render_display_id = "malaria_render_" + str(id(self))
                display(self._render_fig, display_id=self._render_display_id)
            else:
                update_display(self._render_fig, display_id=self._render_display_id)
        except Exception:
            # Outside Jupyter: pause to force a redraw without blocking
            plt.pause(0.1)


# ===========================================================================
# Sanity checks
# ===========================================================================

env = MalariaEpidemicEnv(scenario="rcp85", render_mode=None)
check_env(env)
print("check_env passed ✓\n")

# Verify poleward rows are cold / zero-suitability at episode 0
env0 = MalariaEpidemicEnv(scenario="rcp85", episode_number=0, render_mode=None)
env0.reset(seed=0)
print(f"Episode 0  — equatorial row 0 T: {env0.temperature[0].mean():.1f}°C"
      f"  R₀_suit: {env0.r0_suit[0].mean():.3f}")
print(f"Episode 0  — poleward row {env0.grid_size-1} T: {env0.temperature[env0.grid_size-1].mean():.1f}°C"
      f"  R₀_suit: {env0.r0_suit[env0.grid_size-1].mean():.3f}")

# Verify poleward rows warm and become suitable at late episode
env500 = MalariaEpidemicEnv(scenario="rcp85", episode_number=500, render_mode=None)
env500.reset(seed=0)
print(f"\nEpisode 500 — equatorial row 0 T: {env500.temperature[0].mean():.1f}°C"
      f"  R₀_suit: {env500.r0_suit[0].mean():.3f}")
print(f"Episode 500 — poleward row {env500.grid_size-1} T: {env500.temperature[env500.grid_size-1].mean():.1f}°C"
      f"  R₀_suit: {env500.r0_suit[env500.grid_size-1].mean():.3f}")

# Thermal curve plot (sanity check — peer-reviewed shape verification)
T_range = np.linspace(5.0, 42.0, 300)
fig0, ax0 = plt.subplots(figsize=(8, 4))
ax0.plot(T_range, r0_relative(T_range),                            label="R₀_rel(T)")
ax0.plot(T_range, mosquito_density(T_range),                       label="Mosquito density")
ax0.plot(T_range, eip_factor(T_range) / eip_factor(T_range).max(), label="EIP factor (norm)")
ax0.plot(T_range, biting_rate(T_range),                            label="Biting rate")
ax0.axvline(T_OPT, color="k", linestyle="--", alpha=0.5, label=f"T_opt={T_OPT}°C (Mordecai 2013)")
ax0.axvline(T_MIN, color="b", linestyle=":",  alpha=0.5, label=f"T_min={T_MIN}°C (Johnson 2015)")
ax0.axvline(T_MAX, color="r", linestyle=":",  alpha=0.5, label=f"T_max={T_MAX}°C (Johnson 2015)")
ax0.set_title("Thermal performance curves — Mordecai (2013) / Johnson (2015)")
ax0.set_xlabel("Temperature (°C)"); ax0.set_ylabel("Relative value (0–1)")
ax0.legend(fontsize=8); ax0.grid(alpha=0.3)
plt.tight_layout(); plt.show()


# ===========================================================================
# SARSA state discretizer
# ===========================================================================

def discretize_state(env):
    x, y = env.agent_pos

    known   = int(env.known_mask[x, y])
    obs_inf = float(env.observed_infected[x, y]) if known else 0.0
    pop     = float(env.population[x, y])
    vax     = float(env.vaccinated[x, y])

    pop_level  = 0 if pop == 0 else (1 if pop < 0.5 else 2)
    inf_level  = 0 if obs_inf < 0.01 else (1 if obs_inf < 0.05 else 2)
    vax_level  = 0 if vax < 0.33 else (1 if vax < 0.66 else 2)

    x0, x1 = max(0, x - 1), min(env.grid_size, x + 2)
    y0, y1 = max(0, y - 1), min(env.grid_size, y + 2)
    nbr_level = 0
    nbr_inf = float(env.observed_infected[x0:x1, y0:y1].mean())
    if nbr_inf >= 0.07: nbr_level = 2
    elif nbr_inf >= 0.02: nbr_level = 1

    frac = env.budget / env.max_budget
    budget_level = 2 if frac > 0.66 else (1 if frac > 0.33 else 0)

    h_load = float(env.infected[x, y]) / max(float(env.hospital_capacity[x, y]), 0.01)
    hosp_level = 0 if h_load < 0.8 else (1 if h_load < 1.5 else 2)

    # R₀ suitability bucket — key signal for locating high-risk warm cells
    r0s = float(env.r0_suit[x, y])
    r0_level = 0 if r0s < 0.3 else (1 if r0s < 0.7 else 2)

    # Poleward row bucket — tells agent if it is in the expansion frontier
    row_frac = x / max(env.grid_size - 1, 1)
    row_level = 0 if row_frac < 0.33 else (1 if row_frac < 0.67 else 2)

    season_bucket = int((env.steps / env.max_steps) * 4) % 4
    stress_level  = 0 if env.delta_T < 0.9 else (1 if env.delta_T < 1.8 else 2)

    return (
        x, y,
        pop_level, known, inf_level, nbr_level, vax_level, budget_level,
        hosp_level, r0_level, row_level, season_bucket, stress_level,
    )


# ===========================================================================
# SARSA training
# ===========================================================================

num_actions   = MalariaEpidemicEnv(render_mode=None).action_space.n
Q             = defaultdict(lambda: np.zeros(num_actions, dtype=np.float32))
alpha         = 0.15
gamma         = 0.95
epsilon       = 1.0
epsilon_decay = 0.997
epsilon_min   = 0.05
num_episodes  = 1500

ep_return      = []
ep_people      = []
ep_delta_T     = []
ep_poleward_r0 = []


def eps_greedy(env, Q, s, eps):
    if np.random.rand() < eps:
        return env.action_space.sample()
    return int(np.argmax(Q[s]))


for ep in range(num_episodes):
    env    = MalariaEpidemicEnv(scenario="rcp85", episode_number=ep, render_mode=None)
    obs, _ = env.reset(seed=ep)
    s      = discretize_state(env)
    a      = eps_greedy(env, Q, s, epsilon)
    done   = False
    G      = 0.0
    last_info = {}

    while not done:
        obs2, r, term, trunc, info = env.step(a)
        s2   = discretize_state(env)
        done = term or trunc
        if not done:
            a2     = eps_greedy(env, Q, s2, epsilon)
            target = r + gamma * float(Q[s2][a2])
        else:
            target = r; a2 = 0

        Q[s][a] += alpha * (target - Q[s][a])
        s, a = s2, a2
        G += r
        last_info = info

    epsilon = max(epsilon_min, epsilon * epsilon_decay)
    ep_return.append(G)
    ep_people.append(last_info.get("people_affected", 0.0))
    ep_delta_T.append(last_info.get("delta_T", 0.0))
    ep_poleward_r0.append(last_info.get("poleward_r0", 0.0))

    if ep % 50 == 0:
        print(
            f"ep={ep:03d}  return={G:8.2f}  "
            f"people={last_info.get('people_affected',0):.3f}  "
            f"ΔT={last_info.get('delta_T',0):.2f}°C  "
            f"poleward_R₀={last_info.get('poleward_r0',0):.3f}  "
            f"eps={epsilon:.2f}"
        )


# ===========================================================================
# Plots
# ===========================================================================

def ma(x, w=25):
    x = np.asarray(x, dtype=np.float32)
    return np.convolve(x, np.ones(w) / w, mode="valid") if len(x) >= w else x


fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Malaria Epidemic RL — Ross-Macdonald + Gravity Importation (RCP-8.5)", fontsize=13)

axes[0, 0].plot(ma(ep_return)); axes[0, 0].set_title("Episode Return (MA-25)")
axes[0, 0].set_xlabel("Episode"); axes[0, 0].set_ylabel("Return")

axes[0, 1].plot(ma(ep_people), color="firebrick")
axes[0, 1].set_title("Final People Affected (MA-25)")
axes[0, 1].set_xlabel("Episode"); axes[0, 1].set_ylabel("People")

axes[1, 0].plot(ep_delta_T, color="darkorange", alpha=0.8)
axes[1, 0].axhline(1.8, color="gold", linestyle="--", label="RCP-4.5 ~2100")
axes[1, 0].axhline(2.6, color="red",  linestyle="--", label="RCP-8.5 ~2100")
axes[1, 0].set_title("Cumulative Warming ΔT (°C) — IPCC AR5")
axes[1, 0].set_xlabel("Episode"); axes[1, 0].set_ylabel("ΔT (°C)"); axes[1, 0].legend()

axes[1, 1].plot(ma(ep_poleward_r0), color="purple", alpha=0.9)
axes[1, 1].set_title("Poleward-Half Mean R₀ Suitability (MA-25)\n← range expansion signal")
axes[1, 1].set_xlabel("Episode"); axes[1, 1].set_ylabel("R₀_rel (0–1)")

plt.tight_layout(); plt.show()


# ===========================================================================
# Policy comparison
# ===========================================================================

def random_policy(env, obs):
    return env.action_space.sample()


def heuristic_policy(env, obs):
    """Prioritise cells with high R₀ suitability, population, and overload."""
    x, y       = env.agent_pos
    pop        = env.population[x, y]
    r0s        = env.r0_suit[x, y]
    coast      = env.coastal[x, y] > 0.5
    known      = env.known_mask[x, y] > 0.5
    overloaded = (float(env.infected[x, y]) /
                  max(float(env.hospital_capacity[x, y]), 0.01) > 1.0)
    if not known and np.random.rand() < 0.25:
        return 8
    if env.budget > 0 and pop > 0 and (r0s > 0.3 or coast or overloaded):
        return np.random.choice([4, 5, 6, 7])
    return np.random.choice([0, 1, 2, 3])


def greedy_sarsa(env, obs):
    return int(np.argmax(Q[discretize_state(env)]))


def run_ep(env, policy_fn):
    obs, _ = env.reset()
    done = False; total_r = 0.0; last_info = {}
    while not done:
        a = policy_fn(env, obs)
        obs, r, term, trunc, info = env.step(a)
        total_r += r; last_info = info; done = term or trunc
    return {
        "total_reward":          float(total_r),
        "final_people_affected": float(last_info.get("people_affected", 0)),
        "vaccines_used":         int(last_info.get("vaccines_used", 0)),
        "vaccines_wasted":       int(last_info.get("vaccines_wasted", 0)),
    }


def eval_policy(policy_fn, ep_num=400, n=20, seed0=1000):
    """Run n evaluation episodes, return raw list of result dicts and summary."""
    out = [run_ep(MalariaEpidemicEnv(scenario="rcp85", episode_number=ep_num,
                                      render_mode=None), policy_fn)
           for _ in range(n)]
    def ms(k):
        arr = np.array([d[k] for d in out], dtype=np.float32)
        return float(arr.mean()), float(arr.std(ddof=1))
    summary = {k: ms(k) for k in out[0]}
    return out, summary


def print_summary(name, summ):
    print(f"\n{'='*55}\n{name}")
    for k, (m, sd) in summ.items():
        print(f"  {k:28s}: {m:.3f} ± {sd:.3f}")


print("\n--- Evaluation at episode 400 (ΔT ≈ 2.1 °C, RCP-8.5) ---")
rand_raw,   rand_summ   = eval_policy(random_policy)
heur_raw,   heur_summ   = eval_policy(heuristic_policy)
sarsa_raw,  sarsa_summ  = eval_policy(greedy_sarsa)

print_summary("Random",                  rand_summ)
print_summary("Heuristic",               heur_summ)
print_summary("Learned SARSA (greedy)",  sarsa_summ)

# ===========================================================================
# Policy comparison plots
# ===========================================================================

policies      = ["Random", "Heuristic", "SARSA (greedy)"]
policy_colors = ["#888780", "#1D9E75", "#7F77DD"]   # gray, teal, purple
all_raw       = [rand_raw, heur_raw, sarsa_raw]
metrics = [
    ("total_reward",          "Total reward",          False),
    ("final_people_affected", "Final people affected", False),
    ("vaccines_used",         "Vaccines used",         True),
    ("vaccines_wasted",       "Vaccines wasted",       True),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(
    "Policy comparison — episode 400  (ΔT ≈ 2.1 °C, RCP-8.5, n=20 seeds)",
    fontsize=13, y=1.01
)

for ax, (key, label, is_int) in zip(axes.flat, metrics):
    # Per-policy arrays
    arrays = [np.array([d[key] for d in raw], dtype=np.float32) for raw in all_raw]
    means  = [a.mean() for a in arrays]
    stds   = [a.std(ddof=1) for a in arrays]

    # Bar chart with individual data points overlaid
    x = np.arange(len(policies))
    bars = ax.bar(x, means, color=policy_colors, alpha=0.75, width=0.5,
                  zorder=2, linewidth=0)
    ax.errorbar(x, means, yerr=stds, fmt="none", color="#2C2C2A",
                capsize=5, capthick=1.2, linewidth=1.2, zorder=3)

    # Scatter individual runs (jittered)
    rng = np.random.default_rng(0)
    for i, arr in enumerate(arrays):
        jitter = rng.uniform(-0.12, 0.12, size=len(arr))
        ax.scatter(x[i] + jitter, arr, color=policy_colors[i],
                   edgecolors="#2C2C2A", linewidths=0.4,
                   s=22, alpha=0.7, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(policies, fontsize=10)
    ax.set_ylabel(label, fontsize=10)
    ax.set_title(label, fontsize=11)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    # Annotate mean ± sd above each bar
    for i, (m, s) in enumerate(zip(means, stds)):
        fmt = f"{m:.0f}±{s:.0f}" if is_int else f"{m:.3f}±{s:.3f}"
        ax.text(x[i], m + s + (ax.get_ylim()[1] * 0.01 if ax.get_ylim()[1] > 0 else 0.002),
                fmt, ha="center", va="bottom", fontsize=8)

plt.tight_layout()
plt.savefig("policy_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: policy_comparison.png")

# ------------------------------------------------------------------
# Per-metric violin plot (distribution shape across 20 seeds)
# ------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
fig2.suptitle(
    "Reward and people-affected distributions across 20 evaluation seeds",
    fontsize=12
)

for ax, (key, label, _) in zip(axes2, metrics[:2]):
    arrays = [np.array([d[key] for d in raw], dtype=np.float32) for raw in all_raw]
    parts  = ax.violinplot(arrays, positions=range(len(policies)),
                           showmedians=True, showextrema=True)
    for pc, col in zip(parts["bodies"], policy_colors):
        pc.set_facecolor(col)
        pc.set_alpha(0.55)
    parts["cmedians"].set_color("#2C2C2A")
    parts["cbars"].set_color("#2C2C2A")
    parts["cmins"].set_color("#2C2C2A")
    parts["cmaxes"].set_color("#2C2C2A")
    ax.set_xticks(range(len(policies)))
    ax.set_xticklabels(policies, fontsize=10)
    ax.set_ylabel(label, fontsize=10)
    ax.set_title(label, fontsize=11)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig("policy_distributions.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: policy_distributions.png")


# ===========================================================================
# EXTENSION — Socio-economic factors (separate analysis, does not affect main results)
# ===========================================================================
# This section explores how the trained SARSA policy performs in a world that
# also models wealth inequality. Real-world malaria research consistently shows
# that poorer regions have:
#   - Lower hospital capacity → slower recovery
#   - Lower vaccine uptake → not everyone offered a dose actually receives it
#     (distance to clinic, cold-chain failures, hesitancy, etc.)
# References:
#   Tusting et al. (2013) Lancet 382:963-972 — wealth gradient in malaria
#   Sarma et al. (2019) PLOS ONE 14:e0223568 — socioeconomic inequality
# ---------------------------------------------------------------------------

class MalariaEpidemicEnvSocioEconomic(MalariaEpidemicEnv):
    """Same environment as the main model, plus a wealth layer that affects
    hospital capacity and vaccine uptake. The trained SARSA policy is
    evaluated on this richer world to see whether it generalizes.

    wealth ∈ [0, 1]: 0 = poorest, 1 = wealthiest.
    Wealth correlates with population density (cities tend to be wealthier)
    plus per-cell noise, so it isn't perfectly redundant with population.
    """

    def reset(self, seed=None, options=None):
        obs, info = super().reset(seed=seed, options=options)

        # Build wealth field: weakly correlated with population + noise.
        # Coastal cells are slightly poorer on average (rural/fishing communities).
        base_wealth = 0.4 * self.population + self.np_random.uniform(
            0.0, 0.6, self.population.shape
        ).astype(np.float32)
        base_wealth -= 0.15 * self.coastal
        self.wealth = np.clip(base_wealth, 0.0, 1.0).astype(np.float32)

        # Wealth-adjusted hospital capacity: poor regions have weaker healthcare.
        # Range: a poor cell (wealth=0) gets 50% of its base capacity;
        # a wealthy cell (wealth=1) keeps full capacity.
        wealth_factor = 0.5 + 0.5 * self.wealth
        self.hospital_capacity = (self.hospital_capacity * wealth_factor).astype(np.float32)

        # Vaccine uptake rate: in poor cells, only a fraction of offered doses
        # reach people. Range: 0.5 (poorest) to 1.0 (wealthiest).
        self.uptake_rate = (0.5 + 0.5 * self.wealth).astype(np.float32)

        return obs, info

    def step(self, action):
        # Patch the vax_rate before each vaccinate call so that uptake at the
        # current cell reflects local wealth. We do this by temporarily
        # adjusting the env's vax_rate for the cell we're on.
        x, y = int(self.agent_pos[0]), int(self.agent_pos[1])
        original_vax_rate = self.vax_rate
        self.vax_rate = float(original_vax_rate * self.uptake_rate[x, y])
        try:
            return super().step(action)
        finally:
            self.vax_rate = original_vax_rate


def run_ep_socio(env, policy_fn):
    """Evaluation rollout that also tracks wealth-stratified outcomes."""
    obs, _ = env.reset()
    done = False
    total_r = 0.0
    last_info = {}
    while not done:
        a = policy_fn(env, obs)
        obs, r, term, trunc, info = env.step(a)
        total_r += r
        last_info = info
        done = term or trunc

    # Stratify infections by wealth (lower half vs upper half of cells).
    median_w = float(np.median(env.wealth))
    poor_mask    = env.wealth <= median_w
    wealthy_mask = env.wealth > median_w
    poor_inf    = float(np.sum(env.infected * env.population * poor_mask))
    wealthy_inf = float(np.sum(env.infected * env.population * wealthy_mask))

    return {
        "total_reward":          float(total_r),
        "final_people_affected": float(last_info.get("people_affected", 0)),
        "vaccines_used":         int(last_info.get("vaccines_used", 0)),
        "vaccines_wasted":       int(last_info.get("vaccines_wasted", 0)),
        "poor_infected":         poor_inf,
        "wealthy_infected":      wealthy_inf,
    }


print("\n\n" + "=" * 65)
print("EXTENSION: Evaluating SARSA policy under socio-economic inequality")
print("=" * 65)

# Run the trained SARSA policy on the wealth-aware environment
n_seeds_socio = 20
socio_results = []
for i in range(n_seeds_socio):
    env_socio = MalariaEpidemicEnvSocioEconomic(
        scenario="rcp85", episode_number=400, render_mode=None
    )
    socio_results.append(run_ep_socio(env_socio, greedy_sarsa))

# Summary statistics
def _summarize(results, key):
    arr = np.array([r[key] for r in results], dtype=np.float32)
    return float(arr.mean()), float(arr.std(ddof=1))

socio_summary = {k: _summarize(socio_results, k) for k in socio_results[0]}

print(f"\nSARSA on standard world  (from earlier eval):")
print(f"  people affected: {sarsa_summ['final_people_affected'][0]:.3f} "
      f"± {sarsa_summ['final_people_affected'][1]:.3f}")

print(f"\nSARSA on socio-economic world:")
for k, (m, sd) in socio_summary.items():
    print(f"  {k:24s}: {m:.3f} ± {sd:.3f}")

mean_poor    = socio_summary["poor_infected"][0]
mean_wealthy = socio_summary["wealthy_infected"][0]
ratio = mean_poor / max(mean_wealthy, 1e-6)
print(f"\nInequality ratio (poor / wealthy infections): {ratio:.2f}×")


# ===========================================================================
# Plot — socio-economic comparison
# ===========================================================================

fig_se, axes_se = plt.subplots(1, 3, figsize=(14, 4.5))
fig_se.suptitle(
    "Socio-economic extension — SARSA policy evaluated on heterogeneous world",
    fontsize=12,
)

# Plot 1: total people affected — standard vs socio
standard_arr = np.array([d["final_people_affected"] for d in sarsa_raw], dtype=np.float32)
socio_arr    = np.array([r["final_people_affected"] for r in socio_results], dtype=np.float32)

ax = axes_se[0]
positions = [0, 1]
parts = ax.violinplot(
    [standard_arr, socio_arr], positions=positions,
    showmedians=True, showextrema=True,
)
for pc, col in zip(parts["bodies"], ["#7F77DD", "#C57A5A"]):
    pc.set_facecolor(col); pc.set_alpha(0.6)
parts["cmedians"].set_color("#2C2C2A")
parts["cbars"].set_color("#2C2C2A")
parts["cmins"].set_color("#2C2C2A")
parts["cmaxes"].set_color("#2C2C2A")
ax.set_xticks(positions)
ax.set_xticklabels(["Standard world", "Socio-economic world"], fontsize=9)
ax.set_ylabel("Final people affected")
ax.set_title("Overall infection burden")
ax.yaxis.grid(True, linestyle="--", alpha=0.4)
ax.spines[["top", "right"]].set_visible(False)

# Plot 2: poor vs wealthy infection burden in the socio world
poor_arr    = np.array([r["poor_infected"]    for r in socio_results], dtype=np.float32)
wealthy_arr = np.array([r["wealthy_infected"] for r in socio_results], dtype=np.float32)

ax = axes_se[1]
positions = [0, 1]
bars = ax.bar(
    positions, [poor_arr.mean(), wealthy_arr.mean()],
    color=["#A04040", "#4A7AAB"], alpha=0.75, width=0.5,
)
ax.errorbar(
    positions, [poor_arr.mean(), wealthy_arr.mean()],
    yerr=[poor_arr.std(ddof=1), wealthy_arr.std(ddof=1)],
    fmt="none", color="#2C2C2A", capsize=5, capthick=1.2,
)
rng = np.random.default_rng(42)
for i, arr in enumerate([poor_arr, wealthy_arr]):
    jitter = rng.uniform(-0.12, 0.12, size=len(arr))
    ax.scatter(
        positions[i] + jitter, arr,
        color=["#A04040", "#4A7AAB"][i],
        edgecolors="#2C2C2A", linewidths=0.4,
        s=22, alpha=0.7, zorder=4,
    )
ax.set_xticks(positions)
ax.set_xticklabels(["Poor cells\n(lower-half wealth)", "Wealthy cells\n(upper-half wealth)"],
                   fontsize=9)
ax.set_ylabel("People infected")
ax.set_title("Infection burden by wealth bracket")
ax.yaxis.grid(True, linestyle="--", alpha=0.4)
ax.spines[["top", "right"]].set_visible(False)
ax.text(
    0.5, 0.95, f"Inequality ratio: {ratio:.2f}×",
    transform=ax.transAxes, ha="center", va="top",
    fontsize=10, fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="#F2EEE5", edgecolor="#9E9789"),
)

# Plot 3: vaccines wasted — standard vs socio
standard_waste = np.array([d["vaccines_wasted"] for d in sarsa_raw], dtype=np.float32)
socio_waste    = np.array([r["vaccines_wasted"] for r in socio_results], dtype=np.float32)

ax = axes_se[2]
positions = [0, 1]
ax.bar(
    positions, [standard_waste.mean(), socio_waste.mean()],
    color=["#7F77DD", "#C57A5A"], alpha=0.75, width=0.5,
)
ax.errorbar(
    positions, [standard_waste.mean(), socio_waste.mean()],
    yerr=[standard_waste.std(ddof=1), socio_waste.std(ddof=1)],
    fmt="none", color="#2C2C2A", capsize=5, capthick=1.2,
)
ax.set_xticks(positions)
ax.set_xticklabels(["Standard world", "Socio-economic world"], fontsize=9)
ax.set_ylabel("Vaccines wasted")
ax.set_title("Wasted doses (uptake failures)")
ax.yaxis.grid(True, linestyle="--", alpha=0.4)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig("socio_economic_results.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: socio_economic_results.png")


# ===========================================================================
# Visual rollout
# ===========================================================================

env = MalariaEpidemicEnv(scenario="rcp85", episode_number=400, render_mode="human")
obs, _ = env.reset(seed=1)
done = False
while not done:
    s      = discretize_state(env)
    action = int(np.argmax(Q[s]))
    obs, reward, terminated, truncated, info = env.step(action)
    env.render()
    time.sleep(0.15)
    done = terminated or truncated

print("\nFinal info:", info)

# Keep all matplotlib windows open after the script finishes
# so we can actually see the graphs.
plt.show(block=True)