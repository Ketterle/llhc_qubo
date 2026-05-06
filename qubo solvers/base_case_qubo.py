import json
import time
from itertools import combinations, product
from math import log2
from pathlib import Path

import dimod
import numpy as np


# ============================================================
# Basic utilities
# ============================================================

def all_binary_words(length):
    return ["".join(bits) for bits in product("01", repeat=length)]


def build_dictionary_by_length(L, forbidden_words=None, forbidden_prefixes=None):
    """
    Full binary dictionary up to length L, then remove forbidden words/prefixes.

    Base case:
      forbidden_words = set()
      forbidden_prefixes = tuple()
      -> no restrictions.
    """
    forbidden_words = set(forbidden_words or [])
    forbidden_prefixes = tuple(forbidden_prefixes or [])

    words_by_length = {}
    for length in range(1, L + 1):
        words = all_binary_words(length)
        words_by_length[length] = [
            word for word in words
            if word not in forbidden_words
            and not word.startswith(forbidden_prefixes)
        ]

    return words_by_length


def flatten_words(words_by_length):
    return [
        word
        for length in range(1, max(words_by_length) + 1)
        for word in words_by_length[length]
    ]


def prefix_conflict(u, v):
    return u != v and (u.startswith(v) or v.startswith(u))


def build_prefix_edges(words):
    return [(u, v) for u, v in combinations(words, 2) if prefix_conflict(u, v)]


def harmonic_probabilities(n):
    raw = [1.0 / (i + 1) for i in range(n)]
    total = sum(raw)
    return [x / total for x in raw]


def entropy(probabilities):
    return -sum(p * log2(p) for p in probabilities if p > 0)


def save_json(path, payload):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def instance_metadata(case_label, n, L, probabilities, forbidden_words, forbidden_prefixes, words_by_length):
    admissible_words = flatten_words(words_by_length)
    full_dictionary_size = sum(2 ** length for length in range(1, L + 1))
    edges = build_prefix_edges(admissible_words)
    possible_pairs = len(admissible_words) * (len(admissible_words) - 1) / 2

    return {
        "case_label": case_label,
        "n": n,
        "L": L,
        "probability_model": "harmonic_normalized",
        "entropy": entropy(probabilities),
        "forbidden_words": sorted(forbidden_words),
        "forbidden_prefixes": sorted(forbidden_prefixes),
        "dictionary_size_full": full_dictionary_size,
        "dictionary_size_admissible": len(admissible_words),
        "restriction_severity": 1.0 - len(admissible_words) / full_dictionary_size,
        "prefix_conflicts": len(edges),
        "prefix_conflict_density": 0.0 if possible_pairs == 0 else len(edges) / possible_pairs,
        "binary_variables_length": n * L,
        "binary_variables_selection": len(admissible_words),
        "binary_variables_total": n * L + len(admissible_words),
    }


# ============================================================
# Base Case Configuration
# ============================================================

CASE_LABEL = "base"
N = 8
L = 4

FORBIDDEN_WORDS = set()
FORBIDDEN_PREFIXES = tuple()

LAMBDA_A = 3.0
LAMBDA_C = 3.0
LAMBDA_P = 5.0
LAMBDA_O = 3.0

ROUNDS = 10
TIME_LIMIT_SECONDS = 3
MAX_DWAVE_SECONDS = 30

OUTPUT_RESULTS_JSON = "base_case_qubo_results.json"
OUTPUT_SUMMARY_JSON = "base_case_qubo_summary.json"


# ============================================================
# Variable naming
# ============================================================

def x_name(i, length):
    return f"x_{i}_{length}"


def y_name(word):
    return f"y_{word}"


# ============================================================
# Backend
# ============================================================

class LeapHybridBackend:
    name = "leap_hybrid"

    def __init__(self, token=None, time_limit=None):
        from dwave.system import LeapHybridSampler

        kwargs = {}
        if token is not None:
            kwargs["token"] = token

        self.sampler = LeapHybridSampler(**kwargs)
        self.time_limit = time_limit

    def solve(self, bqm, num_reads):
        sample_kwargs = {}
        if self.time_limit is not None:
            sample_kwargs["time_limit"] = self.time_limit

        results = []
        for _ in range(num_reads):
            sampleset = self.sampler.sample(bqm, **sample_kwargs)
            for rec in sampleset.data(["sample", "energy"]):
                results.append({"sample": dict(rec.sample), "energy": float(rec.energy)})

        return results

    def describe_params(self):
        return {
            "backend": self.name,
            "solver": self.sampler.solver.name,
            "time_limit_s": self.time_limit if self.time_limit is not None else "auto_minimum",
        }


# ============================================================
# Prefix code reconstruction
# ============================================================

def build_symbol_to_code_map(sample, words_by_length, n, L):
    """
    Reconstruct a deterministic symbol-to-codeword assignment.
    """
    lengths_per_symbol = {}
    for i in range(1, n + 1):
        active = [length for length in range(1, L + 1) if sample.get(x_name(i, length), 0) == 1]
        lengths_per_symbol[i] = active[0] if len(active) == 1 else None

    selected_by_length = {}
    for length in range(1, L + 1):
        selected_by_length[length] = sorted(
            word for word in words_by_length[length]
            if sample.get(y_name(word), 0) == 1
        )

    counters = {length: 0 for length in range(1, L + 1)}
    assignment = {}

    for i in range(1, n + 1):
        length = lengths_per_symbol[i]
        if length is None:
            assignment[str(i)] = None
            continue

        idx = counters[length]
        if idx < len(selected_by_length[length]):
            assignment[str(i)] = selected_by_length[length][idx]
            counters[length] += 1
        else:
            assignment[str(i)] = None

    return assignment


# ============================================================
# Sample evaluation
# ============================================================

def evaluate_sample(sample, words_by_length, prefix_edges, n, L, probabilities):
    def get_x(i, length):
        return sample.get(x_name(i, length), 0)

    def get_y(word):
        return sample.get(y_name(word), 0)

    H_cost = 0.0
    for i in range(1, n + 1):
        for length in range(1, L + 1):
            H_cost += probabilities[i - 1] * length * get_x(i, length)

    H_ass = 0.0
    for i in range(1, n + 1):
        s = sum(get_x(i, length) for length in range(1, L + 1))
        H_ass += (1 - s) ** 2

    H_coupl = 0.0
    for length in range(1, L + 1):
        X_l = sum(get_x(i, length) for i in range(1, n + 1))
        Y_l = sum(get_y(word) for word in words_by_length[length])
        H_coupl += (X_l - Y_l) ** 2

    H_pref = 0.0
    active_prefix_conflicts = []
    for u, v in prefix_edges:
        value = get_y(u) * get_y(v)
        H_pref += value
        if value == 1:
            active_prefix_conflicts.append((u, v))

    H_ord = 0.0
    for i in range(1, n):
        for j in range(i + 1, n + 1):
            for length_i in range(1, L + 1):
                for length_j in range(1, length_i):
                    H_ord += get_x(i, length_i) * get_x(j, length_j)

    feasible = bool(H_ass == 0 and H_coupl == 0 and H_pref == 0 and H_ord == 0)

    return {
        "H_cost": H_cost,
        "H_ass": H_ass,
        "H_coupl": H_coupl,
        "H_pref": H_pref,
        "H_ord": H_ord,
        "feasible": feasible,
        "active_prefix_conflicts": active_prefix_conflicts,
    }


# ============================================================
# QUBO construction
# ============================================================

def build_bqm():
    # ------------------------------------------------------------------
    # STEP 2: Build the admissible dictionary and the conflict graph
    # ------------------------------------------------------------------
    probabilities = harmonic_probabilities(N)
    words_by_length = build_dictionary_by_length(L, FORBIDDEN_WORDS, FORBIDDEN_PREFIXES)
    admissible_words = flatten_words(words_by_length)
    prefix_edges = build_prefix_edges(admissible_words)

    var_order = (
        [x_name(i, length) for i in range(1, N + 1) for length in range(1, L + 1)]
        + [y_name(word) for word in admissible_words]
    )
    var_index = {v: idx for idx, v in enumerate(var_order)}

    # ------------------------------------------------------------------
    # STEP 3: Generate the QUBO formulation
    # ------------------------------------------------------------------
    Q = {}

    def add_Q(u, v, coeff):
        key = tuple(sorted((u, v)))
        Q[key] = Q.get(key, 0.0) + coeff

    for i in range(1, N + 1):
        for length in range(1, L + 1):
            add_Q(x_name(i, length), x_name(i, length), probabilities[i - 1] * length)

    offset = 0.0

    for i in range(1, N + 1):
        offset += LAMBDA_A
        for length in range(1, L + 1):
            add_Q(x_name(i, length), x_name(i, length), -LAMBDA_A)
        for l1, l2 in combinations(range(1, L + 1), 2):
            add_Q(x_name(i, l1), x_name(i, l2), 2.0 * LAMBDA_A)

    for length in range(1, L + 1):
        words_l = words_by_length[length]

        for i in range(1, N + 1):
            add_Q(x_name(i, length), x_name(i, length), LAMBDA_C)

        for i, j in combinations(range(1, N + 1), 2):
            add_Q(x_name(i, length), x_name(j, length), 2.0 * LAMBDA_C)

        for word in words_l:
            add_Q(y_name(word), y_name(word), LAMBDA_C)

        for u, v in combinations(words_l, 2):
            add_Q(y_name(u), y_name(v), 2.0 * LAMBDA_C)

        for i in range(1, N + 1):
            for word in words_l:
                add_Q(x_name(i, length), y_name(word), -2.0 * LAMBDA_C)

    for u, v in prefix_edges:
        add_Q(y_name(u), y_name(v), LAMBDA_P)

    for i in range(1, N):
        for j in range(i + 1, N + 1):
            for length_i in range(1, L + 1):
                for length_j in range(1, length_i):
                    add_Q(x_name(i, length_i), x_name(j, length_j), LAMBDA_O)

    Qmat = np.zeros((len(var_order), len(var_order)), dtype=float)
    for (u, v), coeff in Q.items():
        i = var_index[u]
        j = var_index[v]
        if i == j:
            Qmat[i, j] += coeff
        else:
            Qmat[i, j] += coeff / 2.0
            Qmat[j, i] += coeff / 2.0

    bqm = dimod.BinaryQuadraticModel.from_qubo(Q, offset=offset)
    metadata = instance_metadata(
        CASE_LABEL, N, L, probabilities,
        FORBIDDEN_WORDS, FORBIDDEN_PREFIXES, words_by_length,
    )

    return bqm, Qmat, var_order, words_by_length, prefix_edges, probabilities, metadata


# ============================================================
# Solving helpers
# ============================================================

def better(candidate, incumbent):
    if incumbent is None:
        return True

    c = candidate["best_feasible_energy"]
    b = incumbent["best_feasible_energy"]

    if c is not None and b is None:
        return True
    if c is not None and b is not None:
        return c < b
    if c is None and b is None:
        return candidate["energy"] < incumbent["energy"]

    return False


def solve_one_round(round_idx, backend, bqm, Qmat, var_order, words_by_length, prefix_edges, probabilities):
    # ------------------------------------------------------------------
    # STEP 4: Solve the QUBO problem
    # ------------------------------------------------------------------
    start = time.time()
    results = backend.solve(bqm=bqm, num_reads=1)
    elapsed = time.time() - start

    if not results:
        raise RuntimeError("No sample returned by backend.")

    best_entry = min(results, key=lambda r: r["energy"])
    sample = dict(best_entry["sample"])
    energy = float(best_entry["energy"])

    # ------------------------------------------------------------------
    # STEP 5: Post-process the solution
    # ------------------------------------------------------------------
    selected_x = sorted(v for v, val in sample.items() if v.startswith("x_") and val == 1)
    selected_y = sorted(v for v, val in sample.items() if v.startswith("y_") and val == 1)

    z = np.array([sample[v] for v in var_order], dtype=float)
    energy_matrix = float(z @ Qmat @ z + bqm.offset)

    symbol_to_code = build_symbol_to_code_map(sample, words_by_length, N, L)

    # ------------------------------------------------------------------
    # STEP 6: Verify feasibility
    # ------------------------------------------------------------------
    metrics = evaluate_sample(sample, words_by_length, prefix_edges, N, L, probabilities)

    best_feasible_payload = None
    if metrics["feasible"]:
        best_feasible_payload = {
            "energy": energy,
            "symbol_to_code": symbol_to_code,
            "metrics": metrics,
        }

    return {
        "round": round_idx,
        "seed": None,
        "backend_name": "leap_hybrid",
        "case_label": CASE_LABEL,
        "energy": energy,
        "energy_matrix": energy_matrix,
        "feasible": metrics["feasible"],
        "H_cost": metrics["H_cost"],
        "H_ass": metrics["H_ass"],
        "H_coupl": metrics["H_coupl"],
        "H_pref": metrics["H_pref"],
        "H_ord": metrics["H_ord"],
        "elapsed_seconds": elapsed,
        "selected_x": selected_x,
        "selected_y": selected_y,
        "symbol_to_code": symbol_to_code,
        "assigned_lengths": {
            str(i): [length for length in range(1, L + 1) if sample.get(x_name(i, length), 0) == 1]
            for i in range(1, N + 1)
        },
        "active_words_by_length": {
            str(length): [
                word for word in words_by_length[length]
                if sample.get(y_name(word), 0) == 1
            ]
            for length in range(1, L + 1)
        },
        "best_feasible_energy": energy if metrics["feasible"] else None,
        "best_feasible_payload": best_feasible_payload,
    }


def build_summary(payload):
    rounds = payload["all_rounds"]
    feasible = [r for r in rounds if r["best_feasible_energy"] is not None]
    feasible_energies = np.array([r["best_feasible_energy"] for r in feasible])
    hcosts = np.array([r["H_cost"] for r in rounds])
    best = payload["best_overall"]

    best_so_far = []
    current = float("inf")
    for r in feasible:
        current = min(current, r["best_feasible_energy"])
        best_so_far.append(current)

    checkpoints = [1, 5, 10, len(best_so_far)]

    return {
        "experiment": f"qubo_{CASE_LABEL}",
        "case_label": CASE_LABEL,
        "config": payload["parameters"],
        "rounds": len(rounds),
        "planned_dwave_seconds": payload["planned_dwave_seconds"],
        "total_elapsed_s": payload["total_elapsed_seconds"],
        "feasibility_rate": len(feasible) / len(rounds) if rounds else 0.0,
        "feasible_energy_stats": {
            "min": float(feasible_energies.min()) if len(feasible_energies) else None,
            "max": float(feasible_energies.max()) if len(feasible_energies) else None,
            "mean": float(feasible_energies.mean()) if len(feasible_energies) else None,
            "std": float(feasible_energies.std()) if len(feasible_energies) else None,
            "median": float(np.median(feasible_energies)) if len(feasible_energies) else None,
        },
        "hcost_stats": {
            "min": float(hcosts.min()) if len(hcosts) else None,
            "mean": float(hcosts.mean()) if len(hcosts) else None,
            "std": float(hcosts.std()) if len(hcosts) else None,
            "median": float(np.median(hcosts)) if len(hcosts) else None,
        },
        "convergence": {
            str(cp): float(best_so_far[cp - 1])
            for cp in checkpoints
            if 0 < cp <= len(best_so_far)
        },
        "best_solution": {
            "round": best["round"],
            "energy": best["best_feasible_energy"],
            "H_cost": best["H_cost"],
            "selected_y": best["selected_y"],
            "code": best["symbol_to_code"],
        },
        "instance": payload["instance"],
    }


# ============================================================
# Main experiment
# ============================================================

def run_qubo_experiment():
    # ------------------------------------------------------------------
    # STEP 1: Define the instance and execution budget
    # ------------------------------------------------------------------
    probabilities = harmonic_probabilities(N)
    planned_dwave_seconds = ROUNDS * TIME_LIMIT_SECONDS

    if planned_dwave_seconds > MAX_DWAVE_SECONDS:
        raise ValueError(
            f"Planned D-Wave time is {planned_dwave_seconds}s, "
            f"but MAX_DWAVE_SECONDS is {MAX_DWAVE_SECONDS}s."
        )

    print("=" * 80)
    print(f"QUBO SOLVER — BASE CASE (n = {N}, L = {L})")
    print("=" * 80)
    print(f"  Case label:          {CASE_LABEL}")
    print(f"  n:                   {N}")
    print(f"  L:                   {L}")
    print(f"  Forbidden words:     {len(FORBIDDEN_WORDS)}")
    print(f"  Forbidden prefixes:  {FORBIDDEN_PREFIXES if FORBIDDEN_PREFIXES else 'None'}")
    print(f"  Rounds:              {ROUNDS}")
    print(f"  Time limit / round:  {TIME_LIMIT_SECONDS}s")
    print(f"  Planned D-Wave time: {planned_dwave_seconds}s")
    print(f"  lambda_a:            {LAMBDA_A}")
    print(f"  lambda_c:            {LAMBDA_C}")
    print(f"  lambda_p:            {LAMBDA_P}")
    print(f"  lambda_o:            {LAMBDA_O}")

    # ------------------------------------------------------------------
    # STEP 2 + STEP 3: Build BQM (done in build_bqm())
    # ------------------------------------------------------------------
    (
        bqm,
        Qmat,
        var_order,
        words_by_length,
        prefix_edges,
        probabilities,
        metadata,
    ) = build_bqm()

    backend = LeapHybridBackend(time_limit=TIME_LIMIT_SECONDS)

    print(f"\n  Admissible words:    {metadata['dictionary_size_admissible']}")
    print(f"  Prefix conflicts:    {metadata['prefix_conflicts']}")
    print(f"  Binary variables:    {len(bqm.variables)}")
    print(f"  Solver:              {backend.describe_params()['solver']}")

    # ------------------------------------------------------------------
    # STEP 4 to STEP 6: Execute rounds (done in solve_one_round())
    # ------------------------------------------------------------------
    global_start = time.time()
    all_rounds = []
    best_overall = None

    print("\n=== START QUBO EXPERIMENT ===")

    for round_idx in range(1, ROUNDS + 1):
        result = solve_one_round(
            round_idx=round_idx,
            backend=backend,
            bqm=bqm,
            Qmat=Qmat,
            var_order=var_order,
            words_by_length=words_by_length,
            prefix_edges=prefix_edges,
            probabilities=probabilities,
        )

        all_rounds.append(result)

        if better(result, best_overall):
            best_overall = result

        print(
            f"[{round_idx:03d}/{ROUNDS}] "
            f"energy={result['energy']:.6f} "
            f"H_cost={result['H_cost']:.6f} "
            f"feasible={result['feasible']} "
            f"elapsed={result['elapsed_seconds']:.1f}s"
        )

    total_elapsed = time.time() - global_start

    # ------------------------------------------------------------------
    # STEP 7: Evaluate and report results
    # ------------------------------------------------------------------
    payload = {
        "script_type": "plain_qubo_fixed",
        "case_label": CASE_LABEL,
        "rounds": ROUNDS,
        "planned_dwave_seconds": planned_dwave_seconds,
        "parameters": {
            "backend": "leap_hybrid",
            "n": N,
            "L": L,
            "forbidden_words": sorted(FORBIDDEN_WORDS),
            "forbidden_prefixes": sorted(FORBIDDEN_PREFIXES),
            "lambda_a": LAMBDA_A,
            "lambda_c": LAMBDA_C,
            "lambda_p": LAMBDA_P,
            "lambda_o": LAMBDA_O,
            "time_limit_seconds": TIME_LIMIT_SECONDS,
            "rounds": ROUNDS,
            "probabilities": probabilities,
            **backend.describe_params(),
        },
        "instance": metadata,
        "total_elapsed_seconds": total_elapsed,
        "best_overall": best_overall,
        "all_rounds": all_rounds,
    }

    summary = build_summary(payload)

    save_json(OUTPUT_RESULTS_JSON, payload)
    save_json(OUTPUT_SUMMARY_JSON, summary)

    print("\n" + "=" * 80)
    print(f"QUBO EXPERIMENT FINISHED — {CASE_LABEL}")
    print("=" * 80)
    print(f"  Total elapsed time: {total_elapsed:.3f}s")
    print(f"  Feasibility rate:   {summary['feasibility_rate']:.3f}")
    print(f"  Best H_cost:        {summary['best_solution']['H_cost']}")
    print(f"  Results:            {OUTPUT_RESULTS_JSON}")
    print(f"  Summary:            {OUTPUT_SUMMARY_JSON}")


if __name__ == "__main__":
    run_qubo_experiment()