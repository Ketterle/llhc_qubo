import json
import time
from itertools import combinations, product
from math import log2, gcd
from functools import reduce
from pathlib import Path

import numpy as np

from ortools.sat.python import cp_model


# ============================================================
# Basic utilities
# ============================================================

def all_binary_words(length):
    return ["".join(bits) for bits in product("01", repeat=length)]


def build_dictionary_by_length(L, forbidden_words=None, forbidden_prefixes=None):
    """
    Full binary dictionary up to length L, then remove forbidden words/prefixes.
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
# Intermediate Case Configuration
# ============================================================

CASE_LABEL = "intermediate_classical"
N = 8
L = 4

FORBIDDEN_WORDS = set()
FORBIDDEN_PREFIXES = ("11",)

OUTPUT_RESULTS_JSON = "intermediate_cpsat_results.json"
OUTPUT_SUMMARY_JSON = "intermediate_cpsat_summary.json"


# ============================================================
# Helper functions for CP-SAT
# ============================================================

def scale_probabilities_to_ints(probabilities, scale=100000):
    """Convert real-valued probabilities into integer coefficients."""
    raw_weights = [int(round(scale * p)) for p in probabilities]

    if all(weight == 0 for weight in raw_weights):
        raise ValueError("All scaled probabilities are zero. Increase scale.")

    common_divisor = reduce(gcd, raw_weights)
    return [weight // common_divisor for weight in raw_weights]


def reconstruct_assignment(lengths, selected_words):
    """Reconstruct a deterministic symbol-to-codeword assignment."""
    words_by_length = {}

    for word in selected_words:
        words_by_length.setdefault(len(word), []).append(word)

    for length in words_by_length:
        words_by_length[length].sort()

    assignment = {}

    for symbol_index, length in enumerate(lengths, start=1):
        if length not in words_by_length or not words_by_length[length]:
            raise RuntimeError(
                f"No selected word of length {length} for symbol {symbol_index}."
            )

        assignment[str(symbol_index)] = words_by_length[length].pop(0)

    return assignment


def evaluate_solution(lengths, selected_words, probabilities):
    """Evaluate objective value and feasibility of the selected code."""
    expected_length = sum(
        probabilities[i] * lengths[i] for i in range(len(lengths))
    )

    active_conflicts = [
        (u, v)
        for u, v in combinations(selected_words, 2)
        if prefix_conflict(u, v)
    ]

    order_feasible = all(
        lengths[i] <= lengths[i + 1]
        for i in range(len(lengths) - 1)
    )

    return {
        "expected_length": expected_length,
        "prefix_feasible": len(active_conflicts) == 0,
        "active_prefix_conflicts": active_conflicts,
        "order_feasible": order_feasible,
        "feasible": len(active_conflicts) == 0 and order_feasible,
    }


# ============================================================
# Main experiment
# ============================================================

def run_classical_experiment():
    # ------------------------------------------------------------------
    # STEP 1: Define the instance and execution budget
    # ------------------------------------------------------------------
    probabilities = harmonic_probabilities(N)

    print("=" * 80)
    print(f"CLASSICAL SOLVER (CP-SAT) — INTERMEDIATE CASE (n = {N}, L = {L})")
    print("=" * 80)
    print(f"  Case label:          {CASE_LABEL}")
    print(f"  n:                   {N}")
    print(f"  L:                   {L}")
    print(f"  Forbidden words:     {len(FORBIDDEN_WORDS)}")
    print(f"  Forbidden prefixes:  {FORBIDDEN_PREFIXES}")

    # ------------------------------------------------------------------
    # STEP 2: Build the admissible dictionary and the conflict graph
    # ------------------------------------------------------------------
    words_by_length = build_dictionary_by_length(L, FORBIDDEN_WORDS, FORBIDDEN_PREFIXES)
    admissible_words = flatten_words(words_by_length)
    prefix_edges = build_prefix_edges(admissible_words)

    if len(admissible_words) < N:
        raise ValueError(
            f"Only {len(admissible_words)} admissible words available for {N} symbols."
        )

    print(f"  Admissible words:    {len(admissible_words)}")
    print(f"  Prefix conflicts:    {len(prefix_edges)}")

    # ------------------------------------------------------------------
    # STEP 3: Generate the CP-SAT formulation
    # ------------------------------------------------------------------
    model = cp_model.CpModel()

    x = {
        (i, length): model.NewBoolVar(f"x_{i}_{length}")
        for i in range(1, N + 1)
        for length in range(1, L + 1)
    }

    y = {
        word: model.NewBoolVar(f"y_{word}")
        for word in admissible_words
    }

    length_var = {
        i: model.NewIntVar(1, L, f"len_{i}")
        for i in range(1, N + 1)
    }

    # Assignment constraint: exactly one length per symbol
    for i in range(1, N + 1):
        model.AddExactlyOne(x[i, length] for length in range(1, L + 1))
        model.Add(
            length_var[i]
            == sum(length * x[i, length] for length in range(1, L + 1))
        )

    # Coupling constraint: number of symbols per length = selected words
    for length in range(1, L + 1):
        model.Add(
            sum(x[i, length] for i in range(1, N + 1))
            == sum(y[word] for word in words_by_length[length])
        )

    # Prefix-free constraint
    for u, v in prefix_edges:
        model.Add(y[u] + y[v] <= 1)

    # Ordering constraint: probability-length monotonicity
    for i in range(1, N):
        for j in range(i + 1, N + 1):
            model.Add(length_var[i] <= length_var[j])

    # Objective: minimize expected codeword length
    integer_weights = scale_probabilities_to_ints(probabilities, scale=100000)

    model.Minimize(
        sum(integer_weights[i - 1] * length_var[i] for i in range(1, N + 1))
    )

    print(f"  Variables:           {len(x) + len(y) + len(length_var)}")

    # ------------------------------------------------------------------
    # STEP 4: Solve the problem
    # ------------------------------------------------------------------
    print("\n=== START EXPERIMENT ===")
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 60.0
    solver.parameters.num_search_workers = 8
    solver.parameters.log_search_progress = False

    start = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - start

    status_name = solver.StatusName(status)

    print(f"  Solver status:       {status_name}")
    print(f"  Wall time:           {elapsed:.3f}s")
    print(f"  Branches:            {solver.NumBranches()}")
    print(f"  Conflicts:           {solver.NumConflicts()}")

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print(f"\n❌ SOLVER FAILED: {status_name}")
        return

    # ------------------------------------------------------------------
    # STEP 5: Post-process the solution
    # ------------------------------------------------------------------
    lengths = [solver.Value(length_var[i]) for i in range(1, N + 1)]
    selected_words = sorted(
        word for word in admissible_words if solver.Value(y[word]) == 1
    )

    assignment = reconstruct_assignment(lengths, selected_words)

    # ------------------------------------------------------------------
    # STEP 6: Verify feasibility
    # ------------------------------------------------------------------
    metrics = evaluate_solution(lengths, selected_words, probabilities)

    print(f"  Feasible:            {metrics['feasible']}")
    print(f"  Prefix-free:         {metrics['prefix_feasible']}")
    print(f"  Order satisfied:     {metrics['order_feasible']}")

    # ------------------------------------------------------------------
    # STEP 7: Evaluate and report results
    # ------------------------------------------------------------------
    print(f"\n=== RESULTS ===")
    print(f"  Status:              {status_name}")
    print(f"  Expected length:     {metrics['expected_length']:.6f}")
    print(f"  Lengths:             {lengths}")
    print(f"  Selected words:      {selected_words}")

    print(f"\n=== PREFIX CODE ===")
    for symbol, word in sorted(assignment.items(), key=lambda x: int(x[0])):
        print(f"  symbol {symbol}: {word}")

    if metrics["active_prefix_conflicts"]:
        print(f"\n  Active prefix conflicts: {metrics['active_prefix_conflicts']}")

    # Save results
    round_result = {
        "round": 1,
        "seed": None,
        "backend_name": "cp_sat_classical_reference",
        "case_label": CASE_LABEL,
        "status": status_name,
        "feasible": metrics["feasible"],
        "objective_expected_length": metrics["expected_length"],
        "lengths": lengths,
        "selected_words": selected_words,
        "symbol_to_code": assignment,
        "metrics": metrics,
        "elapsed_seconds": elapsed,
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "instance": instance_metadata(
            CASE_LABEL, N, L, probabilities,
            FORBIDDEN_WORDS, FORBIDDEN_PREFIXES, words_by_length,
        ),
    }

    payload = {
        "script_type": "cp_sat_classical_reference",
        "case_label": CASE_LABEL,
        "rounds": 1,
        "parameters": {
            "backend": "cp_sat_classical_reference",
            "n": N,
            "L": L,
            "forbidden_words": sorted(FORBIDDEN_WORDS),
            "forbidden_prefixes": sorted(FORBIDDEN_PREFIXES),
            "probabilities": probabilities,
            "max_time_seconds": 60.0,
            "num_search_workers": 8,
        },
        "total_elapsed_seconds": elapsed,
        "best_overall": round_result,
        "all_rounds": [round_result],
    }

    summary = {
        "case_label": CASE_LABEL,
        "backend": "cp_sat_classical_reference",
        "status": status_name,
        "feasible": metrics["feasible"],
        "objective_expected_length": metrics["expected_length"],
        "best_solution": {
            "lengths": lengths,
            "selected_words": selected_words,
            "code": assignment,
        },
        "instance": round_result["instance"],
    }

    save_json(OUTPUT_RESULTS_JSON, payload)
    save_json(OUTPUT_SUMMARY_JSON, summary)

    print(f"\n{'=' * 80}")
    print("EXPERIMENT FINISHED")
    print(f"{'=' * 80}")
    print(f"  Results: {OUTPUT_RESULTS_JSON}")
    print(f"  Summary: {OUTPUT_SUMMARY_JSON}")


if __name__ == "__main__":
    run_classical_experiment()