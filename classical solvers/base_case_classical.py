import json
import time
from itertools import combinations, product
from math import log2
from pathlib import Path

import numpy as np


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
# Base Case Configuration
# ============================================================

CASE_LABEL = "base"
N = 8
L = 4

FORBIDDEN_WORDS = set()
FORBIDDEN_PREFIXES = tuple()

OUTPUT_RESULTS_JSON = "base_case_package_merge_results.json"
OUTPUT_SUMMARY_JSON = "base_case_package_merge_summary.json"


# ============================================================
# Classical solver: Package-Merge / exhaustive search
# ============================================================

def canonical_code_from_lengths(lengths):
    """
    Build deterministic canonical prefix code from nondecreasing code lengths.
    Symbols are assumed to be ordered by non-increasing probability.
    """
    pairs = sorted((length, symbol) for symbol, length in enumerate(lengths, start=1))

    code = 0
    previous_length = 0
    assignment = {}

    for length, symbol in pairs:
        code <<= length - previous_length
        assignment[str(symbol)] = format(code, f"0{length}b")
        code += 1
        previous_length = length

    return assignment


def evaluate_code(lengths, symbol_to_code, probabilities):
    """Evaluate objective value and feasibility of a prefix code."""
    selected_words = list(symbol_to_code.values())
    active_prefix_conflicts = [
        (u, v)
        for u, v in combinations(selected_words, 2)
        if prefix_conflict(u, v)
    ]
    order_feasible = all(lengths[i] <= lengths[i + 1] for i in range(len(lengths) - 1))
    expected_length = sum(probabilities[i] * lengths[i] for i in range(len(lengths)))

    return {
        "expected_length": expected_length,
        "prefix_feasible": len(active_prefix_conflicts) == 0,
        "active_prefix_conflicts": active_prefix_conflicts,
        "order_feasible": order_feasible,
        "feasible": len(active_prefix_conflicts) == 0 and order_feasible,
    }


def optimal_length_limited_lengths(probabilities, L):
    """
    Exhaustive search over feasible length-count profiles satisfying Kraft's inequality.

    This is the classical reference for the unrestricted case.
    It enumerates all profiles and selects the one minimizing expected codeword length.
    """
    n = len(probabilities)
    weights = sorted(probabilities, reverse=True)

    best_cost = float("inf")
    best_counts = None

    def rec(length, remaining_symbols, kraft_units, counts):
        nonlocal best_cost, best_counts

        if length == L + 1:
            if remaining_symbols == 0 and kraft_units <= 2 ** L:
                lengths = []
                for l in range(1, L + 1):
                    lengths.extend([l] * counts.get(l, 0))

                if len(lengths) != n:
                    return

                cost = sum(w * l for w, l in zip(weights, lengths))
                if cost < best_cost:
                    best_cost = cost
                    best_counts = dict(counts)
            return

        if remaining_symbols < 0 or kraft_units > 2 ** L:
            return

        for c in range(remaining_symbols + 1):
            counts[length] = c
            rec(
                length + 1,
                remaining_symbols - c,
                kraft_units + c * (2 ** (L - length)),
                counts,
            )

        counts.pop(length, None)

    rec(length=1, remaining_symbols=n, kraft_units=0, counts={})

    if best_counts is None:
        raise RuntimeError("No feasible length-limited prefix profile found.")

    lengths = []
    for l in range(1, L + 1):
        lengths.extend([l] * best_counts.get(l, 0))

    return lengths, best_counts, best_cost


# ============================================================
# Main experiment
# ============================================================

def run_classical_experiment():
    # ------------------------------------------------------------------
    # STEP 1: Define the instance and execution budget
    #   Fix n (number of symbols), L (maximum codeword length), and the
    #   probability vector p (normalized, sorted in non-increasing order).
    #   Forbidden words and prefixes define the dictionary restriction.
    # ------------------------------------------------------------------
    probabilities = harmonic_probabilities(N)

    print("=" * 80)
    print(f"CLASSICAL SOLVER (PACKAGE-MERGE) — BASE CASE (n = {N}, L = {L})")
    print("=" * 80)
    print(f"  Case label:          {CASE_LABEL}")
    print(f"  n:                   {N}")
    print(f"  L:                   {L}")
    print(f"  Forbidden words:     {len(FORBIDDEN_WORDS)}")
    print(f"  Forbidden prefixes:  {FORBIDDEN_PREFIXES if FORBIDDEN_PREFIXES else 'None'}")

    # ------------------------------------------------------------------
    # STEP 2: Build the admissible dictionary and the conflict graph
    #   - words_by_length[l] contains the admissible binary words of length l
    #   - prefix edges are computed for metadata / validation
    # ------------------------------------------------------------------
    words_by_length = build_dictionary_by_length(L, FORBIDDEN_WORDS, FORBIDDEN_PREFIXES)
    admissible_words = flatten_words(words_by_length)
    prefix_edges = build_prefix_edges(admissible_words)

    print(f"  Admissible words:    {len(admissible_words)}")
    print(f"  Prefix conflicts:    {len(prefix_edges)}")

    # ------------------------------------------------------------------
    # STEP 3: Generate the optimization formulation
    #   The Package-Merge / exhaustive search formulation enumerates
    #   feasible length-count profiles satisfying Kraft's inequality
    #   and selects the one minimizing expected codeword length.
    # ------------------------------------------------------------------
    # (formulation is implicit in optimal_length_limited_lengths)

    # ------------------------------------------------------------------
    # STEP 4: Solve the problem
    #   Exhaustive search over length-count profiles with Kraft pruning.
    # ------------------------------------------------------------------
    print("\n=== START EXPERIMENT ===")
    start = time.time()
    lengths, length_counts, objective = optimal_length_limited_lengths(probabilities, L)
    elapsed = time.time() - start

    print(f"  Solver completed in {elapsed:.3f}s")
    print(f"  Optimal objective:   {objective:.6f}")

    # ------------------------------------------------------------------
    # STEP 5: Post-process the solution
    #   Build a canonical prefix code from the optimal length profile.
    #   The canonical code assigns codewords deterministically: symbols
    #   are sorted by length (ties broken by symbol index), and codewords
    #   are assigned in lexicographic order.
    # ------------------------------------------------------------------
    symbol_to_code = canonical_code_from_lengths(lengths)

    # ------------------------------------------------------------------
    # STEP 6: Verify feasibility
    #   Check that:
    #     - selected words form a prefix-free set
    #     - lengths are monotonically non-decreasing (ordering constraint)
    # ------------------------------------------------------------------
    metrics = evaluate_code(lengths, symbol_to_code, probabilities)

    print(f"  Feasible:            {metrics['feasible']}")
    print(f"  Prefix-free:         {metrics['prefix_feasible']}")
    print(f"  Order satisfied:     {metrics['order_feasible']}")

    # ------------------------------------------------------------------
    # STEP 7: Evaluate and report results
    #   Report the objective value, length profile, selected words,
    #   feasibility status, and the final prefix code.
    # ------------------------------------------------------------------
    print(f"\n=== RESULTS ===")
    print(f"  Expected length:     {metrics['expected_length']:.6f}")
    print(f"  Lengths:             {lengths}")
    print(f"  Length counts:       {length_counts}")

    print(f"\n=== PREFIX CODE ===")
    for symbol, word in sorted(symbol_to_code.items(), key=lambda x: int(x[0])):
        print(f"  symbol {symbol}: {word}")

    if metrics["active_prefix_conflicts"]:
        print(f"\n  Active prefix conflicts: {metrics['active_prefix_conflicts']}")

    # Save results
    round_result = {
        "round": 1,
        "seed": None,
        "backend_name": "package_merge_classical_reference",
        "case_label": CASE_LABEL,
        "status": "OPTIMAL",
        "feasible": metrics["feasible"],
        "objective_expected_length": metrics["expected_length"],
        "lengths": lengths,
        "length_counts": {str(k): v for k, v in length_counts.items()},
        "selected_words": sorted(symbol_to_code.values()),
        "symbol_to_code": symbol_to_code,
        "metrics": metrics,
        "elapsed_seconds": elapsed,
        "instance": instance_metadata(
            CASE_LABEL, N, L, probabilities,
            FORBIDDEN_WORDS, FORBIDDEN_PREFIXES, words_by_length,
        ),
    }

    payload = {
        "script_type": "package_merge_classical_reference",
        "case_label": CASE_LABEL,
        "rounds": 1,
        "parameters": {
            "backend": "package_merge_classical_reference",
            "n": N,
            "L": L,
            "forbidden_words": sorted(FORBIDDEN_WORDS),
            "forbidden_prefixes": sorted(FORBIDDEN_PREFIXES),
            "probabilities": probabilities,
        },
        "total_elapsed_seconds": elapsed,
        "best_overall": round_result,
        "all_rounds": [round_result],
    }

    summary = {
        "experiment": "package_merge_base",
        "case_label": CASE_LABEL,
        "backend": "package_merge_classical_reference",
        "status": "OPTIMAL",
        "feasible": metrics["feasible"],
        "objective_expected_length": metrics["expected_length"],
        "best_solution": {
            "lengths": lengths,
            "length_counts": {str(k): v for k, v in length_counts.items()},
            "selected_words": sorted(symbol_to_code.values()),
            "code": symbol_to_code,
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