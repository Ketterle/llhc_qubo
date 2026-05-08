import json
import time
from itertools import combinations, product
from math import gcd, log2
from functools import reduce
from pathlib import Path

from ortools.sat.python import cp_model


# ============================================================
# Basic utilities
# ============================================================

def all_binary_words(length):
    return ["".join(bits) for bits in product("01", repeat=length)]


def build_dictionary_by_length(L, forbidden_words=None, forbidden_prefixes=None):
    """
    Full binary dictionary up to length L, then remove forbidden words/prefixes.

    For the advanced case, the same restricted dictionary pattern used
    in the QUBO implementation is applied:

      length = 1  -> keep all 2 words
      length >= 2 -> keep 3 * 2^(length - 2) words

    Afterwards, explicit forbidden words and forbidden prefixes are removed.
    """
    forbidden_words = set(forbidden_words or [])
    forbidden_prefixes = tuple(forbidden_prefixes or [])

    words_by_length = {}

    for length in range(1, L + 1):
        words = all_binary_words(length)

        if length == 1:
            kept_words = words[:]
        else:
            n_keep = 3 * (2 ** (length - 2))
            kept_words = words[:n_keep]

        words_by_length[length] = [
            word for word in kept_words
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


def scale_probabilities_to_ints(probabilities, scale=100000):
    """
    Convert real-valued probabilities into integer coefficients.

    CP-SAT requires integer objective coefficients. Multiplying all
    probabilities by the same factor preserves the relative ordering
    of feasible solutions.
    """
    raw_weights = [int(round(scale * p)) for p in probabilities]

    if all(weight == 0 for weight in raw_weights):
        raise ValueError("All scaled probabilities are zero. Increase scale.")

    common_divisor = reduce(gcd, raw_weights)
    return [weight // common_divisor for weight in raw_weights]


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
# Advanced Case Configuration
# ============================================================

CASE_LABEL = "abecedary_classical"
N = 24
L = 6

FORBIDDEN_WORDS = set()
FORBIDDEN_PREFIXES = tuple()

MAX_TIME_SECONDS = 60.0
NUM_SEARCH_WORKERS = 8
PROBABILITY_SCALE = 100000
ENFORCE_PROBABILITY_LENGTH_ORDER = True
LOG_SEARCH_PROGRESS = True

OUTPUT_RESULTS_JSON = "abecedary_cpsat_results.json"
OUTPUT_SUMMARY_JSON = "abecedary_cpsat_summary.json"


# ============================================================
# Classical solver: CP-SAT
# ============================================================

def reconstruct_assignment(lengths, selected_words):
    """
    Reconstruct a deterministic symbol-to-codeword assignment.

    Symbols are assumed to be sorted by non-increasing probability.
    For each length, selected words are assigned in lexicographic order.
    """
    words_by_length = {}

    for word in selected_words:
        words_by_length.setdefault(len(word), []).append(word)

    for length in words_by_length:
        words_by_length[length].sort()

    assignment = {}

    for symbol_index, length in enumerate(lengths, start=1):
        if length not in words_by_length or not words_by_length[length]:
            raise RuntimeError(
                f"No selected word of length {length} "
                f"for symbol {symbol_index}."
            )

        assignment[str(symbol_index)] = words_by_length[length].pop(0)

    return assignment


def evaluate_code(lengths, symbol_to_code, probabilities):
    """Evaluate objective value and feasibility of a prefix code."""
    selected_words = list(symbol_to_code.values())

    active_prefix_conflicts = [
        (u, v)
        for u, v in combinations(selected_words, 2)
        if prefix_conflict(u, v)
    ]

    order_feasible = all(
        lengths[i] <= lengths[i + 1]
        for i in range(len(lengths) - 1)
    )

    expected_length = sum(
        probabilities[i] * lengths[i]
        for i in range(len(lengths))
    )

    return {
        "expected_length": expected_length,
        "prefix_feasible": len(active_prefix_conflicts) == 0,
        "active_prefix_conflicts": active_prefix_conflicts,
        "order_feasible": order_feasible,
        "feasible": len(active_prefix_conflicts) == 0 and order_feasible,
    }


def solve_cpsat(
    L,
    n,
    probabilities,
    words_by_length,
    enforce_probability_length_order=True,
    max_time_seconds=60.0,
    num_search_workers=8,
    log_search_progress=True,
    probability_scale=100000,
):
    """
    Solve the length-limited prefix code problem with restricted dictionary
    using Google OR-Tools CP-SAT solver.

    Variables:
        x[i,l] : 1 if symbol i is assigned length l
        y[w]   : 1 if word w is selected
        len[i] : selected length for symbol i

    Objective:
        Minimize the expected codeword length.
    """
    if len(probabilities) != n:
        raise ValueError("The probability vector must have length n.")

    if any(probabilities[i] < probabilities[i + 1] for i in range(len(probabilities) - 1)):
        raise ValueError("Probabilities must be sorted in non-increasing order.")

    admissible_words = flatten_words(words_by_length)
    prefix_edges = build_prefix_edges(admissible_words)

    if len(admissible_words) < n:
        raise ValueError(
            f"Only {len(admissible_words)} admissible words "
            f"available for {n} symbols."
        )

    # ------------------------------------------------------------------
    # STEP 3: Generate the CP-SAT formulation
    #   Build the constraint satisfaction model with:
    #     - Assignment constraint: each symbol gets exactly one length
    #     - Coupling constraint: word counts match length assignments
    #     - Prefix constraint: no two selected words are prefixes
    #     - Ordering constraint: monotonicity of lengths w.r.t. probabilities
    #   Objective: minimize expected codeword length.
    # ------------------------------------------------------------------
    model = cp_model.CpModel()

    x = {
        (i, length): model.NewBoolVar(f"x_{i}_{length}")
        for i in range(1, n + 1)
        for length in range(1, L + 1)
    }

    y = {
        word: model.NewBoolVar(f"y_{word}")
        for word in admissible_words
    }

    length_var = {
        i: model.NewIntVar(1, L, f"len_{i}")
        for i in range(1, n + 1)
    }

    # Assignment constraint: exactly one length per symbol
    for i in range(1, n + 1):
        model.AddExactlyOne(x[i, length] for length in range(1, L + 1))
        model.Add(
            length_var[i]
            == sum(length * x[i, length] for length in range(1, L + 1))
        )

    # Coupling constraint: number of symbols per length equals
    # number of selected words of that length
    for length in range(1, L + 1):
        model.Add(
            sum(x[i, length] for i in range(1, n + 1))
            == sum(y[word] for word in words_by_length[length])
        )

    # Prefix constraint: no two selected words may be prefixes of each other
    for u, v in prefix_edges:
        model.Add(y[u] + y[v] <= 1)

    # Ordering constraint: higher-probability symbols should not receive
    # longer codewords than lower-probability symbols
    if enforce_probability_length_order:
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                model.Add(length_var[i] <= length_var[j])

    # Objective: minimize expected codeword length using integer weights
    integer_weights = scale_probabilities_to_ints(
        probabilities,
        scale=probability_scale,
    )

    model.Minimize(
        sum(integer_weights[i - 1] * length_var[i] for i in range(1, n + 1))
    )

    print(f"  Variables:           {len(x) + len(y) + len(length_var)}")
    print(f"  Length variables:    {len(x)}")
    print(f"  Selection variables: {len(y)}")
    print(f"  Auxiliary variables: {len(length_var)}")
    print(f"  Model built successfully")

    # ------------------------------------------------------------------
    # STEP 4: Solve the problem
    #   The CP-SAT model is solved using Google OR-Tools.
    # ------------------------------------------------------------------
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = max_time_seconds
    solver.parameters.num_search_workers = num_search_workers
    solver.parameters.log_search_progress = log_search_progress

    start = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - start

    status_name = solver.StatusName(status)

    print(f"\n  Solver status:       {status_name}")
    print(f"  Wall time:           {solver.WallTime():.3f}s")
    print(f"  Elapsed time:        {elapsed:.3f}s")
    print(f"  Branches:            {solver.NumBranches()}")
    print(f"  Conflicts:           {solver.NumConflicts()}")

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {
            "status": status_name,
            "feasible": False,
            "objective_expected_length": None,
            "lengths": None,
            "length_counts": None,
            "selected_words": None,
            "symbol_to_code": None,
            "metrics": None,
            "branches": solver.NumBranches(),
            "conflicts": solver.NumConflicts(),
            "wall_time": solver.WallTime(),
            "elapsed_seconds": elapsed,
        }

    # ------------------------------------------------------------------
    # STEP 5: Post-process the solution
    #   Extract the assigned lengths from x variables and the selected
    #   words from y variables. Reconstruct the symbol-to-codeword
    #   mapping by assigning words in lexicographic order within each length.
    # ------------------------------------------------------------------
    lengths = [solver.Value(length_var[i]) for i in range(1, n + 1)]

    selected_words = sorted(
        word for word in admissible_words
        if solver.Value(y[word]) == 1
    )

    symbol_to_code = reconstruct_assignment(lengths, selected_words)

    length_counts = {
        length: sum(1 for value in lengths if value == length)
        for length in range(1, L + 1)
        if sum(1 for value in lengths if value == length) > 0
    }

    # ------------------------------------------------------------------
    # STEP 6: Verify feasibility
    #   Check that:
    #     - selected words form a prefix-free set
    #     - lengths are monotonically non-decreasing
    # ------------------------------------------------------------------
    metrics = evaluate_code(lengths, symbol_to_code, probabilities)

    # ------------------------------------------------------------------
    # STEP 7: Evaluate and report results
    #   Report the objective value, length profile, selected words,
    #   feasibility status, and solver statistics.
    # ------------------------------------------------------------------
    print(f"\n=== RESULTS ===")
    print(f"  Status:              {status_name}")
    print(f"  Feasible:            {metrics['feasible']}")
    print(f"  Expected length:     {metrics['expected_length']:.6f}")
    print(f"  Prefix-free:         {metrics['prefix_feasible']}")
    print(f"  Order satisfied:     {metrics['order_feasible']}")
    print(f"  Lengths:             {lengths}")
    print(f"  Length counts:       {length_counts}")
    print(f"  Selected words:      {selected_words}")

    print(f"\n=== PREFIX CODE ===")
    for symbol, word in sorted(symbol_to_code.items(), key=lambda x: int(x[0])):
        print(f"  symbol {symbol}: {word}")

    if metrics["active_prefix_conflicts"]:
        print(f"\n  Active prefix conflicts: {metrics['active_prefix_conflicts']}")

    return {
        "status": status_name,
        "feasible": metrics["feasible"],
        "objective_expected_length": metrics["expected_length"],
        "lengths": lengths,
        "length_counts": length_counts,
        "selected_words": selected_words,
        "symbol_to_code": symbol_to_code,
        "metrics": metrics,
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "wall_time": solver.WallTime(),
        "elapsed_seconds": elapsed,
    }


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
    print(f"CLASSICAL SOLVER (CP-SAT) — ADVANCED CASE (n = {N}, L = {L})")
    print("=" * 80)
    print(f"  Case label:          {CASE_LABEL}")
    print(f"  n:                   {N}")
    print(f"  L:                   {L}")
    print(f"  p:                   [{', '.join(f'{pi:.6f}' for pi in probabilities[:5])}, ...]")
    print(f"  sum(p):              {sum(probabilities):.10f}")
    print(f"  Entropy:             {entropy(probabilities):.6f}")
    print(f"  Forbidden words:     {len(FORBIDDEN_WORDS)}")
    print(f"  Forbidden prefixes:  {FORBIDDEN_PREFIXES if FORBIDDEN_PREFIXES else 'None'}")
    print(f"  Max time:            {MAX_TIME_SECONDS}s")
    print(f"  Workers:             {NUM_SEARCH_WORKERS}")
    print(f"  Probability scale:   {PROBABILITY_SCALE}")

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

    print("\n=== START EXPERIMENT ===")

    result = solve_cpsat(
        L=L,
        n=N,
        probabilities=probabilities,
        words_by_length=words_by_length,
        enforce_probability_length_order=ENFORCE_PROBABILITY_LENGTH_ORDER,
        max_time_seconds=MAX_TIME_SECONDS,
        num_search_workers=NUM_SEARCH_WORKERS,
        log_search_progress=LOG_SEARCH_PROGRESS,
        probability_scale=PROBABILITY_SCALE,
    )

    metadata = instance_metadata(
        CASE_LABEL,
        N,
        L,
        probabilities,
        FORBIDDEN_WORDS,
        FORBIDDEN_PREFIXES,
        words_by_length,
    )

    round_result = {
        "round": 1,
        "seed": None,
        "backend_name": "cpsat_classical_reference",
        "case_label": CASE_LABEL,
        "status": result["status"],
        "feasible": result["feasible"],
        "objective_expected_length": result["objective_expected_length"],
        "lengths": result["lengths"],
        "length_counts": (
            {str(k): v for k, v in result["length_counts"].items()}
            if result["length_counts"] is not None
            else None
        ),
        "selected_words": result["selected_words"],
        "symbol_to_code": result["symbol_to_code"],
        "metrics": result["metrics"],
        "branches": result["branches"],
        "conflicts": result["conflicts"],
        "wall_time": result["wall_time"],
        "elapsed_seconds": result["elapsed_seconds"],
        "instance": metadata,
    }

    payload = {
        "script_type": "cpsat_classical_reference",
        "case_label": CASE_LABEL,
        "rounds": 1,
        "parameters": {
            "backend": "cpsat_classical_reference",
            "n": N,
            "L": L,
            "forbidden_words": sorted(FORBIDDEN_WORDS),
            "forbidden_prefixes": sorted(FORBIDDEN_PREFIXES),
            "probabilities": probabilities,
            "max_time_seconds": MAX_TIME_SECONDS,
            "num_search_workers": NUM_SEARCH_WORKERS,
            "probability_scale": PROBABILITY_SCALE,
            "enforce_probability_length_order": ENFORCE_PROBABILITY_LENGTH_ORDER,
        },
        "total_elapsed_seconds": result["elapsed_seconds"],
        "best_overall": round_result,
        "all_rounds": [round_result],
    }

    summary = {
        "case_label": CASE_LABEL,
        "backend": "cpsat_classical_reference",
        "status": result["status"],
        "feasible": result["feasible"],
        "objective_expected_length": result["objective_expected_length"],
        "best_solution": {
            "lengths": result["lengths"],
            "length_counts": (
                {str(k): v for k, v in result["length_counts"].items()}
                if result["length_counts"] is not None
                else None
            ),
            "selected_words": result["selected_words"],
            "code": result["symbol_to_code"],
        },
        "solver_stats": {
            "wall_time": result["wall_time"],
            "elapsed_seconds": result["elapsed_seconds"],
            "branches": result["branches"],
            "conflicts": result["conflicts"],
        },
        "instance": metadata,
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