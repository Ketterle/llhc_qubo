# QUBO formulation solver for Prefix-Free Code Optimization

## Overview

This project implements and compares **two solver approaches** for the **length-limited prefix-free binary code problem**:

1. **Quantum-Hybrid Solvers** (QUBO formulation + D-Wave Leap Hybrid)
2. **Classical Solvers** (Package-Merge / CP-SAT with Google OR-Tools)

**Each of the three cases is solved using BOTH approaches** to enable direct comparison of quantum vs. classical optimization.

---

## Three Problem Cases

All three cases are solved with **both Classical and QUBO solvers**:

| Case | N | L | Forbidden Prefix | Classical Solver | QUBO Solver |
|------|---|---|------------------|------------------|-------------|
| **Base** | 8 | 4 | None | Package-Merge | D-Wave Hybrid |
| **Intermediate** | 8 | 4 | "11" | CP-SAT | D-Wave Hybrid |
| **Advanced (Abecedary)** | 24 | 6 | "11" | CP-SAT | D-Wave Hybrid |

**Goal:** Compare classical exact/heuristic methods vs. quantum-hybrid approaches across increasing problem complexity.

---

## Problem Statement

**Given:**
- N symbols with non-increasing probability distribution
- L maximum codeword length
- A restricted dictionary of admissible binary words
- Constraints: Prefix-free, length-limited, probability-length monotonicity

**Goal:** Find a prefix-free binary code that **minimizes expected codeword length**.

### Constraints:

1. **Prefix-Free**: No codeword is a prefix of another
2. **Length-Limited**: All codewords have length ≤ L
3. **Assignment**: Each symbol gets exactly one codeword length
4. **Coupling**: Number of selected words per length matches assigned lengths
5. **Probability-Length Ordering**: Higher probability symbols get shorter codewords

### Mathematical Formulation:

**Objective Function:**
```
Minimize: H_code = Σ(i=1 to N) p_i · l_i
```

where:
- `p_i` = probability of symbol i
- `l_i` = assigned codeword length for symbol i

**QUBO Hamiltonian:**
```
H = H_cost + λ_A·H_ass + λ_C·H_coupl + λ_P·H_pref + λ_O·H_ord
```

where:
- `H_cost` = Expected codeword length (objective)
- `H_ass` = Assignment constraint penalty
- `H_coupl` = Coupling constraint penalty
- `H_pref` = Prefix-free constraint penalty
- `H_ord` = Probability-length ordering penalty
- `λ_A, λ_C, λ_P, λ_O` = Penalty coefficients

---

## Probability Model

All three cases use **harmonic (Zipf) distribution**:

```
p_i = (1/i) / Σ(j=1 to N)(1/j),  ∀ i ∈ {1, ..., N}
```

### Base & Intermediate Cases (N=8):

**Formula:**
```
p_i = (1/i) / Σ(j=1 to 8)(1/j),  ∀ i ∈ {1, ..., 8}
```

**Normalization constant:** `Σ(j=1 to 8)(1/j) = 2.7178571429`

**Probability values:**
- p₁ = 0.367937
- p₂ = 0.183968
- p₃ = 0.122646
- p₄ = 0.091984
- p₅ = 0.073587
- p₆ = 0.061323
- p₇ = 0.052562
- p₈ = 0.045992

**Entropy:** H(X) = 2.619715 bits/symbol

### Advanced Case (N=24):

**Formula:**
```
p_i = (1/i) / Σ(j=1 to 24)(1/j),  ∀ i ∈ {1, ..., 24}
```

**Normalization constant:** `Σ(j=1 to 24)(1/j) = 3.7759581778`

**Probability values (first 5):**
- p₁ = 0.264833
- p₂ = 0.132417
- p₃ = 0.088278
- p₄ = 0.066208
- p₅ = 0.052967
- ... (continues to p₂₄ = 0.011035)

**Entropy:** H(X) = 3.843676 bits/symbol

---

## Dictionary Structures

### Base & Intermediate Cases: Full Dictionary

**Structure:**
- All `2^l` binary words of length l for `1 ≤ l ≤ L`
- **Total words:** `Σ(l=1 to L) 2^l = 2^(L+1) - 2`

**For L=4:**
- Length 1: 2 words (0, 1)
- Length 2: 4 words (00, 01, 10, 11)
- Length 3: 8 words
- Length 4: 16 words
- **Total: 30 words**

**Filtering (Intermediate only):**
- Remove all words starting with "11"
- Results in **27 admissible words**

### Advanced Case: Abecedary Dictionary

**Structure:**
```
Number of words at length l:
  n_l = { 2,           if l = 1
        { 3·2^(l-2),   if l ≥ 2
```

**Rationale:** Models linguistic distributions where short words are scarce but longer words grow exponentially.

**For L=6 (before filtering):**
- Length 1: 2 words
- Length 2: 3 words
- Length 3: 6 words
- Length 4: 12 words
- Length 5: 24 words
- Length 6: 48 words
- **Total: 95 words**

**After filtering "11" prefix:**
- **95 admissible words** (abecedary structure prevents "11" in first 3 positions)

---

## Project Structure

```
.
├── QUANTUM SOLVERS (QUBO Formulation - D-Wave Hybrid)
│   ├── base_case_uniform.py              # QUBO for Base Case
│   ├── intermediate_case_uniform.py       # QUBO for Intermediate Case
│   └── advanced_case_uniform.py           # QUBO for Advanced Case (Abecedary)
│
├── CLASSICAL SOLVERS (Reference Implementations)
│   ├── base_case_classical.py             # Package-Merge for Base Case
│   ├── intermediate_case_classical.py     # CP-SAT for Intermediate Case
│   └── advanced_case_classical.py         # CP-SAT for Advanced Case (Abecedary)
│
├── DOCUMENTATION
    ├── QUBO_Matrices_All_Cases.pdf        # 118-page QUBO matrix visualizations
```

---

## Installation & Dependencies

### Prerequisites
- **Python 3.8+**

### Required Packages

**For Classical Solvers (Package-Merge & CP-SAT):**
```bash
pip install numpy google-ortools
```

**For QUBO Solvers (D-Wave Leap Hybrid):**
```bash
pip install numpy dimod dwave-system
```

**Install Everything:**
```bash
pip install numpy dimod google-ortools dwave-system
```

### D-Wave Leap Setup (Required for QUBO Solvers)

1. **Get a free D-Wave Leap account:** https://cloud.dwavesys.com/leap/signup/
2. **Get your API token** from the dashboard
3. **Configure your environment:**
   ```bash
   dwave config create
   ```
   Follow the prompts to enter your API token

4. **Test your connection:**
   ```bash
   dwave ping
   ```

---

## How to Run

### Complete Workflow: Run Both Approaches for Each Case

Each case should be run with **both solvers** for comparison:

---

### **BASE CASE (N=8, L=4, Full Dictionary, No Restrictions)**

**1. Classical Solver (Package-Merge - Exact):**
```bash
python base_case_classical.py
```
**Expected time:** < 1 second  
**Output:** 
- `base_case_package_merge_results.json`
- `base_case_package_merge_summary.json`

**2. QUBO Solver (D-Wave Hybrid):**
```bash
python base_case_qubo.py
```
**Expected time:** 15-20 minutes (100 rounds × 10s each)  
**Output:**
- `base_case_qubo_results.json`
- `base_case_qubo_summary.json`

---

### **INTERMEDIATE CASE (N=8, L=4, Full Dictionary, Forbidden "11")**

**1. Classical Solver (CP-SAT):**
```bash
python intermediate_classical.py
```
**Expected time:** 1-5 seconds  
**Output:**
- `intermediate_case_cpsat_results.json`
- `intermediate_case_cpsat_summary.json`

**2. QUBO Solver (D-Wave Hybrid):**
```bash
python intermediate_qubo.py
```
**Expected time:** 15-20 minutes (100 rounds × 10s each)  
**Output:**
- `intermediate_qubo_results.json`
- `intermediate_qubo_summary.json`

---

### **ADVANCED CASE (N=24, L=6, Abecedary Dictionary, Forbidden "11")**

**1. Classical Solver (CP-SAT):**
```bash
python abecedary_classical.py
```
**Expected time:** 10-60 seconds  
**Output:**
- `advanced_case_cpsat_results.json`
- `advanced_case_cpsat_summary.json`

**2. QUBO Solver (D-Wave Hybrid):**
```bash
python abecedary_qubo.py
```
**Expected time:** 15-20 minutes (100 rounds × 10s each)  
**Output:**
- `advanced_qubo_results.json`
- `advanced_qubo_summary.json`

---

## Configuring QUBO Runtime Parameters

For the QUBO solvers, the number of independent hybrid solver calls and the time assigned to each call can be tuned directly in each `*_uniform.py` script through the following configuration variables:

```python
# Configuration section in *_uniform.py scripts
ROUNDS = 100                    # Number of independent solver calls
TIME_LIMIT_SECONDS = 10         # Time per D-Wave Hybrid call (seconds)
MAX_DWAVE_SECONDS = 1000        # Maximum total D-Wave time allowed (seconds)
```

**Example configurations:**

**Quick test (1-2 minutes):**
```python
ROUNDS = 10
TIME_LIMIT_SECONDS = 3
MAX_DWAVE_SECONDS = 30
```

**Standard run (15-20 minutes):**
```python
ROUNDS = 100
TIME_LIMIT_SECONDS = 10
MAX_DWAVE_SECONDS = 1000
```

**Extended run (50-60 minutes):**
```python
ROUNDS = 300
TIME_LIMIT_SECONDS = 10
MAX_DWAVE_SECONDS = 3000
```

**Note:** Actual D-Wave time consumed will be `ROUNDS × TIME_LIMIT_SECONDS` (not including overhead).

---

## Understanding the Output

### JSON Result Structure

All solvers produce **two JSON files**:

#### 1. `*_results.json` - Complete Results

**Example structure:**
```json
{
  "script_type": "plain_qubo_fixed",
  "case_label": "intermediate",
  "rounds": 100,
  "planned_dwave_seconds": 300,
  "parameters": {
    "backend": "leap_hybrid",
    "n": 8,
    "L": 4,
    "forbidden_prefixes": ["11"],
    "lambda_a": 3.0,
    "lambda_c": 3.0,
    "lambda_p": 5.0,
    "lambda_o": 3.0
  },
  "best_overall": {
    "round": 42,
    "energy": 17.234567,
    "feasible": true,
    "H_cost": 2.156234,
    "H_ass": 0.0,
    "H_coupl": 0.0,
    "H_pref": 0.0,
    "H_ord": 0.0,
    "symbol_to_code": {
      "1": "0",
      "2": "10",
      "3": "010",
      ...
    }
  },
  "all_rounds": [...]
}
```

#### 2. `*_summary.json` - Statistics Summary

**Example structure:**
```json
{
  "experiment": "qubo_intermediate",
  "case_label": "intermediate",
  "rounds": 100,
  "total_elapsed_s": 312.45,
  "feasibility_rate": 0.87,
  "feasible_energy_stats": {
    "min": 17.234567,
    "max": 18.456789,
    "mean": 17.891234,
    "std": 0.234567,
    "median": 17.823456
  },
  "best_solution": {
    "round": 42,
    "energy": 17.234567,
    "H_cost": 2.156234,
    "code": {...}
  }
}
```

### Interpretation of Metrics

**For QUBO Solvers:**

| Metric | Meaning | Target | Interpretation |
|--------|---------|--------|----------------|
| **energy** | Total QUBO energy | Minimize | Lower = better solution |
| **H_cost** | Expected codeword length | Minimize | Objective function value |
| **H_ass** | Assignment violations | = 0 | Must be 0 for feasibility |
| **H_coupl** | Coupling violations | = 0 | Must be 0 for feasibility |
| **H_pref** | Prefix conflicts | = 0 | Must be 0 for feasibility |
| **H_ord** | Ordering violations | = 0 | Must be 0 for feasibility |
| **feasible** | All constraints OK | true | Solution is valid |

**For Classical Solvers:**

| Metric | Meaning | Target | Interpretation |
|--------|---------|--------|----------------|
| **status** | Solver status | OPTIMAL | Solution quality |
| **objective_expected_length** | Expected length | Minimize | Same as H_cost |
| **feasible** | All constraints OK | true | Solution is valid |
| **prefix_feasible** | No prefix conflicts | true | Prefix-free property |
| **order_feasible** | Length ordering OK | true | Monotonicity satisfied |

---

## STEP 1-7 Execution Flow

All scripts follow a **uniform 7-step structure**:

### STEP 1: Define the Instance and Execution Budget
- Set N (number of symbols)
- Set L (maximum codeword length)
- Generate probability vector (harmonic distribution)
- Define forbidden words/prefixes
- Set solver parameters (rounds, time limits)

### STEP 2: Build the Admissible Dictionary and Conflict Graph
- Enumerate binary words up to length L
- Apply dictionary restrictions:
  - Base: Keep all words
  - Intermediate: Remove words starting with "11"
  - Advanced: Apply abecedary structure, then remove "11" prefixes
- Compute prefix-conflict pairs

### STEP 3: Generate the Optimization Formulation
- **QUBO solvers:** Build Hamiltonian with penalty terms
  ```
  H = H_cost + λ_A·H_ass + λ_C·H_coupl + λ_P·H_pref + λ_O·H_ord
  ```
- **Classical solvers:** 
  - Package-Merge: Enumerate Kraft-feasible profiles
  - CP-SAT: Build constraint satisfaction model

### STEP 4: Solve the Problem
- **QUBO:** Submit BQM to D-Wave Leap Hybrid (multiple rounds)
- **Classical:** 
  - Package-Merge: Exhaustive search with pruning
  - CP-SAT: Branch-and-bound with constraint propagation
- Track solver statistics

### STEP 5: Post-Process the Solution
- Extract assigned lengths from solution
- Select admissible words
- Reconstruct symbol-to-codeword mapping
- For QUBO: Track best feasible across all rounds

### STEP 6: Verify Feasibility
- Check all constraints satisfied:
  - Assignment: Each symbol has exactly one length
  - Coupling: Word counts match length assignments
  - Prefix-free: No word is prefix of another
  - Ordering: Higher probability → shorter length

### STEP 7: Evaluate and Report Results
- Compute final metrics
- Save JSON results (full + summary)
- Print summary to console

---

## Comparison: Classical vs. QUBO

### Approach Comparison

| Aspect | Classical (Package-Merge / CP-SAT) | QUBO (D-Wave Hybrid) |
|--------|-----------------------------------|----------------------|
| **Methodology** | Exact/heuristic optimization | Quantum-inspired annealing |
| **Hardware** | CPU (multi-threaded) | Hybrid quantum-classical |
| **Formulation** | Constraint satisfaction | Quadratic unconstrained |
| **Guarantees** | Optimal (PM) / Near-optimal (SAT) | Best-effort feasible |
| **Setup** | Standard Python packages | D-Wave Leap account required |
| **Time/Round** | 0.01s - 60s | 10s per round |
| **Typical Rounds** | 1 | 100 |
| **Total Time** | 0.01s - 60s | 15-20 minutes |
| **Cost** | Free | Leap minutes (free tier available) |

### Problem Complexity Comparison

| Metric | Base | Intermediate | Advanced |
|--------|------|--------------|----------|
| **N (Symbols)** | 8 | 8 | 24 |
| **L (Max Length)** | 4 | 4 | 6 |
| **Full Dictionary Size** | 30 | 30 | 126 |
| **Admissible Words** | 30 | 27 | 95 |
| **Prefix Conflicts** | 20 | 15 | 387 |
| **x Variables** | 32 | 32 | 144 |
| **y Variables** | 30 | 27 | 95 |
| **Total Variables** | 62 | 59 | 239 |
| **QUBO Matrix Size** | 62×62 | 59×59 | 239×239 |
| **Classical Solve Time** | <0.01s | 1-5s | 10-60s |
| **QUBO Solve Time** | ~15-20 min | ~15-20 min | ~15-20 min |

---

## Expected Results

### Theoretical Lower Bounds (Shannon Entropy)

The **Shannon entropy** provides a theoretical lower bound on expected codeword length:

| Case | Entropy H(X) | Interpretation |
|------|--------------|----------------|
| **Base** | 2.6197 bits/symbol | Minimum possible average length |
| **Intermediate** | 2.6197 bits/symbol | Same distribution as Base |
| **Advanced** | 3.8437 bits/symbol | More symbols = higher entropy |

### Typical Solution Quality

**Base Case:**
- **Classical (Package-Merge):** ~2.62 bits (OPTIMAL, matches entropy)
- **QUBO (D-Wave):** ~2.62-2.70 bits (high feasibility rate)

**Intermediate Case:**
- **Classical (CP-SAT):** ~2.16-2.30 bits (OPTIMAL)
- **QUBO (D-Wave):** ~2.20-2.40 bits (feasibility rate ~80%)

**Advanced Case:**
- **Classical (CP-SAT):** ~3.90-4.10 bits (OPTIMAL/near-optimal)
- **QUBO (D-Wave):** ~4.00-4.30 bits (feasibility rate ~60-80%)

**Note:** QUBO solvers may find infeasible solutions (constraint violations). The scripts track both:
- Best overall energy (may be infeasible)
- Best feasible energy (all constraints satisfied)

---

## Quick Start Summary

**Minimal setup to run everything:**

```bash
# Install dependencies
pip install numpy dimod google-ortools dwave-system

# Configure D-Wave (for QUBO solvers)
dwave config create

# Run all 6 experiments (Classical + QUBO for each case)
python base_case_classical.py
python base_case_uniform.py
python intermediate_case_classical.py
python intermediate_case_uniform.py
python advanced_case_classical.py
python advanced_case_uniform.py

# Results are saved as JSON files for analysis
```

---

## Summary: Running Complete Experiments

For **each case**, run both solvers to compare:

```bash
# BASE CASE
python base_case_classical.py      # Classical
python base_case_uniform.py        # QUBO

# INTERMEDIATE CASE
python intermediate_case_classical.py    # Classical
python intermediate_case_uniform.py      # QUBO

# ADVANCED CASE
python advanced_case_classical.py   # Classical
python advanced_case_uniform.py     # QUBO
```

**Ready to explore classical vs. quantum optimization for prefix-free codes!** 🚀
