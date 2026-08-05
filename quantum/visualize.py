"""
RNA Secondary Structure Visualization.

Draws arc diagrams of predicted vs. reference secondary structures, in the
same spirit as the "before/after" visual in WISER's own reference deck
(https://alexgalda.github.io/quantum_mRNA_optimization/#/5), but going one
step further: pairs are color-coded by whether they're correct, wrong, or
missed relative to the true ViennaRNA MFE structure, so the figure itself
shows where a solver succeeded or where the QUBO formulation has a gap,
rather than just showing "an answer."
"""

from typing import Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

BASE_COLORS = {"A": "#4C72B0", "U": "#DD8452", "G": "#55A868", "C": "#C44E52"}


def pairs_from_structure(structure: str) -> Set[Tuple[int, int]]:
    """Parses dot-bracket notation into a set of (i, j) index pairs, i < j."""
    stack, pairs = [], set()
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError(f"Unbalanced structure at position {i}: {structure}")
            pairs.add((stack.pop(), i))
    if stack:
        raise ValueError(f"Unbalanced structure, unclosed '(' remaining: {structure}")
    return pairs


def _draw_backbone(ax, sequence: str):
    n = len(sequence)
    ax.set_xlim(-1, n)
    ax.set_ylim(-0.7, max(n / 2.2, 1.5) + 1)
    ax.axis("off")
    ax.plot(range(n), [0] * n, color="#888888", linewidth=1.5, zorder=1)
    for i, base in enumerate(sequence):
        ax.scatter(i, 0, s=140, color=BASE_COLORS.get(base, "#999999"),
                   zorder=3, edgecolor="white", linewidth=0.8)
        ax.text(i, -0.55, base, ha="center", va="center", fontsize=8)


def _draw_arc(ax, i: int, j: int, color: str, style: str = "-", alpha: float = 0.9, lw: float = 1.6):
    """Draws an upward semicircular arc between positions i and j.
    Order-independent: uses abs() so it doesn't matter which of i, j is larger,
    a previous version assumed i > j, which silently flipped arcs below the
    baseline whenever a (smaller, larger) tuple was passed straight through.
    """
    mid = (i + j) / 2
    radius = abs(i - j) / 2
    theta = np.linspace(0, np.pi, 100)
    ax.plot(mid + radius * np.cos(theta), radius * np.sin(theta) * 0.9,
            color=color, linewidth=lw, alpha=alpha, linestyle=style, zorder=2)


def plot_structure(ax, sequence: str, structure: str, title: str = ""):
    """Draws a single arc diagram for one structure, no comparison."""
    _draw_backbone(ax, sequence)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    for (i, j) in pairs_from_structure(structure):
        _draw_arc(ax, i, j, color="#333333")


def plot_comparison(ax, sequence: str, predicted: str, reference: str, title: str = ""):
    """
    Draws one arc diagram with pairs color-coded against a reference structure:
      - green  solid : pair present in both (correct)
      - red    solid : pair only in `predicted` (wrong)
      - gray  dashed : pair only in `reference` (missed)

    Raises ValueError if `predicted`/`reference`/`sequence` lengths disagree.
    """
    if not (len(sequence) == len(predicted) == len(reference)):
        raise ValueError(
            f"Length mismatch: sequence={len(sequence)}, "
            f"predicted={len(predicted)}, reference={len(reference)}"
        )

    _draw_backbone(ax, sequence)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")

    pred_pairs = pairs_from_structure(predicted)
    ref_pairs = pairs_from_structure(reference)

    for (i, j) in pred_pairs & ref_pairs:
        _draw_arc(ax, i, j, color="#2ca02c")
    for (i, j) in pred_pairs - ref_pairs:
        _draw_arc(ax, i, j, color="#d62728")
    for (i, j) in ref_pairs - pred_pairs:
        _draw_arc(ax, i, j, color="#888888", style="--", alpha=0.7, lw=1.3)

    legend = [
        Line2D([0], [0], color="#2ca02c", lw=2, label="Correct pair"),
        Line2D([0], [0], color="#d62728", lw=2, label="Wrong pair (predicted only)"),
        Line2D([0], [0], color="#888888", lw=2, linestyle="--", label="Missed pair (true only)"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=8, frameon=False)


def compare_figure(sequence: str, predicted: str, reference: str,
                    title: Optional[str] = None, figsize=(7, 4.2)):
    """Convenience wrapper: builds a single-panel comparison figure and returns it."""
    fig, ax = plt.subplots(figsize=figsize)
    plot_comparison(ax, sequence, predicted, reference,
                     title or f"{sequence}: predicted vs. reference")
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Quick self-test using the two documented cases from this project:
    # one exact match, one known formulation-gap mismatch.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
    plot_comparison(axes[0], "GCGCAUACGC", "(((....)))", "(((....)))",
                     "GCGCAUACGC -- CVaR-VQE (exact match)")
    plot_comparison(axes[1], "GCAUCGUAGC", "((.(...)))", "..........",
                     "GCAUCGUAGC -- CVaR-VQE (formulation gap)")
    plt.tight_layout()
    plt.savefig("visualize_selftest.png", dpi=150)
    print("Self-test figure saved to visualize_selftest.png")
