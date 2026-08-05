"""
RNA Secondary Structure Visualization.

Draws arc diagrams of predicted vs. reference secondary structures, in the
same spirit as the "before/after" visual in WISER's own reference deck
(https://alexgalda.github.io/quantum_mRNA_optimization/#/5), but going two
steps further.

1. The predicted structure and the ViennaRNA reference are drawn as *mirrored*
   arc fields around one shared sequence backbone: the solver's answer above,
   ground truth below. Agreement therefore reads as visual symmetry, and any
   disagreement shows up immediately as an arc with nothing facing it.
2. Every pair is color-coded by whether it is correct, wrong, or missed, so
   the figure shows where a solver succeeded or where the QUBO formulation has
   a gap, rather than just showing "an answer".

Nothing here invents structures: every function takes the dot-bracket strings
it is given and draws exactly those.
"""

from typing import Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

BASE_COLORS = {"A": "#4C72B0", "U": "#DD8452", "G": "#55A868", "C": "#C44E52"}

MATCH_COLOR = "#2ca02c"    # pair in both predicted and reference
WRONG_COLOR = "#d62728"    # pair predicted but not in the reference
MISSED_COLOR = "#7f7f7f"   # pair in the reference the solver did not predict

# Vertical clearance kept free around the backbone so the bases stay readable.
_BAND = 0.38

# Padding above the upper field / below the lower field. The lower one is
# larger because the dot-bracket readout is anchored there.
_PAD_TOP = 0.80
_PAD_BOTTOM = 1.40


def _field_extents(pred_pairs, ref_pairs):
    """How far the upper and lower arc fields actually reach, so a field with
    no pairs (or only short ones) doesn't reserve empty space."""
    def reach(pairs):
        widest = max((abs(j - i) / 2 for (i, j) in pairs), default=0.0)
        return max(_BAND + widest * 0.9, 1.0)
    return reach(pred_pairs), reach(ref_pairs)


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


def _draw_shared_backbone(ax, sequence: str):
    """Backbone for the mirrored layout: bases sit *on* the centre line, with
    the letter inside the marker, leaving both fields clear for arcs."""
    n = len(sequence)
    ax.plot(range(n), [0] * n, color="#cccccc", linewidth=1.2, zorder=1)
    for i, base in enumerate(sequence):
        ax.scatter(i, 0, s=200, color=BASE_COLORS.get(base, "#999999"),
                   zorder=3, edgecolor="white", linewidth=1.0)
        ax.text(i, 0, base, ha="center", va="center", fontsize=6.5,
                color="white", fontweight="bold", zorder=4)


def _draw_arc(ax, i: int, j: int, color: str, style: str = "-", alpha: float = 0.9,
              lw: float = 1.6, direction: int = 1, offset: float = 0.0):
    """Draws a semicircular arc between positions i and j.

    `direction=+1` draws it above the backbone, `-1` below; `offset` lifts the
    endpoints clear of the bases.

    Order-independent: uses abs() so it doesn't matter which of i, j is larger,
    a previous version assumed i > j, which silently flipped arcs below the
    baseline whenever a (smaller, larger) tuple was passed straight through.
    """
    mid = (i + j) / 2
    radius = abs(i - j) / 2
    theta = np.linspace(0, np.pi, 200)
    ax.plot(mid + radius * np.cos(theta),
            direction * (offset + radius * np.sin(theta) * 0.9),
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
    Single-sided arc diagram with pairs color-coded against a reference:
      - green  solid : pair present in both (correct)
      - red    solid : pair only in `predicted` (wrong)
      - gray  dashed : pair only in `reference` (missed)

    Compact, but it stacks everything into one field. For the side-by-side
    "what the solver said vs. what ViennaRNA says" read, prefer
    `plot_mirrored_comparison`.

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
        _draw_arc(ax, i, j, color=MATCH_COLOR)
    for (i, j) in pred_pairs - ref_pairs:
        _draw_arc(ax, i, j, color=WRONG_COLOR)
    for (i, j) in ref_pairs - pred_pairs:
        _draw_arc(ax, i, j, color="#888888", style="--", alpha=0.7, lw=1.3)

    legend = [
        Line2D([0], [0], color=MATCH_COLOR, lw=2, label="Correct pair"),
        Line2D([0], [0], color=WRONG_COLOR, lw=2, label="Wrong pair (predicted only)"),
        Line2D([0], [0], color="#888888", lw=2, linestyle="--", label="Missed pair (true only)"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=8, frameon=False)


def plot_mirrored_comparison(ax, sequence: str, predicted: str, reference: str,
                             title: str = "", predicted_label: str = "Predicted",
                             reference_label: str = "ViennaRNA MFE"):
    """
    Mirrored "before/after" arc diagram around one shared sequence backbone:

      above the backbone : the pairs `predicted` contains
      below the backbone : the pairs `reference` contains

    Color coding is shared by both fields:
      - green  solid : pair in both, drawn on both sides (so it looks symmetric)
      - red    solid : predicted only, drawn above with nothing facing it
      - gray  dashed : reference only, drawn below with nothing facing it

    A perfect prediction is therefore a perfectly symmetric figure, and every
    asymmetry is exactly one disagreement.

    Raises ValueError if `predicted`/`reference`/`sequence` lengths disagree.
    """
    if not (len(sequence) == len(predicted) == len(reference)):
        raise ValueError(
            f"Length mismatch: sequence={len(sequence)}, "
            f"predicted={len(predicted)}, reference={len(reference)}"
        )

    pred_pairs = pairs_from_structure(predicted)
    ref_pairs = pairs_from_structure(reference)
    matched = pred_pairs & ref_pairs
    wrong = pred_pairs - ref_pairs
    missed = ref_pairs - pred_pairs

    n = len(sequence)
    up, down = _field_extents(pred_pairs, ref_pairs)
    top, bottom = up + _PAD_TOP, down + _PAD_BOTTOM

    # Note the y-limits are deliberately asymmetric: a matched pair is still
    # drawn at +r and -r, so it stays mirror-symmetric about the backbone,
    # but neither field reserves room it isn't using.
    ax.set_xlim(-1.2, n + 0.2)
    ax.set_ylim(-bottom, top)
    ax.axis("off")

    _draw_shared_backbone(ax, sequence)

    # Upper field: what the solver produced.
    for (i, j) in matched:
        _draw_arc(ax, i, j, MATCH_COLOR, direction=1, offset=_BAND, lw=1.9)
    for (i, j) in wrong:
        _draw_arc(ax, i, j, WRONG_COLOR, direction=1, offset=_BAND, lw=1.9)

    # Lower field: the reference structure.
    for (i, j) in matched:
        _draw_arc(ax, i, j, MATCH_COLOR, direction=-1, offset=_BAND, lw=1.9)
    for (i, j) in missed:
        _draw_arc(ax, i, j, MISSED_COLOR, direction=-1, offset=_BAND,
                  style="--", alpha=0.85, lw=1.7)

    # Which field is which.
    ax.text(-1.1, up + 0.38, f"{predicted_label}  (above)", fontsize=8.5,
            fontweight="bold", color="#333333", ha="left", va="center")
    ax.text(-1.1, -(down + 0.38), f"{reference_label}  (below)", fontsize=8.5,
            fontweight="bold", color="#333333", ha="left", va="center")

    # An empty field is a real result, not a rendering failure, so say so.
    if not pred_pairs:
        ax.text((n - 1) / 2, up * 0.55, "no base pairs predicted",
                fontsize=8.5, style="italic", color="#999999", ha="center", va="center")
    if not ref_pairs:
        ax.text((n - 1) / 2, -down * 0.55, "no base pairs in the reference structure",
                fontsize=8.5, style="italic", color="#999999", ha="center", va="center")

    # Dot-bracket strings, so the figure carries the raw answer too.
    width = max(len(predicted_label), len(reference_label))
    ax.text(-1.1, -bottom + 0.52, f"{predicted_label:>{width}} : {predicted}", fontsize=7.5,
            family="monospace", color="#555555", ha="left", va="center")
    ax.text(-1.1, -bottom + 0.22, f"{reference_label:>{width}} : {reference}", fontsize=7.5,
            family="monospace", color="#555555", ha="left", va="center")

    if title:
        ax.set_title(title, fontsize=11.5, fontweight="bold", pad=12)

    legend = [
        Line2D([0], [0], color=MATCH_COLOR, lw=2.2, label="Match (both)"),
        Line2D([0], [0], color=WRONG_COLOR, lw=2.2, label="Predicted only (wrong)"),
        Line2D([0], [0], color=MISSED_COLOR, lw=2.2, linestyle="--", label="Reference only (missed)"),
    ]
    counts = (f"{len(matched)} matched  |  {len(pred_pairs)} predicted  |  "
              f"{len(ref_pairs)} reference")
    ax.legend(handles=legend, loc="upper right", fontsize=8, frameon=False,
              title=counts, title_fontsize=8, alignment="left")


def compare_figure(sequence: str, predicted: str, reference: str,
                   title: Optional[str] = None, figsize=None,
                   predicted_label: str = "Predicted",
                   reference_label: str = "ViennaRNA MFE",
                   mirrored: bool = True):
    """Convenience wrapper: builds a single-panel comparison figure and returns it.

    `mirrored=True` (default) uses the two-field before/after layout;
    `mirrored=False` falls back to the compact single-field diagram.
    """
    if figsize is None:
        n = len(sequence)
        width = max(7.0, 0.40 * n + 4.0)
        if mirrored:
            # Height follows the actual vertical span, so panels whose fields
            # are short (or empty) don't come out mostly white space.
            up, down = _field_extents(pairs_from_structure(predicted),
                                      pairs_from_structure(reference))
            span = (up + _PAD_TOP) + (down + _PAD_BOTTOM)
            figsize = (width, max(3.8, 0.42 * span + 1.0))
        else:
            figsize = (width, 4.2)

    fig, ax = plt.subplots(figsize=figsize)
    if mirrored:
        plot_mirrored_comparison(
            ax, sequence, predicted, reference,
            title or f"{sequence}: {predicted_label} vs. {reference_label}",
            predicted_label=predicted_label, reference_label=reference_label,
        )
    else:
        plot_comparison(ax, sequence, predicted, reference,
                        title or f"{sequence}: predicted vs. reference")
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Quick self-test on the two documented cases from this project: one exact
    # match, one known formulation-gap mismatch. The `predicted` strings are the
    # structures the CVaR-VQE solver actually produces for these sequences (see
    # notebook 04); the reference side is folded live rather than pasted in, so
    # this self-test cannot drift away from ViennaRNA.
    import RNA

    cases = [
        ("GCGCAUACGC", "(((....)))"),
        ("GCAUCGUAGC", "((.(...)))"),
    ]

    fig, axes = plt.subplots(1, len(cases), figsize=(15, 5.4))
    for ax, (seq, predicted) in zip(axes, cases):
        reference, _ = RNA.fold(seq)
        plot_mirrored_comparison(ax, seq, predicted, reference,
                                 title=f"{seq} -- CVaR-VQE vs. ViennaRNA MFE",
                                 predicted_label="CVaR-VQE")
    plt.tight_layout()
    plt.savefig("visualize_selftest.png", dpi=150)
    print("Self-test figure saved to visualize_selftest.png")
