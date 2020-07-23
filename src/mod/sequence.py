import pandas as pd
import numpy as np
from collections import defaultdict
from pyrolite_meltsutil.util.tables import phasename
from pyrolite.util.spatial import levenshtein_distance


def get_appearance_sequence(
    phases,
    ignore=["bulk", "cumulate", "solid"],
    mode="descending",
    variable="temperature",
    flatten=False,
):
    """

    Could use a masking approach as below to get multiple-appearances of phases
    for weird PT paths.
    """
    if "experiment" in phases.columns:
        if len(phases.experiment.unique()) > 1:
            # multi-experiment dataframe, return dictionary indexed by experiment
            return {
                e: get_appearance_sequence(
                    phases.loc[phases.experiment == e],
                    ignore=ignore,
                    mode=mode,
                    variable=variable,
                )
                for e in phases.experiment.unique()
            }
    _phaseIDs = phases.phaseID.unique()
    sequence = defaultdict(list)
    if mode == "descending":
        reverse = True
        first_appearance = lambda phs: phases.loc[phases.phaseID == p, variable].max()
    else:
        reverse = False
        first_appearance = lambda phs: phases.loc[phases.phaseID == p, variable].min()

    for p in _phaseIDs:
        if p not in ignore and not pd.isnull(p):
            sequence[int(first_appearance(p))].append(p)
    appearances = sorted(
        [(k, v) for (k, v) in sequence.items()], key=lambda x: x[0], reverse=reverse,
    )
    if not flatten:
        return appearances
    else:
        return [p for (value, phs) in appearances for p in phs]


def get_assemblage_sequence(
    phases,
    ignore=["bulk", "cumulate", "solid"],
    mode="descending",
    variable="temperature",
):
    """
    Get change points for assemblages, indexed by a specific variable.
    """
    if "experiment" in phases.columns:
        if len(phases.experiment.unique()) > 1:
            # multi-experiment dataframe, return dictionary indexed by experiment
            return {
                e: get_assemblage_sequence(
                    phases.loc[phases.experiment == e],
                    ignore=ignore,
                    mode=mode,
                    variable=variable,
                )
                for e in phases.experiment.unique()
            }

    phase_seq = get_appearance_sequence(
        phases, flatten=True, ignore=ignore, mode=mode, variable=variable
    )

    mask = pd.DataFrame(index=phases.index.unique(), columns=phase_seq, dtype=bool)
    mask.loc[:, :] = False

    for p in mask.columns:
        subdf = phases.loc[phases.phaseID == p, :]
        mask[p] = ~pd.isnull(subdf.reindex(index=mask.index)["phase"])

    assemblage_sequence = mask @ np.array(mask.columns)

    change_indexes = pd.Series(index=assemblage_sequence.index, dtype=bool)
    change_indexes[1:] = ~(
        assemblage_sequence[1:].values == assemblage_sequence[:-1].values
    )
    return assemblage_sequence, change_indexes


def sequence_distance(
    A,
    B=None,
    phasenames=True,
    ignore_trailing=False,
    ignore=["bulk", "cumulate", "solid"],
    mode="descending",
    variable="temperature",
):
    """
    Compute a distance between two sequences.

    Parameters
    -----------
    A : :class:`pandas.DataFrame`
        Experiment DataFrame to extract first sequence from.
    B : :class:`list` | :class:`pandas.DataFrame`
        List of minerals or second experiment DataFrame to extract sequence from.
    phasenames : :class:`bool`
        Whether to compare phase names rather than phase IDs (e.g. clinopyroxene_0).
    ignore_trailing : :class:`bool`
        Whether to clip sequences to the minimum length (i.e. compare only the first
        few items).
    ignore : :class:`list`
        Phases to ignore.
    mode : :class:`str`
        Ascending or descending.
    variable : :class:`str`
        Index varialbe. Typically temperature for melting or crystallisation
        experiements.

    Returns
    -------
    :class:`int`
        Distance between two sequences.
    """

    seq_A = get_appearance_sequence(
        A, ignore=ignore, mode=mode, variable=variable, flatten=True
    )
    if isinstance(B, list):
        seq_B = B
    elif isinstance(B, pd.DataFrame):
        seq_B = get_appearance_sequence(
            B, ignore=ignore, mode=mode, variable=variable, flatten=True
        )
    else:
        raise NotImplementedError

    if phasenames:
        seq_A = [phasename(a) for a in seq_A]
        seq_B = [phasename(b) for b in seq_B]

    if ignore_trailing:
        min_len = min(len(seq_A), len(seq_B))
        return seq_A, seq_B, levenshtein_distance(seq_A[:min_len], seq_B[:min_len])
    else:
        return seq_A, seq_B, levenshtein_distance(seq_A, seq_B)
