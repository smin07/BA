"""Plotting helpers for consistent figure styling."""

from __future__ import annotations


def shared_histogram_bins(*series, bins: int = 16):
    """Build one bin edge array that can be reused across multiple histograms."""

    import numpy as np

    values = [np.asarray(item).ravel() for item in series if np.asarray(item).size]
    if not values:
        return np.array([0.0, 1.0])

    combined = np.concatenate(values)
    minimum = float(np.nanmin(combined))
    maximum = float(np.nanmax(combined))

    if minimum == maximum:
        return np.array([minimum - 0.5, maximum + 0.5])

    return np.linspace(minimum, maximum, bins)
