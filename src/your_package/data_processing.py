"""Data-processing helpers for the analysis pipeline."""

from __future__ import annotations


def drop_constant_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of *frame* without columns that contain one unique value."""

    import pandas as pd

    constant_columns = [
        column for column in frame.columns if frame[column].nunique(dropna=False) <= 1
    ]
    return frame.drop(columns=constant_columns)


def standardize_numeric_columns(
    frame: pd.DataFrame,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Return a standardized copy of *frame* for the selected numeric columns."""

    import pandas as pd

    standardized = frame.copy()
    target_columns = columns or standardized.select_dtypes(include="number").columns.tolist()

    for column in target_columns:
        series = standardized[column]
        deviation = series.std(ddof=0)
        if deviation == 0 or pd.isna(deviation):
            standardized[column] = 0.0
        else:
            standardized[column] = (series - series.mean()) / deviation

    return standardized
