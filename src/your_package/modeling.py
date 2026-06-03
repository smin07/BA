"""Model helpers for regression-based analyses."""

from __future__ import annotations


def fit_ols_model(y, X):
    """Fit an OLS model with an intercept and return the fitted result."""

    from statsmodels.regression.linear_model import OLS
    from statsmodels.tools import add_constant

    design_matrix = add_constant(X, prepend=True)
    return OLS(y, design_matrix).fit()
