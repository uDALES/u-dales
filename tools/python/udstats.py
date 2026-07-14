"""Pure statistical helpers for uDALES field/facet data.

Stateless array-in/array-out functions extracted from UDBase so they can be
used and tested independently of case state. UDBase keeps thin static-method
wrappers for backward compatibility.
"""
from __future__ import annotations

from typing import Optional

import numpy as np


def time_average(var: np.ndarray, other: Optional[np.ndarray] = None):
    """
    Time-average variables.
    
    Parameters
    ----------
    var : ndarray
        Variable to average. Time must be the last dimension.
    other : ndarray, optional
        Second variable (same shape as ``var``) to compute covariance with.
    
    Returns
    -------
    tuple
        - (mean, variance) when only ``var`` is provided
        - (x_mean, y_mean, covariance) when ``other`` is provided
    
    Notes
    -----
    Variance/covariance are computed by delegating to ``merge_stat`` with
    zero instantaneous contributions.
    """
    X = np.asarray(var)
    n = X.shape[-1]
    
    if other is None:
        zeros = np.zeros_like(X)
        return merge_stat(X, zeros, n)
    
    Y = np.asarray(other)
    zeros = np.zeros_like(X)
    return merge_stat(X, Y, zeros, n)

def merge_stat(X: np.ndarray, *args, Y: Optional[np.ndarray] = None,
               XpXp: Optional[np.ndarray] = None,
               XpYp: Optional[np.ndarray] = None):
    """
    Merge short-term statistics into longer windows.

    Supported call patterns:
      - merge_stat(X, n)
      - merge_stat(X, XpXp, n)                  
      - merge_stat(X, Y, XpYp, n)               
      - merge_stat(X, n, XpXp=array)            # single variable mean/variance
      - merge_stat(X, n, Y=array, XpYp=array)   # two variables mean/cov

    Parameters
    ----------
    X : ndarray
        First variable. Final dimension is time or short-term windows.
    args : tuple
        Positional parameters parsed according to the patterns above.
    Y : ndarray, optional
        Second variable (same trailing dimension as ``X``).
    XpXp : ndarray, optional
        Variance contribution for ``X`` aligned with ``X``.
    XpYp : ndarray, optional
        Covariance contribution for ``X`` and ``Y`` aligned with ``X``.

    Returns
    -------
    tuple | ndarray
        ``Xmean`` if only ``X`` provided;
        ``Xmean, var`` if ``XpXp`` given;
        ``Xmean, Ymean, cov`` if ``Y`` provided (and optionally ``XpYp``).

    Notes
    -----
    Discards the oldest samples that do not fill a complete window so the
    most recent data is retained. Variance/covariance combine the mean of
    short-window contributions with the variance/covariance of the short
    means inside each merged window.
    """
    # Parse positional arguments to support both MATLAB and Python styles
    n = None
    X = np.asarray(X)
    if len(args) == 1:
        n = int(args[0])
    elif len(args) == 2 and Y is None:
        # MATLAB style: merge_stat(X, XpXp, n)
        XpXp = np.asarray(args[0])
        n = int(args[1])
    elif len(args) == 3:
        # MATLAB style: merge_stat(X, Y, XpYp, n)
        Y = np.asarray(args[0])
        XpYp = np.asarray(args[1])
        n = int(args[2])
    else:
        raise ValueError("merge_stat expects 1, 2, or 3 positional arguments after X")

    if n <= 0:
        raise ValueError("n must be positive")
    if X.shape[-1] < n:
        raise ValueError("Not enough samples to form a single merged window")

    nwin = X.shape[-1] // n
    start = X.shape[-1] - nwin * n  # discard oldest incomplete window
    X_use = X[..., start:]
    X_group = X_use.reshape(*X.shape[:-1], nwin, n)
    Xmean = X_group.mean(axis=-1)

    if Y is None:
        if XpXp is None:
            return Xmean

        XpXp = np.asarray(XpXp)
        if XpXp.shape[-1] != X.shape[-1]:
            raise ValueError("XpXp must match X shape in the last dimension")
        XpXp_use = XpXp[..., start:]
        XpXp_group = XpXp_use.reshape(*XpXp.shape[:-1], nwin, n)
        within = XpXp_group.mean(axis=-1)
        between = ((X_group - Xmean[..., None]) ** 2).mean(axis=-1)
        return Xmean, within + between

    Y = np.asarray(Y)
    if Y.shape[-1] < n:
        raise ValueError("Y does not have enough samples to form a merged window")
    if Y.shape[-1] != X.shape[-1]:
        raise ValueError("X and Y must share the same length in the last dimension")

    Y_use = Y[..., start:]
    Y_group = Y_use.reshape(*Y.shape[:-1], nwin, n)
    Ymean = Y_group.mean(axis=-1)

    if XpYp is None:
        cov_within = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)
        between = 0.0
    else:
        XpYp = np.asarray(XpYp)
        if XpYp.shape[-1] != X.shape[-1]:
            raise ValueError("XpYp must match X and Y in the last dimension")
        XpYp_use = XpYp[..., start:]
        XpYp_group = XpYp_use.reshape(*XpYp.shape[:-1], nwin, n)
        cov_within = XpYp_group.mean(axis=-1)
        between = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)

    cov = cov_within + between
    return Xmean, Ymean, cov

def coarsegrain_field(var: np.ndarray, Lflt: np.ndarray,
                      xm: np.ndarray, ym: np.ndarray) -> np.ndarray:
    """
    Apply 2D periodic box filters to a 3D field.

    Parameters
    ----------
    var : ndarray, shape (nx, ny, nz)
        Field data (time or height on the last axis).
    Lflt : array-like
        Filter lengths (meters). Multiple lengths are allowed.
    xm, ym : array-like
        Grid coordinates in meters; must be 1D and roughly uniform.

    Returns
    -------
    ndarray
        Filtered field with shape (nx, ny, nz, n_filters).
    """
    var = np.asarray(var)
    if var.ndim != 3:
        raise ValueError("var must be 3D with shape (nx, ny, nz)")

    xm = np.asarray(xm).ravel()
    ym = np.asarray(ym).ravel()
    if xm.size < 2 or ym.size < 2:
        raise ValueError("xm and ym must contain at least two points")

    dx = float(np.mean(np.diff(xm)))
    dy = float(np.mean(np.diff(ym)))
    if dx <= 0 or dy <= 0:
        raise ValueError("Grid spacings must be positive")

    Lflt_arr = np.atleast_1d(Lflt)
    nx, ny, nz = var.shape
    out = np.empty((nx, ny, nz, len(Lflt_arr)))

    # Build filters using half-width convention (matches MATLAB implementation)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    di = np.minimum(ii, nx - ii)  # periodic distance in i
    dj = np.minimum(jj, ny - jj)  # periodic distance in j

    for i, L in enumerate(Lflt_arr):
        ngx = max(int(round((L / dx) / 2.0)), 1)
        ngy = max(int(round((L / dy) / 2.0)), 1)
        mask = (di <= ngx) & (dj <= ngy)
        kernel = mask.astype(float)
        kernel /= kernel.sum()
        k_hat = np.fft.fftn(kernel)

        for k in range(nz):
            v_hat = np.fft.fftn(var[:, :, k])
            out[:, :, k, i] = np.real(np.fft.ifftn(v_hat * k_hat))

    return out
