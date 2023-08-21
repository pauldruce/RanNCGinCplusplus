import numpy as np


def specvar(eigs, xmax, xstep):
    # This is just splitting the x-axis
    # ? (so xfit is the discretisation of t in the heat kernel?).
    xfit = np.linspace(0., xmax, num=xstep)

    # d is an array with entries being lists of all the eigenvalues for each matrix size in the list M
    d = eigs
    da = d  # Don't need to square, as this is for Laplacian
    ds3 = np.zeros(xstep)
    dv3 = np.zeros(xstep)
    for t in np.arange(xstep):
        x = xfit[t]
        # Boltzmann weight for spec dim/var
        sx = np.exp(-x * da)
        # slx= top of spectral dimensions t * \lambda^2 e^{-t*\lambda^2}
        slx = da * sx * x
        # slx2 = top of first term in spectral variances
        slx2 = da * da * sx * x * x
        # summing the top of spectral dimensions
        sxs = sx.sum()
        # sxs = spectral dimensions?
        # (Paul) changes the redefining of sx to just a new variable sd
        sd = slx.sum() / sxs
        # sx2 = first term in spectral variance?
        sx2 = slx2.sum() / sxs
        # dv3 = spectral variance
        dv3[t] = 2 * sx2 - 2 * sd * sd
        # ds3 = spectral dimension
        ds3[t] = (2 * sd)
    return ds3, dv3
