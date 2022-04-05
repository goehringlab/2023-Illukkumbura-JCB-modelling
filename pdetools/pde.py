import numpy as np


def diffusion(concs, dx=1, pad=True):
    """
    Simulate single diffusion time step

    :param concs: 1D array of concentrations across space
    :param dx: spatial distance between points
    :return:
    """

    if pad:
        concs_ = np.r_[concs[0], concs, concs[-1]]
        d = concs_[:-2] - 2 * concs_[1:-1] + concs_[2:]
    else:
        d = concs[:-2] - 2 * concs[1:-1] + concs[2:]
    return d / (dx ** 2)


def pdeRK(dxdt, X0, Tmax, deltat, t_eval, killfunc=None, stabilitycheck=False, maxstep=None):
    """

    Function for solving system of PDEs using adaptive Runge-Kutta method
    Adapted from Hubatsch et al., 2019 (see https://github.com/lhcgeneva/PARmodelling)

    :param dxdt: function, takes list of 1D arrays corresponding to concentrations over space, returns list of gradients
    :param X0: starting conditions, list of 1D arrays
    :param Tmax: maximum time
    :param deltat: initial timestep
    :param t_eval: array of timepoints to save
    :param killfunc: optional func, takes same input as dxdt, integration will terminate when func returns True
    :param stabilitycheck: if True, will terminate integration when system changes by less that 1% per minute

    Returns
    X: final state
    time: final time
    X_stored: saved states according to t_eval
    t_stored: times corresponding to X_stored. Will not be exactly the same as t_eval but should be close

    """

    # Adaptive step size parameters
    atol = 0.000001
    rtol = 0.000001

    # 5TH ORDER RK COEFFICIENTS for Dormand-Prince
    a21, a31, a32, a41, a42, a43 = 1 / 5, 3 / 40, 9 / 40, 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63 = 9017 / 3168, -355 / 33, 46732 / 5247
    a64, a65 = 49 / 176, -5103 / 18656
    a71, a72, a73, a74 = 35 / 384, 0, 500 / 1113, 125 / 192
    a75, a76 = -2187 / 6784, 11 / 84

    b1, b2, b3, b4, b5 = 35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784
    b6, b7 = 11 / 84, 0

    bs1, bs2, bs3, bs4 = 5179 / 57600, 0, 7571 / 16695, 393 / 640
    bs5, bs6, bs7 = -92097 / 339200, 187 / 2100, 1 / 40

    # Set up
    X = X0
    testsoln = dxdt(X)
    nvars = len(testsoln)
    t_stored = []
    X_stored = [[] for _ in range(nvars)]
    n_stored_times = 0
    if maxstep is None:
        maxstep = 10 * deltat

    # Run
    time = 0
    updated = None
    terminate = False
    while not terminate:
        if time > Tmax:
            terminate = True

        # Calculate increments for RK45
        if (time == 0) or not updated:
            X1 = dxdt(X)
        else:
            X1 = X7

        X2 = dxdt([X[i] + deltat * (a21 * X1[i]) for i in range(nvars)])
        X3 = dxdt([X[i] + deltat * (a31 * X1[i] + a32 * X2[i]) for i in range(len(X))])
        X4 = dxdt([X[i] + deltat * (a41 * X1[i] + a42 * X2[i] + a43 * X3[i]) for i in range(nvars)])
        X5 = dxdt([X[i] + deltat * (a51 * X1[i] + a52 * X2[i] + a53 * X3[i] + a54 * X4[i]) for i in range(nvars)])
        X6 = dxdt([X[i] + deltat * (a61 * X1[i] + a62 * X2[i] + a63 * X3[i] + a64 * X4[i] + a65 * X5[i]) for i in
                   range(nvars)])
        X7 = dxdt([X[i] + deltat * (a71 * X1[i] + a73 * X3[i] + a74 * X4[i] + a75 * X5[i] + a76 * X6[i]) for i in
                   range(nvars)])

        # Update concentrations using A1-A6 and P1-P6, coefficient for A7 and P7 is 0.
        Xn_new = [X[i] + deltat * (b1 * X1[i] + b3 * X3[i] + b4 * X4[i] + b5 * X5[i] + b6 * X6[i]) for i in
                  range(nvars)]  # b2/7=0

        # Compute difference between fourth and fifth order
        deltaXnerr = [max(abs(
            (b1 - bs1) * X1[i] + (b3 - bs3) * X3[i] + (b4 - bs4) * X4[i] + (b5 - bs5) * X5[i] + (b6 - bs6) * X6[
                i] - bs7 * X7[i])) for i in range(nvars)]  # b7 is zero

        # Get maximum concentrations for An and Pn
        yXn = [np.maximum(max(abs(Xn_new[i])), max(abs(X[i]))) for i in range(nvars)]

        # Get error scale, combining relative and absolute error
        scaleXn = [atol + yXn[i] * rtol for i in range(nvars)]

        # Compute total error as norm of maximum errors for each species scaled by the error scale
        errs = [(deltaXnerr[i] / scaleXn[i]) ** 2 for i in range(nvars)]
        totalerror = np.sqrt(sum(errs) / nvars)

        # Compute new timestep
        # sometimes see "RuntimeWarning: divide by zero encountered in double_scalars". Need to look into
        dtnew = 0.8 * deltat * abs(1 / totalerror) ** (1 / 5)

        # Upper and lower bound for timestep to avoid changing too fast
        if dtnew > maxstep:
            dtnew = maxstep
        elif dtnew < deltat / 5:
            dtnew = deltat / 5

        # Compute max percentage change
        change = max(max(abs(X[i] - Xn_new[i]) / Xn_new[i]) * (60 / dtnew) for i in range(nvars))

        # Set timestep for next round
        deltat = dtnew

        # Accept step if error is on the order of error scale or below
        if totalerror < 1:
            time += deltat
            X = Xn_new
            updated = True

            # Store
            if sum(time > t_eval) > n_stored_times:
                t_stored.append(time)
                for i in range(nvars):
                    X_stored[i].append(X[i])
                n_stored_times += 1

            # Kill function
            if killfunc is not None:
                brk = killfunc(X)
                if brk:
                    t_stored.append(time)
                    for i in range(nvars):
                        X_stored[i].append(X[i])
                    break

            # Check stability:
            if stabilitycheck:
                if change < 0.001:
                    t_stored.append(time)
                    for i in range(nvars):
                        X_stored[i].append(X[i])
                    break

        else:
            updated = False

    t_stored_ = np.array(t_stored)
    X_stored_ = [np.array(X_stored[i]) for i in range(nvars)]
    return X, time, X_stored_, t_stored_
