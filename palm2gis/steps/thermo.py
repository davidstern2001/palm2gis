"""Potential temperature -> actual air temperature.

PALM's prognostic variable is the potential temperature `theta`; the air
temperature is `T = theta * Pi`, where Pi is the Exner function
`(p/p0)^(Rd/cp)`.

The pressure needed for Pi is NOT PALM's `p` output: that is the
PERTURBATION pressure (deviation from the hydrostatic basic state, order
tens of Pa), not the absolute pressure. Using it as absolute pressure
gives nonsense. Its actual effect on T is negligible anyway - 20 Pa out of
1e5 Pa shifts T by about 0.02 K.

Instead Pi is integrated hydrostatically. Combining hydrostatic balance
with the definition of theta gives an exact relation with no pressure
left in it:

    dPi/dz = -g / (cp * theta)

so Pi is obtained by integrating up the column from its surface value
`(p_surface/p0)^kappa`, using the horizontally averaged theta profile from
the data itself. To leading order this reduces to the familiar
`T ~ theta - (g/cp) * z`, i.e. about -0.98 K per 100 m.

Two caveats:

* The result is only as good as `surface_pressure`. T scales as
  (p/p0)^0.286, so a 10 hPa error moves the derived temperature by about
  0.85 K - a constant offset, harmless for spatial patterns and gradients,
  but it matters if you compare against measurements.
* This reproduces PALM's own `ta` closely near the surface but not
  exactly, because PALM builds its basic-state pressure from the INITIAL
  temperature profile. If you need exact agreement, ask PALM for `ta`
  directly (add it to `data_output`) rather than deriving it here.
"""

import numpy as np

# ------------------------------
# 1. CONSTANTS
# ------------------------------
G = 9.81                 # gravitational acceleration [m/s2]
CP = 1005.0              # specific heat of dry air at constant pressure
RD = 287.0               # gas constant for dry air [J/(kg K)]
KAPPA = RD / CP          # ~0.286
P0 = 1.0e5               # reference pressure [Pa]
# Dry adiabatic lapse rate [K/km]. A constant-theta column IS a dry
# adiabat, so the theta -> ta offset must shrink with height at this rate
# whatever the stratification - a cheap check that the profile is sane.
LAPSE_DRY = G / CP * 1000.0      # ~9.76 K/km


# ------------------------------
# 2. EXNER FUNCTION
# ------------------------------
def exner_profile(z, theta_bar, p_surface=101325.0):
    """Exner function on the levels `z`, integrated from the surface.

    z         : heights [m] above the level where p = p_surface
                (PALM's z = 0, i.e. the domain base at origin_z).
    theta_bar : horizontally averaged potential temperature [K] on those
                levels; gaps (all-terrain levels) are filled from the
                nearest valid level.
    p_surface : pressure at z = 0 [Pa].
    """
    z = np.asarray(z, dtype="f8")
    th = np.asarray(theta_bar, dtype="f8").copy()

    # fill non-finite entries (levels entirely inside terrain)
    good = np.isfinite(th) & (th > 0)
    if not good.any():
        raise ValueError("No valid potential temperature for the Exner "
                         "profile.")
    idx = np.arange(th.size)
    th = np.interp(idx, idx[good], th[good])

    # p_surface applies at z = 0 (the domain base), which is NOT necessarily
    # the lowest OUTPUT level: a bounds z-cut can start well above it. Anchor
    # the profile at z = 0 and descend/ascend to z[0] first, otherwise the
    # whole column is offset by G/CP * z[0] - about 0.04 K for a cut starting
    # at 4 m, but 1 K for one starting at 100 m.
    pi = np.empty_like(z)
    pi[0] = (p_surface / P0) ** KAPPA - G / CP * z[0] / th[0]
    for k in range(1, z.size):
        theta_mid = 0.5 * (th[k] + th[k - 1])
        pi[k] = pi[k - 1] - G / CP * (z[k] - z[k - 1]) / theta_mid
    return pi


def exner_simple(z, theta_ref=290.0, p_surface=101325.0):
    """Exner function assuming a constant reference theta.

    Used where no profile is available (the surface output has no 3D
    temperature field). The error against the integrated profile is well
    under 0.1 K over the first few hundred metres.
    """
    z = np.asarray(z, dtype="f8")
    return (p_surface / P0) ** KAPPA - G / CP * z / float(theta_ref)


# ------------------------------
# 3. CONVERSION
# ------------------------------
def theta_to_ta(theta, pi):
    """Air temperature from potential temperature and Exner function.

    `pi` broadcasts against `theta`: pass a (nz,) profile for a
    (nz, ny, nx) field, or a per-element array for scattered surfaces.
    """
    return np.asarray(theta) * np.asarray(pi)


# ------------------------------
# 4. KELVIN -> CELSIUS
# ------------------------------
# Temperature variables by PALM naming, shared with the raw-input guard in
# georef.py: `theta*` (potential temperature), `tsurf*` / `t_surf*`
# (surface temperature), `ta*` (air temperature).
K0 = 273.15
TEMP_PREFIXES = ("theta", "tsurf", "t_surf", "ta")


def is_temperature(name):
    """True when *name* is a temperature variable by PALM naming."""
    return str(name).lower().startswith(TEMP_PREFIXES)


def kelvin_to_celsius(vals, fill=None):
    """Subtract 273.15 from the VALID entries of *vals* only.

    Fill sentinels (== *fill*) and non-finite entries are left untouched,
    so NoData stays NoData instead of becoming ~-9726 degC. Returns a new
    float array; the input is not modified.

    All physics (the Exner/ta derivation) runs in kelvin - this is applied
    at WRITE time only, as the last operation before a value leaves the
    program.
    """
    out = np.array(vals, dtype="f8", copy=True)
    good = np.isfinite(out)
    if fill is not None:
        good &= out != fill
    out[good] -= K0
    return out
