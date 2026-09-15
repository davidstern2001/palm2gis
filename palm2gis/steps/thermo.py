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
def fill_profile_gaps(theta_bar):
    """Gap-fill a horizontally averaged theta profile; report the gaps.

    Returns (filled, n_filled). A level lying entirely inside terrain or
    buildings has no air in it, so its horizontal mean is NaN - routine in
    a coarse domain, where the lowest output levels can sit below the
    terrain across the whole cut. Those levels are interpolated from the
    nearest valid ones (flat extrapolation at the ends), because the
    hydrostatic integration has to step through them to reach the levels
    that do carry air.
    """
    th = np.asarray(theta_bar, dtype="f8").copy()
    good = np.isfinite(th) & (th > 0)
    if not good.any():
        raise ValueError("No valid potential temperature for the Exner "
                         "profile: every level of the cut is inside "
                         "terrain. Check the bounds z-range.")
    idx = np.arange(th.size)
    return np.interp(idx, idx[good], th[good]), int((~good).sum())


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
    th, _ = fill_profile_gaps(theta_bar)

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
# 4. TEMPERATURE SCALE
# ------------------------------
# Temperature variables by PALM naming, shared with the raw-input guard in
# georef.py: `theta*` (potential temperature), `tsurf*` / `t_surf*`
# (surface temperature), `ta*` (air temperature).
#
# TEMP_PREFIXES only IDENTIFIES temperature variables. It does NOT say
# what scale they are on. Assuming kelvin here was a real bug: PALM's
# reference tables give `ta` and `ta_2m*` in DEGREES CELSIUS, while
# `theta`, `tsurf*` and `t_surf` are kelvin. Subtracting 273.15 from an
# already-Celsius `ta` produced a voxel layer reading -253 degC. The
# scale now comes from the `units` attribute; see temperature_scale().
#
# Kept deliberately IDENTICAL to palm_postproc/steps/coord.py. The two
# tools are separate packages by design and share no code, so this block
# is the one thing that has to be edited in both places at once.
K0 = 273.15
TEMP_PREFIXES = ("theta", "tsurf", "t_surf", "ta")

# Potential temperature is exempt from `celsius`. It is a kelvin-defined
# quantity and "potential temperature in degC" is a confusing thing to
# hand to anyone; `celsius` governs actual temperatures only.
CELSIUS_EXEMPT_PREFIXES = ("theta",)

# Unit strings, lowercased and stripped. The prefix entries matter: PALM
# writes its units into a fixed-length character buffer and real files
# come out carrying "degree_" - truncated mid-word - so an exact-match
# list rejects the very data this is meant to handle.
KELVIN_UNITS = ("k", "kelvin", "degree_k", "degrees_k", "deg_k")
CELSIUS_UNITS = ("degc", "deg_c", "degree_c", "degrees_c",
                 "celsius", "\N{DEGREE SIGN}c", "c")
CELSIUS_PREFIXES = ("degree", "degrees", "deg")


class TemperatureUnitError(RuntimeError):
    """Raised when a temperature variable's scale cannot be determined."""


def is_temperature(name):
    """True when *name* is a temperature variable by PALM naming."""
    return str(name).lower().startswith(TEMP_PREFIXES)


def is_celsius_exempt(name):
    """True when *name* stays in kelvin regardless of the celsius setting."""
    return str(name).lower().startswith(CELSIUS_EXEMPT_PREFIXES)


def temperature_scale(units, name="", override=None):
    """Return 'K' or 'C' for a temperature variable.

    Resolution order: explicit config override, then the `units` string.
    There is deliberately NO fallback that guesses from the values. A
    temperature silently off by 273.15 still looks like a temperature,
    which is exactly how the -253 degC voxels survived a release; failing
    loudly is cheaper than the alternative.

    Raises TemperatureUnitError when neither source resolves.
    """
    if override:
        for key, scale in override.items():
            if str(name).lower() == str(key).lower():
                s = str(scale).strip().upper()[:1]
                if s in ("K", "C"):
                    return s
                raise TemperatureUnitError(
                    f"temperature_units override for '{name}' is {scale!r}; "
                    f"expected 'K' or 'C'")

    u = str(units or "").strip().lower().replace(" ", "")
    if u in KELVIN_UNITS:
        return "K"
    if u in CELSIUS_UNITS or u.startswith(CELSIUS_PREFIXES):
        return "C"

    raise TemperatureUnitError(
        f"cannot determine the temperature scale of '{name}': units="
        f"{units!r}. Add an entry to `temperature_units` in the config "
        f"(e.g. {{{name or 'var'}: C}}) to state it explicitly.")


def convert_temperature(vals, source, target, fill=None):
    """Put *vals* on the *target* scale, converting only if needed.

    `source` and `target` are 'K' or 'C'. Returns the input unchanged
    (as a float array) when they match, so an already-Celsius PALM
    variable is never touched.
    """
    if source == target:
        return np.array(vals, dtype="f8", copy=True)
    if source == "K":
        return kelvin_to_celsius(vals, fill=fill)
    return celsius_to_kelvin(vals, fill=fill)


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


def celsius_to_kelvin(vals, fill=None):
    """Add 273.15 to the VALID entries of *vals* only.

    The mirror of kelvin_to_celsius, needed now that `celsius: false`
    means "emit kelvin" rather than "do nothing": PALM's own `ta`
    arrives in degC, so producing a kelvin output requires converting it.
    """
    out = np.array(vals, dtype="f8", copy=True)
    good = np.isfinite(out)
    if fill is not None:
        good &= out != fill
    out[good] += K0
    return out


# ------------------------------
# 5. ABSOLUTE HUMIDITY
# ------------------------------
# Virtual temperature factor Rv/Rd - 1. Moist air is lighter than dry air
# at the same T and p; ignoring it overstates the density by ~0.6 % per
# 10 g/kg of vapour.
EPS_V = 0.608


def humidity_scale(units):
    """Factor that turns PALM's `q` into kg/kg, from its `units` string.

    PALM writes q in kg/kg. A g/kg label is honoured; anything unreadable
    is taken as kg/kg, which the caller sanity-checks against the values.
    """
    u = str(units or "").lower().replace(" ", "")
    if u.startswith(("g/kg", "gkg-1", "gkg^-1")):
        return 1.0e-3
    return 1.0


def pressure_from_exner(pi):
    """Absolute pressure [Pa] from the Exner function."""
    return P0 * np.asarray(pi, dtype="f8") ** (1.0 / KAPPA)


def absolute_humidity(q, ta_k, p):
    """Water vapour density [g/m3] from specific humidity [kg/kg], air
    temperature [K] and absolute pressure [Pa]:

        rho_v = q * rho,  rho = p / (Rd * Tv),  Tv = T * (1 + 0.608 q)

    Arrays broadcast, so `p` may be a (nz, 1, 1) profile against
    (nz, ny, nx) fields. Masked inputs stay masked.
    """
    tv = ta_k * (1.0 + EPS_V * q)
    return 1.0e3 * q * p / (RD * tv)

