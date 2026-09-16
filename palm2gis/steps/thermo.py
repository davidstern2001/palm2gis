"""Potential temperature -> actual air temperature.

PALM's prognostic variable is the potential temperature `theta`; the air
temperature is `T = theta * Pi`, where Pi is the Exner function
`(p/p0)^(Rd/cp)`.

The Pi that matters is PALM's OWN, because PALM converted between T and
theta with it (theta_surface = t_surface / Pi). PALM builds it once, at
initialisation, as a hydrostatic basic state of a CONSTANT-theta column:

    Pi(z) = (p_s / p0)^kappa - g * z / (cp * theta_ref)

p_s     - `surface_pressure` at z = 0 (the domain base, at origin_z)
theta_ref - the namelist `pt_surface`

It is linear in z and does not change during the run. Verified on real
output (holesovice_stern_adapt): ta / theta fits this line with a residual
below 0.0001 K over 1 km and every timestep, while the surface pressure in
the dynamic driver moved by several hPa. Integrating the CURRENT theta
profile instead, as earlier versions did, departed from PALM's Pi by about
0.3 K per km.

PALM's `p` output is the PERTURBATION pressure and is not used.

Where p_s and theta_ref come from is decided in georef.resolve_basic_state.
"""

import numpy as np

# ------------------------------
# 1. CONSTANTS
# ------------------------------
# PALM's own values (basic_constants_and_equations_mod), so the inversion
# of PALM's T <-> theta conversion is exact.
G = 9.81                 # gravitational acceleration [m/s2]
CP = 1005.0              # specific heat of dry air at constant pressure
RD = 287.0               # gas constant for dry air [J/(kg K)]
KAPPA = RD / CP          # ~0.286
P0 = 1.0e5               # reference pressure [Pa]
# PALM namelist defaults, used only when nothing else is known.
P_SURFACE_DEFAULT = 101325.0     # surface_pressure = 1013.25 hPa
PT_SURFACE_DEFAULT = 300.0       # pt_surface [K]


# ------------------------------
# 2. EXNER FUNCTION (PALM BASIC STATE)
# ------------------------------
def exner_basic_state(z, p_surface, theta_ref):
    """PALM's Exner function on the heights `z`.

    z         : heights [m] above the domain base (PALM's z = 0, at
                origin_z) - NOT absolute heights.
    p_surface : basic-state pressure at z = 0 [Pa].
    theta_ref : basic-state potential temperature (`pt_surface`) [K].
    """
    z = np.asarray(z, dtype="f8")
    return ((float(p_surface) / P0) ** KAPPA
            - G / CP * z / float(theta_ref))


def fit_basic_state(z, pi_bar):
    """Recover (p_surface [Pa], theta_ref [K], max residual [K]) from a
    horizontally averaged Exner profile `pi_bar` = mean(ta / theta) on `z`.

    A linear fit, because PALM's basic state is linear in z. The residual
    is expressed as a temperature error (|Pi - fit| * theta_ref), so a
    value near zero confirms that the file follows the model. Returns None
    when fewer than two levels are usable or Pi does not fall with z.
    """
    z = np.asarray(z, dtype="f8")
    pi_bar = np.asarray(pi_bar, dtype="f8")
    ok = np.isfinite(z) & np.isfinite(pi_bar) & (pi_bar > 0)
    if ok.sum() < 2:
        return None
    slope, pi_s = np.polyfit(z[ok], pi_bar[ok], 1)
    # A real basic state falls by g/(cp theta_ref) ~ 3e-5 per metre; a
    # slope that is zero to rounding means ta / theta is not one.
    if not slope < -1.0e-9 or not pi_s > 0:
        return None
    theta_ref = G / CP / -slope
    resid = float(np.max(np.abs(pi_bar[ok] - (pi_s + slope * z[ok]))))
    return P0 * pi_s ** (1.0 / KAPPA), theta_ref, resid * theta_ref


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
# (surface temperature - the prefix also covers palm2gis' own
# `t_surface`), `ta*` (air temperature).
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

