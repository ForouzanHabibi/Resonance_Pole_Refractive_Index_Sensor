"""2D TM_z refractive-index sensor based on exact dielectric-cylinder Mie scattering.

An infinitely long, lossless, homogeneous dielectric cylinder of radius R sits
in an infinite lossless liquid (n_env). Nonmagnetic media have equal mu.
The incident field is E_z; the outgoing radial wave is H_m^(1).
Time convention: exp(-i omega t). All poles are in Im(k0 R) < 0.

Dimensionless complex vacuum wave number x = k0 R = (2 pi / lambda0) R.
The external and interior Bessel arguments are n_env*x and n_core*x.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import root, minimize_scalar
from scipy.special import jv, jvp, hankel1, h1vp


@dataclass(frozen=True)
class SensorGeometry:
    radius_nm: float = 500.0
    core_index: float = 3.0
    def __post_init__(self):
        if self.radius_nm <= 0 or self.core_index <= 0:
            raise ValueError('Radius and core refractive index must be positive')


def denominator(m: int, x, n_core: float, n_env: float):
    """Outgoing quasi-normal-mode condition for E_z cylinder scattering.

    D = n_env J_m(n_core x) H_m^(1)'(n_env x)
        - n_core J_m'(n_core x) H_m^(1)(n_env x).
    """
    if m < 0 or n_core <= 0 or n_env <= 0:
        raise ValueError('Invalid angular order or refractive index')
    z = np.asarray(x, dtype=np.complex128)
    c, e = n_core*z, n_env*z
    return n_env*jv(m,c)*h1vp(m,e) - n_core*jvp(m,c)*hankel1(m,e)


def coefficient(m: int, x, n_core: float, n_env: float):
    """Exact partial-wave a_m multiplying the outgoing H_m^(1) field.

    For a real-frequency lossless cylinder, |1+2*a_m| = 1, checking flux.
    Setting n_core=n_env makes a_m vanish (homogeneous space).
    """
    z = np.asarray(x, dtype=np.complex128)
    c, e = n_core*z, n_env*z
    N = n_core*jvp(m,c)*jv(m,e) - n_env*jvp(m,e)*jv(m,c)
    return N/denominator(m,z,n_core,n_env)


def pole(m: int, n_env: float, n_core: float, guess: complex) -> complex:
    """Continue a selected resonance with a local complex root search."""
    def f(v):
        d=denominator(m,complex(v[0],v[1]),n_core,n_env)
        return [float(np.real(d)),float(np.imag(d))]
    fit=root(f,[guess.real,guess.imag],tol=1e-12)
    result=complex(*fit.x)
    if (not np.isfinite(result) or result.real <= 0 or result.imag >= 0
        or abs(denominator(m,result,n_core,n_env)) > 1e-8):
        raise RuntimeError(f'Root search failed for m={m}, n_env={n_env}, seed={guess}')
    return result


def track_mode(m: int, n_env_values, n_core: float, first_guess: complex):
    values=np.asarray(n_env_values,dtype=float)
    if values.ndim!=1 or len(values)==0 or np.any(values<=0):
        raise ValueError('n_env_values must be a nonempty 1D positive array')
    guesses=complex(first_guess)
    out=[]
    for ne in values:
        guesses=pole(m,float(ne),n_core,guesses)
        out.append(guesses)
    out=np.array(out,dtype=np.complex128)
    if np.any(np.abs(np.diff(out))>0.3):
        raise RuntimeError('Mode tracking detected a large pole jump')
    return out


def pole_wavelength_nm(z, radius_nm: float):
    """Resonance wavelength defined from Re(omega), not Re(1/omega)."""
    return 2*np.pi*radius_nm/np.real(z)


def quality_factor(z):
    """Radiative Q for nondispersive, lossless cylinder with exp(-i omega t)."""
    return -np.real(z)/(2*np.imag(z))


def pole_linewidth_nm(z, radius_nm: float):
    """Narrow-line approximation Delta_lambda = lambda_res/Q (not fitted FWHM)."""
    return pole_wavelength_nm(z,radius_nm)/quality_factor(z)


def implicit_sensitivity_nm_per_riu(m:int, z:complex, n_core:float,
                                     n_env:float, radius_nm:float, h:float=1e-5):
    """Pole-based d(lambda_res)/dn_env via D(x,n)=0 and finite derivatives.

    dx/dn = -(partial D / partial n)/(partial D / partial x).
    Derivatives are evaluated on the complex analytically continued denominator.
    """
    Dne=(denominator(m,z,n_core,n_env+h)-denominator(m,z,n_core,n_env-h))/(2*h)
    Dx=(denominator(m,z+h,n_core,n_env)-denominator(m,z-h,n_core,n_env))/(2*h)
    if abs(Dx)<1e-11:raise RuntimeError('Near-degenerate mode: derivative unstable')
    dzdn=-Dne/Dx
    return float(-2*np.pi*radius_nm*np.real(dzdn)/np.real(z)**2)


def spectral_peak(m:int, n_env:float, geometry:SensorGeometry,
                  approx_center_nm:float, span_nm:float=100.0):
    """Locate local maximum of |a_m|^2 near a given tracked resonance."""
    if span_nm <=0:raise ValueError('span_nm must be positive')
    def score(wavelength_nm):
        x=2*np.pi*geometry.radius_nm/wavelength_nm
        return -float(abs(coefficient(m,x,geometry.core_index,n_env))**2)
    fit=minimize_scalar(score,bounds=(approx_center_nm-span_nm/2,
                                      approx_center_nm+span_nm/2),method='bounded',
                         options={'xatol':1e-9})
    if not fit.success:raise RuntimeError('Spectral peak optimization failed')
    if abs(fit.x-approx_center_nm)>0.47*span_nm:
        raise RuntimeError('Peak is too close to the window edge')
    return float(fit.x),float(-fit.fun)


def spectral_fwhm(m:int, n_env:float, geometry:SensorGeometry,
                  center_nm:float, window_nm:float=160.0, points:int=10001):
    """FWHM of partial-wave scattering |a_m|^2 at half the local peak maximum.

    The linewidth is measured on the actual spectrum, rather than equated to
    lambda/Q. The baseline is not subtracted; a crossing on each side is required.
    """
    wl=np.linspace(center_nm-window_nm/2,center_nm+window_nm/2,points)
    x=2*np.pi*geometry.radius_nm/wl
    y=np.abs(coefficient(m,x,geometry.core_index,n_env))**2
    peak_idx=int(np.argmax(y)); halfway=float(y[peak_idx])/2
    left=np.flatnonzero(y[:peak_idx] < halfway)
    right=np.flatnonzero(y[peak_idx+1:] < halfway)
    if not len(left) or not len(right):
        raise RuntimeError('No complete half-maximum width in spectral window')
    li=left[-1]; ri=peak_idx+right[0]+1
    left_cross=np.interp(halfway,[y[li],y[li+1]],[wl[li],wl[li+1]])
    right_cross=np.interp(halfway,[y[ri],y[ri-1]],[wl[ri],wl[ri-1]])
    return float(right_cross-left_cross)
