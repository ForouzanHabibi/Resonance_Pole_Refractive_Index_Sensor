"""Validated mode-order and radius design study for an idealized 2D RI sensor.

Reuse the exact lossless TM_z dielectric-cylinder model in mie_sensor.py.
Each azimuthal m=3..7 is assigned a radius such that its *spectral* peak is
1550 nm at n_env=1.33. This permits comparison in the same operating band.
For fixed mode/index, varying radius scales wavelengths but leaves Q and
spectral FOM invariant in a nondispersive Maxwell model.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from mie_sensor import (SensorGeometry, pole, track_mode, pole_wavelength_nm,
                        quality_factor, pole_linewidth_nm, spectral_peak,
                        spectral_fwhm, implicit_sensitivity_nm_per_riu)

BASE_ENV=1.33
TARGET_PEAK_NM=1550.0
ENVIRONMENTS=np.round(np.linspace(1.33,1.39,7), 2)
MODE_SEEDS={3:1.61-0.054j,4:2.028-0.031j,5:2.438-0.017j,
            6:2.84-0.009j,7:3.2355-0.00425j}


@dataclass(frozen=True)
class ModeReport:
    m:int
    radius_nm:float
    n_core:float
    envs:np.ndarray
    poles:np.ndarray
    peak_nm:np.ndarray
    fwhm_nm:np.ndarray
    pole_wavelength_nm:np.ndarray
    q:np.ndarray
    peak_slope_nm_per_riu:float
    pole_slope_nm_per_riu:float
    peak_fom_per_riu:float
    pole_local_sensitivity_nm_per_riu:float


def spectral_window_nm(z, radius_nm):
    """Adapt spectral search/integration windows to a tracked mode's linewidth."""
    width=pole_linewidth_nm(z, radius_nm)
    return max(25., min(230.,2.*width)), max(45., min(330.,3.2*width))


def report_mode(m:int, radius_nm:float, n_core:float=3.0,
                environments=ENVIRONMENTS, seed:complex|None=None) -> ModeReport:
    """Independent real-axis peak fit and complex-pole fit for selected branch."""
    envs=np.asarray(environments, dtype=float)
    if len(envs)<3 or np.any(np.diff(envs)<=0) or np.any(envs<=0):
        raise ValueError('Need at least three ordered, distinct positive ambient indices')
    if m not in MODE_SEEDS and seed is None:
        raise ValueError('Provide a pole seed for this angular order')
    geom=SensorGeometry(radius_nm=float(radius_nm),core_index=float(n_core))
    poles=track_mode(m,envs,n_core, MODE_SEEDS[m] if seed is None else seed)
    pole_lam=pole_wavelength_nm(poles,geom.radius_nm)
    peaks=[]; widths=[]
    for p,ne,wl0 in zip(poles,envs,pole_lam):
        peak_span,fwhm_span=spectral_window_nm(p,geom.radius_nm)
        peak,_=spectral_peak(m,float(ne),geom,float(wl0),span_nm=peak_span)
        width=spectral_fwhm(m,float(ne),geom,peak,window_nm=fwhm_span)
        peaks.append(peak);widths.append(width)
    peaks=np.asarray(peaks);widths=np.asarray(widths)
    s_peak=float(np.polyfit(envs,peaks,1)[0]);s_pole=float(np.polyfit(envs,pole_lam,1)[0])
    return ModeReport(m,geom.radius_nm,n_core,envs,poles,peaks,widths,pole_lam,
                      quality_factor(poles),s_peak,s_pole,s_peak/widths[0],
                      implicit_sensitivity_nm_per_riu(m,poles[0],n_core,envs[0],geom.radius_nm))


def radius_for_target_spectral_peak(m:int,target_peak_nm:float=TARGET_PEAK_NM,
                                   n_core:float=3.,n_env:float=BASE_ENV)->float:
    """Choose R so the real-axis |a_m|^2 peak is at target_peak_nm.

    Exact scale symmetry means one direct peak evaluation suffices: lambda_peak
    is proportional to R when material indices are wavelength independent.
    """
    if target_peak_nm <= 0:
        raise ValueError('Target wavelength must be positive')
    seed=MODE_SEEDS[m]
    z=pole(m,n_env,n_core,seed)
    r_provisional=target_peak_nm*z.real/(2*np.pi)
    g=SensorGeometry(r_provisional,n_core)
    span,_=spectral_window_nm(z,r_provisional)
    peak,_=spectral_peak(m,n_env,g,target_peak_nm,span_nm=span)
    return float(r_provisional*target_peak_nm/peak)


def compare_modes(modes=(3,4,5,6,7),target_peak_nm=TARGET_PEAK_NM):
    reports=[]
    for m in modes:
        radius=radius_for_target_spectral_peak(m,target_peak_nm)
        report=report_mode(m,radius)
        if abs(report.peak_nm[0]-target_peak_nm)>2e-5:
            raise RuntimeError('Failed to align spectral peak to target')
        reports.append(report)
    return reports


def compare_radii(m=4,radii_nm=(350,425,500,575,650)):
    return [report_mode(m,float(r)) for r in radii_nm]
