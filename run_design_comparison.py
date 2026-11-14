"""Run sensor comparative design: m=3..7 near 1550 nm; m=4 radius scaling."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mie_sensor import coefficient
from sensor_design import compare_modes, compare_radii

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True)


def save_csv(name, rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def run():
    modes=compare_modes()
    radii=compare_radii()
    rows_mode=[];rows_cal=[];rows_radius=[]
    for a in modes:
        rows_mode.append({
            'm':a.m,'radius_nm_for_peak_1550':f'{a.radius_nm:.8f}',
            'n_env_reference':f'{a.envs[0]:.2f}',
            'peak_reference_nm':f'{a.peak_nm[0]:.7f}',
            'Q_reference':f'{a.q[0]:.7f}',
            'peak_sensitivity_fit_nm_per_RIU':f'{a.peak_slope_nm_per_riu:.7f}',
            'pole_sensitivity_fit_nm_per_RIU':f'{a.pole_slope_nm_per_riu:.7f}',
            'peak_FWHM_reference_nm':f'{a.fwhm_nm[0]:.7f}',
            'peak_FOM_reference_per_RIU':f'{a.peak_fom_per_riu:.7f}',
            'local_pole_sensitivity_reference_nm_per_RIU':f'{a.pole_local_sensitivity_nm_per_riu:.7f}'})
        for i,ne in enumerate(a.envs):
            rows_cal.append({'m':a.m,'radius_nm':f'{a.radius_nm:.6f}',
                             'n_env':f'{ne:.2f}','spectral_peak_nm':f'{a.peak_nm[i]:.8f}',
                             'spectral_FWHM_nm':f'{a.fwhm_nm[i]:.8f}',
                             'pole_wavelength_nm':f'{a.pole_wavelength_nm[i]:.8f}',
                             'pole_Q':f'{a.q[i]:.8f}'})
    for a in radii:
        rows_radius.append({'m':a.m,'radius_nm':f'{a.radius_nm:.1f}',
                            'n_env_reference':f'{a.envs[0]:.2f}',
                            'peak_reference_nm':f'{a.peak_nm[0]:.8f}',
                            'Q_reference':f'{a.q[0]:.8f}',
                            'peak_sensitivity_nm_per_RIU':f'{a.peak_slope_nm_per_riu:.8f}',
                            'peak_FWHM_nm':f'{a.fwhm_nm[0]:.8f}',
                            'peak_FOM_per_RIU':f'{a.peak_fom_per_riu:.8f}'})
    save_csv('mode_design_metrics.csv',rows_mode)
    save_csv('mode_design_calibration.csv',rows_cal)
    save_csv('radius_scaling_metrics.csv',rows_radius)

    # Mode-order tradeoff: Q, sensitivity, and FOM, same approx operating wavelength.
    q=np.array([a.q[0] for a in modes]);S=np.array([a.peak_slope_nm_per_riu for a in modes]);
    fom=np.array([a.peak_fom_per_riu for a in modes]);orders=np.array([a.m for a in modes])
    fig,axs=plt.subplots(1,2,figsize=(11.2,4.5))
    axs[0].plot(q,S,'o-',lw=1.7)
    axs[1].plot(q,fom,'s-',lw=1.7)
    for i,m in enumerate(orders):
        axs[0].annotate(f'm={m}',(q[i],S[i]),xytext=(4,6),textcoords='offset points',fontsize=9)
        axs[1].annotate(f'm={m}',(q[i],fom[i]),xytext=(4,6),textcoords='offset points',fontsize=9)
    axs[0].set(xlabel='Radiative Q (pole)',ylabel='Measured-style peak slope (nm/RIU)',title='Higher Q does not guarantee larger sensitivity')
    axs[1].set(xlabel='Radiative Q (pole)',ylabel=r'Spectral FOM (RIU$^{-1}$)',title='Narrower linewidth improves spectral FOM')
    for ax in axs: ax.grid(alpha=.25)
    fig.suptitle('Mode-order comparison: spectral peaks aligned at 1550 nm (n_env=1.33)',fontsize=12)
    fig.tight_layout();fig.savefig(OUT/'06_mode_order_q_sensitivity_fom.png',dpi=180);plt.close(fig)

    # All spectral peaks aligned initially; compare index-calibration slopes.
    fig,ax=plt.subplots(figsize=(8.6,5))
    for a in modes:
        ax.plot(a.envs,a.peak_nm,'o-',label=f'm={a.m}; R={a.radius_nm:.0f} nm; S={a.peak_slope_nm_per_riu:.1f}')
    ax.set(xlabel='Surrounding refractive index (RIU)',ylabel='Spectral peak wavelength (nm)',
           title='Mode-dependent index calibration (same initial peak: 1550 nm)')
    ax.grid(alpha=.25);ax.legend(fontsize=8.5)
    fig.tight_layout();fig.savefig(OUT/'07_mode_order_calibration.png',dpi=180);plt.close(fig)

    # Radius effect: S and FWHM scale together; exact-scale FOM stays constant.
    r=np.array([a.radius_nm for a in radii]);S_r=np.array([a.peak_slope_nm_per_riu for a in radii]);
    W_r=np.array([a.fwhm_nm[0] for a in radii]);F_r=np.array([a.peak_fom_per_riu for a in radii]);
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.5))
    r_ref=2  # index of R=500 nm in radius-sweep array
    axs[0].plot(r,S_r/S_r[r_ref],'o-',label='Peak sensitivity / S(500 nm)')
    axs[0].plot(r,W_r/W_r[r_ref],'s--',label='FWHM / FWHM(500 nm)')
    axs[0].set(xlabel='Cylinder radius (nm)',ylabel='Ratio to value at R=500 nm',
               title='Sensitivity and linewidth scale together')
    axs[0].legend(fontsize=8)
    axs[1].plot(r,F_r,'o-',label='Computed peak FOM')
    axs[1].axhline(F_r[r_ref],ls='--',lw=1.25,label='Scale-invariant FOM')
    axs[1].set(xlabel='Cylinder radius (nm)',ylabel=r'Spectral FOM (RIU$^{-1}$)',
               title='FOM unchanged by pure radius scaling',ylim=(0,2.4))
    axs[1].legend(fontsize=9)
    for ax in axs:ax.grid(alpha=.25)
    fig.suptitle('Radius scaling at m=4: Q and spectral FOM unchanged (nondispersive 2D model)',fontsize=11.3)
    fig.tight_layout();fig.savefig(OUT/'08_radius_scaling_invariance.png',dpi=180);plt.close(fig)

    # Spectral shape at nenv=1.33: compare FWHM on a common wavelength scale.
    wl=np.linspace(1400,1700,3001)
    fig,ax=plt.subplots(figsize=(9,5))
    for a in modes:
        xx=2*np.pi*a.radius_nm/wl
        y=np.abs(coefficient(a.m,xx,a.n_core,a.envs[0]))**2
        ax.plot(wl,y,lw=1.6,label=f'm={a.m}, Q={a.q[0]:.0f}, FWHM={a.fwhm_nm[0]:.1f} nm')
    ax.set(xlabel='Vacuum wavelength (nm)',ylabel=r'Partial-wave intensity $|a_m|^2$',
           title='Comparable operating wavelength, differing radiative linewidths',xlim=(1400,1700))
    ax.legend(fontsize=8.5);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'09_mode_order_spectra.png',dpi=180);plt.close(fig)

    print('Mode-order sensor design: all reference peaks aligned at 1550 nm')
    print('m | Radius (nm) | Q | Peak S (nm/RIU) | FWHM (nm) | Peak FOM (1/RIU)')
    for a in modes:
        print(f'{a.m} | {a.radius_nm:10.2f} | {a.q[0]:7.2f} | {a.peak_slope_nm_per_riu:15.2f} | {a.fwhm_nm[0]:9.2f} | {a.peak_fom_per_riu:17.2f}')
    print('\nRadius study, m=4:')
    for a in radii:
        print(f'R={a.radius_nm:6.0f} nm; S={a.peak_slope_nm_per_riu:8.3f}; FWHM={a.fwhm_nm[0]:8.3f}; Q={a.q[0]:8.3f}; FOM={a.peak_fom_per_riu:8.4f}')
    print('Saved outputs/06-09 PNG, mode_design_metrics.csv, mode_design_calibration.csv, radius_scaling_metrics.csv')

if __name__=='__main__':run()
