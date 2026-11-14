"""Run validated Mie resonance-index-sensing benchmark and save 5 figures/CSVs."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mie_sensor import (SensorGeometry, coefficient, denominator, track_mode,
                        pole_wavelength_nm, quality_factor, pole_linewidth_nm,
                        implicit_sensitivity_nm_per_riu, spectral_peak,
                        spectral_fwhm)

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True)
GEOM=SensorGeometry(radius_nm=500.0, core_index=3.0)
ENVS=np.round(np.linspace(1.33,1.39,7),2)
SEEDS={3:1.61-0.054j,4:2.028-0.031j,5:2.438-0.017j}


def write_csv(name, rows, columns):
    with (OUT/name).open('w',newline='',encoding='utf-8') as file:
        w=csv.DictWriter(file,fieldnames=columns);w.writeheader();w.writerows(rows)


def run():
    scans={}
    for m,seed in SEEDS.items():
        z=track_mode(m,ENVS,GEOM.core_index,seed)
        lam=pole_wavelength_nm(z,GEOM.radius_nm)
        Q=quality_factor(z)
        linewidth=pole_linewidth_nm(z,GEOM.radius_nm)
        S_local=np.array([implicit_sensitivity_nm_per_riu(m,p,GEOM.core_index,ne,GEOM.radius_nm)
                          for p,ne in zip(z,ENVS)])
        scans[m]={'z':z,'lambda':lam,'Q':Q,'linewidth':linewidth,'S':S_local,
                  'FOM':S_local/linewidth}
    base=scans[4]
    peaks=np.array([spectral_peak(4,ne,GEOM,lam,span_nm=120)[0]
                    for ne,lam in zip(ENVS,base['lambda'])])
    fwhms=np.array([spectral_fwhm(4,ne,GEOM,peak,window_nm=180)
                     for ne,peak in zip(ENVS,peaks)])
    fit_pole=np.polyfit(ENVS,base['lambda'],1)
    fit_peak=np.polyfit(ENVS,peaks,1)

    rows=[]
    for m,a in scans.items():
        for i,ne in enumerate(ENVS):
            d={'m':m,'n_env':f'{ne:.2f}',
               'pole_real_x':f'{a["z"][i].real:.12g}',
               'pole_imag_x':f'{a["z"][i].imag:.12g}',
               'lambda_pole_nm':f'{a["lambda"][i]:.8f}',
               'Q_rad':f'{a["Q"][i]:.8f}',
               'pole_linewidth_nm':f'{a["linewidth"][i]:.8f}',
               'pole_sensitivity_nm_RIU':f'{a["S"][i]:.8f}',
               'pole_FOM_RIU_inverse':f'{a["FOM"][i]:.8f}',
               'spectral_peak_nm':f'{peaks[i]:.8f}' if m==4 else '',
               'spectral_FWHM_nm':f'{fwhms[i]:.8f}' if m==4 else ''}
            rows.append(d)
    write_csv('sensor_pole_metrics.csv',rows,list(rows[0]))
    wl=np.linspace(1420,1690,2701)
    spectra=[]
    for ne in ENVS:
        y=np.abs(coefficient(4,2*np.pi*GEOM.radius_nm/wl,GEOM.core_index,ne))**2
        spectra.extend([{'n_env':f'{ne:.2f}','wavelength_nm':f'{w:.5f}',
                         'partial_scattering_intensity':f'{v:.9g}'} for w,v in zip(wl,y)])
    write_csv('sensor_spectra.csv',spectra,['n_env','wavelength_nm','partial_scattering_intensity'])

    # Figure 1: sensor signal under index changes
    fig,ax=plt.subplots(figsize=(9,5))
    for ne in [1.33,1.35,1.37,1.39]:
        y=np.abs(coefficient(4,2*np.pi*GEOM.radius_nm/wl,GEOM.core_index,ne))**2
        ax.plot(wl,y,lw=1.75,label=f'n_env={ne:.2f}')
    ax.set(xlabel='Vacuum wavelength (nm)',ylabel=r'Partial-wave intensity $|a_4|^2$',
           title='Refractive-index sensing: shift of the m=4 scattering peak',xlim=(1470,1640))
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'01_index_dependent_spectra.png',dpi=180);plt.close(fig)

    # Figure 2: extracted peak and pole wavelength versus ambient RI
    fig,ax=plt.subplots(figsize=(8.6,5))
    ax.plot(ENVS,base['lambda'],'o-',label='Complex-pole resonance wavelength')
    ax.plot(ENVS,peaks,'s-',label='Spectral |a4|² maximum')
    ax.plot(ENVS,np.polyval(fit_pole,ENVS),'--',lw=1,alpha=.75,
            label=f'Pole slope: {fit_pole[0]:.1f} nm/RIU')
    ax.plot(ENVS,np.polyval(fit_peak,ENVS),'--',lw=1,alpha=.75,
            label=f'Peak slope: {fit_peak[0]:.1f} nm/RIU')
    ax.set(xlabel='Surrounding refractive index (RIU)',ylabel='Resonance wavelength (nm)',
           title='Calibration: resonance pole versus observable spectral peak')
    ax.grid(alpha=.25);ax.legend(fontsize=9);fig.tight_layout()
    fig.savefig(OUT/'02_sensor_calibration.png',dpi=180);plt.close(fig)

    # Figure 3: poles in complex frequency plane for m3,m4,m5
    fig,ax=plt.subplots(figsize=(8.7,5))
    for m,a in scans.items():
        line=ax.plot(a['z'].real,a['z'].imag,'o-',label=f'm={m}')[0]
        ax.annotate(f'{ENVS[0]:.2f}',(a['z'][0].real,a['z'][0].imag),
                    xytext=(4,5),textcoords='offset points',fontsize=8,color=line.get_color())
        ax.annotate(f'{ENVS[-1]:.2f}',(a['z'][-1].real,a['z'][-1].imag),
                    xytext=(4,-12),textcoords='offset points',fontsize=8,color=line.get_color())
    ax.axhline(0,color='grey',ls='--',lw=.9)
    ax.set(xlabel=r'Re$(k_0R)$',ylabel=r'Im$(k_0R)$',
           title='Complex-frequency pole trajectories as surrounding index increases')
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'03_resonance_pole_trajectories.png',dpi=180);plt.close(fig)

    # Figure 4: Q vs sensitivity; m differs with same physical core and R
    fig,ax=plt.subplots(figsize=(8.4,5))
    for m,a in scans.items():
        ax.plot(a['Q'],a['S'],'o-',label=f'm={m}')
        ax.annotate(f'{ENVS[0]:.2f}',(a['Q'][0],a['S'][0]),
                    xytext=(4,5),textcoords='offset points',fontsize=8)
    ax.set(xlabel='Radiative Q from complex pole',ylabel='Local pole sensitivity (nm/RIU)',
           title='Resonance linewidth versus refractive-index sensitivity')
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'04_q_sensitivity_tradeoff.png',dpi=180);plt.close(fig)

    # Figure 5: spectral linewidth and FOM from actual intensity plus pole metrics
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5))
    axes[0].plot(ENVS,fwhms,'o-',label='Spectrum FWHM')
    axes[0].plot(ENVS,base['linewidth'],'s--',label='Pole λ/Q approximation')
    axes[0].set(xlabel='Surrounding index',ylabel='Linewidth (nm)',title='Two linewidth definitions')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.25)
    axes[1].plot(ENVS,fit_peak[0]/fwhms,'o-',label='Peak slope / spectral FWHM')
    axes[1].plot(ENVS,base['FOM'],'s--',label='Local pole S / (λ/Q)')
    axes[1].set(xlabel='Surrounding index',ylabel=r'FOM (RIU$^{-1}$)',title='Refractive-index sensing FOM')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.25)
    fig.suptitle('Sensor metrics: experimental-style readout versus pole diagnostics',fontsize=12)
    fig.tight_layout();fig.savefig(OUT/'05_linewidth_and_fom.png',dpi=180);plt.close(fig)

    summary=[
       {'metric':'core_index','value':f'{GEOM.core_index:.2f}','unit':'RIU'},
       {'metric':'radius','value':f'{GEOM.radius_nm:.1f}','unit':'nm'},
       {'metric':'ambient_RI_min','value':f'{ENVS.min():.2f}','unit':'RIU'},
       {'metric':'ambient_RI_max','value':f'{ENVS.max():.2f}','unit':'RIU'},
       {'metric':'channel','value':'4','unit':'azimuthal m'},
       {'metric':'pole_sensitivity_linear','value':f'{fit_pole[0]:.8f}','unit':'nm/RIU'},
       {'metric':'peak_sensitivity_linear','value':f'{fit_peak[0]:.8f}','unit':'nm/RIU'},
       {'metric':'pole_sensitivity_implicit_at_start','value':f'{base["S"][0]:.8f}','unit':'nm/RIU'},
       {'metric':'Q_rad_at_start','value':f'{base["Q"][0]:.8f}','unit':'-'},
       {'metric':'FWHM_at_start','value':f'{fwhms[0]:.8f}','unit':'nm'},
       {'metric':'FOM_peak_at_start','value':f'{fit_peak[0]/fwhms[0]:.8f}','unit':'RIU^-1'},
       {'metric':'max_abs_D_of_all_poles','value':f'{max(abs(denominator(m,z,GEOM.core_index,ne)) for m,a in scans.items() for z,ne in zip(a["z"],ENVS)):.4g}','unit':'-'},
    ]
    write_csv('sensor_summary.csv',summary,['metric','value','unit'])
    print('Optical resonance refractive-index sensor benchmark')
    print(f'Core n={GEOM.core_index:.2f}, radius={GEOM.radius_nm:.1f} nm, ambient n={ENVS.min():.2f}..{ENVS.max():.2f}')
    print(f'm=4 pole calibration sensitivity:  {fit_pole[0]:.3f} nm/RIU')
    print(f'm=4 spectral-peak sensitivity:    {fit_peak[0]:.3f} nm/RIU')
    print(f'm=4 Q at n_env={ENVS[0]:.2f}:             {base["Q"][0]:.3f}')
    print(f'm=4 spectral FWHM at start:      {fwhms[0]:.3f} nm')
    print(f'm=4 spectral peak FOM at start:  {fit_peak[0]/fwhms[0]:.3f} RIU^-1')
    print('Wrote outputs/01..05 PNG, sensor_pole_metrics.csv, sensor_spectra.csv, sensor_summary.csv')

if __name__=='__main__':
    run()
