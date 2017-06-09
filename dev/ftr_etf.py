#!/usr/bin/env python
# -*- coding: utf-8 -*-
import click
import numpy as np
from tqdm import tqdm

def simulate_slopes(filt, n, shape):
    """Simulate slopes"""
    phase_ft = np.random.randn(n, *shape) + 1j * np.random.randn(n, *shape)
    phase_p = np.empty_like(phase_ft)
    for i,phi in enumerate(tqdm(phase_ft)):
        phase_p[i,...] = filt(phase_ft[i,...])
    return phase_ft, phase_p

def plot_phases(axes, phase_ft, phase_ip, length=256, rate=1e3, label=None):
    """Plot phases."""
    from controltheory.periodogram import periodogram
    from controltheory.fourier import frequencies
    psd_ft = periodogram(phase_ft, length=256)
    psd_ip = periodogram(phase_ip, length=256)
    freq = frequencies(256, rate=1e3) # Simulate at 1kHz
    
    ax_psd, ax_etf = axes
    line, = ax_psd.plot(freq, psd_ip[:,10,10], ls="-", label=label)
    ax_psd.plot(freq, psd_ft[:,10,10], ls="--", color=line.get_color())
    ax_etf.plot(freq, psd_ip[:,10,10] / psd_ft[:,10,10], ls='-', color=line.get_color(), label=label)
    

@click.command()
def main():
    """Run a simulation of an ETF."""
    global plt
    import matplotlib.pyplot as plt
    from FTR.lqg import LQGFilter
    
    n = 1024 * 10
    shape = (32, 32)
    
    int_filt = LQGFilter.generate_integrator(0.2, 0.995, shape)
    kal_filt = LQGFilter.generate_kalman_filter(0.6, 0.995, shape)
    
    fig, axes = plt.subplots(2, 1, sharex=True)
    phase_ft, phase_ip = simulate_slopes(int_filt, n, shape)
    plot_phases(axes, phase_ft, phase_ip, label="Integrator")
    
    phase_ft, phase_ip = simulate_slopes(kal_filt, n, shape)
    plot_phases(axes, phase_ft, phase_ip, label="Kalman")
    
    ax_psd, ax_etf = axes
    ax_psd.set_xscale('log')
    ax_psd.set_yscale('log')
    ax_psd.set_xlim(1, 500)
    ax_psd.grid(True)
    
    ax_etf.set_yscale('log')
    ax_etf.grid(True)
    ax_etf.set_xlabel("Frequency (Hz)")
    ax_etf.legend()
    plt.show()

if __name__ == '__main__':
    main()