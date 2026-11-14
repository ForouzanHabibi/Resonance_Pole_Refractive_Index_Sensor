"""Regression checks for near-1550-nm mode tradeoffs and exact radius scaling."""
import unittest
import numpy as np
from mie_sensor import SensorGeometry, coefficient, quality_factor, pole
from sensor_design import (BASE_ENV, TARGET_PEAK_NM, compare_modes,compare_radii,
                           report_mode, radius_for_target_spectral_peak)


class DesignComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modes=compare_modes()
        cls.radii=compare_radii()

    def test_reference_peaks_all_aligned(self):
        for r in self.modes:
            self.assertLess(abs(r.peak_nm[0]-TARGET_PEAK_NM),2e-5)

    def test_mode_residuals_are_small(self):
        from mie_sensor import denominator
        for report in self.modes:
            for p, ne in zip(report.poles,report.envs):
                self.assertLess(abs(denominator(report.m,p,report.n_core,ne)),1e-8)

    def test_all_tracked_modes_have_positive_Q_and_sensitivity(self):
        for a in self.modes:
            self.assertTrue(np.all(a.q>0))
            self.assertTrue(a.peak_slope_nm_per_riu>0)
            self.assertTrue(a.peak_fom_per_riu>0)

    def test_higher_q_is_not_highest_raw_sensitivity(self):
        a,b=self.modes[2],self.modes[-1]
        self.assertGreater(b.q[0],a.q[0])
        self.assertLess(b.peak_slope_nm_per_riu,a.peak_slope_nm_per_riu)
        self.assertGreater(b.peak_fom_per_riu,a.peak_fom_per_riu)

    def test_radius_scaling_keeps_Q_and_FOM(self):
        reference=self.radii[2]
        for rep in self.radii:
            self.assertAlmostEqual(rep.q[0],reference.q[0],delta=1e-8)
            self.assertAlmostEqual(rep.peak_fom_per_riu,reference.peak_fom_per_riu,delta=5e-4)

    def test_radius_scaling_wavelength_sensitivity_linewidth(self):
        reference=self.radii[2]
        for rep in self.radii:
            scale=rep.radius_nm/reference.radius_nm
            self.assertAlmostEqual(rep.peak_nm[0]/reference.peak_nm[0],scale,delta=1e-5)
            self.assertAlmostEqual(rep.peak_slope_nm_per_riu/reference.peak_slope_nm_per_riu,scale,delta=1e-4)
            self.assertAlmostEqual(rep.fwhm_nm[0]/reference.fwhm_nm[0],scale,delta=5e-4)

    def test_dimensionless_scattering_is_scale_invariant(self):
        # Same x = 2*pi*R/lambda, so exact coefficient is identical.
        r1,r2=350.,650.
        wl=np.array([1300.,1550.,1800.])
        c1=coefficient(4,2*np.pi*r1/wl,3.,1.33)
        c2=coefficient(4,2*np.pi*r2/(wl*r2/r1),3.,1.33)
        np.testing.assert_allclose(c1,c2,rtol=1e-12,atol=1e-12)

    def test_invalid_target_or_ambient_grid_rejected(self):
        with self.assertRaises(ValueError):radius_for_target_spectral_peak(4,-100.)
        with self.assertRaises(ValueError):report_mode(4,500,environments=[1.33,1.33,1.39])

if __name__=='__main__':unittest.main()
