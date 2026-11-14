"""Physics checks, consistency checks, and sensing metric regressions."""
import unittest
import numpy as np
from mie_sensor import (SensorGeometry, denominator, coefficient, pole,
                        track_mode, pole_wavelength_nm, quality_factor,
                        pole_linewidth_nm, implicit_sensitivity_nm_per_riu,
                        spectral_peak, spectral_fwhm)

class MieSensorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geom=SensorGeometry(radius_nm=500,core_index=3)
        cls.p4=pole(4,1.33,3.0,2.028-0.031j)

    def test_index_matched_cylinder_does_not_scatter(self):
        x=np.linspace(.4,5,50)
        self.assertLess(float(np.max(np.abs(coefficient(3,x,1.33,1.33)))),1e-12)

    def test_lossless_partial_wave_unitarity(self):
        x=np.linspace(.4,5,100)
        for ne in (1.0,1.33,1.39):
            for m in (0,2,3,4,5):
                S=1+2*coefficient(m,x,3.0,ne)
                self.assertLess(float(np.max(np.abs(abs(S)-1))),3e-12)

    def test_poles_in_lower_half_plane(self):
        self.assertLess(self.p4.imag,0)
        self.assertGreater(self.p4.real,0)
        self.assertLess(abs(denominator(4,self.p4,3.0,1.33)),1e-9)

    def test_mode_tracking(self):
        env=np.arange(1.33,1.4,.01)
        z=track_mode(4,env,3.0,self.p4)
        self.assertEqual(len(z),len(env))
        self.assertTrue(np.all(np.diff(z.real)<0))
        for ne,p in zip(env,z):
            self.assertLess(abs(denominator(4,p,3.0,ne)),1e-8)

    def test_sensitivity_matches_finite_difference_of_tracked_pole(self):
        h=1e-4
        plus=pole(4,1.33+h,3.,self.p4)
        minus=pole(4,1.33-h,3.,self.p4)
        fd=(pole_wavelength_nm(plus,500)-pole_wavelength_nm(minus,500))/(2*h)
        predicted=implicit_sensitivity_nm_per_riu(4,self.p4,3,1.33,500)
        self.assertAlmostEqual(predicted,fd,delta=2e-3)

    def test_sensitivity_positive(self):
        sens=implicit_sensitivity_nm_per_riu(4,self.p4,3,1.33,500)
        self.assertGreater(sens,0)
        self.assertLess(sens,1000)

    def test_pole_q_positive(self):
        q=quality_factor(self.p4)
        self.assertGreater(q,10)
        self.assertAlmostEqual(q,self.p4.real/(-2*self.p4.imag))

    def test_pole_linewidth(self):
        linewidth=pole_linewidth_nm(self.p4,500)
        self.assertAlmostEqual(linewidth,pole_wavelength_nm(self.p4,500)/quality_factor(self.p4))

    def test_peak_within_telecom_band(self):
        center=pole_wavelength_nm(self.p4,500)
        peak,intensity=spectral_peak(4,1.33,self.geom,center)
        self.assertTrue(1510 < peak < 1590)
        self.assertAlmostEqual(intensity,1.0,delta=1e-8)

    def test_peak_shift_increases_with_ambient_index(self):
        z1=pole(4,1.39,3.,self.p4)
        p0,_=spectral_peak(4,1.33,self.geom,pole_wavelength_nm(self.p4,500))
        p1,_=spectral_peak(4,1.39,self.geom,pole_wavelength_nm(z1,500))
        self.assertGreater(p1-p0,4)
        self.assertLess(p1-p0,15)

    def test_spectral_fwhm_positive(self):
        pk,_=spectral_peak(4,1.33,self.geom,pole_wavelength_nm(self.p4,500))
        fwhm=spectral_fwhm(4,1.33,self.geom,pk)
        self.assertTrue(30<fwhm<80)

    def test_spectral_and_pole_width_not_assumed_equal(self):
        pk,_=spectral_peak(4,1.33,self.geom,pole_wavelength_nm(self.p4,500))
        fwhm=spectral_fwhm(4,1.33,self.geom,pk)
        pole_width=pole_linewidth_nm(self.p4,500)
        self.assertGreater(abs(fwhm-pole_width),.5)

    def test_nonpositive_index_rejected(self):
        with self.assertRaises(ValueError):denominator(4,1.0,3,0)
        with self.assertRaises(ValueError):SensorGeometry(radius_nm=-1)
        with self.assertRaises(ValueError):track_mode(4,[-1],3,self.p4)

    def test_peak_invalid_window_rejected(self):
        with self.assertRaises(ValueError):spectral_peak(4,1.33,self.geom,1550,-10)

if __name__=='__main__':unittest.main()
