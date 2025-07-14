#!/usr/bin/env python3
"""
Siril Log Analyzer with Tilt Estimation (based on solve-field coordinate hints)
"""

import re
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
from astropy import units as u
from scipy.optimize import minimize

class SirilLogAnalyzer:
    def __init__(self):
        self.results = []

    def parse_log_file(self, log_file):
        print(f"\U0001F4CA Analyzing Siril log: {log_file}")
        print("=" * 60)

        current_image = None
        current_stellina_ra = None
        current_stellina_dec = None

        with open(log_file, 'r') as f:
            for line in f:
                line = line.strip()

                # Extract image filename
                if 'dark_calibrated_img-' in line:
                    match = re.search(r'dark_calibrated_(img-\d+\.fits)', line)
                    if match:
                        current_image = match.group(1)

                # Extract Stellina RA/DEC hint passed to solve-field
                if 'solve-field with coordinate hints' in line:
                    match = re.search(r'RA=([\d.]+)°, Dec=([\d.]+)°', line)
                    if match:
                        current_stellina_ra = float(match.group(1))
                        current_stellina_dec = float(match.group(2))

                # Extract solved RA/DEC after plate solve
                elif 'solve-field found coordinates' in line:
                    match = re.search(r'RA=([\d.]+)°, Dec=([\d.]+)°', line)
                    if match:
                        solved_ra = float(match.group(1))
                        solved_dec = float(match.group(2))
                        if current_image:
                            self.results.append({
                                'image': current_image,
                                'status': 'SUCCESS',
                                'ra': solved_ra,
                                'dec': solved_dec,
                                'stellina_ra': current_stellina_ra,
                                'stellina_dec': current_stellina_dec
                            })
                            current_stellina_ra = None
                            current_stellina_dec = None

        self.fit_tilt_model()

    def fit_tilt_model(self):
        print("\n\U0001F4CB Fitting tilt model...")
        location = EarthLocation(lat=52.25*u.deg, lon=0.07*u.deg, height=0*u.m)

        true_alt, true_az = [], []
        est_alt, est_az = [], []
        now = Time(datetime.now())

        for r in self.results:
            try:
                if 'ra' in r and 'dec' in r and 'stellina_ra' in r and 'stellina_dec' in r:
                    solved = SkyCoord(ra=r['ra']*u.deg, dec=r['dec']*u.deg, frame='icrs')
                    stl = SkyCoord(ra=r['stellina_ra']*u.deg, dec=r['stellina_dec']*u.deg, frame='icrs')

                    solved_altaz = solved.transform_to(AltAz(obstime=now, location=location))
                    stl_altaz = stl.transform_to(AltAz(obstime=now, location=location))

                    true_alt.append(solved_altaz.alt.deg)
                    true_az.append(solved_altaz.az.deg)
                    est_alt.append(stl_altaz.alt.deg)
                    est_az.append(stl_altaz.az.deg)
            except Exception:
                continue

        true_alt = np.array(true_alt)
        true_az = np.radians(np.array(true_az))
        est_alt = np.array(est_alt)
        est_az = np.radians(np.array(est_az))

        delta_alt = true_alt - est_alt
        delta_az = np.degrees(true_az - est_az)

        def model(params):
            theta_n, theta_e = params
            dalt = np.degrees(theta_n * np.cos(est_az) + theta_e * np.sin(est_az))
            daz = np.degrees((theta_e * np.cos(est_az) - theta_n * np.sin(est_az)) / np.tan(np.radians(est_alt)))
            return np.mean((delta_alt - dalt)**2 + (delta_az - daz)**2)

        res = minimize(model, [0, 0], method='Nelder-Mead')
        theta_n, theta_e = res.x

        print(f"  • North tilt θ_N = {np.degrees(theta_n):.4f}°")
        print(f"  • East tilt  θ_E = {np.degrees(theta_e):.4f}°")


def main():
    if len(sys.argv) != 2:
        print("Usage: python siril_log_analyzer.py <log_file>")
        sys.exit(1)

    log_file = Path(sys.argv[1])
    if not log_file.exists():
        print(f"❌ Error: Log file '{log_file}' not found")
        sys.exit(1)

    analyzer = SirilLogAnalyzer()
    analyzer.parse_log_file(log_file)

if __name__ == "__main__":
    main()
