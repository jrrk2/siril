#!/usr/bin/env python3
"""
Siril Log Analyzer
Summarizes siril-cli processing logs to extract key statistics and results
"""

import re
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict

class SirilLogAnalyzer:
    def __init__(self):
        self.stats = {
            'total_images': 0,
            'successful_solves': 0,
            'failed_solves': 0,
            'commands_processed': 0,
            'start_time': None,
            'end_time': None
        }
        self.results = []
        self.failures = []
        self.commands = []
        
    def parse_log_file(self, log_file):
        """Parse siril-cli log file and extract key information"""
        print(f"📊 Analyzing Siril log: {log_file}")
        print("=" * 60)
        
        current_image = None
        current_command = None
        in_wcs_block = False
        wcs_data = {}
        
        with open(log_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Track timestamps
                if 'log: Welcome to siril' in line:
                    self.stats['start_time'] = datetime.now()
                
                # Track commands
                if 'log: Running command:' in line:
                    cmd_match = re.search(r'log: Running command: (\w+)', line)
                    if cmd_match:
                        current_command = cmd_match.group(1)
                        self.stats['commands_processed'] += 1
                
                # Track image loading
                if 'log: Reading FITS: file' in line:
                    img_match = re.search(r'file (img-\d+\.fits)', line)
                    if img_match:
                        current_image = img_match.group(1)
                        self.stats['total_images'] += 1
                
                # Track plate solving results
                if 'log: Siril solve succeeded' in line:
                    self.stats['successful_solves'] += 1
                    if current_image:
                        # Look for preceding WCS data
                        result = {
                            'image': current_image,
                            'status': 'SUCCESS',
                            'ra': wcs_data.get('crval1', 'N/A'),
                            'dec': wcs_data.get('crval2', 'N/A'),
                            'focal_length': None,
                            'pixel_size': None,
                            'field_of_view': None,
                            'stars_matched': wcs_data.get('nmatched', 'N/A')
                        }
                        self.results.append(result)
                        wcs_data = {}  # Reset for next image
                
                elif 'log: Plate solving failed' in line:
                    self.stats['failed_solves'] += 1
                    if current_image:
                        failure = {
                            'image': current_image,
                            'status': 'FAILED',
                            'reason': line.split('failed: ')[-1] if 'failed: ' in line else 'Unknown'
                        }
                        self.failures.append(failure)
                
                # Extract WCS data
                if '****Current WCS data*************' in line:
                    in_wcs_block = True
                elif in_wcs_block and line.startswith('crval1'):
                    match = re.search(r'crval1\s*=\s*([\d\.e\+\-]+)', line)
                    if match:
                        wcs_data['crval1'] = float(match.group(1))
                elif in_wcs_block and line.startswith('crval2'):
                    match = re.search(r'crval2\s*=\s*([\d\.e\+\-]+)', line)
                    if match:
                        wcs_data['crval2'] = float(match.group(1))
                elif in_wcs_block and 'nmatched:' in line:
                    match = re.search(r'nmatched:\s*(\d+)', line)
                    if match:
                        wcs_data['nmatched'] = int(match.group(1))
                elif in_wcs_block and '*' in line and 'Current WCS' not in line:
                    in_wcs_block = False
                
                # Extract additional metadata
                if current_image and 'log: Focal length:' in line:
                    match = re.search(r'Focal length:\s*([\d\.]+)\s*mm', line)
                    if match and self.results:
                        self.results[-1]['focal_length'] = float(match.group(1))
                
                if current_image and 'log: Pixel size:' in line:
                    match = re.search(r'Pixel size:\s*([\d\.]+)', line)
                    if match and self.results:
                        self.results[-1]['pixel_size'] = float(match.group(1))
                
                if current_image and 'log: Field of view:' in line:
                    match = re.search(r'Field of view:\s*(.+)', line)
                    if match and self.results:
                        self.results[-1]['field_of_view'] = match.group(1).strip()
        
        self.stats['end_time'] = datetime.now()
    
    def print_summary(self):
        """Print comprehensive summary of processing results"""
        
        # Overall Statistics
        print("📈 PROCESSING STATISTICS")
        print("-" * 30)
        print(f"Total Images Processed: {self.stats['total_images']}")
        print(f"Successful Plate Solves: {self.stats['successful_solves']}")
        print(f"Failed Plate Solves: {self.stats['failed_solves']}")
        
        if self.stats['total_images'] > 0:
            success_rate = (self.stats['successful_solves'] / self.stats['total_images']) * 100
            print(f"Success Rate: {success_rate:.1f}%")
        
        print(f"Total Commands: {self.stats['commands_processed']}")
        print()
        
        # Successful Results Summary
        if self.results:
            print("✅ SUCCESSFUL PLATE SOLVES")
            print("-" * 40)
            print(f"{'Image':<15} {'RA (deg)':<12} {'Dec (deg)':<12} {'Stars':<8} {'Focal':<8}")
            print("-" * 40)
            
            ra_values = []
            dec_values = []
            
            for result in self.results:
                ra = result['ra']
                dec = result['dec']
                stars = result['stars_matched']
                focal = result['focal_length']
                
                if isinstance(ra, float):
                    ra_values.append(ra)
                    ra_str = f"{ra:.4f}"
                else:
                    ra_str = str(ra)
                
                if isinstance(dec, float):
                    dec_values.append(dec)
                    dec_str = f"{dec:.4f}"
                else:
                    dec_str = str(dec)
                
                focal_str = f"{focal:.1f}" if focal else "N/A"
                
                print(f"{result['image']:<15} {ra_str:<12} {dec_str:<12} {stars:<8} {focal_str:<8}")
            
            # Statistics for successful solves
            if ra_values and dec_values:
                print()
                print("📊 COORDINATE STATISTICS")
                print("-" * 25)
                print(f"RA Range: {min(ra_values):.4f}° to {max(ra_values):.4f}°")
                print(f"Dec Range: {min(dec_values):.4f}° to {max(dec_values):.4f}°")
                print(f"RA Mean: {sum(ra_values)/len(ra_values):.4f}°")
                print(f"Dec Mean: {sum(dec_values)/len(dec_values):.4f}°")
                
                # Convert to HMS/DMS for first successful solve
                if ra_values and dec_values:
                    ra_deg = ra_values[0]
                    dec_deg = dec_values[0]
                    
                    # RA to HMS
                    ra_hours = ra_deg / 15.0
                    h = int(ra_hours)
                    m = int((ra_hours - h) * 60)
                    s = ((ra_hours - h) * 60 - m) * 60
                    
                    # Dec to DMS
                    d = int(dec_deg)
                    am = int(abs(dec_deg - d) * 60)
                    asec = (abs(dec_deg - d) * 60 - am) * 60
                    
                    print(f"First Solve (HMS/DMS): {h:02d}h {m:02d}m {s:05.2f}s, {d:+03d}° {am:02d}' {asec:04.1f}\"")
            print()
        
        # Failed Results Summary
        if self.failures:
            print("❌ FAILED PLATE SOLVES")
            print("-" * 40)
            
            failure_reasons = defaultdict(int)
            failed_images = []
            
            for failure in self.failures:
                failed_images.append(failure['image'])
                reason = failure['reason']
                failure_reasons[reason] += 1
            
            print(f"Failed Images: {', '.join(failed_images)}")
            print()
            print("Failure Reasons:")
            for reason, count in failure_reasons.items():
                print(f"  • {reason}: {count} image(s)")
            print()
        
        # Processing Quality Assessment
        print("🎯 QUALITY ASSESSMENT")
        print("-" * 25)
        
        if self.stats['total_images'] > 0:
            success_rate = (self.stats['successful_solves'] / self.stats['total_images']) * 100
            
            if success_rate >= 80:
                quality = "EXCELLENT 🏆"
            elif success_rate >= 60:
                quality = "GOOD 👍" 
            elif success_rate >= 40:
                quality = "FAIR 👌"
            else:
                quality = "POOR 👎"
            
            print(f"Overall Quality: {quality}")
            print(f"Success Rate: {success_rate:.1f}%")
            
            if success_rate < 70:
                print("\n💡 RECOMMENDATIONS:")
                if "could not be aligned" in str(self.failures):
                    print("  • Some images have poor star alignment - check tracking quality")
                if len(self.failures) > 3:
                    print("  • Consider enabling more aggressive quality filtering")
                    print("  • Check telescope focus and tracking during observation")
        
        print()
        print("=" * 60)
        print("🚀 Analysis Complete!")

def main():
    if len(sys.argv) != 2:
        print("Usage: python siril_log_analyzer.py <log_file>")
        print("Example: python siril_log_analyzer.py siril_processing.log")
        sys.exit(1)
    
    log_file = Path(sys.argv[1])
    
    if not log_file.exists():
        print(f"❌ Error: Log file '{log_file}' not found")
        sys.exit(1)
    
    analyzer = SirilLogAnalyzer()
    try:
        analyzer.parse_log_file(log_file)
        analyzer.print_summary()
    except Exception as e:
        print(f"❌ Error analyzing log file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
