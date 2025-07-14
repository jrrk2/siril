#!/usr/bin/env python3
"""
Stellina Coordinate Error Analysis Script

This script analyzes the coordinate error trends between Stellina mount coordinates
and actual solve-field plate solving results from processing logs.

Usage:
    python stellina_coordinate_analysis.py <log_file_path> [output_dir]
    
Examples:
    python stellina_coordinate_analysis.py ~/siril.txt
    python stellina_coordinate_analysis.py stellina_log.txt ./analysis_results
"""

import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import seaborn as sns
import sys
import os
import argparse

def parse_stellina_log(log_file_path):
    """
    Parse Stellina processing log to extract coordinate data
    """
    data = []
    
    with open(log_file_path, 'r') as f:
        content = f.read()
    
    # Pattern to match processing entries
    processing_pattern = r'Processing \d+ of \d+: dark_calibrated_img-(\d+)\.fits'
    stellina_coord_pattern = r'Using coordinates from calibrated file: RA=([\d.]+)°, Dec=([\d.]+)°'
    solve_field_pattern = r'solve-field found coordinates: RA=([\d.]+)°, Dec=([\d.]+)°'
    time_pattern = r'Extracted DATE-OBS: \'([^\']+)\''
    
    # Find all processing blocks
    lines = content.split('\n')
    
    current_image = None
    current_time = None
    stellina_ra = None
    stellina_dec = None
    solved_ra = None
    solved_dec = None
    
    for line in lines:
        # Check for image processing start
        proc_match = re.search(processing_pattern, line)
        if proc_match:
            # Save previous entry if complete
            if (current_image and stellina_ra and stellina_dec and 
                solved_ra and solved_dec and current_time):
                data.append({
                    'image_num': int(current_image),
                    'time': current_time,
                    'stellina_ra': float(stellina_ra),
                    'stellina_dec': float(stellina_dec),
                    'solved_ra': float(solved_ra),
                    'solved_dec': float(solved_dec),
                    'ra_error': float(stellina_ra) - float(solved_ra),
                    'dec_error': float(stellina_dec) - float(solved_dec)
                })
            
            # Reset for new image
            current_image = proc_match.group(1)
            stellina_ra = stellina_dec = solved_ra = solved_dec = None
            continue
        
        # Extract time
        time_match = re.search(time_pattern, line)
        if time_match:
            current_time = time_match.group(1)
            continue
            
        # Extract Stellina coordinates
        stellina_match = re.search(stellina_coord_pattern, line)
        if stellina_match:
            stellina_ra = stellina_match.group(1)
            stellina_dec = stellina_match.group(2)
            continue
            
        # Extract solved coordinates
        solved_match = re.search(solve_field_pattern, line)
        if solved_match:
            solved_ra = solved_match.group(1)
            solved_dec = solved_match.group(2)
            continue
    
    # Don't forget the last entry
    if (current_image and stellina_ra and stellina_dec and 
        solved_ra and solved_dec and current_time):
        data.append({
            'image_num': int(current_image),
            'time': current_time,
            'stellina_ra': float(stellina_ra),
            'stellina_dec': float(stellina_dec),
            'solved_ra': float(solved_ra),
            'solved_dec': float(solved_dec),
            'ra_error': float(stellina_ra) - float(solved_ra),
            'dec_error': float(stellina_dec) - float(solved_dec)
        })
    
    return pd.DataFrame(data)

def calculate_error_statistics(df):
    """Calculate error statistics"""
    stats = {
        'ra_error_mean': df['ra_error'].mean(),
        'ra_error_std': df['ra_error'].std(),
        'ra_error_rms': np.sqrt(np.mean(df['ra_error']**2)),
        'dec_error_mean': df['dec_error'].mean(),
        'dec_error_std': df['dec_error'].std(),
        'dec_error_rms': np.sqrt(np.mean(df['dec_error']**2)),
        'total_error_rms': np.sqrt(np.mean(df['ra_error']**2 + df['dec_error']**2)),
        'max_ra_error': df['ra_error'].abs().max(),
        'max_dec_error': df['dec_error'].abs().max(),
        'count': len(df)
    }
    return stats

def plot_coordinate_analysis(df, output_dir='.'):
    """Create comprehensive coordinate error analysis plots"""
    
    # Set style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Convert time strings to datetime for better plotting
    df['datetime'] = pd.to_datetime(df['time'])
    df['minutes_elapsed'] = (df['datetime'] - df['datetime'].min()).dt.total_seconds() / 60
    
    # 1. RA/Dec Error vs Time
    plt.subplot(3, 3, 1)
    plt.plot(df['minutes_elapsed'], df['ra_error'], 'b-', label='RA Error', alpha=0.7)
    plt.plot(df['minutes_elapsed'], df['dec_error'], 'r-', label='Dec Error', alpha=0.7)
    plt.xlabel('Time (minutes from start)')
    plt.ylabel('Error (degrees)')
    plt.title('Coordinate Errors vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Error magnitude vs time
    plt.subplot(3, 3, 2)
    error_magnitude = np.sqrt(df['ra_error']**2 + df['dec_error']**2)
    plt.plot(df['minutes_elapsed'], error_magnitude, 'g-', alpha=0.7)
    plt.xlabel('Time (minutes from start)')
    plt.ylabel('Total Error Magnitude (degrees)')
    plt.title('Total Coordinate Error Magnitude')
    plt.grid(True, alpha=0.3)
    
    # 3. RA Error histogram
    plt.subplot(3, 3, 3)
    plt.hist(df['ra_error'], bins=30, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel('RA Error (degrees)')
    plt.ylabel('Frequency')
    plt.title('RA Error Distribution')
    plt.axvline(df['ra_error'].mean(), color='red', linestyle='--', label='Mean')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 4. Dec Error histogram
    plt.subplot(3, 3, 4)
    plt.hist(df['dec_error'], bins=30, alpha=0.7, color='red', edgecolor='black')
    plt.xlabel('Dec Error (degrees)')
    plt.ylabel('Frequency')
    plt.title('Dec Error Distribution')
    plt.axvline(df['dec_error'].mean(), color='blue', linestyle='--', label='Mean')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 5. 2D Error scatter plot
    plt.subplot(3, 3, 5)
    scatter = plt.scatter(df['ra_error'], df['dec_error'], 
                         c=df['minutes_elapsed'], cmap='viridis', alpha=0.6)
    plt.xlabel('RA Error (degrees)')
    plt.ylabel('Dec Error (degrees)')
    plt.title('RA vs Dec Error (colored by time)')
    plt.colorbar(scatter, label='Minutes from start')
    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.axvline(0, color='black', linestyle='-', alpha=0.3)
    plt.grid(True, alpha=0.3)
    
    # 6. Stellina vs Solved coordinates
    plt.subplot(3, 3, 6)
    plt.scatter(df['stellina_ra'], df['solved_ra'], alpha=0.6, color='blue', label='RA')
    plt.scatter(df['stellina_dec'], df['solved_dec'], alpha=0.6, color='red', label='Dec')
    plt.xlabel('Stellina Coordinates (degrees)')
    plt.ylabel('Solve-field Coordinates (degrees)')
    plt.title('Stellina vs Solve-field Coordinates')
    
    # Add perfect correlation line
    all_coords = np.concatenate([df['stellina_ra'], df['stellina_dec'], 
                                df['solved_ra'], df['solved_dec']])
    min_coord, max_coord = all_coords.min(), all_coords.max()
    plt.plot([min_coord, max_coord], [min_coord, max_coord], 'k--', alpha=0.5, label='Perfect correlation')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 7. Error vs image number
    plt.subplot(3, 3, 7)
    plt.plot(df['image_num'], df['ra_error'], 'b.', alpha=0.6, label='RA Error')
    plt.plot(df['image_num'], df['dec_error'], 'r.', alpha=0.6, label='Dec Error')
    plt.xlabel('Image Number')
    plt.ylabel('Error (degrees)')
    plt.title('Coordinate Errors vs Image Number')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 8. Error in arcminutes
    plt.subplot(3, 3, 8)
    ra_error_arcmin = df['ra_error'] * 60
    dec_error_arcmin = df['dec_error'] * 60
    plt.plot(df['minutes_elapsed'], ra_error_arcmin, 'b-', alpha=0.7, label='RA Error')
    plt.plot(df['minutes_elapsed'], dec_error_arcmin, 'r-', alpha=0.7, label='Dec Error')
    plt.xlabel('Time (minutes from start)')
    plt.ylabel('Error (arcminutes)')
    plt.title('Coordinate Errors in Arcminutes')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 9. Statistics summary text
    plt.subplot(3, 3, 9)
    plt.axis('off')
    stats = calculate_error_statistics(df)
    
    stats_text = f"""
    Coordinate Error Statistics:
    
    RA Error:
    • Mean: {stats['ra_error_mean']:.4f}°
    • Std: {stats['ra_error_std']:.4f}°
    • RMS: {stats['ra_error_rms']:.4f}°
    • Max: {stats['max_ra_error']:.4f}°
    
    Dec Error:
    • Mean: {stats['dec_error_mean']:.4f}°
    • Std: {stats['dec_error_std']:.4f}°
    • RMS: {stats['dec_error_rms']:.4f}°
    • Max: {stats['max_dec_error']:.4f}°
    
    Total RMS Error: {stats['total_error_rms']:.4f}°
    ({stats['total_error_rms']*60:.1f} arcminutes)
    
    Images Analyzed: {stats['count']}
    """
    
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
             verticalalignment='top', fontfamily='monospace', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/stellina_coordinate_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return stats

def create_detailed_error_table(df, output_dir='.'):
    """Create a detailed table of errors for specific images"""
    
    # Sort by error magnitude
    df['total_error'] = np.sqrt(df['ra_error']**2 + df['dec_error']**2)
    df_sorted = df.sort_values('total_error', ascending=False)
    
    # Create summary table
    print("\n" + "="*80)
    print("TOP 10 LARGEST COORDINATE ERRORS")
    print("="*80)
    print(f"{'Image':<8} {'Time':<20} {'RA_Err':<8} {'Dec_Err':<8} {'Total':<8} {'RA_Err(\")':<10} {'Dec_Err(\")':<10}")
    print("-"*80)
    
    for i, row in df_sorted.head(10).iterrows():
        ra_err_arcsec = row['ra_error'] * 3600
        dec_err_arcsec = row['dec_error'] * 3600
        print(f"img-{row['image_num']:04d} {row['time']:<20} {row['ra_error']:7.3f}° {row['dec_error']:7.3f}° {row['total_error']:7.3f}° {ra_err_arcsec:8.1f}\" {dec_err_arcsec:8.1f}\"")
    
    # Save to CSV
    output_file = f'{output_dir}/stellina_coordinate_errors.csv'
    df_sorted.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")
    
    return df_sorted

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze Stellina coordinate errors from processing logs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s ~/siril.txt
  %(prog)s stellina_log.txt
  %(prog)s /path/to/log.txt ./analysis_results
  %(prog)s siril.txt --output-dir ./plots --verbose
        """
    )
    
    parser.add_argument('log_file', 
                       help='Path to the Stellina processing log file')
    
    parser.add_argument('--output-dir', '-o',
                       default='.',
                       help='Output directory for plots and CSV files (default: current directory)')
    
    parser.add_argument('--verbose', '-v',
                       action='store_true',
                       help='Enable verbose output')
    
    parser.add_argument('--no-plots',
                       action='store_true',
                       help='Skip generating plots (useful for headless environments)')
    
    parser.add_argument('--csv-only',
                       action='store_true',
                       help='Only generate CSV output, no plots or statistics')
    
    return parser.parse_args()

def main():
    """Main analysis function"""
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Expand user path (handles ~/path)
    log_file_path = os.path.expanduser(args.log_file)
    output_dir = os.path.expanduser(args.output_dir)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        if args.verbose:
            print(f"Created output directory: {output_dir}")
    
    print("Stellina Coordinate Error Analysis")
    print("=" * 50)
    
    # Validate input file
    if not os.path.exists(log_file_path):
        print(f"Error: Could not find log file at '{log_file_path}'")
        print("\nUsage examples:")
        print("  python stellina_coordinate_analysis.py ~/siril.txt")
        print("  python stellina_coordinate_analysis.py /path/to/stellina.log")
        sys.exit(1)
    
    try:
        # Parse the log file
        print(f"Parsing log file: {log_file_path}")
        if args.verbose:
            file_size = os.path.getsize(log_file_path)
            print(f"Log file size: {file_size:,} bytes")
            
        df = parse_stellina_log(log_file_path)
        
        if df.empty:
            print("No coordinate data found in log file!")
            print("\nThe log file should contain lines with:")
            print("- 'Processing X of Y: dark_calibrated_img-XXXX.fits'")
            print("- 'Using coordinates from calibrated file: RA=X.X°, Dec=Y.Y°'")
            print("- 'solve-field found coordinates: RA=X.X°, Dec=Y.Y°'")
            sys.exit(1)
        
        print(f"Found {len(df)} processed images with coordinate data")
        
        if args.csv_only:
            # Only generate CSV output
            output_file = os.path.join(output_dir, 'stellina_coordinate_errors.csv')
            df['total_error'] = np.sqrt(df['ra_error']**2 + df['dec_error']**2)
            df_sorted = df.sort_values('total_error', ascending=False)
            df_sorted.to_csv(output_file, index=False)
            print(f"CSV output saved to: {output_file}")
            return
        
        # Calculate statistics
        stats = calculate_error_statistics(df)
        
        print(f"\nKey Findings:")
        print(f"• Average RA error: {stats['ra_error_mean']:+.4f}° ({stats['ra_error_mean']*60:+.1f} arcmin)")
        print(f"• Average Dec error: {stats['dec_error_mean']:+.4f}° ({stats['dec_error_mean']*60:+.1f} arcmin)")
        print(f"• Total RMS error: {stats['total_error_rms']:.4f}° ({stats['total_error_rms']*60:.1f} arcmin)")
        print(f"• Maximum RA error: {stats['max_ra_error']:.4f}° ({stats['max_ra_error']*60:.1f} arcmin)")
        print(f"• Maximum Dec error: {stats['max_dec_error']:.4f}° ({stats['max_dec_error']*60:.1f} arcmin)")
        
        if args.verbose:
            print(f"• RA error std dev: {stats['ra_error_std']:.4f}° ({stats['ra_error_std']*60:.1f} arcmin)")
            print(f"• Dec error std dev: {stats['dec_error_std']:.4f}° ({stats['dec_error_std']*60:.1f} arcmin)")
            print(f"• RA RMS error: {stats['ra_error_rms']:.4f}° ({stats['ra_error_rms']*60:.1f} arcmin)")
            print(f"• Dec RMS error: {stats['dec_error_rms']:.4f}° ({stats['dec_error_rms']*60:.1f} arcmin)")
        
        # Create plots (unless disabled)
        if not args.no_plots:
            print(f"\nGenerating analysis plots in: {output_dir}")
            plot_stats = plot_coordinate_analysis(df, output_dir)
        
        # Create detailed error table
        print(f"\nGenerating detailed error analysis in: {output_dir}")
        df_detailed = create_detailed_error_table(df, output_dir)
        
        print("\nAnalysis complete!")
        print(f"Results saved to: {output_dir}")
        
    except FileNotFoundError:
        print(f"Error: Could not find log file at '{log_file_path}'")
        print("\nPlease check the file path and try again.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
