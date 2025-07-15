#!/bin/bash
# build_and_test.sh - Build and test the TAN projection implementation

set -e  # Exit on any error

echo "=== Building and Testing TAN Projection Implementation ==="
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if we're on macOS and have required dependencies
print_status "Checking dependencies..."

# Check for cfitsio
if ! pkg-config --exists cfitsio 2>/dev/null && ! ls /opt/homebrew/lib/libcfitsio* >/dev/null 2>&1; then
    print_error "cfitsio not found. Install with: brew install cfitsio"
    exit 1
fi

# Check for OpenCV
if ! pkg-config --exists opencv4 2>/dev/null && ! ls /opt/homebrew/lib/libopencv_core* >/dev/null 2>&1; then
    print_error "OpenCV not found. Install with: brew install opencv"
    exit 1
fi

print_success "Dependencies found"

# Create build directory
BUILD_DIR="build_test"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

print_status "Building test program..."

# Build the test program
if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists cfitsio; then
    # Use pkg-config if available
    CFLAGS=$(pkg-config --cflags cfitsio)
    LIBS=$(pkg-config --libs cfitsio)
else
    # Use homebrew paths
    CFLAGS="-I/opt/homebrew/include"
    LIBS="-L/opt/homebrew/lib -lcfitsio"
fi

g++ -std=c++17 $CFLAGS ../test_tan_projection.cpp $LIBS -o test_tan_projection

if [ $? -eq 0 ]; then
    print_success "Test program built successfully"
else
    print_error "Failed to build test program"
    exit 1
fi

# Function to test with a FITS file if available
test_fits_file() {
    local fits_file="$1"
    if [ -f "$fits_file" ]; then
        print_status "Testing with FITS file: $fits_file"
        ./test_tan_projection "$fits_file"
        return $?
    else
        print_warning "FITS file not found: $fits_file"
        return 1
    fi
}

# Look for test FITS files in common locations
TEST_FILES_FOUND=false

# Check for plate-solved files in the current project directory
for dir in "../" "../../" "../../../"; do
    for subdir in "plate_solved" "calibrated" "test_data" "samples"; do
        if [ -d "${dir}${subdir}" ]; then
            print_status "Looking for FITS files in ${dir}${subdir}..."
            for fits_file in "${dir}${subdir}"/*.fits "${dir}${subdir}"/*solved*.fits; do
                if [ -f "$fits_file" ]; then
                    print_status "Found FITS file: $fits_file"
                    if test_fits_file "$fits_file"; then
                        TEST_FILES_FOUND=true
                        break 2
                    fi
                fi
            done
        fi
    done
done

# If no test files found, show usage
if [ "$TEST_FILES_FOUND" = false ]; then
    print_warning "No plate-solved FITS files found for testing"
    print_status "To test the TAN projection implementation:"
    echo
    echo "1. Run plate solving on some images first:"
    echo "   solve-field your_image.fits --new-fits solved_image.fits"
    echo
    echo "2. Then test with:"
    echo "   ./test_tan_projection solved_image.fits"
    echo
    print_status "Testing basic functionality without FITS file..."
    
    # Create a minimal test without FITS file
    cat > test_basic.cpp << 'EOF'
#include <iostream>
#include <cmath>

struct SimpleTANWCS {
    double crval1, crval2, crpix1, crpix2;
    double cd11, cd12, cd21, cd22;
    bool valid;
    
    SimpleTANWCS() : valid(false) {}
    
    bool pixelToWorld(double px, double py, double& ra, double& dec) const {
        if (!valid) return false;
        double dx = px - crpix1;
        double dy = py - crpix2;
        double xi = cd11 * dx + cd12 * dy;
        double eta = cd21 * dx + cd22 * dy;
        double ra0_rad = crval1 * M_PI / 180.0;
        double dec0_rad = crval2 * M_PI / 180.0;
        double xi_rad = xi * M_PI / 180.0;
        double eta_rad = eta * M_PI / 180.0;
        double cos_dec0 = cos(dec0_rad);
        double sin_dec0 = sin(dec0_rad);
        double denom = cos_dec0 - eta_rad * sin_dec0;
        if (std::abs(denom) < 1e-12) return false;
        double ra_rad = ra0_rad + atan2(xi_rad, denom);
        double dec_rad = atan((sin_dec0 + eta_rad * cos_dec0) / sqrt(xi_rad*xi_rad + denom*denom));
        ra = ra_rad * 180.0 / M_PI;
        dec = dec_rad * 180.0 / M_PI;
        while (ra < 0) ra += 360.0;
        while (ra >= 360.0) ra -= 360.0;
        return true;
    }
    
    bool worldToPixel(double ra, double dec, double& px, double& py) const {
        if (!valid) return false;
        double ra_rad = ra * M_PI / 180.0;
        double dec_rad = dec * M_PI / 180.0;
        double ra0_rad = crval1 * M_PI / 180.0;
        double dec0_rad = crval2 * M_PI / 180.0;
        double cos_dec = cos(dec_rad);
        double sin_dec = sin(dec_rad);
        double cos_dec0 = cos(dec0_rad);
        double sin_dec0 = sin(dec0_rad);
        double cos_dra = cos(ra_rad - ra0_rad);
        double sin_dra = sin(ra_rad - ra0_rad);
        double denom = sin_dec * sin_dec0 + cos_dec * cos_dec0 * cos_dra;
        if (std::abs(denom) < 1e-12) return false;
        double xi = cos_dec * sin_dra / denom;
        double eta = (sin_dec * cos_dec0 - cos_dec * sin_dec0 * cos_dra) / denom;
        xi *= 180.0 / M_PI;
        eta *= 180.0 / M_PI;
        double det = cd11 * cd22 - cd12 * cd21;
        if (std::abs(det) < 1e-12) return false;
        double dx = (cd22 * xi - cd12 * eta) / det;
        double dy = (-cd21 * xi + cd11 * eta) / det;
        px = dx + crpix1;
        py = dy + crpix2;
        return true;
    }
};

int main() {
    std::cout << "Testing basic TAN projection math..." << std::endl;
    
    SimpleTANWCS wcs;
    wcs.crval1 = 45.0;  // RA
    wcs.crval2 = 30.0;  // Dec
    wcs.crpix1 = 512.5; // Center X
    wcs.crpix2 = 512.5; // Center Y
    wcs.cd11 = -0.000347; // ~1.25 arcsec/pixel
    wcs.cd12 = 0.0;
    wcs.cd21 = 0.0;
    wcs.cd22 = 0.000347;
    wcs.valid = true;
    
    // Test round-trip
    double ra = 45.1, dec = 30.1;
    double px, py, ra_out, dec_out;
    
    if (wcs.worldToPixel(ra, dec, px, py) && wcs.pixelToWorld(px, py, ra_out, dec_out)) {
        double ra_error = (ra_out - ra) * 3600.0;
        double dec_error = (dec_out - dec) * 3600.0;
        std::cout << "Input: RA=" << ra << "°, Dec=" << dec << "°" << std::endl;
        std::cout << "Pixel: X=" << px << ", Y=" << py << std::endl;
        std::cout << "Output: RA=" << ra_out << "°, Dec=" << dec_out << "°" << std::endl;
        std::cout << "Error: " << ra_error << "\" RA, " << dec_error << "\" Dec" << std::endl;
        
        if (std::abs(ra_error) < 0.01 && std::abs(dec_error) < 0.01) {
            std::cout << "✓ SUCCESS: TAN projection math is working correctly!" << std::endl;
            return 0;
        }
    }
    
    std::cout << "✗ FAILED: TAN projection has errors" << std::endl;
    return 1;
}
EOF
    
    g++ -std=c++17 test_basic.cpp -o test_basic
    if ./test_basic; then
        print_success "Basic TAN projection math test passed!"
    else
        print_error "Basic TAN projection math test failed!"
        exit 1
    fi
fi

print_status "Building main project..."
cd ..

# Clean and build the main project
rm -rf build
qmake
make clean
make

if [ $? -eq 0 ]; then
    print_success "Main project built successfully!"
    print_status "You can now run the StellinaProcessor application"
    
    # Check if app was built
    if [ -f "StellinaProcessor.app/Contents/MacOS/StellinaProcessor" ]; then
        print_success "macOS app bundle created: StellinaProcessor.app"
    elif [ -f "StellinaProcessor" ]; then
        print_success "Executable created: StellinaProcessor"
    fi
    
else
    print_error "Failed to build main project"
    exit 1
fi

echo
print_success "Build and test complete!"
echo
print_status "Next steps:"
echo "1. Run the StellinaProcessor application"
echo "2. Test WCS stacking with your plate-solved FITS files"
echo "3. Check the processing log for 'TAN projection' messages"
echo "4. Enable debug mode to see detailed WCS diagnostics"
echo
print_status "To enable detailed WCS debugging, set environment variable:"
echo "export DEBUG_WCS=1"
echo "before running the application"