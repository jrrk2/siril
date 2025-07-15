// test_tan_projection.cpp
// Simple test program to verify TAN projection implementation
// Compile with: g++ -std=c++17 test_tan_projection.cpp -lcfitsio -o test_tan

#include <iostream>
#include <fitsio.h>
#include <cmath>
#include <string>
#include <iomanip>

struct SimpleTANWCS {
    double crval1, crval2;  // Reference RA, Dec (degrees)
    double crpix1, crpix2;  // Reference pixels (1-indexed)
    double cd11, cd12, cd21, cd22;  // CD matrix (degrees/pixel)
    bool valid;
    
    SimpleTANWCS() : crval1(0), crval2(0), crpix1(0), crpix2(0),
                     cd11(0), cd12(0), cd21(0), cd22(0), valid(false) {}
    
    bool pixelToWorld(double px, double py, double& ra, double& dec) const {
        if (!valid) return false;
        
        // Convert to 0-indexed and apply CD matrix
        double dx = px - crpix1;
        double dy = py - crpix2;
        
        double xi = cd11 * dx + cd12 * dy;
        double eta = cd21 * dx + cd22 * dy;
        
        // TAN projection inverse (standard formulae)
        double ra0_rad = crval1 * M_PI / 180.0;
        double dec0_rad = crval2 * M_PI / 180.0;
        double xi_rad = xi * M_PI / 180.0;
        double eta_rad = eta * M_PI / 180.0;
        
        double cos_dec0 = cos(dec0_rad);
        double sin_dec0 = sin(dec0_rad);
        
        double denom = cos_dec0 - eta_rad * sin_dec0;
        if (std::abs(denom) < 1e-12) return false;
        
        double ra_rad = ra0_rad + atan2(xi_rad, denom);
        double dec_rad = atan((sin_dec0 + eta_rad * cos_dec0) / (sqrt(xi_rad*xi_rad + denom*denom)));
        
        ra = ra_rad * 180.0 / M_PI;
        dec = dec_rad * 180.0 / M_PI;
        
        // Normalize RA
        while (ra < 0) ra += 360.0;
        while (ra >= 360.0) ra -= 360.0;
        
        return true;
    }
    
    bool worldToPixel(double ra, double dec, double& px, double& py) const {
        if (!valid) return false;
        
        // TAN projection forward
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
        
        // Convert to degrees
        xi *= 180.0 / M_PI;
        eta *= 180.0 / M_PI;
        
        // Apply inverse CD matrix
        double det = cd11 * cd22 - cd12 * cd21;
        if (std::abs(det) < 1e-12) return false;
        
        double dx = (cd22 * xi - cd12 * eta) / det;
        double dy = (-cd21 * xi + cd11 * eta) / det;
        
        px = dx + crpix1;
        py = dy + crpix2;
        
        return true;
    }
    
    double getPixelScale() const {
        if (!valid) return 0.0;
        double det = std::abs(cd11 * cd22 - cd12 * cd21);
        double scale_deg = sqrt(det);
        return scale_deg * 3600.0; // Convert to arcsec/pixel
    }
    
    void printDiagnostics() const {
        std::cout << "=== SimpleTAN WCS Diagnostics ===" << std::endl;
        std::cout << "Valid: " << (valid ? "Yes" : "No") << std::endl;
        if (!valid) return;
        
        std::cout << "Reference point: RA = " << std::fixed << std::setprecision(6) 
                  << crval1 << "°, Dec = " << crval2 << "°" << std::endl;
        std::cout << "Reference pixel: X = " << std::setprecision(2) 
                  << crpix1 << ", Y = " << crpix2 << std::endl;
        std::cout << "CD matrix:" << std::endl;
        std::cout << "  " << std::scientific << std::setprecision(6) 
                  << cd11 << "  " << cd12 << std::endl;
        std::cout << "  " << cd21 << "  " << cd22 << std::endl;
        
        double det = cd11 * cd22 - cd12 * cd21;
        std::cout << "Matrix determinant: " << det << std::endl;
        std::cout << "Pixel scale: " << std::fixed << std::setprecision(3) 
                  << getPixelScale() << " arcsec/pixel" << std::endl;
        
        // Test reference pixel
        double ra, dec;
        if (pixelToWorld(crpix1, crpix2, ra, dec)) {
            std::cout << "Reference pixel test:" << std::endl;
            std::cout << "  Maps to: RA = " << std::setprecision(6) 
                      << ra << "°, Dec = " << dec << "°" << std::endl;
            std::cout << "  Error: RA = " << (ra - crval1) 
                      << "°, Dec = " << (dec - crval2) << "°" << std::endl;
        }
    }
};

bool loadWCSFromFITS(const std::string& fits_file, SimpleTANWCS& wcs) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    if (fits_open_file(&fptr, fits_file.c_str(), READONLY, &status)) {
        std::cerr << "Failed to open FITS file: " << fits_file 
                  << " (status: " << status << ")" << std::endl;
        return false;
    }
    
    // Read essential WCS keywords with error checking
    auto readKey = [&](const char* keyword, double& value) -> bool {
        int local_status = 0;
        return fits_read_key(fptr, TDOUBLE, keyword, &value, nullptr, &local_status) == 0;
    };
    
    bool success = true;
    success &= readKey("CRVAL1", wcs.crval1);
    success &= readKey("CRVAL2", wcs.crval2);
    success &= readKey("CRPIX1", wcs.crpix1);
    success &= readKey("CRPIX2", wcs.crpix2);
    
    // Try CD matrix first, fall back to CDELT
    if (readKey("CD1_1", wcs.cd11) && readKey("CD1_2", wcs.cd12) &&
        readKey("CD2_1", wcs.cd21) && readKey("CD2_2", wcs.cd22)) {
        std::cout << "Using CD matrix for WCS" << std::endl;
    } else {
        // Use CDELT and assume no rotation
        double cdelt1, cdelt2;
        if (readKey("CDELT1", cdelt1) && readKey("CDELT2", cdelt2)) {
            wcs.cd11 = cdelt1; 
            wcs.cd12 = 0.0;
            wcs.cd21 = 0.0; 
            wcs.cd22 = cdelt2;
            std::cout << "Using CDELT for WCS (no rotation)" << std::endl;
        } else {
            success = false;
            std::cout << "No CD matrix or CDELT found" << std::endl;
        }
    }
    
    // Check coordinate types
    char ctype1[FLEN_VALUE], ctype2[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "CTYPE1", ctype1, nullptr, &status) == 0 &&
        fits_read_key(fptr, TSTRING, "CTYPE2", ctype2, nullptr, &status) == 0) {
        
        std::string ctype1_str(ctype1);
        std::string ctype2_str(ctype2);
        
        // Remove quotes and whitespace
        ctype1_str.erase(std::remove(ctype1_str.begin(), ctype1_str.end(), '\''), ctype1_str.end());
        ctype2_str.erase(std::remove(ctype2_str.begin(), ctype2_str.end(), '\''), ctype2_str.end());
        
        std::cout << "Coordinate types: " << ctype1_str << ", " << ctype2_str << std::endl;
        
        if (ctype1_str.find("TAN") == std::string::npos || 
            ctype2_str.find("TAN") == std::string::npos) {
            std::cout << "Warning: Non-TAN projection detected!" << std::endl;
        }
    }
    
    fits_close_file(fptr, &status);
    
    wcs.valid = success;
    return success;
}

void testRoundTrip(const SimpleTANWCS& wcs, double test_ra, double test_dec) {
    std::cout << "\n=== Round-trip test ===" << std::endl;
    std::cout << "Input: RA = " << std::fixed << std::setprecision(6) 
              << test_ra << "°, Dec = " << test_dec << "°" << std::endl;
    
    // World to pixel
    double px, py;
    if (!wcs.worldToPixel(test_ra, test_dec, px, py)) {
        std::cout << "World to pixel conversion failed!" << std::endl;
        return;
    }
    
    std::cout << "Pixel: X = " << std::setprecision(2) << px 
              << ", Y = " << py << std::endl;
    
    // Pixel back to world
    double ra_out, dec_out;
    if (!wcs.pixelToWorld(px, py, ra_out, dec_out)) {
        std::cout << "Pixel to world conversion failed!" << std::endl;
        return;
    }
    
    std::cout << "Output: RA = " << std::setprecision(6) 
              << ra_out << "°, Dec = " << dec_out << "°" << std::endl;
    
    // Calculate errors
    double ra_error = ra_out - test_ra;
    double dec_error = dec_out - test_dec;
    
    // Handle RA wraparound
    if (ra_error > 180.0) ra_error -= 360.0;
    if (ra_error < -180.0) ra_error += 360.0;
    
    std::cout << "Errors: RA = " << std::scientific << std::setprecision(3) 
              << ra_error << "°, Dec = " << dec_error << "°" << std::endl;
    
    // Convert errors to arcseconds
    double ra_error_arcsec = ra_error * 3600.0;
    double dec_error_arcsec = dec_error * 3600.0;
    
    std::cout << "Errors (arcsec): RA = " << std::fixed << std::setprecision(3) 
              << ra_error_arcsec << "\", Dec = " << dec_error_arcsec << "\"" << std::endl;
    
    if (std::abs(ra_error_arcsec) < 0.01 && std::abs(dec_error_arcsec) < 0.01) {
        std::cout << "✓ EXCELLENT: Round-trip error < 0.01 arcsec" << std::endl;
    } else if (std::abs(ra_error_arcsec) < 0.1 && std::abs(dec_error_arcsec) < 0.1) {
        std::cout << "✓ GOOD: Round-trip error < 0.1 arcsec" << std::endl;
    } else if (std::abs(ra_error_arcsec) < 1.0 && std::abs(dec_error_arcsec) < 1.0) {
        std::cout << "⚠ ACCEPTABLE: Round-trip error < 1 arcsec" << std::endl;
    } else {
        std::cout << "✗ POOR: Large round-trip error detected!" << std::endl;
    }
}

void testCorners(const SimpleTANWCS& wcs, int image_width, int image_height) {
    std::cout << "\n=== Corner pixel test ===" << std::endl;
    std::cout << "Image dimensions: " << image_width << " x " << image_height << std::endl;
    
    struct Corner {
        std::string name;
        double x, y;
    };
    
    Corner corners[] = {
        {"Bottom-left", 1.0, 1.0},
        {"Bottom-right", double(image_width), 1.0},
        {"Top-left", 1.0, double(image_height)},
        {"Top-right", double(image_width), double(image_height)},
        {"Center", double(image_width)/2.0 + 0.5, double(image_height)/2.0 + 0.5}
    };
    
    for (const auto& corner : corners) {
        double ra, dec;
        if (wcs.pixelToWorld(corner.x, corner.y, ra, dec)) {
            std::cout << corner.name << " (" << std::fixed << std::setprecision(1) 
                      << corner.x << ", " << corner.y << "): "
                      << "RA = " << std::setprecision(6) << ra 
                      << "°, Dec = " << dec << "°" << std::endl;
        } else {
            std::cout << corner.name << ": Conversion failed!" << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <plate_solved_fits_file>" << std::endl;
        std::cout << "This program tests the TAN projection implementation by:" << std::endl;
        std::cout << "1. Loading WCS from a plate-solved FITS file" << std::endl;
        std::cout << "2. Testing coordinate transformations" << std::endl;
        std::cout << "3. Verifying round-trip accuracy" << std::endl;
        return 1;
    }
    
    std::string fits_file = argv[1];
    std::cout << "Testing TAN projection with file: " << fits_file << std::endl;
    std::cout << "=============================================" << std::endl;
    
    // Load WCS from FITS file
    SimpleTANWCS wcs;
    if (!loadWCSFromFITS(fits_file, wcs)) {
        std::cerr << "Failed to load WCS from FITS file!" << std::endl;
        return 1;
    }
    
    std::cout << "\nWCS loaded successfully!" << std::endl;
    
    // Print diagnostics
    wcs.printDiagnostics();
    
    // Get image dimensions
    fitsfile *fptr = nullptr;
    int status = 0;
    long naxes[2];
    
    if (fits_open_file(&fptr, fits_file.c_str(), READONLY, &status) == 0) {
        if (fits_get_img_size(fptr, 2, naxes, &status) == 0) {
            std::cout << "\nImage dimensions: " << naxes[0] << " x " << naxes[1] << " pixels" << std::endl;
        }
        fits_close_file(fptr, &status);
    }
    
    // Test corner pixels
    if (status == 0) {
        testCorners(wcs, naxes[0], naxes[1]);
    }
    
    // Test round-trip accuracy with reference point
    testRoundTrip(wcs, wcs.crval1, wcs.crval2);
    
    // Test round-trip with nearby points
    testRoundTrip(wcs, wcs.crval1 + 0.1, wcs.crval2 + 0.1);
    testRoundTrip(wcs, wcs.crval1 - 0.1, wcs.crval2 - 0.1);
    
    // Test with field center if different from reference
    if (naxes[0] > 0 && naxes[1] > 0) {
        double center_x = naxes[0] / 2.0 + 0.5;
        double center_y = naxes[1] / 2.0 + 0.5;
        
        double center_ra, center_dec;
        if (wcs.pixelToWorld(center_x, center_y, center_ra, center_dec)) {
            std::cout << "\n=== Field center test ===" << std::endl;
            testRoundTrip(wcs, center_ra, center_dec);
        }
    }
    
    std::cout << "\n=============================================" << std::endl;
    std::cout << "TAN projection test complete!" << std::endl;
    std::cout << "If all round-trip errors are < 1 arcsec, the implementation is working correctly." << std::endl;
    
    return 0;
}