// Updated WcsAstrometricStacker.cpp with TAN projection implementation
// Replace the WCS-related functions with these implementations

#include "WcsAstrometricStacker.h"
#include <QFileInfo>
#include <QTimer>
#include <QDateTime>
#include <QApplication>
#include <QTextStream>
#include <QRegularExpression>
#include <QDebug>
#include <cmath>
#include <algorithm>

// SimpleTANWCS implementation
bool SimpleTANWCS::pixelToWorld(double px, double py, double& ra, double& dec) const {
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

bool SimpleTANWCS::worldToPixel(double ra, double dec, double& px, double& py) const {
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

void SimpleTANWCS::printDiagnostics() const {
    if (!valid) {
        qDebug() << "WCS is invalid";
        return;
    }
    
    qDebug() << "=== SimpleTAN WCS Diagnostics ===";
    qDebug() << "Reference point: RA =" << crval1 << "°, Dec =" << crval2 << "°";
    qDebug() << "Reference pixel: X =" << crpix1 << ", Y =" << crpix2;
    qDebug() << "CD matrix:";
    qDebug() << "  " << cd11 << cd12;
    qDebug() << "  " << cd21 << cd22;
    
    double det = cd11 * cd22 - cd12 * cd21;
    qDebug() << "Matrix determinant:" << det;
    qDebug() << "Pixel scale:" << getPixelScale() << "arcsec/pixel";
    
    // Test reference pixel
    double ra, dec;
    if (pixelToWorld(crpix1, crpix2, ra, dec)) {
        qDebug() << "Reference pixel maps to: RA =" << ra << "°, Dec =" << dec << "°";
        qDebug() << "Should equal CRVAL1/CRVAL2 - Error: RA =" 
                 << (ra - crval1) << "°, Dec =" << (dec - crval2) << "°";
    }
}

double SimpleTANWCS::getPixelScale() const {
    if (!valid) return 0.0;
    
    // Calculate pixel scale from CD matrix determinant
    double det = std::abs(cd11 * cd22 - cd12 * cd21);
    double scale_deg = sqrt(det);
    return scale_deg * 3600.0; // Convert to arcsec/pixel
}

// Updated WCSAstrometricStacker constructor
WCSAstrometricStacker::WCSAstrometricStacker(QObject *parent)
    : QObject(parent)
    , m_progress_bar(nullptr)
    , m_status_label(nullptr)
    , m_stacking_active(false)
    , m_current_image_index(0)
    , m_processing_timer(new QTimer(this))
    , m_total_processing_time(0.0)
    , m_pixels_processed(0)
    , m_stacked_image()
    , m_pixels_rejected(0)
{
    // No WCSLIB initialization needed
    
    // Set up processing timer for non-blocking operation
    m_processing_timer->setSingleShot(true);
    //    connect(m_processing_timer, &QTimer::timeout, this, &WCSAstrometricStacker::processNextImage);
}

WCSAstrometricStacker::~WCSAstrometricStacker() {
    // No WCSLIB cleanup needed
}

// Updated loadWCSFromFITS function
bool WCSAstrometricStacker::loadWCSFromFITS(const QString &fits_file, WCSImageData &img_data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fits_file.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        emit errorOccurred(QString("Failed to open FITS file: %1 (status: %2)").arg(fits_file).arg(status));
        return false;
    }
    
    img_data.filename = fits_file;
    img_data.solved_filename = fits_file;
    
    // Read essential WCS keywords with error checking
    auto readKey = [&](const char* keyword, double& value) -> bool {
        int local_status = 0;
        return fits_read_key(fptr, TDOUBLE, keyword, &value, nullptr, &local_status) == 0;
    };
    
    bool success = true;
    success &= readKey("CRVAL1", img_data.wcs.crval1);
    success &= readKey("CRVAL2", img_data.wcs.crval2);
    success &= readKey("CRPIX1", img_data.wcs.crpix1);
    success &= readKey("CRPIX2", img_data.wcs.crpix2);
    
    // Try CD matrix first, fall back to CDELT
    if (readKey("CD1_1", img_data.wcs.cd11) && readKey("CD1_2", img_data.wcs.cd12) &&
        readKey("CD2_1", img_data.wcs.cd21) && readKey("CD2_2", img_data.wcs.cd22)) {
        // CD matrix available
        logProcessing("Using CD matrix for WCS");
    } else {
        // Use CDELT and assume no rotation
        double cdelt1, cdelt2;
        if (readKey("CDELT1", cdelt1) && readKey("CDELT2", cdelt2)) {
            img_data.wcs.cd11 = cdelt1; 
            img_data.wcs.cd12 = 0.0;
            img_data.wcs.cd21 = 0.0; 
            img_data.wcs.cd22 = cdelt2;
            logProcessing("Using CDELT for WCS (no rotation)");
        } else {
            success = false;
            logProcessing("No CD matrix or CDELT found");
        }
    }
    
    // Verify coordinate types are TAN projection
    char ctype1[FLEN_VALUE], ctype2[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "CTYPE1", ctype1, nullptr, &status) == 0 &&
        fits_read_key(fptr, TSTRING, "CTYPE2", ctype2, nullptr, &status) == 0) {
        
        QString ctype1_str = QString::fromLatin1(ctype1).trimmed().remove('\'').remove('"');
        QString ctype2_str = QString::fromLatin1(ctype2).trimmed().remove('\'').remove('"');
        
        if (!ctype1_str.contains("TAN") || !ctype2_str.contains("TAN")) {
            logProcessing(QString("Warning: Non-TAN projection detected: %1, %2").arg(ctype1_str).arg(ctype2_str));
            logProcessing("Proceeding anyway - results may be inaccurate for non-TAN projections");
        }
    }
    
    img_data.wcs.valid = success;
    img_data.wcs_valid = success;
    
    if (success) {
        logProcessing(QString("WCS loaded successfully from: %1").arg(QFileInfo(fits_file).fileName()));
        logProcessing(QString("  CRVAL: %1, %2").arg(img_data.wcs.crval1, 0, 'f', 6).arg(img_data.wcs.crval2, 0, 'f', 6));
        logProcessing(QString("  CRPIX: %1, %2").arg(img_data.wcs.crpix1, 0, 'f', 2).arg(img_data.wcs.crpix2, 0, 'f', 2));
        logProcessing(QString("  Pixel scale: %1 arcsec/pixel").arg(img_data.wcs.getPixelScale(), 0, 'f', 3));
        
        // Print full diagnostics in debug mode
        if (qEnvironmentVariableIsSet("DEBUG_WCS")) {
            img_data.wcs.printDiagnostics();
        }
    }
    
    // Read image dimensions and data
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, &status)) {
        emit errorOccurred(QString("Failed to get image size: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    long totalPixels = naxes[0] * naxes[1];
    std::vector<float> pixels(totalPixels);
    
    if (fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status)) {
        emit errorOccurred(QString("Failed to read image data: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Convert to OpenCV Mat (FITS is row-major, OpenCV expects this)
    img_data.image = cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
    
    fits_close_file(fptr, &status);
    
    logProcessing(QString("Successfully loaded WCS and image data: %1x%2 pixels")
                 .arg(naxes[0]).arg(naxes[1]));
    
    return success;
}

// Updated computeOptimalWCS function
bool WCSAstrometricStacker::computeOptimalWCS() {
    if (m_images.empty()) {
        emit errorOccurred("No images available for WCS calculation");
        return false;
    }
    
    logProcessing("Computing optimal output WCS using TAN projection...");
    
    // Calculate bounding box in world coordinates
    double ra_min = 360.0, ra_max = 0.0;
    double dec_min = 90.0, dec_max = -90.0;
    
    // Use the first valid WCS as reference
    const SimpleTANWCS* ref_wcs = nullptr;
    for (const auto& img : m_images) {
        if (img->wcs_valid) {
            ref_wcs = &img->wcs;
            break;
        }
    }
    
    if (!ref_wcs) {
        emit errorOccurred("No valid WCS found in any image");
        return false;
    }
    
    logProcessing(QString("Using reference WCS from first image: center RA=%1°, Dec=%2°")
                 .arg(ref_wcs->crval1, 0, 'f', 6).arg(ref_wcs->crval2, 0, 'f', 6));
    
    // For each image, find the corner points in world coordinates
    for (const auto& img : m_images) {
        if (!img->wcs_valid) continue;
        
        // Get image corners in pixel coordinates (1-indexed for FITS convention)
        double corners_pix[4][2] = {
            {1.0, 1.0},  // Bottom left
            {double(img->image.cols), 1.0},  // Bottom right
            {1.0, double(img->image.rows)},  // Top left
            {double(img->image.cols), double(img->image.rows)}  // Top right
        };
        
        // Convert each corner to world coordinates
        for (int i = 0; i < 4; ++i) {
            double ra, dec;
            if (img->wcs.pixelToWorld(corners_pix[i][0], corners_pix[i][1], ra, dec)) {
                // Handle RA wrapping at 0/360 degrees
                if (i == 0) {
                    ra_min = ra_max = ra;
                } else {
                    // Check if we cross the RA=0/360 boundary
                    if (std::abs(ra - ra_min) > 180.0) {
                        if (ra < ra_min) ra += 360.0;
                        else ra_min += 360.0;
                    }
                    if (std::abs(ra - ra_max) > 180.0) {
                        if (ra < ra_max) ra += 360.0;
                        else ra_max += 360.0;
                    }
                    
                    ra_min = std::min(ra_min, ra);
                    ra_max = std::max(ra_max, ra);
                }
                
                dec_min = std::min(dec_min, dec);
                dec_max = std::max(dec_max, dec);
            }
        }
    }
    
    // Normalize RA range
    while (ra_min >= 360.0) ra_min -= 360.0;
    while (ra_max >= 360.0) ra_max -= 360.0;
    
    // Add margin
    double ra_margin = (ra_max - ra_min) * 0.05;
    double dec_margin = (dec_max - dec_min) * 0.05;
    
    ra_min = std::max(0.0, ra_min - ra_margin);
    ra_max = std::min(360.0, ra_max + ra_margin);
    dec_min = std::max(-90.0, dec_min - dec_margin);
    dec_max = std::min(90.0, dec_max + dec_margin);
    
    logProcessing(QString("Field bounds: RA=%1° to %2°, Dec=%3° to %4°")
                 .arg(ra_min, 0, 'f', 4).arg(ra_max, 0, 'f', 4)
                 .arg(dec_min, 0, 'f', 4).arg(dec_max, 0, 'f', 4));
    
    // Determine pixel scale - use the first image as reference if not overridden
    if (m_params.output_pixel_scale <= 0.0) {
        m_output_pixel_scale = ref_wcs->getPixelScale();
    } else {
        m_output_pixel_scale = m_params.output_pixel_scale;
    }
    
    logProcessing(QString("Using output pixel scale: %1 arcsec/pixel").arg(m_output_pixel_scale, 0, 'f', 3));
    
    // Calculate output dimensions based on sky coverage and pixel scale
    double ra_span_deg = ra_max - ra_min;
    if (ra_span_deg < 0) ra_span_deg += 360.0;
    
    double dec_span_deg = dec_max - dec_min;
    
    // Convert to pixels using pixel scale
    int width, height;
    
    // Adjust RA span for cos(Dec) factor at center declination
    double dec_center = (dec_min + dec_max) / 2.0;
    double cos_dec_factor = std::cos(dec_center * M_PI / 180.0);
    if (cos_dec_factor < 0.01) cos_dec_factor = 0.01; // Avoid division by very small values
    
    width = int(ra_span_deg * 3600.0 / (m_output_pixel_scale / cos_dec_factor));
    height = int(dec_span_deg * 3600.0 / m_output_pixel_scale);
    
    // Override dimensions if specified
    if (m_params.output_width > 0) width = m_params.output_width;
    if (m_params.output_height > 0) height = m_params.output_height;
    
    // Limit to reasonable size
    width = std::min(std::max(width, 100), 10000);
    height = std::min(std::max(height, 100), 10000);
    
    m_output_size = cv::Size(width, height);
    
    logProcessing(QString("Output dimensions: %1 x %2 pixels").arg(width).arg(height));
    
    // Initialize output WCS as TAN projection centered on field
    double ra_center = (ra_min + ra_max) / 2.0;
    if (ra_max < ra_min) ra_center = fmod(ra_center + 180.0, 360.0); // Handle RA wrap
    
    m_output_wcs.crpix1 = width / 2.0 + 0.5;   // Center pixel X (1-indexed)
    m_output_wcs.crpix2 = height / 2.0 + 0.5;  // Center pixel Y (1-indexed)
    m_output_wcs.crval1 = ra_center;           // RA at reference pixel
    m_output_wcs.crval2 = dec_center;          // Dec at reference pixel
    m_output_wcs.cd11 = -m_output_pixel_scale / 3600.0; // RA step (degrees, negative for normal orientation)
    m_output_wcs.cd12 = 0.0;                            // No rotation
    m_output_wcs.cd21 = 0.0;                            // No rotation
    m_output_wcs.cd22 = m_output_pixel_scale / 3600.0;  // Dec step (degrees)
    m_output_wcs.valid = true;
    
    logProcessing(QString("Output WCS computed: center RA=%1°, Dec=%2°")
                 .arg(ra_center, 0, 'f', 6)
                 .arg(dec_center, 0, 'f', 6));
    
    // Test the output WCS
    double test_ra, test_dec;
    if (m_output_wcs.pixelToWorld(m_output_wcs.crpix1, m_output_wcs.crpix2, test_ra, test_dec)) {
        logProcessing(QString("WCS test: center pixel maps to RA=%1°, Dec=%2° (should match center)")
                     .arg(test_ra, 0, 'f', 6).arg(test_dec, 0, 'f', 6));
    }
    
    return true;
}

// SOLUTION: Two-pass stacking with overlap-aware brightness compensation
// Replace the main stacking loop in stackImages() with this improved version

bool WCSAstrometricStacker::stackImages() {
    updateProgress(0, "Starting TAN projection astrometric stacking...");
    
    if (m_images.empty()) {
        emit errorOccurred("No images loaded for stacking");
        return false;
    }
    
    // Step 1: Compute optimal output WCS and dimensions
    updateProgress(5, "Computing optimal output coordinate system...");
    if (!computeOptimalWCS()) {
        emit errorOccurred("Failed to compute optimal WCS");
        return false;
    }
    
    // Initialize output images
    m_stacked_image = cv::Mat::zeros(m_output_size, CV_32F);
    m_weight_map = cv::Mat::zeros(m_output_size, CV_32F);
    m_overlap_map = cv::Mat::zeros(m_output_size, CV_8U);
    
    // NEW: Create overlap count map first (Pass 1)
    cv::Mat overlap_count = cv::Mat::zeros(m_output_size, CV_32F);
    
    updateProgress(10, "Pass 1: Computing overlap regions...");
    
    // First pass: count overlaps for each output pixel
    for (size_t i = 0; i < m_images.size(); ++i) {
        const auto& img = m_images[i];
        if (!img->wcs_valid || img->image.empty()) continue;
        
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                double ra, dec;
                if (m_output_wcs.pixelToWorld(x + 1.0, y + 1.0, ra, dec)) {
                    double imgpx, imgpy;
                    if (img->wcs.worldToPixel(ra, dec, imgpx, imgpy)) {
                        float srcX = imgpx - 1.0f;
                        float srcY = imgpy - 1.0f;
                        
                        if (srcX >= 0.0f && srcX < img->image.cols - 1 &&
                            srcY >= 0.0f && srcY < img->image.rows - 1) {
                            overlap_count.at<float>(y, x) += 1.0f;
                        }
                    }
                }
            }
        }
    }
    
    updateProgress(30, "Pass 2: Stacking with overlap-aware weights...");
    
    // Prepare for pixel rejection if using sigma clipping
    std::vector<std::vector<float>> pixelStacks;
    std::vector<std::vector<float>> weightStacks;
    bool useSigmaClipping = (m_params.rejection == StackingParams::SIGMA_CLIPPING);
    
    if (useSigmaClipping) {
        pixelStacks.resize(m_output_size.height * m_output_size.width);
        weightStacks.resize(m_output_size.height * m_output_size.width);
    }
    
    // Step 2: Second pass - reproject and combine with overlap compensation
    m_pixels_processed = 0;
    m_pixels_rejected = 0;
    
    for (size_t i = 0; i < m_images.size(); ++i) {
        updateProgress(30 + (i * 50) / m_images.size(), 
                       QString("Processing image %1 of %2: %3")
                       .arg(i+1)
                       .arg(m_images.size())
                       .arg(QFileInfo(m_images[i]->filename).fileName()));
        
        const auto& img = m_images[i];
        if (!img->wcs_valid || img->image.empty()) {
            logProcessing(QString("Skipping invalid image: %1").arg(QFileInfo(img->filename).fileName()));
            continue;
        }
        
        // Base image weight (quality * exposure)
        float baseImageWeight = img->quality_score;
        if (m_params.normalize_exposure && img->exposure_time > 0) {
            baseImageWeight *= img->exposure_time;
        }
        
        int pixelsProcessedThisImage = 0;
        
        // Loop through output pixels with overlap-aware weighting
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                double ra, dec;
                if (m_output_wcs.pixelToWorld(x + 1.0, y + 1.0, ra, dec)) {
                    double imgpx, imgpy;
                    if (img->wcs.worldToPixel(ra, dec, imgpx, imgpy)) {
                        float srcX = imgpx - 1.0f;
                        float srcY = imgpy - 1.0f;
                        
                        if (srcX >= 0.0f && srcX < img->image.cols - 1 &&
                            srcY >= 0.0f && srcY < img->image.rows - 1) {
                            
                            // Bilinear interpolation
                            int x0 = int(floor(srcX));
                            int y0 = int(floor(srcY));
                            int x1 = x0 + 1;
                            int y1 = y0 + 1;
                            
                            if (x0 >= 0 && x1 < img->image.cols && y0 >= 0 && y1 < img->image.rows) {
                                float dx = srcX - x0;
                                float dy = srcY - y0;
                                
                                float v00 = img->image.at<float>(y0, x0);
                                float v01 = img->image.at<float>(y0, x1);
                                float v10 = img->image.at<float>(y1, x0);
                                float v11 = img->image.at<float>(y1, x1);
                                
                                float v0 = v00 * (1 - dx) + v01 * dx;
                                float v1 = v10 * (1 - dx) + v11 * dx;
                                float pixelValue = v0 * (1 - dy) + v1 * dy;
                                
                                if (std::isfinite(pixelValue)) {
                                    // CRITICAL FIX: Overlap-aware weight calculation
                                    float overlapCount = overlap_count.at<float>(y, x);
                                    float pixelWeight = baseImageWeight;
                                    
                                    // Apply overlap compensation factor
                                    if (overlapCount > 1.0f) {
                                        // Reduce weight for high-overlap regions to prevent over-brightening
                                        pixelWeight /= sqrt(overlapCount);
                                    }
                                    
                                    // Alternative: Normalize by expected total weight
                                    // This ensures consistent brightness regardless of overlap
                                    float normalizedWeight = baseImageWeight / overlapCount;
                                    
                                    if (useSigmaClipping) {
                                        int idx = y * m_output_size.width + x;
                                        pixelStacks[idx].push_back(pixelValue);
                                        weightStacks[idx].push_back(normalizedWeight);
                                    } else {
                                        m_stacked_image.at<float>(y, x) += pixelValue * normalizedWeight;
                                        m_weight_map.at<float>(y, x) += normalizedWeight;
                                        m_overlap_map.at<uchar>(y, x)++;
                                    }
                                    
                                    pixelsProcessedThisImage++;
                                    m_pixels_processed++;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        logProcessing(QString("Reprojected image %1: %2 pixels processed")
                     .arg(i+1)
                     .arg(pixelsProcessedThisImage));
    }
    
    // Step 3: Handle sigma clipping or normalize
    if (useSigmaClipping) {
        updateProgress(80, "Performing sigma clipping rejection...");
        
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                int idx = y * m_output_size.width + x;
                std::vector<float>& pixelStack = pixelStacks[idx];
                std::vector<float>& weightStack = weightStacks[idx];
                
                if (!pixelStack.empty()) {
                    // Compute statistics
                    float sum = 0.0f, sumSq = 0.0f;
                    for (float val : pixelStack) {
                        sum += val;
                        sumSq += val * val;
                    }
                    
                    float mean = sum / pixelStack.size();
                    float variance = (sumSq / pixelStack.size()) - (mean * mean);
                    float stddev = sqrt(std::max(0.0f, variance));
                    
                    // Apply sigma clipping
                    float lowThreshold = mean - m_params.sigma_low * stddev;
                    float highThreshold = mean + m_params.sigma_high * stddev;
                    
                    float pixelSum = 0.0f;
                    float weightSum = 0.0f;
                    int validPixels = 0;
                    
                    for (size_t i = 0; i < pixelStack.size(); ++i) {
                        float pixelValue = pixelStack[i];
                        float weight = weightStack[i];
                        
                        if (pixelValue >= lowThreshold && pixelValue <= highThreshold) {
                            pixelSum += pixelValue * weight;
                            weightSum += weight;
                            validPixels++;
                        } else {
                            m_pixels_rejected++;
                        }
                    }
                    
                    if (weightSum > 0) {
                        m_stacked_image.at<float>(y, x) = pixelSum / weightSum;
                        m_weight_map.at<float>(y, x) = weightSum;
                        m_overlap_map.at<uchar>(y, x) = validPixels;
                    }
                }
            }
        }
    } else {
        updateProgress(80, "Normalizing stacked image...");
        
        // Final normalization - this should now produce uniform brightness
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                float totalWeight = m_weight_map.at<float>(y, x);
                if (totalWeight > 0) {
                    m_stacked_image.at<float>(y, x) /= totalWeight;
                }
            }
        }
    }
    
    updateProgress(90, "Applying final brightness compensation...");
    
    // OPTIONAL: Additional global brightness equalization
    // Calculate median brightness in high-overlap vs low-overlap regions
    cv::Mat highOverlapMask = (overlap_count >= 3.0f);
    cv::Mat lowOverlapMask = (overlap_count >= 1.0f) & (overlap_count < 3.0f);
    
    if (cv::countNonZero(highOverlapMask) > 100 && cv::countNonZero(lowOverlapMask) > 100) {
        cv::Scalar highOverlapMean = cv::mean(m_stacked_image, highOverlapMask);
        cv::Scalar lowOverlapMean = cv::mean(m_stacked_image, lowOverlapMask);
        
        if (highOverlapMean[0] > 0 && lowOverlapMean[0] > 0) {
            float brightness_ratio = lowOverlapMean[0] / highOverlapMean[0];
            logProcessing(QString("Brightness compensation ratio: %1").arg(brightness_ratio, 0, 'f', 3));
            
            // Apply gentle global correction if needed
            if (brightness_ratio < 0.8 || brightness_ratio > 1.2) {
                for (int y = 0; y < m_output_size.height; ++y) {
                    for (int x = 0; x < m_output_size.width; ++x) {
                        float overlap = overlap_count.at<float>(y, x);
                        if (overlap >= 3.0f) {
                            // Slightly reduce brightness in high-overlap regions
                            m_stacked_image.at<float>(y, x) *= 0.95f;
                        } else if (overlap >= 1.0f && overlap < 3.0f) {
                            // Slightly boost brightness in low-overlap regions
                            m_stacked_image.at<float>(y, x) *= 1.05f;
                        }
                    }
                }
                logProcessing("Applied gentle brightness equalization");
            }
        }
    }
    
    updateProgress(100, "TAN projection stacking successfully completed");
    emit stackingComplete(true);
    
    return true;
}

// Updated saveResult function
bool WCSAstrometricStacker::saveResult(const QString &output_path) {
    if (m_stacked_image.empty()) {
        return false;
    }
    
    updateProgress(0, "Saving TAN projection result...");
    
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = QString("!%1").arg(output_path).toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status)) {
        return false;
    }
    
    long naxes[2] = {m_stacked_image.cols, m_stacked_image.rows};
    if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status)) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Convert OpenCV Mat to FITS format
    long totalPixels = m_stacked_image.rows * m_stacked_image.cols;
    std::vector<float> pixels(totalPixels);
    
    for (int y = 0; y < m_stacked_image.rows; ++y) {
        for (int x = 0; x < m_stacked_image.cols; ++x) {
            pixels[y * m_stacked_image.cols + x] = m_stacked_image.at<float>(y, x);
        }
    }
    
    if (fits_write_img(fptr, TFLOAT, 1, totalPixels, pixels.data(), &status)) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Write complete WCS keywords for TAN projection
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &m_output_wcs.crval1, "Reference RA (degrees)", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &m_output_wcs.crval2, "Reference Dec (degrees)", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &m_output_wcs.crpix1, "Reference pixel X", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &m_output_wcs.crpix2, "Reference pixel Y", &status);
    
    // Write CD matrix
    fits_write_key(fptr, TDOUBLE, "CD1_1", &m_output_wcs.cd11, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD1_2", &m_output_wcs.cd12, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD2_1", &m_output_wcs.cd21, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD2_2", &m_output_wcs.cd22, "Coordinate matrix element", &status);
    
    // Write coordinate types
    char ctype1[] = "RA---TAN";
    char ctype2[] = "DEC--TAN";
    char* ctype1_ptr = ctype1;
    char* ctype2_ptr = ctype2;
    fits_write_key(fptr, TSTRING, "CTYPE1", &ctype1_ptr, "Coordinate type", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", &ctype2_ptr, "Coordinate type", &status);
    
    // Write coordinate units
    char cunit1[] = "deg";
    char cunit2[] = "deg";
    char* cunit1_ptr = cunit1;
    char* cunit2_ptr = cunit2;
    fits_write_key(fptr, TSTRING, "CUNIT1", &cunit1_ptr, "Coordinate unit", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", &cunit2_ptr, "Coordinate unit", &status);
    
    // Add processing information
    int nimages = m_images.size();
    fits_write_key(fptr, TINT, "NSTACKED", &nimages, "Number of stacked images", &status);
    
    double pixelScale = m_output_wcs.getPixelScale();
    fits_write_key(fptr, TDOUBLE, "PIXSCALE", &pixelScale, "Pixel scale (arcsec/pixel)", &status);
    
    // Add processing method
    QString method = QString("TAN_PROJECTION_%1").arg(static_cast<int>(m_params.combination));
    QByteArray methodBytes = method.toLocal8Bit();
    char* methodPtr = methodBytes.data();
    fits_write_key(fptr, TSTRING, "STACKMET", &methodPtr, "Stacking method", &status);
    
    // Add history
    QString history = QString("Stacked %1 images using TAN projection astrometric alignment").arg(nimages);
    QByteArray historyBytes = history.toLocal8Bit();
    fits_write_history(fptr, historyBytes.data(), &status);
    
    fits_close_file(fptr, &status);
    
    updateProgress(100, "Save complete");
    
    return (status == 0);
}

// Missing member function implementations

void WCSAstrometricStacker::setProgressWidgets(QProgressBar *progress, QLabel *status) {
    m_progress_bar = progress;
    m_status_label = status;
}

bool WCSAstrometricStacker::addImage(const QString &fits_file, const QString &solved_fits_file) {
    // Use the solved_fits_file if provided, otherwise use the fits_file
    QString fileToUse = solved_fits_file.isEmpty() ? fits_file : solved_fits_file;
    return addPlatesolveDFITSFile(fileToUse);
}

bool WCSAstrometricStacker::addImageWithMetadata(const QString &fits_file, const StellinaImageData &stellina_data) {
    return addImageFromStellinaData(fits_file, stellina_data);
}

void WCSAstrometricStacker::setStackingParameters(const StackingParams &params) {
    m_params = params;
    logProcessing(QString("Updated stacking parameters: method=%1, rejection=%2")
                 .arg(static_cast<int>(params.combination))
                 .arg(static_cast<int>(params.rejection)));
}

bool WCSAstrometricStacker::addPlatesolveDFITSFile(const QString &solved_fits_file) {
    auto img_data = std::make_unique<WCSImageData>();
    
    // Load FITS image and WCS from plate-solved file
    if (!loadWCSFromFITS(solved_fits_file, *img_data)) {
        emit errorOccurred(QString("Failed to load WCS from: %1").arg(solved_fits_file));
        return false;
    }
    
    // Extract basic metadata from FITS headers
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = solved_fits_file.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status) == 0) {
        // Try to read exposure time
        double exptime = 10.0;
        if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exptime, nullptr, &status) == 0) {
            img_data->exposure_time = exptime;
        }
        status = 0;
        
        // Try to read observation time
        char dateobs[FLEN_VALUE];
        if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) == 0) {
            QString dateStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
            img_data->obs_time = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss");
            if (!img_data->obs_time.isValid()) {
                // Try alternative formats
                img_data->obs_time = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss.zzz");
            }
        }
        
        fits_close_file(fptr, &status);
    }
    
    // Extract image statistics and compute quality score
    extractImageStatistics(*img_data);
    computeImageQualityScore(*img_data);
    
    m_images.push_back(std::move(img_data));
    
    emit imageProcessed(solved_fits_file);
    logProcessing(QString("Added FITS file: %1 (quality: %2)")
                 .arg(QFileInfo(solved_fits_file).fileName())
                 .arg(m_images.back()->quality_score, 0, 'f', 3));
    
    return true;
}

bool WCSAstrometricStacker::addImageFromStellinaData(const QString &fits_file, 
                                                    const StellinaImageData &stellina_data) {
    auto img_data = std::make_unique<WCSImageData>();
    
    // Load FITS image and WCS from plate-solved file
    if (!loadWCSFromFITS(fits_file, *img_data)) {
        emit errorOccurred(QString("Failed to load WCS from: %1").arg(fits_file));
        return false;
    }
    
    // Integrate Stellina-specific data
    img_data->exposure_time = stellina_data.exposureSeconds;
    img_data->obs_time = QDateTime::fromString(stellina_data.dateObs, "yyyy-MM-ddThh:mm:ss");
    if (!img_data->obs_time.isValid()) {
        img_data->obs_time = QDateTime::fromString(stellina_data.dateObs, "yyyy-MM-ddThh:mm:ss.zzz");
    }
    
    // Use Stellina quality metrics if available
    if (stellina_data.hasValidCoordinates) {
        img_data->stellina_correction_magnitude = sqrt(pow(stellina_data.altitude, 2) + 
                                                      pow(stellina_data.azimuth, 2));
    }
    
    // Extract image statistics and compute quality score
    extractImageStatistics(*img_data);
    computeImageQualityScore(*img_data);
    
    m_images.push_back(std::move(img_data));
    
    emit imageProcessed(fits_file);
    logProcessing(QString("Added Stellina image: %1 (quality: %2)")
                 .arg(QFileInfo(fits_file).fileName())
                 .arg(m_images.back()->quality_score, 0, 'f', 3));
    
    return true;
}

bool WCSAstrometricStacker::extractImageStatistics(WCSImageData &img_data) {
    if (img_data.image.empty()) return false;
    
    // Calculate basic statistics
    cv::Scalar mean, stddev;
    cv::meanStdDev(img_data.image, mean, stddev);
    
    img_data.background_level = mean[0];
    img_data.noise_level = stddev[0];
    
    // Estimate star count by counting bright pixels
    cv::Mat mask;
    double threshold = img_data.background_level + 3 * img_data.noise_level;
    cv::threshold(img_data.image, mask, threshold, 255, cv::THRESH_BINARY);
    
    // Convert to 8-bit for contour detection
    cv::Mat mask8;
    mask.convertTo(mask8, CV_8U);
    
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(mask8, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    
    img_data.star_count = 0;
    for (const auto &contour : contours) {
        double area = cv::contourArea(contour);
        if (area > 10 && area < 1000) {  // Reasonable star size range
            img_data.star_count++;
        }
    }
    
    return true;
}

bool WCSAstrometricStacker::computeImageQualityScore(WCSImageData &img_data) {
    double quality = 1.0;
    
    // Factor in noise level (lower noise = higher quality)
    if (img_data.noise_level > 0) {
        double snr = img_data.background_level / img_data.noise_level;
        quality *= std::min(1.0, snr / 50.0);  // Normalize to reasonable SNR range
    }
    
    // Factor in star count (more stars = better registration)
    if (img_data.star_count > 0) {
        quality *= std::min(1.0, img_data.star_count / 200.0);
    }
    
    // Factor in Stellina-specific metrics
    if (img_data.stellina_stars_used > 0) {
        quality *= std::min(1.0, img_data.stellina_stars_used / 100.0);
    }
    
    // Reduce quality for images with large correction magnitudes
    if (img_data.stellina_correction_magnitude > 100.0) {
        quality *= 0.5;  // Penalize images that needed large corrections
    }
    
    img_data.quality_score = std::max(0.1, std::min(1.0, quality));  // Clamp to [0.1, 1.0]
    
    return true;
}

void WCSAstrometricStacker::startStacking() {
    if (m_images.empty()) {
        emit errorOccurred("No images loaded for stacking");
        return;
    }
    
    m_stacking_active = true;
    m_current_image_index = 0;
    
    updateProgress(0, "Starting TAN projection stacking...");
    
    if (!stackImages()) {
        emit errorOccurred("Stacking process failed");
        return;
    }
    
    analyzeImageQuality();
    emit qualityAnalysisComplete();
}

void WCSAstrometricStacker::cancelStacking() {
    m_stacking_active = false;
    m_processing_timer->stop();
    emit statusUpdated("Stacking cancelled by user");
}

void WCSAstrometricStacker::analyzeImageQuality() {
    // Quality analysis is already done during image loading
    // This could be expanded to do more detailed analysis
    logProcessing("Image quality analysis complete");
    
    if (m_images.empty()) return;
    
    // Calculate aggregate statistics
    double totalQuality = 0.0;
    double totalExposure = 0.0;
    int totalStars = 0;
    
    for (const auto& img : m_images) {
        totalQuality += img->quality_score;
        totalExposure += img->exposure_time;
        totalStars += img->star_count;
    }
    
    double avgQuality = totalQuality / m_images.size();
    double avgStars = double(totalStars) / m_images.size();
    
    logProcessing(QString("Quality analysis: %1 images, avg quality: %2, avg stars: %3, total exposure: %4s")
                 .arg(m_images.size())
                 .arg(avgQuality, 0, 'f', 3)
                 .arg(avgStars, 0, 'f', 1)
                 .arg(totalExposure, 0, 'f', 1));
}

QString WCSAstrometricStacker::getQualityReport() const {
    QStringList report;
    
    report << "=== TAN Projection Astrometric Stacking Quality Report ===";
    report << "";
    report << QString("Total images processed: %1").arg(m_images.size());
    
    if (!m_stacked_image.empty()) {
        report << QString("Output dimensions: %1 x %2 pixels").arg(m_output_size.width).arg(m_output_size.height);
    }
    
    if (m_output_wcs.valid) {
        report << QString("Output pixel scale: %1 arcsec/pixel").arg(m_output_wcs.getPixelScale(), 0, 'f', 2);
        report << QString("Field center: RA=%1°, Dec=%2°").arg(m_output_wcs.crval1, 0, 'f', 6).arg(m_output_wcs.crval2, 0, 'f', 6);
    }
    
    report << QString("Total exposure time: %1 seconds").arg(getTotalExposureTime(), 0, 'f', 1);
    report << QString("Average image quality: %1").arg(getAverageQuality(), 0, 'f', 3);
    
    if (m_pixels_rejected > 0) {
        report << QString("Pixels rejected by sigma clipping: %1").arg(m_pixels_rejected);
    }
    
    report << "";
    report << "Individual Image Quality:";
    
    for (size_t i = 0; i < m_images.size(); ++i) {
        const auto &img = m_images[i];
        report << QString("  %1: Quality=%2, Stars=%3, Exposure=%4s")
                     .arg(QFileInfo(img->filename).fileName())
                     .arg(img->quality_score, 0, 'f', 3)
                     .arg(img->star_count)
                     .arg(img->exposure_time, 0, 'f', 1);
    }
    
    report << "";
    report << "Processing Log:";
    for (const QString &entry : m_processing_log) {
        report << QString("  %1").arg(entry);
    }
    
    return report.join("\n");
}

double WCSAstrometricStacker::getTotalExposureTime() const {
    double total = 0.0;
    for (const auto &img : m_images) {
        total += img->exposure_time;
    }
    return total;
}

double WCSAstrometricStacker::getAverageQuality() const {
    if (m_images.empty()) return 0.0;
    
    double total = 0.0;
    for (const auto &img : m_images) {
        total += img->quality_score;
    }
    return total / m_images.size();
}

void WCSAstrometricStacker::updateProgress(int percentage, const QString &message) {
    if (m_progress_bar) {
        m_progress_bar->setValue(percentage);
    }
    if (m_status_label) {
        m_status_label->setText(message);
    }
    
    emit progressUpdated(percentage);
    emit statusUpdated(message);
    logProcessing(message);
}

void WCSAstrometricStacker::logProcessing(const QString &message) {
    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    m_processing_log.append(QString("[%1] %2").arg(timestamp).arg(message));
}

void WCSAstrometricStacker::finishStacking() {
    m_stacking_active = false;
    emit stackingComplete(true);
    emit statusUpdated("TAN projection stacking completed successfully");
}

// Add this diagnostic function to WCSAstrometricStacker class
// Call it after the first pass to understand overlap patterns

void WCSAstrometricStacker::analyzeOverlapDistribution(const cv::Mat& overlap_count) {
    // Calculate overlap statistics
    double minOverlap, maxOverlap;
    cv::minMaxLoc(overlap_count, &minOverlap, &maxOverlap);
    
    cv::Scalar meanOverlap = cv::mean(overlap_count, overlap_count > 0);
    
    // Count pixels by overlap level
    int noOverlap = cv::countNonZero(overlap_count == 0);
    int singleOverlap = cv::countNonZero(overlap_count == 1);
    int doubleOverlap = cv::countNonZero(overlap_count == 2);
    int tripleOverlap = cv::countNonZero(overlap_count >= 3);
    int totalPixels = m_output_size.area();
    
    logProcessing("=== OVERLAP ANALYSIS ===");
    logProcessing(QString("Overlap range: %1 to %2 images per pixel")
                 .arg(minOverlap, 0, 'f', 1)
                 .arg(maxOverlap, 0, 'f', 1));
    logProcessing(QString("Mean overlap: %1 images per pixel")
                 .arg(meanOverlap[0], 0, 'f', 2));
    
    logProcessing(QString("Pixel distribution:"));
    logProcessing(QString("  No coverage: %1 pixels (%2%)")
                 .arg(noOverlap)
                 .arg(noOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Single image: %1 pixels (%2%)")
                 .arg(singleOverlap)
                 .arg(singleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Two images: %1 pixels (%2%)")
                 .arg(doubleOverlap)
                 .arg(doubleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Three+ images: %1 pixels (%2%)")
                 .arg(tripleOverlap)
                 .arg(tripleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    
    // This helps identify if your brightness problem correlates with overlap
    if (singleOverlap > totalPixels * 0.3) {
        logProcessing("WARNING: High percentage of single-overlap pixels may cause brightness variations");
    }
    
    if (tripleOverlap > totalPixels * 0.5) {
        logProcessing("INFO: Good overlap coverage - applying compensation for center over-brightening");
    }
}

// Also add this method to save overlap map for visual inspection
void WCSAstrometricStacker::saveOverlapMap(const cv::Mat& overlap_count, const QString& output_path) {
    // Convert overlap count to 8-bit for visualization
    cv::Mat overlap_vis;
    overlap_count.convertTo(overlap_vis, CV_8U, 255.0 / 5.0); // Scale assuming max 5 overlaps
    
    // Save as FITS for quantitative analysis
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = QString("!%1").arg("_overlap.fits").toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status) == 0) {
        long naxes[2] = {overlap_count.cols, overlap_count.rows};
        if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status) == 0) {
            
            std::vector<float> pixels(overlap_count.rows * overlap_count.cols);
            for (int y = 0; y < overlap_count.rows; ++y) {
                for (int x = 0; x < overlap_count.cols; ++x) {
                    pixels[y * overlap_count.cols + x] = overlap_count.at<float>(y, x);
                }
            }
            
            if (fits_write_img(fptr, TFLOAT, 1, pixels.size(), pixels.data(), &status) == 0) {
                logProcessing(QString("Saved overlap map"));
            }
        }
        fits_close_file(fptr, &status);
    }
}
