// WcsAstrometricStacker.cpp
// Implementation of the WCS Astrometric Stacker

#include "WcsAstrometricStacker.h"
#include <QFileInfo>
#include <QTimer>
#include <QDateTime>
#include <QApplication>
#include <QTextStream>
#include <QRegularExpression>
#include <cmath>
#include <algorithm>

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
    // wcsini(1, 2, &m_output_wcs);
    memset(&m_output_wcs, 0, sizeof(m_output_wcs));
    
    // Set up processing timer for non-blocking operation
    m_processing_timer->setSingleShot(true);
    connect(m_processing_timer, &QTimer::timeout, this, &WCSAstrometricStacker::processNextImage);
}

WCSAstrometricStacker::~WCSAstrometricStacker() {
    wcsfree(&m_output_wcs);
}

bool WCSAstrometricStacker::addImage(const QString &fits_file, const QString &solved_fits_file) {
    // Use the solved_fits_file if provided, otherwise use the fits_file
    QString fileToUse = solved_fits_file.isEmpty() ? fits_file : solved_fits_file;
    return addPlatesolveDFITSFile(fileToUse);
}

bool WCSAstrometricStacker::addImageWithMetadata(const QString &fits_file, const StellinaImageData &stellina_data) {
    return addImageFromStellinaData(fits_file, stellina_data);
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

void WCSAstrometricStacker::setStackingParameters(const StackingParams &params) {
    m_params = params;
    logProcessing(QString("Updated stacking parameters: method=%1, rejection=%2")
                 .arg(static_cast<int>(params.combination))
                 .arg(static_cast<int>(params.rejection)));
}

void WCSAstrometricStacker::setProgressWidgets(QProgressBar *progress, QLabel *status) {
    m_progress_bar = progress;
    m_status_label = status;
}

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
    
    // First, get the number of header keywords properly
    int nkeys = 0;
    if (fits_get_hdrspace(fptr, &nkeys, nullptr, &status)) {
        emit errorOccurred(QString("Failed to get header space: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    if (nkeys == 0) {
        emit errorOccurred(QString("FITS file has no header keywords: %1").arg(fits_file));
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read WCS header with proper parameters
    char *header = nullptr;
    int nkeysret = 0;
    
    // Use fits_hdr2str with correct parameters:
    // - excludecomm = 1 (exclude comment and history)
    // - nkeywords = 0 (return all keywords)
    if (fits_hdr2str(fptr, 1, nullptr, 0, &header, &nkeysret, &status)) {
        emit errorOccurred(QString("Failed to read FITS header: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    if (!header || nkeysret == 0) {
        emit errorOccurred(QString("Empty FITS header returned: %1 (nkeys: %2)").arg(fits_file).arg(nkeysret));
        if (header) free(header);
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Debug output to understand what we got
    logProcessing(QString("FITS header read: %1 keywords from %2").arg(nkeysret).arg(QFileInfo(fits_file).fileName()));
    
    // Parse WCS using wcspih
    int nreject = 0, nwcs = 0;
    struct wcsprm *wcs = nullptr;
    
    // Use wcspih to parse WCS from header string
    int wcspih_status = wcspih(header, nkeysret, WCSHDR_all, 2, &nreject, &nwcs, &wcs);
    
    if (wcspih_status) {
        emit errorOccurred(QString("WCS parsing failed: %1 (wcspih status: %2, nwcs: %3, nreject: %4)")
                          .arg(fits_file).arg(wcspih_status).arg(nwcs).arg(nreject));
        free(header);
        fits_close_file(fptr, &status);
        return false;
    }
    
    if (nwcs == 0) {
        emit errorOccurred(QString("No WCS found in FITS file: %1 (rejected: %2)").arg(fits_file).arg(nreject));
        free(header);
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Copy the first WCS structure
    img_data.wcs = wcs[0];
    if (wcsset(&img_data.wcs)) {
        emit errorOccurred(QString("WCS setup failed: %1").arg(fits_file));
        free(header);
        wcsvfree(&nwcs, &wcs);
        fits_close_file(fptr, &status);
        return false;
    }
    
    img_data.wcs_valid = true;
    
    // Clean up WCS memory
    wcsvfree(&nwcs, &wcs);
    free(header);
    
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
    
    // Convert to OpenCV Mat (note: FITS is row-major, OpenCV expects this)
    img_data.image = cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
    
    fits_close_file(fptr, &status);
    
    logProcessing(QString("Successfully loaded WCS and image data from: %1").arg(QFileInfo(fits_file).fileName()));
    return (status == 0);
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
    cv::threshold(img_data.image, mask, img_data.background_level + 3 * img_data.noise_level, 
                  255, cv::THRESH_BINARY);
    
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
    quality *= std::min(1.0, img_data.star_count / 200.0);
    
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
}

QString WCSAstrometricStacker::getQualityReport() const {
    QStringList report;
    
    report << "=== WCS Astrometric Stacking Quality Report ===";
    report << "";
    report << QString("Total images processed: %1").arg(m_images.size());
    report << QString("Output dimensions: %1 x %2 pixels").arg(m_output_size.width).arg(m_output_size.height);
    report << QString("Output pixel scale: %1 arcsec/pixel").arg(m_output_pixel_scale, 0, 'f', 2);
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

/*
* wcsp2s() - Pixel-to-world transformation
* ----------------------------------------
* wcsp2s() transforms pixel coordinates to world coordinates.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length
*                       nelem but containing wcs.naxis coordinate elements.
*                       Thus nelem must equal or exceed the value of the
*                       NAXIS keyword unless ncoord == 1, in which case nelem
*                       is not used.
*
*   pixcrd    const double[ncoord][nelem]
*                       Array of pixel coordinates.
*
* Returned:
*   imgcrd    double[ncoord][nelem]
*                       Array of intermediate world coordinates.  For
*                       celestial axes, imgcrd[][wcs.lng] and
*                       imgcrd[][wcs.lat] are the projected x-, and
*                       y-coordinates in pseudo "degrees".  For spectral
*                       axes, imgcrd[][wcs.spec] is the intermediate spectral
*                       coordinate, in SI units.  For time axes,
*                       imgcrd[][wcs.time] is the intermediate time
*                       coordinate.
*
*   phi,theta double[ncoord]
*                       Longitude and latitude in the native coordinate system
*                       of the projection [deg].
*
*   world     double[ncoord][nelem]
*                       Array of world coordinates.  For celestial axes,
*                       world[][wcs.lng] and world[][wcs.lat] are the
*                       celestial longitude and latitude [deg].  For spectral
*                       axes, world[][wcs.spec] is the spectral coordinate, in
*                       SI units.  For time axes, world[][wcs.time] is the
*                       time coordinate.
*
*   stat      int[ncoord]
*                       Status return value for each coordinate:
*                         0: Success.
*                        1+: A bit mask indicating invalid pixel coordinate
*                            element(s).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Linear transformation matrix is singular.
*                         4: Inconsistent or unrecognized coordinate axis
*                            types.
*                         5: Invalid parameter value.
*                         6: Invalid coordinate transformation parameters.
*                         7: Ill-conditioned coordinate transformation
*                            parameters.
*                         8: One or more of the pixel coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*/

int wcsp2s(struct wcsprm *wcs, int ncoord, int nelem, const double pixcrd[],
           double imgcrd[], double phi[], double theta[], double world[],
           int stat[]);

bool WCSAstrometricStacker::stackImages() {
    updateProgress(0, "Starting WCS astrometric stacking...");
    
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
    
    // Prepare for pixel rejection if using sigma clipping
    std::vector<cv::Mat> pixelStacks;
    std::vector<cv::Mat> weightStacks;
    bool useSigmaClipping = (m_params.rejection == StackingParams::SIGMA_CLIPPING);
    
    if (useSigmaClipping) {
        pixelStacks.resize(m_output_size.height * m_output_size.width);
        weightStacks.resize(m_output_size.height * m_output_size.width);
    }
    
    // Step 2: Reproject and combine each image
    m_pixels_processed = 0;
    m_pixels_rejected = 0;
    
    for (size_t i = 0; i < m_images.size(); ++i) {
        updateProgress(10 + (i * 80) / m_images.size(), 
                       QString("Processing image %1 of %2: %3")
                       .arg(i+1)
                       .arg(m_images.size())
                       .arg(QFileInfo(m_images[i]->filename).fileName()));
        
        // Get current image data
        const auto& img = m_images[i];
        if (!img->wcs_valid || img->image.empty()) {
            logProcessing(QString("Skipping invalid image: %1").arg(QFileInfo(img->filename).fileName()));
            continue;
        }
        
        // Create temporary reprojected image and weight map
        cv::Mat reprojected = cv::Mat::zeros(m_output_size, CV_32F);
        cv::Mat weights = cv::Mat::zeros(m_output_size, CV_32F);
        
        // Calculate image weight based on quality and exposure
        float imageWeight = img->quality_score;
        if (m_params.normalize_exposure && img->exposure_time > 0) {
            imageWeight *= img->exposure_time;
        }
        
        // Loop through output pixels and find corresponding input pixels
        #pragma omp parallel for
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                // Convert output pixel (x,y) to world coordinates (RA,Dec)
                double pixcrd[2] = {x + 1.0, y + 1.0};  // FITS is 1-indexed
                double world[2] = {0.0, 0.0};
		double imgcrd[2] = {0.0, 0.0};    // Add missing
                double phi = 0.0, theta = 0.0;
                int stat[1] = {0};                // Add missing
		
                // Use WCS to transform from pixel to world coordinates
                if (wcsp2s(&m_output_wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat) == 0) {
                    // Now convert world coordinates to input image pixel coordinates
                    double imgpix[2] = {0.0, 0.0};
                    double imgcrd2[2] = {0.0, 0.0}; // Add missing
		    int stat2[1] = {0};             // Add missing
    
                    if (wcss2p(&img->wcs, 1, 2, world, &phi, &theta, imgcrd2, imgpix, stat2) == 0) {
                        // Check if the pixel is within the input image bounds
                        if (imgpix[0] >= 1.0 && imgpix[0] <= img->image.cols &&
                            imgpix[1] >= 1.0 && imgpix[1] <= img->image.rows) {
                            
                            // Convert to 0-indexed for OpenCV
                            float srcX = imgpix[0] - 1.0;
                            float srcY = imgpix[1] - 1.0;
                            
                            // Perform bilinear interpolation to get pixel value
                            float pixelValue = 0.0f;
                            
                            // Bilinear interpolation
                            int x0 = floor(srcX);
                            int y0 = floor(srcY);
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
                                
                                pixelValue = v0 * (1 - dy) + v1 * dy;
                                
                                // Skip if pixel is invalid (NaN or Inf)
                                if (std::isfinite(pixelValue)) {
                                    if (useSigmaClipping) {
                                        // For sigma clipping, store all values for later processing
                                        int idx = y * m_output_size.width + x;
                                        pixelStacks[idx].push_back(pixelValue);
                                        weightStacks[idx].push_back(imageWeight);
                                    } else {
                                        // For direct stacking methods
                                        #pragma omp critical
                                        {
                                            reprojected.at<float>(y, x) = pixelValue;
                                            weights.at<float>(y, x) = imageWeight;
                                        }
                                    }
                                    
                                    #pragma omp atomic
                                    m_pixels_processed++;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // If not using sigma clipping, add this image's contribution to the stack
        if (!useSigmaClipping) {
            // Combine this image with the stack
            for (int y = 0; y < m_output_size.height; ++y) {
                for (int x = 0; x < m_output_size.width; ++x) {
                    float w = weights.at<float>(y, x);
                    if (w > 0) {
                        m_stacked_image.at<float>(y, x) += reprojected.at<float>(y, x) * w;
                        m_weight_map.at<float>(y, x) += w;
                        m_overlap_map.at<uchar>(y, x)++;
                    }
                }
            }
        }
        
        logProcessing(QString("Reprojected image %1: %2 pixels processed")
                     .arg(i+1)
                     .arg(m_pixels_processed));
    }
    
    // Step 3: Perform pixel rejection if using sigma clipping
    if (useSigmaClipping) {
        updateProgress(90, "Performing sigma clipping rejection...");
        
        // Process each pixel stack
        #pragma omp parallel for
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                int idx = y * m_output_size.width + x;
                cv::Mat& pixelStack = pixelStacks[idx];
                cv::Mat& weightStack = weightStacks[idx];
                
                if (!pixelStack.empty()) {
                    // Compute mean and standard deviation
                    cv::Scalar mean, stddev;
                    cv::meanStdDev(pixelStack, mean, stddev);
                    
                    float pixelMean = mean[0];
                    float pixelStd = stddev[0];
                    
                    // Apply sigma clipping
                    float lowThreshold = pixelMean - m_params.sigma_low * pixelStd;
                    float highThreshold = pixelMean + m_params.sigma_high * pixelStd;
                    
                    float pixelSum = 0.0f;
                    float weightSum = 0.0f;
                    int validPixels = 0;
                    
                    for (int i = 0; i < pixelStack.rows; ++i) {
                        float pixelValue = pixelStack.at<float>(i);
                        float weight = weightStack.at<float>(i);
                        
                        if (pixelValue >= lowThreshold && pixelValue <= highThreshold) {
                            pixelSum += pixelValue * weight;
                            weightSum += weight;
                            validPixels++;
                        } else {
                            #pragma omp atomic
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
        
        logProcessing(QString("Sigma clipping completed: rejected %1 of %2 pixels (%3%)")
                     .arg(m_pixels_rejected)
                     .arg(m_pixels_processed)
                     .arg(m_pixels_processed > 0 ? (m_pixels_rejected * 100.0 / m_pixels_processed) : 0.0, 0, 'f', 1));
    } else {
        // Step 3 (alternative): Normalize by weights for non-sigma-clipping methods
        updateProgress(90, "Normalizing stacked image...");
        
        // Combine based on selected method
        switch (m_params.combination) {
            case StackingParams::WEIGHTED_MEAN:
                // Already handled during stacking (weighted accumulation)
                // Just need to normalize by total weight
                for (int y = 0; y < m_output_size.height; ++y) {
                    for (int x = 0; x < m_output_size.width; ++x) {
                        if (m_weight_map.at<float>(y, x) > 0) {
                            m_stacked_image.at<float>(y, x) /= m_weight_map.at<float>(y, x);
                        }
                    }
                }
                break;
                
            case StackingParams::MEDIAN:
                // This would be handled differently - not implemented in this simplified version
                // Would require storing all pixel values for each output pixel
                logProcessing("Warning: Median stacking not fully implemented in this version");
                break;
                
            default:
                // Default to weighted mean
                for (int y = 0; y < m_output_size.height; ++y) {
                    for (int x = 0; x < m_output_size.width; ++x) {
                        if (m_weight_map.at<float>(y, x) > 0) {
                            m_stacked_image.at<float>(y, x) /= m_weight_map.at<float>(y, x);
                        }
                    }
                }
                break;
        }
    }
    
    // Step 4: Final post-processing
    updateProgress(95, "Performing final image adjustments...");
    
    // Create a mask for pixels with no data
    cv::Mat mask = (m_weight_map > 0);
    
    // Normalize image to [0, 1] range for display
    double minVal, maxVal;
    cv::minMaxLoc(m_stacked_image, &minVal, &maxVal, nullptr, nullptr, mask);
    
    // Log statistics
    logProcessing(QString("Stacking complete - Image statistics: min=%1, max=%2, mean=%3")
                 .arg(minVal, 0, 'f', 2)
                 .arg(maxVal, 0, 'f', 2)
                 .arg(cv::mean(m_stacked_image, mask)[0], 0, 'f', 2));
    
    // Optional: Stretch for better contrast
    if (maxVal > minVal) {
        logProcessing("Applying linear stretch for improved contrast");
        cv::Mat temp;
        m_stacked_image.copyTo(temp, mask);
        temp = (temp - minVal) / (maxVal - minVal);
        temp.copyTo(m_stacked_image, mask);
    }
    
    // Create noise map if requested
    if (m_params.create_weight_map) {
        // Convert weight map to noise estimate (approximate)
        m_noise_map = cv::Mat::zeros(m_output_size, CV_32F);
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                float weight = m_weight_map.at<float>(y, x);
                if (weight > 0) {
                    // Noise decreases with sqrt(N) where N is effective number of images
                    m_noise_map.at<float>(y, x) = 1.0f / sqrt(weight);
                }
            }
        }
    }
    
    updateProgress(100, "Stacking successfully completed");
    emit stackingComplete(true);
    
    return true;
}

/*
* wcsinit() - Default constructor for the wcsprm struct
* -----------------------------------------------------
* wcsinit() optionally allocates memory for arrays in a wcsprm struct and sets
* all members of the struct to default values.
*
* PLEASE NOTE: every wcsprm struct should be initialized by wcsinit(),
* possibly repeatedly.  On the first invokation, and only the first
* invokation, wcsprm::flag must be set to -1 to initialize memory management,
* regardless of whether wcsinit() will actually be used to allocate memory.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for the
*                       crpix, etc. arrays.  Please note that memory is never
*                       allocated by wcsinit() for the auxprm, tabprm, nor
*                       wtbarr structs.
*
*                       If false, it is assumed that pointers to these arrays
*                       have been set by the user except if they are null
*                       pointers in which case memory will be allocated for
*                       them regardless.  (In other words, setting alloc true
*                       saves having to initalize these pointers to zero.)
*
*   naxis     int       The number of world coordinate axes.  This is used to
*                       determine the length of the various wcsprm vectors and
*                       matrices and therefore the amount of memory to
*                       allocate for them.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
*                       Note that, in order to initialize memory management,
*                       wcsprm::flag should be set to -1 when wcs is
*                       initialized for the first time (memory leaks may
*                       result if it had already been initialized).
*
* Given:
*   npvmax    int       The number of PVi_ma keywords to allocate space for.
*                       If set to -1, the value of the global variable NPVMAX
*                       will be used.  This is potentially thread-unsafe if
*                       wcsnpv() is being used dynamically to alter its value.
*
*   npsmax    int       The number of PSi_ma keywords to allocate space for.
*                       If set to -1, the value of the global variable NPSMAX
*                       will be used.  This is potentially thread-unsafe if
*                       wcsnps() is being used dynamically to alter its value.
*
*   ndpmax    int       The number of DPja or DQia keywords to allocate space
*                       for.  If set to -1, the value of the global variable
*                       NDPMAX will be used.  This is potentially
*                       thread-unsafe if disndp() is being used dynamically to
*                       alter its value.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*/
int wcsinit(int alloc, int naxis, struct wcsprm *wcs, int npvmax, int npsmax,
            int ndpmax);

bool WCSAstrometricStacker::computeOptimalWCS() {
    if (m_images.empty()) {
        emit errorOccurred("No images available for WCS calculation");
        return false;
    }
    
    logProcessing("Computing optimal output WCS...");
    
    // Initialize WCS structure
    int naxis = 2;
    int status = wcsinit(true, naxis, &m_output_wcs, -1, -1, -1);
    if (status) {
        emit errorOccurred(QString("Failed to initialize WCS structure (error: %1)").arg(status));
        return false;
    }
    
    // Calculate bounding box in world coordinates
    double ra_min = 360.0, ra_max = 0.0;
    double dec_min = 90.0, dec_max = -90.0;
    
    // Use the first valid WCS as reference
    struct wcsprm* ref_wcs = nullptr;
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
    
    // For each image, find the corner points in world coordinates
    for (const auto& img : m_images) {
        if (!img->wcs_valid) continue;
        wcsprm *wcs = &(img->wcs);

        // Get image corners in pixel coordinates (1-indexed for WCSLIB)
        double corners_pix[4][2] = {
            {1.0, 1.0},  // Bottom left
            {double(img->image.cols), 1.0},  // Bottom right
            {1.0, double(img->image.rows)},  // Top left
            {double(img->image.cols), double(img->image.rows)}  // Top right
        };
        
        // Convert each corner to world coordinates
        for (int i = 0; i < 4; ++i) {
            double world[2] = {0.0, 0.0};
	    double imgcrd[2] = {0.0, 0.0};  // Add missing intermediate coordinates
            double phi = 0.0, theta = 0.0;
            int stat[1] = {0};              // Add missing status array
	    
            if (wcsp2s(wcs, 1, 2, corners_pix[i], imgcrd, &phi, &theta, world, stat) == 0) {
                // Handle RA wrapping at 0/360 degrees
                if (i == 0) {
                    ra_min = ra_max = world[0];
                } else {
                    // Check if we cross the RA=0/360 boundary
                    if (std::abs(world[0] - ra_min) > 180.0) {
                        if (world[0] < ra_min) world[0] += 360.0;
                        else ra_min += 360.0;
                    }
                    if (std::abs(world[0] - ra_max) > 180.0) {
                        if (world[0] < ra_max) world[0] += 360.0;
                        else ra_max += 360.0;
                    }
                    
                    ra_min = std::min(ra_min, world[0]);
                    ra_max = std::max(ra_max, world[0]);
                }
                
                dec_min = std::min(dec_min, world[1]);
                dec_max = std::max(dec_max, world[1]);
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
    
    // Determine pixel scale - use the first image as reference if not overridden
    if (m_params.output_pixel_scale <= 0.0) {
        // Calculate pixel scale from reference WCS
        // For simplicity, estimate from CD matrix if available
        if (ref_wcs->altlin & 2) { // Has PC and CDELT
            m_output_pixel_scale = std::abs(ref_wcs->cdelt[0] * 3600.0); // Convert degrees to arcsec
        } else if (ref_wcs->altlin & 4) { // Has CD
            m_output_pixel_scale = std::sqrt(ref_wcs->cd[0] * ref_wcs->cd[0] + 
                                           ref_wcs->cd[2] * ref_wcs->cd[2]) * 3600.0;
        } else {
            // Default for Stellina
            m_output_pixel_scale = 1.25; // arcseconds per pixel
        }
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
    
    width = int(ra_span_deg * 3600.0 / (m_output_pixel_scale * cos_dec_factor));
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
    
    m_output_wcs.naxis = 2;
    m_output_wcs.crpix[0] = width / 2.0 + 0.5;  // Center pixel X (1-indexed)
    m_output_wcs.crpix[1] = height / 2.0 + 0.5; // Center pixel Y (1-indexed)
    m_output_wcs.crval[0] = ra_center;          // RA at reference pixel
    m_output_wcs.crval[1] = dec_center;         // Dec at reference pixel
    m_output_wcs.cdelt[0] = -m_output_pixel_scale / 3600.0; // RA step (degrees)
    m_output_wcs.cdelt[1] = m_output_pixel_scale / 3600.0;  // Dec step (degrees)
    
    // Set projection to TAN (tangent plane)
    strncpy(m_output_wcs.ctype[0], "RA---TAN", 9);
    strncpy(m_output_wcs.ctype[1], "DEC--TAN", 9);
    
    // Initialize PC matrix as identity
    m_output_wcs.pc[0] = 1.0;
    m_output_wcs.pc[1] = 0.0;
    m_output_wcs.pc[2] = 0.0;
    m_output_wcs.pc[3] = 1.0;
    
    // Set up the rest of the WCS structure
    m_output_wcs.lonpole = 180.0;
    m_output_wcs.latpole = 90.0;
    m_output_wcs.restfrq = 0.0;
    m_output_wcs.restwav = 0.0;
    
    // Set alt lin flags to indicate PC and CDELT are present
    m_output_wcs.altlin = 2;
    
    // Set up the rest of the WCS struct that we need
    m_output_wcs.flag = 0;
    
    // Set units to degrees
    strncpy(m_output_wcs.cunit[0], "deg", 4);
    strncpy(m_output_wcs.cunit[1], "deg", 4);
    
    // Initialize the WCS system
    status = wcsset(&m_output_wcs);
    if (status) {
        emit errorOccurred(QString("Failed to set up WCS structure (error: %1)").arg(status));
        return false;
    }
    
    logProcessing(QString("Output WCS computed: center RA=%1°, Dec=%2°")
                 .arg(ra_center, 0, 'f', 6)
                 .arg(dec_center, 0, 'f', 6));
    
    return true;
}

bool WCSAstrometricStacker::saveResult(const QString &output_path) {
    if (m_stacked_image.empty()) {
        return false;
    }
    
    updateProgress(0, "Saving result...");
    
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
    
    // Write basic WCS keywords
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &m_output_wcs.crval[0], "Reference RA", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &m_output_wcs.crval[1], "Reference Dec", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &m_output_wcs.crpix[0], "Reference pixel X", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &m_output_wcs.crpix[1], "Reference pixel Y", &status);
    
    char ctype1[] = "RA---TAN";
    char ctype2[] = "DEC--TAN";
    char* ctype1_ptr = ctype1;
    char* ctype2_ptr = ctype2;
    fits_write_key(fptr, TSTRING, "CTYPE1", &ctype1_ptr, "Coordinate type", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", &ctype2_ptr, "Coordinate type", &status);
    
    // Add processing information
    int nimages = m_images.size();
    fits_write_key(fptr, TINT, "NSTACKED", &nimages, "Number of stacked images", &status);
    
    fits_close_file(fptr, &status);
    
    updateProgress(100, "Save complete");
    
    return (status == 0);
}

void WCSAstrometricStacker::processNextImage() {
    // This could be used for progressive processing if needed
    // Currently the main stacking is done synchronously in stackImages()
}

void WCSAstrometricStacker::finishStacking() {
    emit stackingComplete(true);
    emit statusUpdated("Stacking completed successfully");
}

// MOC file will be automatically generated and linked by qmake
