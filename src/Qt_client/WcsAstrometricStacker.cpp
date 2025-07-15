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

// Simplified implementation of key methods for demonstration
// In a full implementation, these would contain the complete WCS reprojection logic

bool WCSAstrometricStacker::computeOptimalWCS() {
    // Simplified - just use the first image's WCS for now
    if (m_images.empty()) return false;
    
    m_output_wcs = m_images[0]->wcs;
    m_output_size = cv::Size(m_images[0]->image.cols, m_images[0]->image.rows);
    m_output_pixel_scale = 1.25; // Default pixel scale for Stellina
    
    return true;
}

bool WCSAstrometricStacker::stackImages() {
    updateProgress(0, "Starting WCS stacking...");
    
    if (!computeOptimalWCS()) {
        return false;
    }
    
    updateProgress(50, "Combining images...");
    
    // Simplified stacking - just average the images for demonstration
    m_stacked_image = cv::Mat::zeros(m_output_size, CV_32F);
    m_weight_map = cv::Mat::zeros(m_output_size, CV_32F);
    
    for (const auto &img : m_images) {
        m_stacked_image += img->image * img->quality_score;
        m_weight_map += img->quality_score;
    }
    
    // Normalize by total weight
    cv::divide(m_stacked_image, m_weight_map, m_stacked_image);
    
    updateProgress(100, "Stacking complete");
    emit stackingComplete(true);
    
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
