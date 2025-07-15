// Complete WCS-based Astrometric Stacker with Overlap Intensity Correction
// Add this to StellinaProcessor project

#ifndef WCS_ASTROMETRIC_STACKER_H
#define WCS_ASTROMETRIC_STACKER_H

#include <QObject>
#include <QString>
#include <QList>
#include <QDateTime>
#include <QProgressBar>
#include <QLabel>
#include <opencv2/opencv.hpp>
#include <wcslib/wcs.h>
#include <fitsio.h>
#include <vector>
#include <memory>
#include "StellinaProcessor.h"

// Enhanced image data structure for WCS stacking
struct WCSImageData {
    cv::Mat image;                    // Image pixel data (32-bit float)
    cv::Mat weight_map;               // Quality weight map
    struct wcsprm wcs;                // WCS coordinate system
    QString filename;                 // Original filename
    QString solved_filename;          // Plate-solved FITS file
    double quality_score;             // Overall quality (0-1)
    double exposure_time;             // Exposure in seconds
    int star_count;                   // Number of detected stars
    double background_level;          // Background ADU level
    double noise_level;               // Image noise estimate
    QDateTime obs_time;               // Observation time
    bool wcs_valid;                   // WCS successfully loaded
    
    // Stellina-specific quality metrics
    double stellina_correction_magnitude; // From stacking JSON
    int stellina_stars_used;             // Stars used in Stellina registration
    
    WCSImageData() : quality_score(1.0), exposure_time(10.0), star_count(0),
                    background_level(0.0), noise_level(0.0), wcs_valid(false),
                    stellina_correction_magnitude(0.0), stellina_stars_used(0) {
        wcsini(1, 2, &wcs);  // Initialize WCS structure
    }
    
    ~WCSImageData() {
        wcsfree(&wcs);  // Clean up WCS memory
    }
};

class WCSAstrometricStacker : public QObject {
    Q_OBJECT

public:
    explicit WCSAstrometricStacker(QObject *parent = nullptr);
    ~WCSAstrometricStacker();
    
    // Main interface
    bool addImage(const QString &fits_file, const QString &solved_fits_file = "");
    bool addImageWithMetadata(const QString &fits_file, const struct StellinaImageData &stellina_data);
    void setStackingParameters(const StackingParams &params);
    void setProgressWidgets(QProgressBar *progress, QLabel *status);
    
    // Core stacking process
    bool computeOptimalWCS();
    bool stackImages();
    bool saveResult(const QString &output_path);
    
    // Quality analysis
    void analyzeImageQuality();
    void generateQualityReport();
    QString getQualityReport() const { return m_quality_report; }
    
    // Access results
    cv::Mat getStackedImage() const { return m_stacked_image; }
    cv::Mat getWeightMap() const { return m_weight_map; }
    cv::Mat getOverlapMap() const { return m_overlap_map; }
    struct wcsprm getOutputWCS() const { return m_output_wcs; }
    
    // Statistics
    int getImageCount() const { return m_images.size(); }
    double getTotalExposureTime() const;
    double getAverageQuality() const;

signals:
    void progressUpdated(int percentage);
    void statusUpdated(const QString &message);
    void imageProcessed(const QString &filename);
    void stackingComplete(bool success);
    void errorOccurred(const QString &error);

private:
    // Core processing functions
    bool loadWCSImage(const QString &fits_file, WCSImageData &img_data);
    bool extractImageStatistics(WCSImageData &img_data);
    bool computeImageQualityScore(WCSImageData &img_data);
    
    // WCS and coordinate functions
    bool findOptimalOutputWCS();
    bool reprojectImage(const WCSImageData &input, cv::Mat &output, cv::Mat &weight_output);
    bool worldToPixel(const struct wcsprm &wcs, double ra, double dec, double &x, double &y);
    bool pixelToWorld(const struct wcsprm &wcs, double x, double y, double &ra, double &dec);
    
    // Stacking algorithms with overlap correction
    cv::Mat performWeightedMean(const std::vector<cv::Mat> &images, 
                               const std::vector<cv::Mat> &weights,
                               const cv::Mat &overlap_count);
    cv::Mat performSigmaClippedMean(const std::vector<cv::Mat> &images,
                                   const std::vector<cv::Mat> &weights,
                                   const cv::Mat &overlap_count);
    cv::Mat performMedianStack(const std::vector<cv::Mat> &images,
                              const cv::Mat &overlap_count);
    
    // Overlap intensity correction - THE KEY FEATURE
    void computeOverlapCounts(const std::vector<cv::Mat> &weight_maps, cv::Mat &overlap_count);
    void applyOverlapIntensityCorrection(cv::Mat &stacked_image, const cv::Mat &overlap_count);
    void normalizeByOverlapDepth(cv::Mat &image, const cv::Mat &overlap_count);
    
    // Quality weighting
    void computeAdaptiveWeights(std::vector<cv::Mat> &weight_maps);
    void applyExposureNormalization(std::vector<cv::Mat> &images);
    void applyQualityWeighting(std::vector<cv::Mat> &weight_maps);
    
    // Rejection algorithms
    cv::Mat applySigmaClipping(const std::vector<cv::Mat> &images,
                              const std::vector<cv::Mat> &weights,
                              double sigma_low, double sigma_high);
    cv::Mat applyPercentileClipping(const std::vector<cv::Mat> &images,
                                   double percentile_low, double percentile_high);
    
    // Utility functions
    bool saveWCSFITS(const QString &filename, const cv::Mat &image, 
                    const struct wcsprm &wcs, const cv::Mat &weight_map = cv::Mat());
    void updateProgress(int percentage, const QString &message);
    QString formatTime(double seconds);
    
    // Member variables
    std::vector<std::unique_ptr<WCSImageData>> m_images;
    struct wcsprm m_output_wcs;       // Target WCS for output
    cv::Size m_output_size;           // Output image dimensions
    double m_output_pixel_scale;      // Arcseconds per pixel
    
    // Results
    cv::Mat m_stacked_image;          // Final stacked result
    cv::Mat m_weight_map;             // Combined weight map
    cv::Mat m_overlap_map;            // Number of overlapping images per pixel
    cv::Mat m_noise_map;              // Noise estimate map
    
    // Processing parameters
    StackingParams m_params;
    QString m_quality_report;
    
    // UI progress tracking
    QProgressBar *m_progress_bar;
    QLabel *m_status_label;
    
    // Statistics
    double m_total_processing_time;
    int m_pixels_processed;
    int m_pixels_rejected;
};

// Implementation of key overlap correction functions

WCSAstrometricStacker::WCSAstrometricStacker(QObject *parent)
    : QObject(parent)
    , m_progress_bar(nullptr)
    , m_status_label(nullptr)
    , m_total_processing_time(0.0)
    , m_pixels_processed(0)
    , m_pixels_rejected(0)
{
    wcsini(1, 2, &m_output_wcs);
}

WCSAstrometricStacker::~WCSAstrometricStacker() {
    wcsfree(&m_output_wcs);
}

bool WCSAstrometricStacker::addImageWithMetadata(const QString &fits_file, 
                                                const struct StellinaImageData &stellina_data) {
    auto img_data = std::make_unique<WCSImageData>();
    
    // Load FITS image and WCS
    if (!loadWCSImage(fits_file, *img_data)) {
        emit errorOccurred(QString("Failed to load WCS from: %1").arg(fits_file));
        return false;
    }
    
    // Add Stellina-specific quality metrics
    img_data->stellina_correction_magnitude = stellina_data.hasValidCoordinates ? 
        sqrt(pow(stellina_data.altitude, 2) + pow(stellina_data.azimuth, 2)) : 0.0;
    img_data->exposure_time = stellina_data.exposureSeconds;
    img_data->obs_time = QDateTime::fromString(stellina_data.dateObs, "yyyy-MM-ddThh:mm:ss");
    
    // Extract image statistics and compute quality score
    extractImageStatistics(*img_data);
    computeImageQualityScore(*img_data);
    
    m_images.push_back(std::move(img_data));
    
    emit imageProcessed(fits_file);
    return true;
}

bool WCSAstrometricStacker::computeOptimalWCS() {
    if (m_images.empty()) {
        emit errorOccurred("No images loaded for WCS computation");
        return false;
    }
    
    updateProgress(0, "Computing optimal output WCS...");
    
    // Find the coverage area of all images
    double min_ra = 1e10, max_ra = -1e10;
    double min_dec = 1e10, max_dec = -1e10;
    double total_pixel_scale = 0.0;
    int valid_wcs_count = 0;
    
    for (const auto &img : m_images) {
        if (!img->wcs_valid) continue;
        
        // Get corner coordinates of image
        std::vector<std::pair<double, double>> corners = {
            {0, 0}, 
            {img->image.cols-1, 0},
            {0, img->image.rows-1}, 
            {img->image.cols-1, img->image.rows-1}
        };
        
        for (const auto &corner : corners) {
            double ra, dec;
            if (pixelToWorld(img->wcs, corner.first, corner.second, ra, dec)) {
                min_ra = std::min(min_ra, ra);
                max_ra = std::max(max_ra, ra);
                min_dec = std::min(min_dec, dec);
                max_dec = std::max(max_dec, dec);
            }
        }
        
        // Calculate pixel scale from CD matrix
        double pixel_scale = sqrt(img->wcs.cd[0]*img->wcs.cd[0] + img->wcs.cd[1]*img->wcs.cd[1]) * 3600.0;
        total_pixel_scale += pixel_scale;
        valid_wcs_count++;
    }
    
    if (valid_wcs_count == 0) {
        emit errorOccurred("No valid WCS found in any image");
        return false;
    }
    
    // Set up optimal output WCS
    m_output_wcs = m_images[0]->wcs;  // Start with reference WCS
    
    // Center the output WCS on the coverage area
    m_output_wcs.crval[0] = (min_ra + max_ra) / 2.0;   // Center RA
    m_output_wcs.crval[1] = (min_dec + max_dec) / 2.0; // Center Dec
    
    // Use average pixel scale unless overridden
    if (m_params.output_pixel_scale <= 0.0) {
        m_output_pixel_scale = total_pixel_scale / valid_wcs_count;
    } else {
        m_output_pixel_scale = m_params.output_pixel_scale;
    }
    
    // Update CD matrix for desired pixel scale
    double scale_factor = m_output_pixel_scale / 3600.0; // Convert to degrees
    m_output_wcs.cd[0] = -scale_factor;  // RA decreases with increasing X
    m_output_wcs.cd[1] = 0.0;
    m_output_wcs.cd[2] = 0.0;
    m_output_wcs.cd[3] = scale_factor;   // Dec increases with increasing Y
    
    // Calculate output image size
    if (m_params.output_width > 0 && m_params.output_height > 0) {
        m_output_size = cv::Size(m_params.output_width, m_params.output_height);
    } else {
        // Calculate size to encompass all images
        double ra_span = max_ra - min_ra;
        double dec_span = max_dec - min_dec;
        
        // Add 10% padding
        ra_span *= 1.1;
        dec_span *= 1.1;
        
        int width = static_cast<int>((ra_span * cos(m_output_wcs.crval[1] * M_PI / 180.0)) / scale_factor);
        int height = static_cast<int>(dec_span / scale_factor);
        
        m_output_size = cv::Size(width, height);
    }
    
    // Set reference pixel to center
    m_output_wcs.crpix[0] = m_output_size.width / 2.0;
    m_output_wcs.crpix[1] = m_output_size.height / 2.0;
    
    // Finalize WCS
    int status = wcsset(&m_output_wcs);
    if (status != 0) {
        emit errorOccurred("Failed to initialize output WCS");
        return false;
    }
    
    updateProgress(100, QString("Output WCS computed: %1x%2 pixels, %.2f\"/pixel")
                   .arg(m_output_size.width).arg(m_output_size.height).arg(m_output_pixel_scale));
    
    return true;
}

bool WCSAstrometricStacker::stackImages() {
    if (m_images.empty()) {
        emit errorOccurred("No images to stack");
        return false;
    }
    
    if (!computeOptimalWCS()) {
        return false;
    }
    
    updateProgress(0, "Starting astrometric stacking...");
    
    // Reproject all images to common grid
    std::vector<cv::Mat> reprojected_images;
    std::vector<cv::Mat> weight_maps;
    
    for (size_t i = 0; i < m_images.size(); ++i) {
        updateProgress(static_cast<int>(i * 50 / m_images.size()), 
                      QString("Reprojecting image %1/%2").arg(i+1).arg(m_images.size()));
        
        cv::Mat reproj_img, weight_img;
        if (reprojectImage(*m_images[i], reproj_img, weight_img)) {
            reprojected_images.push_back(reproj_img);
            weight_maps.push_back(weight_img);
        } else {
            emit errorOccurred(QString("Failed to reproject: %1").arg(m_images[i]->filename));
        }
    }
    
    if (reprojected_images.empty()) {
        emit errorOccurred("No images successfully reprojected");
        return false;
    }
    
    updateProgress(50, "Computing overlap counts...");
    
    // CRITICAL: Compute overlap counts for intensity correction
    cv::Mat overlap_count;
    computeOverlapCounts(weight_maps, overlap_count);
    m_overlap_map = overlap_count.clone();
    
    updateProgress(60, "Applying quality weighting...");
    
    // Apply exposure normalization if requested
    if (m_params.normalize_exposure) {
        applyExposureNormalization(reprojected_images);
    }
    
    // Apply quality-based weighting
    applyQualityWeighting(weight_maps);
    
    updateProgress(70, "Performing stacking combination...");
    
    // Perform the actual stacking with overlap correction
    cv::Mat result;
    switch (m_params.combination) {
    case StackingParams::WEIGHTED_MEAN:
        result = performWeightedMean(reprojected_images, weight_maps, overlap_count);
        break;
    case StackingParams::SIGMA_CLIPPED_MEAN:
        result = performSigmaClippedMean(reprojected_images, weight_maps, overlap_count);
        break;
    case StackingParams::MEDIAN:
        result = performMedianStack(reprojected_images, overlap_count);
        break;
    case StackingParams::MEAN:
        // Simple mean with overlap correction
        result = performWeightedMean(reprojected_images, weight_maps, overlap_count);
        break;
    default:
        result = performWeightedMean(reprojected_images, weight_maps, overlap_count);
        break;
    }
    
    updateProgress(90, "Applying final intensity correction...");
    
    // CRITICAL: Apply overlap-based intensity correction
    applyOverlapIntensityCorrection(result, overlap_count);
    
    // Create combined weight map
    cv::Mat combined_weights = cv::Mat::zeros(m_output_size, CV_32F);
    for (const auto &weight : weight_maps) {
        combined_weights += weight;
    }
    
    // Store results
    m_stacked_image = result;
    m_weight_map = combined_weights;
    
    updateProgress(100, QString("Stacking complete! %1 images combined.")
                   .arg(reprojected_images.size()));
    
    emit stackingComplete(true);
    return true;
}

// KEY FUNCTION: Compute overlap counts for each pixel
void WCSAstrometricStacker::computeOverlapCounts(const std::vector<cv::Mat> &weight_maps, 
                                                cv::Mat &overlap_count) {
    overlap_count = cv::Mat::zeros(m_output_size, CV_32S);
    
    // Count how many images contribute to each pixel
    for (const auto &weight_map : weight_maps) {
        cv::Mat mask = weight_map > 0;  // Non-zero weights indicate valid pixels
        cv::Mat mask_int;
        mask.convertTo(mask_int, CV_32S);
        overlap_count += mask_int;
    }
    
    // Debug: log overlap statistics
    double min_overlap, max_overlap;
    cv::minMaxLoc(overlap_count, &min_overlap, &max_overlap);
    
    // Count pixels by overlap depth
    std::vector<int> overlap_histogram(static_cast<int>(max_overlap) + 1, 0);
    for (int y = 0; y < overlap_count.rows; ++y) {
        for (int x = 0; x < overlap_count.cols; ++x) {
            int count = overlap_count.at<int>(y, x);
            if (count > 0) {
                overlap_histogram[count]++;
            }
        }
    }
    
    QString overlap_info = "Overlap distribution: ";
    for (size_t i = 1; i < overlap_histogram.size(); ++i) {
        if (overlap_histogram[i] > 0) {
            overlap_info += QString("%1-deep: %2px ").arg(i).arg(overlap_histogram[i]);
        }
    }
    
    updateProgress(-1, overlap_info);  // -1 means don't update percentage
}

// KEY FUNCTION: Apply overlap-based intensity correction
void WCSAstrometricStacker::applyOverlapIntensityCorrection(cv::Mat &stacked_image, 
                                                           const cv::Mat &overlap_count) {
    // This is the crucial step that adjusts intensity based on overlap depth
    // Without this, areas with more overlapping images will appear brighter
    
    for (int y = 0; y < stacked_image.rows; ++y) {
        for (int x = 0; x < stacked_image.cols; ++x) {
            int count = overlap_count.at<int>(y, x);
            
            if (count > 1) {  // Only correct pixels with multiple contributions
                float &pixel = stacked_image.at<float>(y, x);
                
                // Apply correction factor based on overlap depth
                // The intensity should be normalized by the number of contributing images
                // This ensures uniform brightness across the field
                
                if (pixel > 0) {  // Only process valid pixels
                    // Simple approach: divide by overlap count
                    // This assumes each image contributes equally to the intensity
                    pixel = pixel / static_cast<float>(count);
                    
                    // Alternative: More sophisticated weighting could consider
                    // - Quality of each contributing image
                    // - Exposure time differences
                    // - Distance from image centers
                }
            }
        }
    }
}

// Weighted mean stacking with overlap correction
cv::Mat WCSAstrometricStacker::performWeightedMean(const std::vector<cv::Mat> &images, 
                                                   const std::vector<cv::Mat> &weights,
                                                   const cv::Mat &overlap_count) {
    cv::Mat result = cv::Mat::zeros(m_output_size, CV_32F);
    cv::Mat total_weight = cv::Mat::zeros(m_output_size, CV_32F);
    
    // Accumulate weighted pixel values
    for (size_t i = 0; i < images.size(); ++i) {
        cv::Mat weighted_image;
        cv::multiply(images[i], weights[i], weighted_image);
        
        result += weighted_image;
        total_weight += weights[i];
    }
    
    // Normalize by total weight (avoid division by zero)
    cv::Mat mask = total_weight > 1e-6;
    cv::divide(result, total_weight, result, 1.0, -1);  // Only divide where mask is true
    
    // The overlap correction will be applied separately in applyOverlapIntensityCorrection
    return result;
}

// Sigma-clipped mean with overlap correction
cv::Mat WCSAstrometricStacker::performSigmaClippedMean(const std::vector<cv::Mat> &images,
                                                      const std::vector<cv::Mat> &weights,
                                                      const cv::Mat &overlap_count) {
    cv::Mat result = cv::Mat::zeros(m_output_size, CV_32F);
    cv::Mat rejection_map = cv::Mat::zeros(m_output_size, CV_32S);
    
    // For each pixel, collect values and apply sigma clipping
    for (int y = 0; y < m_output_size.height; ++y) {
        for (int x = 0; x < m_output_size.width; ++x) {
            std::vector<float> pixel_values;
            std::vector<float> pixel_weights;
            
            // Collect values from all images at this pixel
            for (size_t i = 0; i < images.size(); ++i) {
                float weight = weights[i].at<float>(y, x);
                if (weight > 0) {
                    pixel_values.push_back(images[i].at<float>(y, x));
                    pixel_weights.push_back(weight);
                }
            }
            
            if (pixel_values.size() >= 3) {  // Need at least 3 values for sigma clipping
                // Calculate mean and standard deviation
                double sum = 0, sum_weights = 0;
                for (size_t i = 0; i < pixel_values.size(); ++i) {
                    sum += pixel_values[i] * pixel_weights[i];
                    sum_weights += pixel_weights[i];
                }
                double mean = sum / sum_weights;
                
                double variance = 0;
                for (size_t i = 0; i < pixel_values.size(); ++i) {
                    double diff = pixel_values[i] - mean;
                    variance += pixel_weights[i] * diff * diff;
                }
                double std_dev = sqrt(variance / sum_weights);
                
                // Apply sigma clipping
                double final_sum = 0, final_weights = 0;
                int rejected = 0;
                
                for (size_t i = 0; i < pixel_values.size(); ++i) {
                    double diff = fabs(pixel_values[i] - mean);
                    double sigma_level = diff / std_dev;
                    
                    bool reject_low = (pixel_values[i] < mean) && (sigma_level > m_params.sigma_low);
                    bool reject_high = (pixel_values[i] > mean) && (sigma_level > m_params.sigma_high);
                    
                    if (!reject_low && !reject_high) {
                        final_sum += pixel_values[i] * pixel_weights[i];
                        final_weights += pixel_weights[i];
                    } else {
                        rejected++;
                    }
                }
                
                if (final_weights > 0) {
                    result.at<float>(y, x) = static_cast<float>(final_sum / final_weights);
                }
                
                rejection_map.at<int>(y, x) = rejected;
                m_pixels_rejected += rejected;
            } else if (!pixel_values.empty()) {
                // Not enough for sigma clipping, use weighted mean
                double sum = 0, sum_weights = 0;
                for (size_t i = 0; i < pixel_values.size(); ++i) {
                    sum += pixel_values[i] * pixel_weights[i];
                    sum_weights += pixel_weights[i];
                }
                result.at<float>(y, x) = static_cast<float>(sum / sum_weights);
            }
            
            m_pixels_processed++;
        }
    }
    
    updateProgress(-1, QString("Sigma clipping: %1 pixels rejected").arg(m_pixels_rejected));
    
    return result;
}

// Apply quality-based weighting using Stellina metrics
void WCSAstrometricStacker::applyQualityWeighting(std::vector<cv::Mat> &weight_maps) {
    for (size_t i = 0; i < weight_maps.size(); ++i) {
        if (i < m_images.size()) {
            double quality = m_images[i]->quality_score;
            
            // Apply quality scaling to weight map
            weight_maps[i] *= static_cast<float>(quality);
            
            // Additional weighting based on Stellina star count
            if (m_images[i]->stellina_stars_used > 0) {
                // More stars = better registration = higher weight
                double star_weight = std::min(1.0, m_images[i]->stellina_stars_used / 100.0);
                weight_maps[i] *= static_cast<float>(star_weight);
            }
        }
    }
}

#endif // WCS_ASTROMETRIC_STACKER_H
