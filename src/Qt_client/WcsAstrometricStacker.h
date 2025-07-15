#ifndef WCS_ASTROMETRIC_STACKER_H
#define WCS_ASTROMETRIC_STACKER_H

#include <QObject>
#include <QString>
#include <QList>
#include <QDateTime>
#include <QProgressBar>
#include <QLabel>
#include <QTimer>
#include <QStringList>
#include <opencv2/opencv.hpp>
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <fitsio.h>
#include <vector>
#include <memory>
#include "StellinaProcessor.h"

// Forward declarations
struct StellinaImageData;
struct StackingParams;

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
        // wcsini(1, 2, &wcs);  // Initialize WCS structure
						memset(&wcs, 0, sizeof(wcs));
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
    QString getQualityReport() const;
    
    // Access results
    cv::Mat getStackedImage() const { return m_stacked_image; }
    cv::Mat getWeightMap() const { return m_weight_map; }
    cv::Mat getOverlapMap() const { return m_overlap_map; }
    struct wcsprm getOutputWCS() const { return m_output_wcs; }
    
    // Statistics
    int getImageCount() const { return m_images.size(); }
    double getTotalExposureTime() const;
    double getAverageQuality() const;
    
    // Processing control
    void startStacking();
    void cancelStacking();

signals:
    void progressUpdated(int percentage);
    void statusUpdated(const QString &message);
    void imageProcessed(const QString &filename);
    void stackingComplete(bool success);
    void errorOccurred(const QString &error);
    void qualityAnalysisComplete();

private slots:
    void processNextImage();

private:
    // Core processing functions
    bool loadWCSFromFITS(const QString &fits_file, WCSImageData &img_data);
    bool extractImageStatistics(WCSImageData &img_data);
    bool computeImageQualityScore(WCSImageData &img_data);
    bool addPlatesolveDFITSFile(const QString &solved_fits_file);
    bool addImageFromStellinaData(const QString &fits_file, const StellinaImageData &stellina_data);
    
    // Utility functions
    void updateProgress(int percentage, const QString &message);
    void logProcessing(const QString &message);
    void finishStacking();
    
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
    QStringList m_processing_log;
    
    // UI progress tracking
    QProgressBar *m_progress_bar;
    QLabel *m_status_label;
    
    // Processing state
    bool m_stacking_active;
    int m_current_image_index;
    QTimer *m_processing_timer;
    
    // Statistics
    double m_total_processing_time;
    int m_pixels_processed;
    int m_pixels_rejected;
};

#endif // WCS_ASTROMETRIC_STACKER_H
