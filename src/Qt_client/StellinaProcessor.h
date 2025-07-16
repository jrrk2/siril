#ifndef STELLINAPROCESSOR_H
#define STELLINAPROCESSOR_H

#include <QMainWindow>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QTextEdit>
#include <QProgressBar>
#include <QLabel>
#include <QGroupBox>
#include <QCheckBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QTimer>
#include <QStringList>
#include <QStandardPaths>
#include <QJsonObject>
#include <QJsonDocument>
#include <QComboBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QSplitter>
#include <QElapsedTimer>
#include <QProcessEnvironment>
#include <fitsio.h>

//============================================================================
// 2. ADD THESE FORWARD DECLARATIONS (if not already present)
//============================================================================

class WCSAstrometricStacker;

// Processing modes
enum ProcessingMode {
    MODE_BASIC_PLATESOLVE = 0,
    MODE_DARK_CALIBRATION = 1,
    MODE_ASTROMETRIC_STACKING = 2,
    MODE_FULL_PIPELINE = 3
};

// Stacking parameters - unified structure for both WCS and traditional stacking
struct StackingParams {
    enum CombinationMethod {
        MEAN,
        MEDIAN, 
        WEIGHTED_MEAN,
        SIGMA_CLIPPED_MEAN,
        MINIMUM,
        MAXIMUM
    };
    
    enum RejectionMethod {
        NO_REJECTION,
        SIGMA_CLIPPING,
        PERCENTILE_CLIPPING,
        LINEAR_FIT_CLIPPING
    };
    
    CombinationMethod combination = WEIGHTED_MEAN;
    RejectionMethod rejection = SIGMA_CLIPPING;
    double sigma_low = 3.0;           // Low sigma clipping threshold
    double sigma_high = 3.0;          // High sigma clipping threshold
    double percentile_low = 5.0;      // Low percentile (%)
    double percentile_high = 95.0;    // High percentile (%)
    bool normalize_exposure = true;   // Normalize by exposure time
    bool apply_flat_correction = false; // Apply flat field correction
    bool apply_brightness_normalization = false; // Apply brightness normalisation
    double output_pixel_scale = 0.0;  // Override pixel scale (0 = auto)
    int output_width = 0;             // Override width (0 = auto)
    int output_height = 0;            // Override height (0 = auto)
    bool create_weight_map = true;    // Generate output weight map
    bool save_intermediate = false;   // Save reprojected images
    QString output_format = "fits";   // Output format
};

// Dark frame information
struct DarkFrame {
    QString filepath;
    int exposure;     // exposure time in seconds
    int temperature;  // sensor temperature in degrees C
    QString binning;  // binning mode (e.g., "1x1", "2x2")
    QString bayerPattern;      // NEW: bayer pattern for this dark frame
    
    DarkFrame() : exposure(0), temperature(0), bayerPattern("RGGB") {}

};

// Update the StellinaImageData structure
struct StellinaImageData {
    QString originalFitsPath;     // Original raw FITS file path
    QString originalJsonPath;     // Original JSON metadata file path  
    QString currentFitsPath;      // Current FITS file path (updated through pipeline)
    QJsonObject metadata;         // Complete JSON metadata
    double altitude;              // Stellina altitude (degrees)
    double azimuth;               // Stellina azimuth (degrees)
    QString dateObs;              // DATE-OBS from FITS header
    bool hasValidCoordinates;     // Whether Alt/Az coordinates are valid
    int exposureSeconds;          // Exposure time in seconds
    int temperatureKelvin;        // Sensor temperature in Kelvin
    QString binning;              // Binning mode (e.g., "1x1", "2x2")
    
    // NEW: Pre-calculated RA/DEC coordinates from coordinate conversion
    double calculatedRA;          // Calculated RA from Alt/Az conversion (degrees)
    double calculatedDec;         // Calculated Dec from Alt/Az conversion (degrees)
    bool hasCalculatedCoords;     // Whether calculated coordinates are available
    
    // NEW: Reversed image support
    bool isReversedImage;         // True if this is a reversed stellina image (img-0001r.fits)
    QString bayerPattern;         // Detected bayer pattern (e.g., "RGGB", "GRBG", "GBRG", "BGGR")
    QString baseName;             // Base name without 'r' suffix (e.g., "img-0001")
    
    StellinaImageData() : altitude(0), azimuth(0), hasValidCoordinates(false), 
                         exposureSeconds(0), temperatureKelvin(284), binning("1x1"),
			  calculatedRA(0), calculatedDec(0), hasCalculatedCoords(false), isReversedImage(false), bayerPattern("RGGB")  {}
    
    // Convenience function to check if pre-calculated coordinates are available
    bool hasPreCalculatedCoords() const {
        return hasCalculatedCoords && calculatedRA != 0.0 && calculatedDec != 0.0;
    }
};

struct StackingCorrectionData {
    QString imageFilename;
    int imageNumber;
    qint64 acqTime;
    double minutesFromStart;
    
    // Mount position
    double stellinaAlt, stellinaAz;
    
    // Stacking corrections (pixels)
    double correctionX, correctionY, correctionRot;
    
    // Registration quality metrics
    int starsUsed;
    QString statusMessage;
    double distanceToCenter;
    
    bool isValid;
    
    StackingCorrectionData() : imageNumber(0), minutesFromStart(0.0), 
                              stellinaAlt(0), stellinaAz(0),
                              correctionX(0), correctionY(0), correctionRot(0),
                              starsUsed(0), distanceToCenter(0), isValid(false) {}
};

struct ProcessedImageData {
    QString filename;
    int imageNumber;
    double stellinaAlt, stellinaAz;
    double predictedRA, predictedDec;
    double solvedRA, solvedDec;
    QString dateObs;
    QDateTime obsTime;
    double minutesFromStart;
    bool isValid;
};


struct CoordinateTestCase {
    QString name;
    double lat;           // Observer latitude (degrees)
    double lon;           // Observer longitude (degrees) 
    QString date;         // UTC date/time string
    double julianDay;     // Julian Day
    double raOfDate;      // RA of date (degrees)
    double decOfDate;     // Dec of date (degrees)
    double expectedAz;    // Expected azimuth (degrees)
    double expectedAlt;   // Expected altitude (degrees)
    double siderealTime;  // Local sidereal time (hours)
    double hourAngle;     // Hour angle (hours)
    QString description;  // Test description
};

class StellinaProcessor : public QMainWindow {
    Q_OBJECT

public:
    explicit StellinaProcessor(QWidget *parent = nullptr);
    ~StellinaProcessor();
    void diagnoseSiderealTimeIssues();
    void testFixedCoordinateConversion();
    void diagnoseStellinaProcessing();
    void verifyCorrectProcessing(); 
    void analyzeRealCoordinateErrors();
    bool processImagePlatesolving_Fixed(const QString &calibratedFitsPath);
    void diagnoseTrackingIssue();
    void debugCoordinateSystem();				
    void testStellinaAzimuthConvention();
    double calculateLST_HighPrecision(double JD, double longitude);
    void diagnoseLSTAccuracy();
    void testTimeDriftFix();
    void analyzeRealStellinaIssue();
    void testRealisticAccuracy();
    void verifyPlatesolvingHints();
    void dumpCoordinateData();
    void dumpCoordinateDataToCSV();
    void analyzeCoordinateDrift();
    void calibrateFromProcessedFiles();
    bool readStellinaDataFromSolvedFits(const QString &fitsPath, ProcessedImageData &data);
    bool readSolveFieldResults(const QString &fitsPath, ProcessedImageData &data);
    void analyzeAndCalibrateFromData(const QList<ProcessedImageData> &imageData, const double &sessionStart);
    void testSystematicOffsetCorrection();
    void verifySystematicOffsetsInUse();
    void plotMountErrors();
    bool readAcqTimeFromFits(const QString &fitsPath, qint64 &acqTime, qint64 &refAcqTime, QDateTime &refUTCTime);
    QDateTime reconstructUTCFromAcqTime(qint64 acqTime, qint64 refAcqTime, const QDateTime &refUTCTime);
    bool convertAltAzToRaDecFromCalibratedFits(const QString &calibratedFitsPath, double &ra, double &dec);
    void compareTimingAccuracy() ;
    void validateAcqTimePreservation() ;
				       
private slots:
    // Debug slot functions
    void onTestConversion();
    void onLoadImageData();
    void onTestBatch();
    void onRunPresetTest();
    void onTestSkyRegion();
  
    // UI slots
    void onSelectSourceDirectory();
    void onSelectDarkDirectory();
    void onSelectCalibratedDirectory();
    void onSelectPlateSolvedDirectory();
    void onSelectStackedDirectory();
    void onStartProcessing();
    void onStopProcessing();
    void onClearLog();
    void onProcessingModeChanged();
    void onRefreshDarkFrames();
//    void onStackingParametersChanged();
//    void onPreviewStacking();
    
    // Siril connection slots
    void onSirilConnected();
    void onSirilDisconnected();
    void onSirilCommandExecuted(const QString &command, bool success);
    void onSirilError(const QString &error);
    
    // Processing slots
    void onProcessingTimer();
    void calibrateFromStackingJSON();
    bool parseStackingJSON(const QString &jsonPath, StackingCorrectionData &data);
    void analyzeStackingCorrections(const QList<StackingCorrectionData> &stackingData, const double &sessionStart);
    void analyzeMosaicCorrections(const QList<StackingCorrectionData> &stackingData,
                                  const QMap<QString, int> &patternCount);
    void onWCSParametersChanged();
    void onStartWCSStacking();
    void onWCSStackingComplete(bool success);
    void onWCSProgressUpdated(int percentage);
    void onWCSStatusUpdated(const QString &message);
    void onSaveWCSResult();
  
private:
    // Session timing state
    static QDateTime s_sessionReferenceTime;
    static qint64 s_sessionReferenceAcqTime;
    static bool s_sessionTimingInitialized;
    // Method declarations only
    void resetSessionTiming();
    bool initializeSessionTiming(const QString &sourceDirectory);
    QDateTime convertAcqTimeToUTC(qint64 acqTime);
    bool convertAltAzToRaDecWithPreciseTiming(double alt, double az, const QString &jsonPath, double &ra, double &dec);
    void validateTimingOffset();
    void compareDriftWithDynamicOffset();
    void updateProcessingToDynamicOffset();
  
    // Mount tilt correction parameters
    struct MountTiltParams {
        double northTilt;           // Static north tilt θ_N in degrees (deprecated)
        double eastTilt;            // Static east tilt θ_E in degrees (deprecated)
        double driftRA;             // RA drift rate in degrees per hour
        double driftDec;            // Dec drift rate in degrees per hour
        double systematicRAOffset;  // Systematic RA offset correction (deprecated)
        double systematicDecOffset; // Systematic Dec offset correction (deprecated)
        double initialRAOffset;     // RA error at session start (t=0)
        double initialDecOffset;    // Dec error at session start (t=0)
        double sessionStart;        // Start time of observing session
        bool enableCorrection;      // Whether to apply correction
        bool enableDriftCorrection; // Whether to apply time-dependent drift correction
        
        MountTiltParams() : northTilt(0.0), eastTilt(0.0), 
                           driftRA(0.0), driftDec(0.0),
                           systematicRAOffset(0.0), systematicDecOffset(0.0),
                           initialRAOffset(0.0), initialDecOffset(0.0),
                           enableCorrection(false), enableDriftCorrection(false) {}
    };
    MountTiltParams m_mountTilt;

    void setupUI();
    void setupBasicTab();
    void setupDarkTab();
    void setupStackingTab();
    void setupLogTab();
    void setupMenu();
    void connectSignals();
    void updateUI();
    void updateConnectionStatus();
    void updateProcessingStatus();
    void logMessage(const QString &message, const QString &color = "black");
    void loadSettings();
    void saveSettings();
    
    // Processing functions
    void startStellinaProcessing();
    void processNextImage();
    bool processImageDarkCalibration(const QString &lightFrame);
    bool processImagePlatesolving(const QString &fitsPath);
    void finishProcessing();
    bool findStellinaImages();
    QJsonObject loadStellinaJson(const QString &jsonPath);
    bool extractCoordinates(const QJsonObject &json, double &alt, double &az);
    bool convertAltAzToRaDec(double alt, double az, const QString &dateObs, double &ra, double &dec);
    bool checkStellinaQuality(const QJsonObject &json);
    QString getStageDescription() const;
    void runCoordinateTestSuite();
    void testSingleCoordinate(const CoordinateTestCase &testCase);
    QList<CoordinateTestCase> getCoordinateTestCases();
    void runRandomTestSubset(int numTests = 20);
    void runAccuracyAnalysis();
    // NEW: Reversed image support methods
    bool isReversedStellinaImage(const QString &filename);
    QString detectBayerPattern(const QString &fitsPath);
    QString getBaseName(const QString &filename);
    bool needsDarkRotation(const QString &lightBayerPattern, const QString &darkBayerPattern);
    bool rotateDarkFrame(const QString &inputDarkPath, const QString &outputDarkPath, 
                        const QString &fromPattern, const QString &toPattern);
    QPoint getBayerPatternOffset(const QString &pattern);
    QString getMatchingDarkKey(const StellinaImageData &imageData);
  QString getRotationMethodDescription(const QString &fromPattern, const QString &toPattern);
  bool performBayerPatternRotation(const std::vector<float> &inputData, 
                                                   std::vector<float> &outputData,
                                                   long width, long height,
                                                   const QString &fromPattern, 
						      const QString &toPattern) ;
  int getBayerRotationType(const QString &fromPattern, const QString &toPattern) ;
  bool rotateBayerImage180(const std::vector<float> &inputData, 
                                           std::vector<float> &outputData,
					      long width, long height);
  bool rotateBayerImage90CW(const std::vector<float> &inputData, 
                                            std::vector<float> &outputData,
					       long width, long height) ;
  bool rotateBayerImage270CW(const std::vector<float> &inputData, 
                                             std::vector<float> &outputData,
						long width, long height);
  bool flipBayerImageHorizontal(const std::vector<float> &inputData, 
                                                std::vector<float> &outputData,
						   long width, long height);
  bool flipBayerImageVertical(const std::vector<float> &inputData, 
                                              std::vector<float> &outputData,
						 long width, long height);
  QStringList findAllMatchingDarkFrames(int targetExposure, int targetTemperature, 
                                                        const QString &targetBinning, 
							   const QString &targetBayerPattern);
  void loadDarkFrames() ;
  
    // Mount tilt correction functions
    void applyMountTiltCorrection(double &alt, double &az, double inputAlt, double inputAz);
    void calibrateMountTilt();
    void testMountTiltCorrection();
    bool loadMountTiltFromSettings();
    void saveMountTiltToSettings();
    void updateTiltUI();
    void performDriftAnalysis(const QList<StackingCorrectionData> &stackingData,
                             const QList<double> &timePoints,
                             const QList<double> &xCorrections,
                             const QList<double> &yCorrections,
                             const double &sessionStart);
    
    // Dark calibration functions
    void scanDarkFrames();
    bool findMatchingDarkFrame(const QString &lightFrame, DarkFrame &darkFrame);
    QStringList findAllMatchingDarkFrames(int targetExposure, int targetTemperature, const QString &targetBinning);
    bool createMasterDark(const QStringList &darkFrames, const QString &outputPath);
    bool createMasterDarkDirect(const QStringList &darkFrames, const QString &outputPath);
    bool applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame);
    int extractExposureTime(const QString &fitsFile);
    int extractTemperature(const QString &fitsFile);
    QString extractBinning(const QString &fitsFile);
    
    // Astrometric stacking functions
    bool performAstrometricStacking();
    bool registerImages(const QStringList &imageList, const QString &referenceImage);
    bool stackRegisteredImages(const QStringList &registeredImages, const QString &outputStack);
    bool createSequence(const QStringList &imageList, const QString &sequenceName);
    bool performGlobalRegistration(const QString &sequenceName);
    bool performStacking(const QString &sequenceName, const StackingParams &params);
    
    // Utility functions
    QString formatProcessingTime(qint64 milliseconds);
    bool validateProcessingInputs();
    void saveProcessingReport();
    QString getOutputDirectoryForCurrentStage() const;
    bool applyMasterDarkDirect(const QString &lightFrame, const QString &masterDark, const QString &outputFrame);
    bool parseObserverLocation(const QString &location, double &lat, double &lon, double &elevation);
    QString extractDateObs(const QString &fitsFile);
    void testLibnovaConversion();
    void testSingleConversion(const QString &testName,
                             double alt, double az, 
                             const QString &dateObs,
                             double expectedRA = 0.0, double expectedDec = 0.0,
                             double currentRA = 0.0, double currentDec = 0.0,
                             double testLat = 0.0, double testLon = 0.0);
    bool runSolveField(const QString &fitsPath, const QString &outputPath, double ra, double dec);
    bool checkSolveFieldInstalled();
    bool createBinnedImageForPlatesolving(const QString &inputPath, const QString &binnedPath);
    bool performCFABinning(const std::vector<float> &inputPixels, std::vector<float> &binnedPixels, 
                          long width, long height, long &binnedWidth, long &binnedHeight);
    QProcessEnvironment createSolveFieldEnvironment();
    QStringList getAstrometryPaths();
    void initializeWCSStacker();
    void setupWCSStackingUI();
    void addWCSMenuItems();
    void updateWCSUI();
    void loadWCSSettings();
    void saveWCSSettings();
    
    // Enhanced processing methods
    bool performAstrometricStackingEnhanced();
    // NEW: Debug tab and components
    QWidget *m_debugTab;
    void setupDebugTab();
    void debugLog(const QString &message) ;
    // Debug UI components
    QGroupBox *m_coordDebugGroup;
    QDoubleSpinBox *m_debugAltSpin;
    QDoubleSpinBox *m_debugAzSpin;
    QLineEdit *m_debugTimeEdit;
    QDoubleSpinBox *m_debugLatSpin;
    QDoubleSpinBox *m_debugLonSpin;
    QPushButton *m_testConversionButton;
    QPushButton *m_loadImageDataButton;
    QPushButton *m_testBatchButton;
    QTextEdit *m_debugResultsEdit;
    
    QGroupBox *m_presetTestsGroup;
    QComboBox *m_presetTestCombo;
    QPushButton *m_runPresetButton;
    
    QGroupBox *m_skyRegionGroup;
    QComboBox *m_skyRegionCombo;
    QPushButton *m_testSkyRegionButton;

    // UI components - Main tabs
    QTabWidget *m_tabWidget;
    QWidget *m_basicTab;
    QWidget *m_darkTab;
    QWidget *m_stackingTab;
    QWidget *m_logTab;
    
    // Connection group
//    QGroupBox *m_connectionGroup;
    QPushButton *m_testConnectionButton;
    QLabel *m_connectionStatus;
    
    // Processing mode group
    QGroupBox *m_modeGroup;
    QComboBox *m_processingModeCombo;
    QLabel *m_modeDescription;
    
    // Input group
    QGroupBox *m_inputGroup;
    QLineEdit *m_sourceDirectoryEdit;
    QPushButton *m_selectSourceButton;
    QLineEdit *m_darkDirectoryEdit;
    QPushButton *m_selectDarkButton;
    QLineEdit *m_calibratedDirectoryEdit;
    QPushButton *m_selectCalibratedButton;
    QLineEdit *m_plateSolvedDirectoryEdit;
    QPushButton *m_selectPlateSolvedButton;
    QLineEdit *m_stackedDirectoryEdit;
    QPushButton *m_selectStackedButton;
    QLabel *m_darkFramesCount;
    QPushButton *m_refreshDarkButton;
    
    // Basic options group
    QGroupBox *m_basicOptionsGroup;
    QCheckBox *m_qualityFilterCheck;
    QCheckBox *m_debugModeCheck;
    QDoubleSpinBox *m_focalLengthSpin;
    QDoubleSpinBox *m_pixelSizeSpin;
    QLineEdit *m_observerLocationEdit;
    
    // Dark calibration options
    QGroupBox *m_darkOptionsGroup;
    QCheckBox *m_autoMatchDarksCheck;
    QSpinBox *m_temperatureToleranceSpin;
    QSpinBox *m_exposureToleranceSpin;
    QTableWidget *m_darkFramesTable;
    // Processing group
    QGroupBox *m_processingGroup;
    QPushButton *m_startButton;
    QPushButton *m_stopButton;
    QProgressBar *m_progressBar;
    QLabel *m_progressLabel;
    QLabel *m_timeEstimateLabel;
    QLabel *m_currentTaskLabel;
    
    // Advanced processing info
    QGroupBox *m_advancedInfoGroup;
    QLabel *m_registrationStatusLabel;
    QLabel *m_stackingStatusLabel;
    QLabel *m_darkCalibrationStatusLabel;
    QProgressBar *m_subTaskProgressBar;
    
    // Log group
    QGroupBox *m_logGroup;
    QTextEdit *m_logTextEdit;
    QPushButton *m_clearLogButton;
    QPushButton *m_saveLogButton;

    // Mount tilt correction UI
    QGroupBox *m_mountTiltGroup;
    QCheckBox *m_enableTiltCorrectionCheck;
    QDoubleSpinBox *m_northTiltSpin;
    QDoubleSpinBox *m_eastTiltSpin;
    QPushButton *m_calibrateTiltButton;
    QPushButton *m_testTiltButton;
    QLabel *m_tiltStatusLabel;
    QCheckBox *m_enableDriftCorrectionCheck;
    QDoubleSpinBox *m_driftRASpin;
    QDoubleSpinBox *m_driftDecSpin;
    QLabel *m_driftStatusLabel;
    
    // WCS Astrometric Stacker
    WCSAstrometricStacker *m_wcsStacker;
    StackingParams m_wcsStackingParams;  // Changed from WCSStackingParams to StackingParams
    
    // WCS Stacking UI elements
    QGroupBox *m_wcsStackingGroup;
    QComboBox *m_wcsCombinationMethodCombo;
    QComboBox *m_wcsRejectionMethodCombo;
    QDoubleSpinBox *m_wcsSigmaLowSpin;
    QDoubleSpinBox *m_wcsSigmaHighSpin;
    QCheckBox *m_wcsNormalizeExposureCheck;
    QCheckBox *m_wcsCreateWeightMapCheck;
    QSpinBox *m_wcsOutputWidthSpin;
    QSpinBox *m_wcsOutputHeightSpin;
    QDoubleSpinBox *m_wcsOutputPixelScaleSpin;
    QPushButton *m_startWCSStackingButton;
    QPushButton *m_saveWCSResultButton;
  
    // Status bar
    QLabel *m_statusLabel;
    QLabel *m_memoryUsageLabel;
    
    // Core components
    QTimer *m_processingTimer;
    
    // Processing state
    bool m_processing;
    ProcessingMode m_processingMode;
    QStringList m_imagesToProcess;
    QList<DarkFrame> m_darkFrames;
    int m_currentImageIndex;
    int m_processedCount;
    int m_errorCount;
    int m_skippedCount;
    int m_darkCalibratedCount;
    int m_registeredCount;
    qint64 m_processingStartTime;
    
    // Processing stages for full pipeline
    enum ProcessingStage {
        STAGE_DARK_CALIBRATION,
        STAGE_PLATE_SOLVING,
        STAGE_REGISTRATION,
        STAGE_STACKING,
        STAGE_COMPLETE
    };
    ProcessingStage m_currentStage;
    
    // Settings - Updated to use multiple directories
    QString m_sourceDirectory;        // Raw light frames
    QString m_darkDirectory;         // Dark frames
    QString m_calibratedDirectory;   // Dark-calibrated light frames
    QString m_plateSolvedDirectory;  // Plate-solved light frames
    QString m_stackedDirectory;      // Final stacked images
    bool m_qualityFilter;
    bool m_debugMode;
    double m_focalLength;
    double m_pixelSize;
    QString m_observerLocation;
    
    // Dark calibration settings
    bool m_autoMatchDarks;
    int m_temperatureTolerance;
    int m_exposureTolerance;
    
    // Stacking settings
    StackingParams m_stackingParams;  // Legacy params for backward compatibility
    
    // File tracking for pipeline
    QStringList m_darkCalibratedFiles;
    QStringList m_plateSolvedFiles;
    QStringList m_registeredFiles;
    QString m_finalStackedImage;
    QString m_sequenceName;

    // Add to private member variables
    QList<StellinaImageData> m_stellinaImageData;  // Tracks metadata through pipeline

    // Stage setup functions for clean pipeline flow
    bool setupDarkCalibrationStage();
    bool setupPlatesolvingStage();
    bool setupStackingStage();
    void handlePipelineStageTransition();
    
    // FITS metadata functions - enhanced version that includes coordinate conversion
    bool writeStellinaMetadataWithCoordinates(const QString &fitsPath, const StellinaImageData &imageData);
    bool readStellinaMetadataFromFits(const QString &fitsPath, StellinaImageData &imageData);
    bool updateProcessingStage(const QString &fitsPath, const QString &stage);
    bool cleanExistingStellinaKeywords(const QString &fitsPath);
    StellinaImageData* findImageDataByPath(const QString &path);
};

#endif // STELLINAPROCESSOR_H
