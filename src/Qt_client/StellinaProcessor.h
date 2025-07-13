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

#include "SirilClient.h"

// cfitsio is required for Siril, so we can assume it's available
#include <fitsio.h>

// Processing modes
enum ProcessingMode {
    MODE_BASIC_PLATESOLVE = 0,
    MODE_DARK_CALIBRATION = 1,
    MODE_ASTROMETRIC_STACKING = 2,
    MODE_FULL_PIPELINE = 3
};

// Dark frame information
struct DarkFrame {
    QString filepath;
    int exposure;     // exposure time in seconds
    int temperature;  // sensor temperature in degrees C
    QString binning;  // binning mode (e.g., "1x1", "2x2")
};

// Stacking parameters
struct StackingParams {
    QString method;           // "sum", "median", "mean", "sigma_clipping"
    QString rejection;        // "none", "sigma", "linear", "percentile"
    double rejectionLow;      // lower rejection threshold
    double rejectionHigh;     // upper rejection threshold
    bool normalizeImages;     // normalize before stacking
    bool applyDrizzle;        // apply drizzle enhancement
    double drizzleScale;      // drizzle scale factor
    QString outputFormat;     // "fits", "tiff", "png"
};

class StellinaProcessor : public QMainWindow {
    Q_OBJECT

public:
    explicit StellinaProcessor(QWidget *parent = nullptr);
    ~StellinaProcessor();

private slots:
    // UI slots
    void onSelectSourceDirectory();
    void onSelectDarkDirectory();
    void onSelectCalibratedDirectory();
    void onSelectPlateSolvedDirectory();
    void onSelectStackedDirectory();
    void onStartProcessing();
    void onStopProcessing();
    void onTestConnection();
    void onClearLog();
    void onProcessingModeChanged();
    void onRefreshDarkFrames();
    void onStackingParametersChanged();
    void onPreviewStacking();
    
    // Siril connection slots
    void onSirilConnected();
    void onSirilDisconnected();
    void onSirilCommandExecuted(const QString &command, bool success);
    void onSirilError(const QString &error);
    
    // Processing slots
    void onProcessingTimer();

private:
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
    
    // Dark calibration functions
    void scanDarkFrames();
    bool findMatchingDarkFrame(const QString &lightFrame, DarkFrame &darkFrame); // Deprecated
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
    double calculateJD(int year, int month, int day, int hour, int minute, int second);
    double calculateLST(double JD, double longitude);
    void altAzToRaDec(double alt, double az, double lat, double lst, double &ra, double &dec);
    bool runSolveField(const QString &fitsPath, const QString &outputPath, double ra, double dec);
    bool checkSolveFieldInstalled();
    bool createBinnedImageForPlatesolving(const QString &inputPath, const QString &binnedPath);
    bool performCFABinning(const std::vector<float> &inputPixels, std::vector<float> &binnedPixels, 
			   long width, long height, long &binnedWidth, long &binnedHeight);

    // UI components - Main tabs
    QTabWidget *m_tabWidget;
    QWidget *m_basicTab;
    QWidget *m_darkTab;
    QWidget *m_stackingTab;
    QWidget *m_logTab;
    
    // Connection group
    QGroupBox *m_connectionGroup;
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
    
    // Basic options group (original options)
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
    
    // Stacking options group
    QGroupBox *m_stackingOptionsGroup;
    QComboBox *m_stackingMethodCombo;
    QComboBox *m_rejectionMethodCombo;
    QDoubleSpinBox *m_rejectionLowSpin;
    QDoubleSpinBox *m_rejectionHighSpin;
    QCheckBox *m_normalizeCheck;
    QCheckBox *m_drizzleCheck;
    QDoubleSpinBox *m_drizzleScaleSpin;
    QComboBox *m_outputFormatCombo;
    QPushButton *m_previewStackButton;
    
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
    
    // Status bar
    QLabel *m_statusLabel;
    QLabel *m_memoryUsageLabel;
    
    // Core components
    SirilClient *m_sirilClient;
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
    StackingParams m_stackingParams;
    
    // File tracking for pipeline
    QStringList m_darkCalibratedFiles;
    QStringList m_plateSolvedFiles;
    QStringList m_registeredFiles;
    QString m_finalStackedImage;
    QString m_sequenceName;
};

#endif // STELLINAPROCESSOR_H
