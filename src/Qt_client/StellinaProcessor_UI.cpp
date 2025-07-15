#include "StellinaProcessor.h"
#include <QStatusBar>
#include <QMenuBar>
#include <QMenu>
#include <QHeaderView>

void StellinaProcessor::setupUI() {
    // Create central widget with tabs
    m_tabWidget = new QTabWidget;
    setCentralWidget(m_tabWidget);
    
    // Create tab widgets
    m_basicTab = new QWidget;
    m_darkTab = new QWidget;
    m_stackingTab = new QWidget;
    m_logTab = new QWidget;
    
    m_tabWidget->addTab(m_basicTab, "Basic Processing");
    m_tabWidget->addTab(m_darkTab, "Dark Calibration");
    m_tabWidget->addTab(m_stackingTab, "Astrometric Stacking");
    m_tabWidget->addTab(m_logTab, "Processing Log");
    
    setupBasicTab();
    setupDarkTab();
    setupStackingTab();
    setupLogTab();
    
    // Status bar with additional info
    m_statusLabel = new QLabel("Ready");
    m_memoryUsageLabel = new QLabel("Memory: 0 MB");
    statusBar()->addWidget(m_statusLabel);
    statusBar()->addPermanentWidget(m_memoryUsageLabel);
}

void StellinaProcessor::setupBasicTab() {
    QVBoxLayout *layout = new QVBoxLayout(m_basicTab);
    
    // Processing mode group
    m_modeGroup = new QGroupBox("Processing Mode");
    QVBoxLayout *modeLayout = new QVBoxLayout(m_modeGroup);
    
    m_processingModeCombo = new QComboBox;
    m_processingModeCombo->addItem("Basic Plate Solving", MODE_BASIC_PLATESOLVE);
    m_processingModeCombo->addItem("Dark Calibration Only", MODE_DARK_CALIBRATION);
    m_processingModeCombo->addItem("Astrometric Stacking", MODE_ASTROMETRIC_STACKING);
    m_processingModeCombo->addItem("Full Pipeline", MODE_FULL_PIPELINE);
    
    m_modeDescription = new QLabel("Basic plate solving with coordinate annotation");
    m_modeDescription->setWordWrap(true);
    m_modeDescription->setStyleSheet("color: gray; font-style: italic;");
    
    modeLayout->addWidget(m_processingModeCombo);
    modeLayout->addWidget(m_modeDescription);
    
    // Input group - Processing Directories
    m_inputGroup = new QGroupBox("Processing Directories");
    QGridLayout *inputLayout = new QGridLayout(m_inputGroup);
    
    // Raw light frames directory
    inputLayout->addWidget(new QLabel("Raw Light Frames:"), 0, 0);
    m_sourceDirectoryEdit = new QLineEdit;
    m_selectSourceButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_sourceDirectoryEdit, 0, 1);
    inputLayout->addWidget(m_selectSourceButton, 0, 2);
    
    // Dark frames directory
    inputLayout->addWidget(new QLabel("Dark Frames:"), 1, 0);
    m_darkDirectoryEdit = new QLineEdit;
    m_selectDarkButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_darkDirectoryEdit, 1, 1);
    inputLayout->addWidget(m_selectDarkButton, 1, 2);
    
    QHBoxLayout *darkInfoLayout = new QHBoxLayout;
    m_darkFramesCount = new QLabel("No dark frames loaded");
    m_refreshDarkButton = new QPushButton("Refresh");
    darkInfoLayout->addWidget(m_darkFramesCount);
    darkInfoLayout->addWidget(m_refreshDarkButton);
    darkInfoLayout->addStretch();
    inputLayout->addLayout(darkInfoLayout, 2, 1, 1, 2);
    
    // Calibrated light frames directory
    inputLayout->addWidget(new QLabel("Calibrated Lights:"), 3, 0);
    m_calibratedDirectoryEdit = new QLineEdit;
    m_selectCalibratedButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_calibratedDirectoryEdit, 3, 1);
    inputLayout->addWidget(m_selectCalibratedButton, 3, 2);
    
    // Plate-solved light frames directory
    inputLayout->addWidget(new QLabel("Plate-Solved Lights:"), 4, 0);
    m_plateSolvedDirectoryEdit = new QLineEdit;
    m_selectPlateSolvedButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_plateSolvedDirectoryEdit, 4, 1);
    inputLayout->addWidget(m_selectPlateSolvedButton, 4, 2);
    
    // Final stacked images directory
    inputLayout->addWidget(new QLabel("Stacked Images:"), 5, 0);
    m_stackedDirectoryEdit = new QLineEdit;
    m_selectStackedButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_stackedDirectoryEdit, 5, 1);
    inputLayout->addWidget(m_selectStackedButton, 5, 2);
    
    // Basic options group
    m_basicOptionsGroup = new QGroupBox("Basic Processing Options");
    QGridLayout *basicOptionsLayout = new QGridLayout(m_basicOptionsGroup);
    
    m_qualityFilterCheck = new QCheckBox("Enable Stellina quality filtering");
    m_qualityFilterCheck->setChecked(true);
    m_debugModeCheck = new QCheckBox("Debug mode (verbose logging)");
    
    basicOptionsLayout->addWidget(m_qualityFilterCheck, 0, 0, 1, 2);
    basicOptionsLayout->addWidget(m_debugModeCheck, 1, 0, 1, 2);
    
    basicOptionsLayout->addWidget(new QLabel("Focal Length (mm):"), 2, 0);
    m_focalLengthSpin = new QDoubleSpinBox;
    m_focalLengthSpin->setRange(50.0, 5000.0);
    m_focalLengthSpin->setValue(400.0);
    m_focalLengthSpin->setSuffix(" mm");
    basicOptionsLayout->addWidget(m_focalLengthSpin, 2, 1);
    
    basicOptionsLayout->addWidget(new QLabel("Pixel Size (μm):"), 3, 0);
    m_pixelSizeSpin = new QDoubleSpinBox;
    m_pixelSizeSpin->setRange(1.0, 20.0);
    m_pixelSizeSpin->setValue(2.40);
    m_pixelSizeSpin->setDecimals(2);
    m_pixelSizeSpin->setSuffix(" μm");
    basicOptionsLayout->addWidget(m_pixelSizeSpin, 3, 1);
    
    basicOptionsLayout->addWidget(new QLabel("Observer Location:"), 4, 0);
    m_observerLocationEdit = new QLineEdit("London");
    basicOptionsLayout->addWidget(m_observerLocationEdit, 4, 1);

    // Mount tilt correction group
    m_mountTiltGroup = new QGroupBox("Mount Tilt Correction");
    QGridLayout *tiltLayout = new QGridLayout(m_mountTiltGroup);

    m_enableTiltCorrectionCheck = new QCheckBox("Enable mount tilt correction");
    m_enableTiltCorrectionCheck->setChecked(false);
    tiltLayout->addWidget(m_enableTiltCorrectionCheck, 0, 0, 1, 3);

    tiltLayout->addWidget(new QLabel("North Tilt θ_N (°):"), 1, 0);
    m_northTiltSpin = new QDoubleSpinBox;
    m_northTiltSpin->setRange(-10.0, 10.0);
    m_northTiltSpin->setValue(1.0832);  // Default from your analysis
    m_northTiltSpin->setDecimals(4);
    m_northTiltSpin->setSuffix("°");
    tiltLayout->addWidget(m_northTiltSpin, 1, 1);

    tiltLayout->addWidget(new QLabel("East Tilt θ_E (°):"), 2, 0);
    m_eastTiltSpin = new QDoubleSpinBox;
    m_eastTiltSpin->setRange(-10.0, 10.0);
    m_eastTiltSpin->setValue(2.4314);   // Default from your analysis
    m_eastTiltSpin->setDecimals(4);
    m_eastTiltSpin->setSuffix("°");
    tiltLayout->addWidget(m_eastTiltSpin, 2, 1);

    QHBoxLayout *tiltButtonLayout = new QHBoxLayout;
    m_calibrateTiltButton = new QPushButton("Calibrate Tilt");
    m_testTiltButton = new QPushButton("Test Correction");
    tiltButtonLayout->addWidget(m_calibrateTiltButton);
    tiltButtonLayout->addWidget(m_testTiltButton);
    tiltButtonLayout->addStretch();
    tiltLayout->addLayout(tiltButtonLayout, 3, 0, 1, 3);

    m_tiltStatusLabel = new QLabel("Tilt correction disabled");
    m_tiltStatusLabel->setStyleSheet("color: gray; font-style: italic;");
    tiltLayout->addWidget(m_tiltStatusLabel, 4, 0, 1, 3);

    m_enableDriftCorrectionCheck = new QCheckBox("Enable drift correction");
    m_enableDriftCorrectionCheck->setChecked(false);
    tiltLayout->addWidget(m_enableDriftCorrectionCheck, 5, 0, 1, 3);

    tiltLayout->addWidget(new QLabel("RA Drift (°/h):"), 6, 0);
    m_driftRASpin = new QDoubleSpinBox;
    m_driftRASpin->setRange(-10.0, 10.0);
    m_driftRASpin->setValue(0.0);
    m_driftRASpin->setDecimals(3);
    m_driftRASpin->setSuffix("°/h");
    tiltLayout->addWidget(m_driftRASpin, 6, 1);

    tiltLayout->addWidget(new QLabel("Dec Drift (°/h):"), 7, 0);
    m_driftDecSpin = new QDoubleSpinBox;
    m_driftDecSpin->setRange(-10.0, 10.0);
    m_driftDecSpin->setValue(0.0);
    m_driftDecSpin->setDecimals(3);
    m_driftDecSpin->setSuffix("°/h");
    tiltLayout->addWidget(m_driftDecSpin, 7, 1);

    m_driftStatusLabel = new QLabel("Drift correction disabled");
    m_driftStatusLabel->setStyleSheet("color: gray; font-style: italic;");
    tiltLayout->addWidget(m_driftStatusLabel, 8, 0, 1, 3);
    
    // Processing group
    m_processingGroup = new QGroupBox("Processing Control");
    QVBoxLayout *processingLayout = new QVBoxLayout(m_processingGroup);
    
    QHBoxLayout *buttonLayout = new QHBoxLayout;
    m_startButton = new QPushButton("Start Processing");
    m_stopButton = new QPushButton("Stop Processing");
    m_stopButton->setEnabled(false);
    
    buttonLayout->addWidget(m_startButton);
    buttonLayout->addWidget(m_stopButton);
    buttonLayout->addStretch();
    
    m_progressBar = new QProgressBar;
    m_progressLabel = new QLabel("Ready");
    m_timeEstimateLabel = new QLabel("");
    m_currentTaskLabel = new QLabel("");
    
    m_subTaskProgressBar = new QProgressBar;
    m_subTaskProgressBar->setVisible(false);
    
    processingLayout->addLayout(buttonLayout);
    processingLayout->addWidget(m_progressBar);
    processingLayout->addWidget(m_subTaskProgressBar);
    processingLayout->addWidget(m_progressLabel);
    processingLayout->addWidget(m_timeEstimateLabel);
    processingLayout->addWidget(m_currentTaskLabel);
    
    // Advanced processing info
    m_advancedInfoGroup = new QGroupBox("Processing Status");
    QGridLayout *advancedLayout = new QGridLayout(m_advancedInfoGroup);
    
    m_darkCalibrationStatusLabel = new QLabel("Dark Calibration: Ready");
    m_registrationStatusLabel = new QLabel("Registration: Ready");
    m_stackingStatusLabel = new QLabel("Stacking: Ready");
    
    advancedLayout->addWidget(m_darkCalibrationStatusLabel, 0, 0);
    advancedLayout->addWidget(m_registrationStatusLabel, 0, 1);
    advancedLayout->addWidget(m_stackingStatusLabel, 0, 2);
    
    // Add all groups to main layout
//    layout->addWidget(m_connectionGroup);
    layout->addWidget(m_modeGroup);
    layout->addWidget(m_inputGroup);
    layout->addWidget(m_basicOptionsGroup);
    layout->addWidget(m_mountTiltGroup);
    layout->addWidget(m_processingGroup);
    layout->addWidget(m_advancedInfoGroup);
    layout->addStretch();
}

void StellinaProcessor::setupDarkTab() {
    QVBoxLayout *layout = new QVBoxLayout(m_darkTab);
    
    // Dark calibration options
    m_darkOptionsGroup = new QGroupBox("Dark Calibration Settings");
    QGridLayout *darkOptionsLayout = new QGridLayout(m_darkOptionsGroup);
    
    m_autoMatchDarksCheck = new QCheckBox("Automatically match dark frames to light frames");
    m_autoMatchDarksCheck->setChecked(true);
    darkOptionsLayout->addWidget(m_autoMatchDarksCheck, 0, 0, 1, 2);
    
    darkOptionsLayout->addWidget(new QLabel("Temperature Tolerance (°C):"), 1, 0);
    m_temperatureToleranceSpin = new QSpinBox;
    m_temperatureToleranceSpin->setRange(1, 20);
    m_temperatureToleranceSpin->setValue(5);
    darkOptionsLayout->addWidget(m_temperatureToleranceSpin, 1, 1);
    
    darkOptionsLayout->addWidget(new QLabel("Exposure Tolerance (%):"), 2, 0);
    m_exposureToleranceSpin = new QSpinBox;
    m_exposureToleranceSpin->setRange(1, 50);
    m_exposureToleranceSpin->setValue(10);
    darkOptionsLayout->addWidget(m_exposureToleranceSpin, 2, 1);
    
    // Add info label about master dark creation
    QLabel *masterDarkInfo = new QLabel("Master dark frames are automatically created from matching dark frames");
    masterDarkInfo->setStyleSheet("color: gray; font-style: italic;");
    masterDarkInfo->setWordWrap(true);
    darkOptionsLayout->addWidget(masterDarkInfo, 3, 0, 1, 2);
    
    // Dark frames table
    QGroupBox *darkTableGroup = new QGroupBox("Available Dark Frames");
    QVBoxLayout *darkTableLayout = new QVBoxLayout(darkTableGroup);
    
    m_darkFramesTable = new QTableWidget;
    m_darkFramesTable->setColumnCount(5);
    QStringList headers = {"Filename", "Exposure (s)", "Temperature (°C)", "Binning", "Count"};
    m_darkFramesTable->setHorizontalHeaderLabels(headers);
    m_darkFramesTable->horizontalHeader()->setStretchLastSection(true);
    m_darkFramesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_darkFramesTable->setAlternatingRowColors(true);
    
    darkTableLayout->addWidget(m_darkFramesTable);
    
    layout->addWidget(m_darkOptionsGroup);
    layout->addWidget(darkTableGroup);
    layout->addStretch();
}

void StellinaProcessor::setupStackingTab() {
    QVBoxLayout *layout = new QVBoxLayout(m_stackingTab);
/*
    // Stacking options group
    m_stackingOptionsGroup = new QGroupBox("Astrometric Stacking Settings");
    QGridLayout *stackingLayout = new QGridLayout(m_stackingOptionsGroup);
    
    stackingLayout->addWidget(new QLabel("Stacking Method:"), 0, 0);
    m_stackingMethodCombo = new QComboBox;
    m_stackingMethodCombo->addItems({"median", "mean", "sum", "sigma_clipping"});
    stackingLayout->addWidget(m_stackingMethodCombo, 0, 1);
    
    stackingLayout->addWidget(new QLabel("Rejection Method:"), 1, 0);
    m_rejectionMethodCombo = new QComboBox;
    m_rejectionMethodCombo->addItems({"none", "sigma", "linear", "percentile"});
    m_rejectionMethodCombo->setCurrentText("sigma");
    stackingLayout->addWidget(m_rejectionMethodCombo, 1, 1);
    
    stackingLayout->addWidget(new QLabel("Low Rejection:"), 2, 0);
    m_rejectionLowSpin = new QDoubleSpinBox;
    m_rejectionLowSpin->setRange(0.1, 10.0);
    m_rejectionLowSpin->setValue(3.0);
    m_rejectionLowSpin->setDecimals(1);
    stackingLayout->addWidget(m_rejectionLowSpin, 2, 1);
    
    stackingLayout->addWidget(new QLabel("High Rejection:"), 3, 0);
    m_rejectionHighSpin = new QDoubleSpinBox;
    m_rejectionHighSpin->setRange(0.1, 10.0);
    m_rejectionHighSpin->setValue(3.0);
    m_rejectionHighSpin->setDecimals(1);
    stackingLayout->addWidget(m_rejectionHighSpin, 3, 1);
    
    m_normalizeCheck = new QCheckBox("Normalize images before stacking");
    m_normalizeCheck->setChecked(true);
    stackingLayout->addWidget(m_normalizeCheck, 4, 0, 1, 2);
    
    m_drizzleCheck = new QCheckBox("Apply drizzle enhancement");
    stackingLayout->addWidget(m_drizzleCheck, 5, 0, 1, 2);
    
    stackingLayout->addWidget(new QLabel("Drizzle Scale:"), 6, 0);
    m_drizzleScaleSpin = new QDoubleSpinBox;
    m_drizzleScaleSpin->setRange(1.0, 3.0);
    m_drizzleScaleSpin->setValue(1.5);
    m_drizzleScaleSpin->setDecimals(1);
    m_drizzleScaleSpin->setEnabled(false);
    stackingLayout->addWidget(m_drizzleScaleSpin, 6, 1);
    
    stackingLayout->addWidget(new QLabel("Output Format:"), 7, 0);
    m_outputFormatCombo = new QComboBox;
    m_outputFormatCombo->addItems({"fits", "tiff", "png"});
    stackingLayout->addWidget(m_outputFormatCombo, 7, 1);
    
    QHBoxLayout *previewLayout = new QHBoxLayout;
    m_previewStackButton = new QPushButton("Preview Stack (First 5 Images)");
    previewLayout->addWidget(m_previewStackButton);
    previewLayout->addStretch();
    stackingLayout->addLayout(previewLayout, 8, 0, 1, 2);
    
    layout->addWidget(m_stackingOptionsGroup);
 */
    setupWCSStackingUI();
    layout->addStretch();
}

void StellinaProcessor::setupLogTab() {
    QVBoxLayout *layout = new QVBoxLayout(m_logTab);
    
    // Log group
    m_logGroup = new QGroupBox("Processing Log");
    QVBoxLayout *logLayout = new QVBoxLayout(m_logGroup);
    
    m_logTextEdit = new QTextEdit;
    
    // Use a monospace font that exists on macOS
    QFont monoFont;
    QStringList monoFonts = {"SF Mono", "Monaco", "Menlo", "Courier New", "monospace"};
    for (const QString &fontName : monoFonts) {
        QFont testFont(fontName);
        if (QFontInfo(testFont).family() == fontName) {
            monoFont = testFont;
            break;
        }
    }
    monoFont.setPointSize(9);
    m_logTextEdit->setFont(monoFont);
    
    QHBoxLayout *logButtonLayout = new QHBoxLayout;
    m_clearLogButton = new QPushButton("Clear Log");
    m_saveLogButton = new QPushButton("Save Log...");
    logButtonLayout->addStretch();
    logButtonLayout->addWidget(m_clearLogButton);
    logButtonLayout->addWidget(m_saveLogButton);
    
    logLayout->addWidget(m_logTextEdit);
    logLayout->addLayout(logButtonLayout);
    
    layout->addWidget(m_logGroup);
}

void StellinaProcessor::setupMenu() {
    QMenuBar *menuBar = this->menuBar();
    
    // File menu
    QMenu *fileMenu = menuBar->addMenu("&File");
    fileMenu->addAction("&Save Log...", this, &StellinaProcessor::onClearLog);
    fileMenu->addSeparator();
    fileMenu->addAction("&Exit", this, &QWidget::close);
    
    // Tools menu
    QMenu *toolsMenu = menuBar->addMenu("&Tools");
    toolsMenu->addAction("&Refresh Dark Frames", this, &StellinaProcessor::onRefreshDarkFrames);
    addWCSMenuItems();
    toolsMenu->addSeparator();

    // Enhanced Mount tilt submenu
    QMenu *tiltMenu = toolsMenu->addMenu("Mount &Tilt");
    /*
    tiltMenu->addAction("&Enhanced Calibration (with Drift)", this, &StellinaProcessor::calibrateMountWithDrift);
    tiltMenu->addSeparator();
    tiltMenu->addAction("&Basic Auto-Calibrate", this, &StellinaProcessor::autoCalibrateTiltFromSolveResults);
    tiltMenu->addAction("&Fine-Tune Parameters", this, &StellinaProcessor::fineTuneMountTilt);
     */
    tiltMenu->addSeparator();
    tiltMenu->addAction("&Test Tilt Correction", this, &StellinaProcessor::testMountTiltCorrection);
    tiltMenu->addAction("&Manual Calibrate", this, &StellinaProcessor::calibrateMountTilt);
    tiltMenu->addSeparator();
    tiltMenu->addAction("Enable/Disable Tilt Correction", [this]() {
        m_mountTilt.enableCorrection = !m_mountTilt.enableCorrection;
        updateTiltUI();
        saveMountTiltToSettings();
        logMessage(QString("Mount tilt correction %1").arg(m_mountTilt.enableCorrection ? "enabled" : "disabled"), "blue");
    });
    tiltMenu->addAction("Enable/Disable Drift Correction", [this]() {
        m_mountTilt.enableDriftCorrection = !m_mountTilt.enableDriftCorrection;
        saveMountTiltToSettings();
        logMessage(QString("Mount drift correction %1").arg(m_mountTilt.enableDriftCorrection ? "enabled" : "disabled"), "blue");
    });
    tiltMenu->addAction("&Auto-Calibrate from Processed Files", this, &StellinaProcessor::calibrateFromProcessedFiles);
    tiltMenu->addAction("&Test Systematic Correction", this, &StellinaProcessor::testSystematicOffsetCorrection);
    tiltMenu->addAction("&Verify Offsets in Use", this, &StellinaProcessor::verifySystematicOffsetsInUse);
    tiltMenu->addAction("&Fast Calibrate from Stacking JSON", this, &StellinaProcessor::calibrateFromStackingJSON);
 
    toolsMenu->addSeparator();

    toolsMenu->addAction("&Clear Log", this, &StellinaProcessor::onClearLog);

    toolsMenu->addAction("Diagnose Sidereal Time", this, &StellinaProcessor::diagnoseSiderealTimeIssues);
    toolsMenu->addAction("Diagnose Coordinate conversion", this, &StellinaProcessor::testFixedCoordinateConversion);
    toolsMenu->addAction("Diagnose Tracking Issue", this, &StellinaProcessor::diagnoseTrackingIssue);
    toolsMenu->addAction("Analyze Coordinate errors", this, &StellinaProcessor::analyzeRealCoordinateErrors);
    toolsMenu->addAction("Test All Coordinate Variations", this, &StellinaProcessor::testAllCoordinateVariations);
    toolsMenu->addAction("Test Time Drift Fix", this, &StellinaProcessor::testTimeDriftFix);
    toolsMenu->addAction("analyzeRealStellinaIssue", this, &StellinaProcessor::analyzeRealStellinaIssue);
    toolsMenu->addAction("testRealisticAccuracy", this, &StellinaProcessor::testRealisticAccuracy);
    toolsMenu->addAction("verifyPlatesolvingHints", this, &StellinaProcessor::verifyPlatesolvingHints);
    // Add to your Tools menu:
    toolsMenu->addAction("Dump Coordinate Data", this, &StellinaProcessor::dumpCoordinateData);
    toolsMenu->addAction("Export Coordinates to CSV", this, &StellinaProcessor::dumpCoordinateDataToCSV);
    toolsMenu->addAction("Analyze Coordinate Drift", this, &StellinaProcessor::analyzeCoordinateDrift);
    toolsMenu->addAction("&Plot Mount Errors", this, &StellinaProcessor::plotMountErrors);    
    // Help menu
    QMenu *helpMenu = menuBar->addMenu("&Help");
    helpMenu->addAction("&About", [this]() {
        QMessageBox::about(this, "About Enhanced Stellina Processor",
                          "Enhanced Stellina Processor v2.0\n\n"
                          "A Qt application for advanced processing of Stellina telescope images\n"
                          "New Features:\n"
                          "• Dark frame calibration with automatic matching\n"
                          "• Astrometric registration and stacking\n"
                          "• Full processing pipeline\n"
                          "• Advanced rejection algorithms\n"
                          "• Drizzle enhancement support\n"
                          "• Comprehensive processing reports");
    });
}
