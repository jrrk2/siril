#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QStatusBar>
#include <QMenuBar>
#include <QMenu>
#include <QRandomGenerator>
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
    m_debugTab = new QWidget;  // ADD THIS LINE
    
    m_tabWidget->addTab(m_basicTab, "Basic Processing");
    m_tabWidget->addTab(m_darkTab, "Dark Calibration");
    m_tabWidget->addTab(m_stackingTab, "Astrometric Stacking");
    m_tabWidget->addTab(m_logTab, "Processing Log");
    m_tabWidget->addTab(m_debugTab, "Debug");
    
    setupBasicTab();
    setupDarkTab();
    setupStackingTab();
    setupLogTab();
    setupDebugTab();  // ADD THIS LINE
    
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
    toolsMenu->addAction("Test Time Drift Fix", this, &StellinaProcessor::testTimeDriftFix);
    toolsMenu->addAction("analyzeRealStellinaIssue", this, &StellinaProcessor::analyzeRealStellinaIssue);
    toolsMenu->addAction("testRealisticAccuracy", this, &StellinaProcessor::testRealisticAccuracy);
    toolsMenu->addAction("verifyPlatesolvingHints", this, &StellinaProcessor::verifyPlatesolvingHints);
    // Add to your Tools menu:
    toolsMenu->addAction("Dump Coordinate Data", this, &StellinaProcessor::dumpCoordinateData);
    toolsMenu->addAction("Export Coordinates to CSV", this, &StellinaProcessor::dumpCoordinateDataToCSV);
    toolsMenu->addAction("Analyze Coordinate Drift", this, &StellinaProcessor::analyzeCoordinateDrift);
    toolsMenu->addAction("&Plot Mount Errors", this, &StellinaProcessor::plotMountErrors);
    toolsMenu->addAction("Validate Timing Offset", this, &StellinaProcessor::validateTimingOffset);
    toolsMenu->addAction("Compare Drift (Dynamic Offset)", this, &StellinaProcessor::compareDriftWithDynamicOffset);
    toolsMenu->addAction("Update Processing (Dynamic Offset)", this, &StellinaProcessor::updateProcessingToDynamicOffset);

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

// =============================================================================
// 3. IMPLEMENT setupDebugTab()
// =============================================================================

void StellinaProcessor::setupDebugTab() {
    QVBoxLayout *layout = new QVBoxLayout(m_debugTab);
    
    // Main coordinate debug group
    m_coordDebugGroup = new QGroupBox("Manual Coordinate Conversion Test");
    QGridLayout *coordLayout = new QGridLayout(m_coordDebugGroup);
    
    // Input controls
    coordLayout->addWidget(new QLabel("Altitude (°):"), 0, 0);
    m_debugAltSpin = new QDoubleSpinBox;
    m_debugAltSpin->setRange(-90.0, 90.0);
    m_debugAltSpin->setValue(42.0410);  // Default from your test data
    m_debugAltSpin->setDecimals(4);
    m_debugAltSpin->setSuffix("°");
    coordLayout->addWidget(m_debugAltSpin, 0, 1);
    
    coordLayout->addWidget(new QLabel("Azimuth (°):"), 0, 2);
    m_debugAzSpin = new QDoubleSpinBox;
    m_debugAzSpin->setRange(0.0, 360.0);
    m_debugAzSpin->setValue(286.8526);  // Default from your test data
    m_debugAzSpin->setDecimals(4);
    m_debugAzSpin->setSuffix("°");
    coordLayout->addWidget(m_debugAzSpin, 0, 3);
    
    coordLayout->addWidget(new QLabel("Observation Time:"), 1, 0);
    m_debugTimeEdit = new QLineEdit("2024-01-09T22:13:29");
    m_debugTimeEdit->setPlaceholderText("YYYY-MM-DDTHH:MM:SS");
    coordLayout->addWidget(m_debugTimeEdit, 1, 1, 1, 3);
    
    coordLayout->addWidget(new QLabel("Observer Latitude (°):"), 2, 0);
    m_debugLatSpin = new QDoubleSpinBox;
    m_debugLatSpin->setRange(-90.0, 90.0);
    m_debugLatSpin->setValue(51.5074);  // London default
    m_debugLatSpin->setDecimals(4);
    m_debugLatSpin->setSuffix("°");
    coordLayout->addWidget(m_debugLatSpin, 2, 1);
    
    coordLayout->addWidget(new QLabel("Observer Longitude (°):"), 2, 2);
    m_debugLonSpin = new QDoubleSpinBox;
    m_debugLonSpin->setRange(-180.0, 180.0);
    m_debugLonSpin->setValue(-0.1278);  // London default
    m_debugLonSpin->setDecimals(4);
    m_debugLonSpin->setSuffix("°");
    coordLayout->addWidget(m_debugLonSpin, 2, 3);
    
    // Buttons
    QHBoxLayout *buttonLayout = new QHBoxLayout;
    m_testConversionButton = new QPushButton("Test Single Conversion");
    m_testConversionButton->setStyleSheet("QPushButton { background-color: #2196F3; color: white; font-weight: bold; }");
    m_loadImageDataButton = new QPushButton("Load from Current Images");
    m_testBatchButton = new QPushButton("Test Batch Processing");
    
    buttonLayout->addWidget(m_testConversionButton);
    buttonLayout->addWidget(m_loadImageDataButton);
    buttonLayout->addWidget(m_testBatchButton);
    buttonLayout->addStretch();
    
    coordLayout->addLayout(buttonLayout, 3, 0, 1, 4);
    
    // Preset tests group
    m_presetTestsGroup = new QGroupBox("Preset Test Cases");
    QHBoxLayout *presetLayout = new QHBoxLayout(m_presetTestsGroup);
    
    m_presetTestCombo = new QComboBox;
    m_presetTestCombo->addItem("Your Original Test Data", 0);
    m_presetTestCombo->addItem("Zenith Test (Alt=90°)", 1);
    m_presetTestCombo->addItem("Horizon Test (Alt=0°)", 2);
    m_presetTestCombo->addItem("North Pole Test", 3);
    m_presetTestCombo->addItem("Celestial Equator Test", 4);
    m_presetTestCombo->addItem("Known Star Positions", 5);
    
    m_runPresetButton = new QPushButton("Run Preset Test");
    m_runPresetButton->setStyleSheet("QPushButton { background-color: #FF9800; color: white; font-weight: bold; }");
    
    presetLayout->addWidget(new QLabel("Test Case:"));
    presetLayout->addWidget(m_presetTestCombo);
    presetLayout->addWidget(m_runPresetButton);
    presetLayout->addStretch();
    
    // Sky region tests
    m_skyRegionGroup = new QGroupBox("Sky Region Tests");
    QHBoxLayout *skyLayout = new QHBoxLayout(m_skyRegionGroup);
    
    m_skyRegionCombo = new QComboBox;
    m_skyRegionCombo->addItem("Orion (Winter Evening)", 0);
    m_skyRegionCombo->addItem("Ursa Major (Spring Evening)", 1);
    m_skyRegionCombo->addItem("Cygnus (Summer Evening)", 2);
    m_skyRegionCombo->addItem("Cassiopeia (Autumn Evening)", 3);
    m_skyRegionCombo->addItem("Southern Cross (Southern Sky)", 4);
    m_skyRegionCombo->addItem("Near Zenith", 5);
    m_skyRegionCombo->addItem("Near Horizon", 6);
    
    m_testSkyRegionButton = new QPushButton("Test Sky Region");
    m_testSkyRegionButton->setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; }");
    
    skyLayout->addWidget(new QLabel("Region:"));
    skyLayout->addWidget(m_skyRegionCombo);
    skyLayout->addWidget(m_testSkyRegionButton);
    skyLayout->addStretch();
    
    // Results display
    QGroupBox *resultsGroup = new QGroupBox("Debug Results");
    QVBoxLayout *resultsLayout = new QVBoxLayout(resultsGroup);
    
    m_debugResultsEdit = new QTextEdit;
    m_debugResultsEdit->setMinimumHeight(300);
    
    // Use monospace font for results
    QFont monoFont;
    QStringList monoFonts = {"SF Mono", "Monaco", "Menlo", "Courier New", "monospace"};
    for (const QString &fontName : monoFonts) {
        QFont testFont(fontName);
        if (QFontInfo(testFont).family() == fontName) {
            monoFont = testFont;
            break;
        }
    }
    monoFont.setPointSize(10);
    m_debugResultsEdit->setFont(monoFont);
    
    QPushButton *clearResultsButton = new QPushButton("Clear Results");
    QPushButton *saveResultsButton = new QPushButton("Save Results...");
    
    QHBoxLayout *resultsButtonLayout = new QHBoxLayout;
    resultsButtonLayout->addStretch();
    resultsButtonLayout->addWidget(clearResultsButton);
    resultsButtonLayout->addWidget(saveResultsButton);
    
    resultsLayout->addWidget(m_debugResultsEdit);
    resultsLayout->addLayout(resultsButtonLayout);

    // Add comprehensive test suite buttons
    QGroupBox *testSuiteGroup = new QGroupBox("Comprehensive Test Suite");
    QVBoxLayout *testSuiteLayout = new QVBoxLayout(testSuiteGroup);
    
    QHBoxLayout *suiteButtonLayout = new QHBoxLayout;
    
    QPushButton *runFullSuiteButton = new QPushButton("Run Full Test Suite (20 tests)");
    runFullSuiteButton->setStyleSheet("QPushButton { background-color: #E91E63; color: white; font-weight: bold; }");
    
    QPushButton *runRandomSubsetButton = new QPushButton("Run Random Subset (10 tests)");
    runRandomSubsetButton->setStyleSheet("QPushButton { background-color: #9C27B0; color: white; font-weight: bold; }");
    
    QPushButton *runAccuracyAnalysisButton = new QPushButton("Accuracy Analysis");
    runAccuracyAnalysisButton->setStyleSheet("QPushButton { background-color: #673AB7; color: white; font-weight: bold; }");
    
    suiteButtonLayout->addWidget(runFullSuiteButton);
    suiteButtonLayout->addWidget(runRandomSubsetButton);
    suiteButtonLayout->addWidget(runAccuracyAnalysisButton);
    
    QLabel *testSuiteInfo = new QLabel("NASA JPL Horizons validation tests with known precise coordinates");
    testSuiteInfo->setStyleSheet("color: gray; font-style: italic;");
    testSuiteInfo->setWordWrap(true);
    
    testSuiteLayout->addLayout(suiteButtonLayout);
    testSuiteLayout->addWidget(testSuiteInfo);
    
    // Add all groups to main layout
    layout->addWidget(m_coordDebugGroup);
    layout->addWidget(m_presetTestsGroup);
    layout->addWidget(m_skyRegionGroup);
    layout->addWidget(resultsGroup);
    layout->addStretch();
    layout->addWidget(testSuiteGroup);
    
    // Connect signals
    connect(runFullSuiteButton, &QPushButton::clicked, this, &StellinaProcessor::runCoordinateTestSuite);
    connect(runRandomSubsetButton, &QPushButton::clicked, [this]() { runRandomTestSubset(10); });
    connect(runAccuracyAnalysisButton, &QPushButton::clicked, this, &StellinaProcessor::runAccuracyAnalysis);
    connect(m_testConversionButton, &QPushButton::clicked, this, &StellinaProcessor::onTestConversion);
    connect(m_loadImageDataButton, &QPushButton::clicked, this, &StellinaProcessor::onLoadImageData);
    connect(m_testBatchButton, &QPushButton::clicked, this, &StellinaProcessor::onTestBatch);
    connect(m_runPresetButton, &QPushButton::clicked, this, &StellinaProcessor::onRunPresetTest);
    connect(m_testSkyRegionButton, &QPushButton::clicked, this, &StellinaProcessor::onTestSkyRegion);
    
    connect(clearResultsButton, &QPushButton::clicked, m_debugResultsEdit, &QTextEdit::clear);
    connect(saveResultsButton, &QPushButton::clicked, [this]() {
        QString fileName = QFileDialog::getSaveFileName(this, "Save Debug Results", 
                                                       "coordinate_debug_results.txt", 
                                                       "Text Files (*.txt)");
        if (!fileName.isEmpty()) {
            QFile file(fileName);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << m_debugResultsEdit->toPlainText();
                debugLog("Debug results saved to: " + fileName);
            }
        }
    });
    
}

// =============================================================================
// 4. IMPLEMENT SLOT FUNCTIONS
// =============================================================================

void StellinaProcessor::onTestConversion() {
    debugLog("=== MANUAL COORDINATE CONVERSION TEST ===");
    
    double alt = m_debugAltSpin->value();
    double az = m_debugAzSpin->value();
    QString timeStr = m_debugTimeEdit->text();
    double lat = m_debugLatSpin->value();
    double lon = m_debugLonSpin->value();
    
    debugLog(QString("Input Parameters:"));
    debugLog(QString("  Altitude:  %1°").arg(alt, 0, 'f', 4));
    debugLog(QString("  Azimuth:   %1°").arg(az, 0, 'f', 4));
    debugLog(QString("  Time:      %1").arg(timeStr));
    debugLog(QString("  Observer:  %1°N, %2°E").arg(lat, 0, 'f', 4).arg(lon, 0, 'f', 4));
    debugLog("");
    
    // Save current observer location and set test location
    QString savedLocation = m_observerLocation;
    m_observerLocation = QString("%1,%2").arg(lat, 0, 'f', 4).arg(lon, 0, 'f', 4);
    
    // Test the conversion
    double ra, dec;
    bool success = convertAltAzToRaDec(alt, az, timeStr, ra, dec);
    
    if (success) {
        debugLog("Conversion Result:");
        debugLog(QString("  RA:  %1° (%2h %3m %4s)")
                    .arg(ra, 0, 'f', 6)
                    .arg(static_cast<int>(ra / 15.0))
                    .arg(static_cast<int>(fmod(ra / 15.0, 1.0) * 60))
                    .arg(fmod(fmod(ra / 15.0, 1.0) * 60, 1.0) * 60, 0, 'f', 1));
        debugLog(QString("  Dec: %1° (%2° %3' %4\")")
                    .arg(dec, 0, 'f', 6)
                    .arg(static_cast<int>(dec))
                    .arg(static_cast<int>(fabs(dec - static_cast<int>(dec)) * 60))
                    .arg(fmod(fabs(dec - static_cast<int>(dec)) * 60, 1.0) * 60, 0, 'f', 1));
        
        // Calculate LST for additional debugging info
        QDateTime obsTime = QDateTime::fromString(timeStr, "yyyy-MM-ddThh:mm:ss");
        obsTime.setTimeSpec(Qt::UTC);
        
        if (obsTime.isValid()) {
            double jd = CoordinateUtils::computeJulianDay(obsTime.date().year(),
                                                          obsTime.date().month(),
                                                          obsTime.date().day(),
                                                          obsTime.time().hour(),
                                                          obsTime.time().minute(),
                                                          obsTime.time().second());
            double lst = calculateLST_HighPrecision(jd, lon);
            
            debugLog("");
            debugLog("Additional Info:");
            debugLog(QString("  Julian Day: %1").arg(jd, 0, 'f', 6));
            debugLog(QString("  LST:        %1h (%2°)").arg(lst, 0, 'f', 4).arg(lst * 15.0, 0, 'f', 2));
            
            // Calculate hour angle
            double hourAngle = lst * 15.0 - ra;
            while (hourAngle < -180.0) hourAngle += 360.0;
            while (hourAngle > 180.0) hourAngle -= 360.0;
            debugLog(QString("  Hour Angle: %1°").arg(hourAngle, 0, 'f', 4));
        }
        
    } else {
        debugLog("❌ CONVERSION FAILED");
    }
    
    // Restore original observer location
    m_observerLocation = savedLocation;
    
    debugLog("=====================================");
    debugLog("");
}

void StellinaProcessor::onLoadImageData() {
    debugLog("=== LOADING IMAGE DATA ===");
    
    if (m_stellinaImageData.isEmpty()) {
        debugLog("❌ No Stellina image data loaded");
        debugLog("Please load images first from the Basic Processing tab");
        debugLog("");
        return;
    }
    
    // Use the first image's data
    const StellinaImageData &firstImage = m_stellinaImageData.first();
    
    m_debugAltSpin->setValue(firstImage.altitude);
    m_debugAzSpin->setValue(firstImage.azimuth);
    m_debugTimeEdit->setText(firstImage.dateObs);
    
    debugLog(QString("✓ Loaded data from: %1").arg(QFileInfo(firstImage.originalFitsPath).fileName()));
    debugLog(QString("  Altitude:  %1°").arg(firstImage.altitude, 0, 'f', 4));
    debugLog(QString("  Azimuth:   %1°").arg(firstImage.azimuth, 0, 'f', 4));
    debugLog(QString("  Time:      %1").arg(firstImage.dateObs));
    
    // If the image has calculated coordinates, show them for comparison
    if (firstImage.hasCalculatedCoords) {
        debugLog(QString("  Has pre-calculated coordinates:"));
        debugLog(QString("    RA:  %1°").arg(firstImage.calculatedRA, 0, 'f', 6));
        debugLog(QString("    Dec: %1°").arg(firstImage.calculatedDec, 0, 'f', 6));
    }
    
    debugLog("");
}

void StellinaProcessor::onTestBatch() {
    debugLog("=== BATCH COORDINATE CONVERSION TEST ===");
    
    if (m_stellinaImageData.isEmpty()) {
        debugLog("❌ No Stellina image data loaded");
        debugLog("Please load images first from the Basic Processing tab");
        debugLog("");
        return;
    }
    
    debugLog(QString("Testing conversion for %1 images:").arg(m_stellinaImageData.size()));
    debugLog("");
    
    // Save current observer location
    QString savedLocation = m_observerLocation;
    double lat = m_debugLatSpin->value();
    double lon = m_debugLonSpin->value();
    m_observerLocation = QString("%1,%2").arg(lat, 0, 'f', 4).arg(lon, 0, 'f', 4);
    
    int successCount = 0;
    double totalRAError = 0.0;
    double totalDecError = 0.0;
    int errorCount = 0;
    
    for (int i = 0; i < qMin(10, m_stellinaImageData.size()); ++i) {  // Limit to first 10 images
        const StellinaImageData &imageData = m_stellinaImageData[i];
        
        double ra, dec;
        bool success = convertAltAzToRaDec(imageData.altitude, imageData.azimuth, imageData.dateObs, ra, dec);
        
        QString fileName = QFileInfo(imageData.originalFitsPath).fileName();
        
        if (success) {
            successCount++;
            debugLog(QString("%1: Alt=%2° Az=%3° → RA=%4° Dec=%5°")
                        .arg(fileName, -15)
                        .arg(imageData.altitude, 0, 'f', 3)
                        .arg(imageData.azimuth, 0, 'f', 3)
                        .arg(ra, 0, 'f', 4)
                        .arg(dec, 0, 'f', 4));
            
            // If we have expected values, calculate errors
            if (imageData.hasCalculatedCoords) {
                double raError = qAbs(ra - imageData.calculatedRA);
                double decError = qAbs(dec - imageData.calculatedDec);
                if (raError > 180.0) raError = 360.0 - raError;  // Handle RA wrap-around
                
                totalRAError += raError;
                totalDecError += decError;
                errorCount++;
                
                debugLog(QString("    Expected: RA=%1° Dec=%2° (errors: RA=%3° Dec=%4°)")
                            .arg(imageData.calculatedRA, 0, 'f', 4)
                            .arg(imageData.calculatedDec, 0, 'f', 4)
                            .arg(raError, 0, 'f', 4)
                            .arg(decError, 0, 'f', 4));
            }
        } else {
            debugLog(QString("%1: ❌ CONVERSION FAILED").arg(fileName));
        }
    }
    
    debugLog("");
    debugLog("Summary:");
    debugLog(QString("  Successful conversions: %1/%2").arg(successCount).arg(qMin(10, m_stellinaImageData.size())));
    
    if (errorCount > 0) {
        double avgRAError = totalRAError / errorCount;
        double avgDecError = totalDecError / errorCount;
        debugLog(QString("  Average errors: RA=%1° (%2 arcmin), Dec=%3° (%4 arcmin)")
                    .arg(avgRAError, 0, 'f', 4)
                    .arg(avgRAError * 60.0, 0, 'f', 1)
                    .arg(avgDecError, 0, 'f', 4)
                    .arg(avgDecError * 60.0, 0, 'f', 1));
        
        if (avgRAError > 1.0 || avgDecError > 1.0) {
            debugLog("  ⚠️  Large systematic errors detected - check conversion algorithm");
        } else {
            debugLog("  ✅ Errors within reasonable range");
        }
    }
    
    // Restore original observer location
    m_observerLocation = savedLocation;
    
    debugLog("=====================================");
    debugLog("");
}

void StellinaProcessor::onRunPresetTest() {
    int testIndex = m_presetTestCombo->currentData().toInt();
    
    debugLog(QString("=== PRESET TEST: %1 ===").arg(m_presetTestCombo->currentText()));
    
    switch (testIndex) {
        case 0: {  // Your Original Test Data
            m_debugAltSpin->setValue(42.0410);
            m_debugAzSpin->setValue(286.8526);
            m_debugTimeEdit->setText("2024-01-09T22:13:29");
            debugLog("Using your original problematic data point");
            debugLog("Expected result: RA ≈ 10.68°, Dec ≈ 41.27°");
            break;
        }
        case 1: {  // Zenith Test
            m_debugAltSpin->setValue(90.0);
            m_debugAzSpin->setValue(0.0);  // Azimuth doesn't matter at zenith
            debugLog("Testing zenith position (Alt=90°)");
            debugLog("Should give Dec = Observer Latitude");
            break;
        }
        case 2: {  // Horizon Test
            m_debugAltSpin->setValue(0.0);
            m_debugAzSpin->setValue(0.0);  // Due North
            debugLog("Testing horizon position due North");
            debugLog("Should give Dec = 90° - Observer Latitude");
            break;
        }
        case 3: {  // North Pole Test
            m_debugAltSpin->setValue(51.5074);  // Observer latitude
            m_debugAzSpin->setValue(0.0);       // Due North
            debugLog("Testing North Celestial Pole");
            debugLog("Should give RA ≈ any value, Dec ≈ 90°");
            break;
        }
        case 4: {  // Celestial Equator Test
            m_debugAltSpin->setValue(90.0 - 51.5074);  // 90° - latitude
            m_debugAzSpin->setValue(90.0);              // Due East
            debugLog("Testing Celestial Equator due East");
            debugLog("Should give Dec ≈ 0°");
            break;
        }
        case 5: {  // Known Star - Polaris
            // Polaris coordinates for the test time
            m_debugAltSpin->setValue(51.5);  // Approximately latitude for London
            m_debugAzSpin->setValue(359.5);  // Very close to North
            debugLog("Testing known star position (Polaris approximation)");
            debugLog("Should give RA ≈ 37.95°, Dec ≈ 89.26°");
            break;
        }
    }
    
    debugLog("");
    onTestConversion();  // Run the actual test
}

void StellinaProcessor::onTestSkyRegion() {
    int regionIndex = m_skyRegionCombo->currentData().toInt();
    
    debugLog(QString("=== SKY REGION TEST: %1 ===").arg(m_skyRegionCombo->currentText()));
    
    // Use winter evening time for most tests
    QString testTime = "2024-01-15T20:00:00";
    m_debugTimeEdit->setText(testTime);
    
    struct TestPosition {
        double alt;
        double az;
        QString description;
        QString expectedRA;
        QString expectedDec;
    };
    
    QList<TestPosition> positions;
    
    switch (regionIndex) {
        case 0:  // Orion (Winter Evening)
            positions = {
                {45.0, 180.0, "Orion Belt region", "~84°", "~0°"},
                {35.0, 190.0, "Orion Nebula region", "~84°", "~-5°"}
            };
            break;
        case 1:  // Ursa Major (Spring Evening)
            testTime = "2024-04-15T20:00:00";
            positions = {
                {70.0, 30.0, "Big Dipper region", "~165°", "~55°"},
                {75.0, 45.0, "Ursa Major", "~170°", "~60°"}
            };
            break;
        case 2:  // Cygnus (Summer Evening)
            testTime = "2024-07-15T20:00:00";
            positions = {
                {85.0, 90.0, "Cygnus/Deneb region", "~310°", "~45°"},
                {60.0, 80.0, "Cygnus Wing", "~300°", "~35°"}
            };
            break;
        case 3:  // Cassiopeia (Autumn Evening)
            testTime = "2024-10-15T20:00:00";
            positions = {
                {75.0, 15.0, "Cassiopeia W", "~15°", "~60°"},
                {70.0, 10.0, "Near Cassiopeia", "~10°", "~55°"}
            };
            break;
        case 4:  // Southern Cross (Southern Sky) - simulate from different latitude
            m_debugLatSpin->setValue(-35.0);  // Southern hemisphere
            positions = {
                {45.0, 180.0, "Southern Cross", "~185°", "-60°"},
                {35.0, 170.0, "Centaurus", "~210°", "-50°"}
            };
            break;
        case 5:  // Near Zenith
            positions = {
                {85.0, 0.0, "Near zenith North", "Variable", "~85°"},
                {85.0, 90.0, "Near zenith East", "Variable", "~85°"},
                {85.0, 180.0, "Near zenith South", "Variable", "~85°"},
                {85.0, 270.0, "Near zenith West", "Variable", "~85°"}
            };
            break;
        case 6:  // Near Horizon
            positions = {
                {5.0, 0.0, "Low North", "Variable", "~40°"},
                {5.0, 90.0, "Low East", "Variable", "~0°"},
                {5.0, 180.0, "Low South", "Variable", "-40°"},
                {5.0, 270.0, "Low West", "Variable", "~0°"}
            };
            break;
    }
    
    m_debugTimeEdit->setText(testTime);
    debugLog(QString("Testing %1 positions in this sky region:").arg(positions.size()));
    debugLog("");
    
    for (const TestPosition &pos : positions) {
        m_debugAltSpin->setValue(pos.alt);
        m_debugAzSpin->setValue(pos.az);
        
        debugLog(QString("Position: %1").arg(pos.description));
        debugLog(QString("Expected: RA %1, Dec %2").arg(pos.expectedRA).arg(pos.expectedDec));
        onTestConversion();
    }
    
    // Reset to London if we changed it
    if (regionIndex == 4) {
        m_debugLatSpin->setValue(51.5074);
    }
}

// =============================================================================
// 5. ADD HELPER FUNCTION FOR DEBUG LOGGING
// =============================================================================

void StellinaProcessor::debugLog(const QString &message) {
    // Add timestamp
    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    QString formattedMessage = QString("[%1] %2").arg(timestamp).arg(message);
    
    // Add to debug results
    m_debugResultsEdit->append(formattedMessage);
    
    // Also add to main log with special formatting
    logMessage(QString("DEBUG: %1").arg(message), "purple");
    
    // Auto-scroll to bottom
    QTextCursor cursor = m_debugResultsEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    m_debugResultsEdit->setTextCursor(cursor);
    }


QList<CoordinateTestCase> StellinaProcessor::getCoordinateTestCases() {
    QList<CoordinateTestCase> testCases;
    
    // Selected high-precision test cases from NASA JPL Horizons
    // Format: {name, lat, lon, date, julianDay, raOfDate, decOfDate, expectedAz, expectedAlt, siderealTime, hourAngle, description}
    
    testCases << CoordinateTestCase{
        "London_Winter_Evening", 51.5074, -0.1278, "2024-01-09T22:13:29", 2460317.426262,
        10.6760, 41.2734, 286.8526, 42.0410, 23.1069, -11.3133,
        "London winter evening - your original problem case"
    };
    
    testCases << CoordinateTestCase{
        "Zenith_Test", 51.5074, -0.1278, "2024-01-09T12:00:00", 2460316.0,
        180.0, 51.5074, 0.0, 90.0, 12.0, 0.0,
        "Object at zenith - should give Dec = latitude"
    };
    
    testCases << CoordinateTestCase{
        "North_Horizon", 51.5074, -0.1278, "2024-01-09T12:00:00", 2460316.0,
        180.0, 38.4926, 0.0, 0.0, 12.0, 0.0,
        "North horizon - should give Dec = 90° - latitude"
    };
    
    // Diverse test cases from different locations and times
    testCases << CoordinateTestCase{
        "Arctic_Summer", 77.2069702520977, 118.639627806683, "1472-08-18T16:47:23", 2258936.199571759,
        156.3041080, 9.9129372, 10.1558, -2.6785, 23.1069400896, -11.313333776,
        "Arctic location, medieval time"
    };
    
    testCases << CoordinateTestCase{
        "Antarctic_Modern", -45.0396581387729, 52.6267552330626, "1527-10-21T22:24:53", 2279088.433946759,
        198.2076359, -5.3974508, 124.4669, -22.8271, 4.5307070225, -8.683135370,
        "Antarctic conditions"
    };
    
    testCases << CoordinateTestCase{
        "Tropical_Future", 75.8072259286161, 308.733126034821, "2549-01-01T20:31:11", 2652064.354988426,
        303.4519425, -21.2960037, 231.4286, -12.7578, 23.8918905112, 3.661761012,
        "Far future tropical observation"
    };
    
    testCases << CoordinateTestCase{
        "Equatorial_Modern", -73.0670209146797, 289.84443213466, "2031-10-13T00:15:55", 2463152.511053241,
        159.4272905, 4.2167361, 205.7468, -19.5860, 21.0132760718, 10.384790039,
        "Near equator, modern era"
    };
    
    testCases << CoordinateTestCase{
        "High_North_Modern", 81.9121742352522, 248.896246076003, "1955-08-26T17:03:14", 2435346.210578704,
        152.0646154, 12.7006088, 145.5122, 19.4223, 7.9250404113, -2.212600613,
        "High northern latitude"
    };
    
    testCases << CoordinateTestCase{
        "Mid_Latitude_Medieval", 68.4688217628369, 60.2632226629549, "1658-07-25T05:53:38", 2326837.745578703,
        116.6004466, 21.5438150, 149.0091, 40.6459, 6.1176176006, -1.655745503,
        "Medieval observation, mid-latitude"
    };
    
    testCases << CoordinateTestCase{
        "Southern_Extreme", -33.2644558278182, 0.545545956075841, "1905-03-24T20:19:06", 2416929.346597222,
        330.6995282, -13.2817491, 209.9314, -38.2334, 8.4634774969, 10.416842280,
        "Southern hemisphere extreme"
    };
    
    testCases << CoordinateTestCase{
        "Moderate_North_Future", 27.0113601589112, 328.857064625173, "2104-06-15T08:30:28", 2489695.854490741,
        38.0922383, 14.5011793, 101.7903, 52.6233, 0.0151652940, -2.524317262,
        "Future moderate northern latitude"
    };
    
    testCases << CoordinateTestCase{
        "Tropical_Modern", 19.0977681214935, 325.206309923673, "1961-10-05T08:10:23", 2437577.840543982,
        218.1331775, -13.1749970, 94.7853, -29.0830, 6.7710808583, -7.771130976,
        "Tropical observation, modern era"
    };
    
    testCases << CoordinateTestCase{
        "Ancient_High_South", -75.6965454927354, 163.712986244958, "1125-06-14T01:12:44", 2132128.550509259,
        88.3563845, 23.5451200, 358.1011, -9.2500, 6.0267231267, 0.136297493,
        "Ancient era, high southern latitude"
    };
    
    testCases << CoordinateTestCase{
        "Medieval_Mid_North", 26.7967632675912, 217.814997034872, "1111-10-25T15:55:01", 2127148.163206019,
        215.7391762, -14.2874437, 107.4392, 2.5997, 9.0780394226, -5.304572322,
        "Medieval era, mid-northern latitude"
    };
    
    testCases << CoordinateTestCase{
        "Ancient_Europe", -71.5030354802627, 322.670493082992, "1063-03-19T18:38:43", 2109396.276886574,
        4.0935808, 1.7854219, 298.4607, 7.1989, 4.3240156388, 4.051110250,
        "Ancient European coordinates"
    };
    
    testCases << CoordinateTestCase{
        "Far_Future_High_North", 40.5685793309708, 179.240599988614, "2633-04-22T05:33:46", 2682854.731782407,
        29.9967331, 12.1818592, 274.8692, 13.1751, 7.5320745447, 5.532292336,
        "Far future, high northern observation"
    };
    
    testCases << CoordinateTestCase{
        "Extreme_South_Future", -89.163023764716, 66.9975288348826, "2491-09-06T15:58:33", 2631128.165659722,
        269.8757219, -23.5979494, 336.9535, 24.3685, 19.5184614332, 1.526746639,
        "Extreme southern latitude, future"
    };
    
    testCases << CoordinateTestCase{
        "High_South_Future", -86.173169090456, 61.5160879461091, "2335-09-07T21:11:15", 2574150.382812500,
        321.2821694, -15.9125458, 314.9858, 18.6381, 0.3640299479, 2.945218653,
        "High southern latitude, mid-future"
    };
    
    testCases << CoordinateTestCase{
        "Pacific_Modern", 30.7931630605358, 243.850042925984, "2145-10-09T00:10:26", 2504786.507245370,
        154.2404887, 11.4925881, 290.4320, -10.7385, 17.6180334349, 7.335334187,
        "Pacific region, modern era"
    };
    
    testCases << CoordinateTestCase{
        "Northern_Historical", 74.5187820270132, 26.0066337748148, "1889-10-17T03:13:28", 2411292.634351852,
        166.5470123, 7.2023871, 111.0623, 13.1225, 6.6794544915, -4.423679661,
        "Northern region, historical period"
    };
    
    testCases << CoordinateTestCase{
        "Southern_Modern", -79.1099454095946, 55.9116751156616, "2000-06-07T12:38:35", 2451703.026793981,
        308.6635246, -18.4762477, 167.4575, 7.8562, 9.4515532358, -11.126015073,
        "Southern region, Y2K era"
    };
    
    testCases << CoordinateTestCase{
        "Equatorial_Future", -7.15816700199404, 105.760433303697, "2440-05-15T05:27:25", 2612386.727372685,
        84.1673057, 22.9939312, 35.9368, 52.4645, 4.0873715721, -1.523782144,
        "Near equator, future observation"
    };
    
    testCases << CoordinateTestCase{
        "Arctic_Far_Future", 82.3817832640091, 46.2141135463267, "3682-04-11T01:06:32", 3065983.546203704,
        268.1839880, -29.2215956, 175.2772, -21.6309, 17.5435237090, -0.335408821,
        "Arctic region, far future"
    };
    
    return testCases;
}

void StellinaProcessor::runCoordinateTestSuite() {
    debugLog("=== COMPREHENSIVE COORDINATE TEST SUITE ===");
    debugLog("Running NASA JPL Horizons validation tests...");
    debugLog("");
    
    QList<CoordinateTestCase> testCases = getCoordinateTestCases();
    
    int passCount = 0;
    int failCount = 0;
    double totalAltError = 0.0;
    double totalAzError = 0.0;
    double maxAltError = 0.0;
    double maxAzError = 0.0;
    
    for (const CoordinateTestCase &testCase : testCases) {
        debugLog(QString("Test: %1").arg(testCase.description));
        debugLog(QString("Location: %1°N, %2°E, Time: %3")
                    .arg(testCase.lat, 0, 'f', 4)
                    .arg(testCase.lon, 0, 'f', 4)
                    .arg(testCase.date));
        
        // Save current settings
        QString savedLocation = m_observerLocation;
        m_observerLocation = QString("%1,%2").arg(testCase.lat, 0, 'f', 6).arg(testCase.lon, 0, 'f', 6);
        
        // Run INVERSE conversion: RA/Dec -> Alt/Az to test our algorithm
        // We have the expected Alt/Az from JPL Horizons, so we convert the given RA/Dec
        // and compare with the expected Alt/Az
        
        double calculatedRA, calculatedDec;
        bool success = convertAltAzToRaDec(testCase.expectedAlt, testCase.expectedAz, testCase.date, calculatedRA, calculatedDec);
        
        if (success) {
            // Calculate errors
            double raError = qAbs(calculatedRA - testCase.raOfDate);
            if (raError > 180.0) raError = 360.0 - raError;  // Handle wrap-around
            double decError = qAbs(calculatedDec - testCase.decOfDate);
            
            // Convert to more meaningful units
            double raErrorArcmin = raError * 60.0;
            double decErrorArcmin = decError * 60.0;
            
            bool testPassed = (raError < 1.0 && decError < 1.0);  // 1 degree tolerance
            
            if (testPassed) {
                passCount++;
                debugLog(QString("✓ PASS: RA error = %1° (%2'), Dec error = %3° (%4')")
                            .arg(raError, 0, 'f', 4)
                            .arg(raErrorArcmin, 0, 'f', 1)
                            .arg(decError, 0, 'f', 4)
                            .arg(decErrorArcmin, 0, 'f', 1));
            } else {
                failCount++;
                debugLog(QString("✗ FAIL: RA error = %1° (%2'), Dec error = %3° (%4')")
                            .arg(raError, 0, 'f', 4)
                            .arg(raErrorArcmin, 0, 'f', 1)
                            .arg(decError, 0, 'f', 4)
                            .arg(decErrorArcmin, 0, 'f', 1));
            }
            
            debugLog(QString("Expected: RA=%1°, Dec=%2°")
                        .arg(testCase.raOfDate, 0, 'f', 4)
                        .arg(testCase.decOfDate, 0, 'f', 4));
            debugLog(QString("Got:      RA=%1°, Dec=%2°")
                        .arg(calculatedRA, 0, 'f', 4)
                        .arg(calculatedDec, 0, 'f', 4));
            
            totalAltError += decError;
            totalAzError += raError;
            maxAltError = qMax(maxAltError, decError);
            maxAzError = qMax(maxAzError, raError);
            
        } else {
            failCount++;
            debugLog("✗ FAIL: Conversion failed");
        }
        
        // Restore settings
        m_observerLocation = savedLocation;
        debugLog("");
    }
    
    // Summary statistics
    debugLog("=== TEST SUITE SUMMARY ===");
    debugLog(QString("Total tests: %1").arg(testCases.size()));
    debugLog(QString("Passed: %1").arg(passCount));
    debugLog(QString("Failed: %1").arg(failCount));
    debugLog(QString("Success rate: %1%").arg(100.0 * passCount / testCases.size(), 0, 'f', 1));
    
    if (passCount > 0) {
        debugLog(QString("Average RA error: %1° (%2 arcmin)")
                    .arg(totalAzError / testCases.size(), 0, 'f', 4)
                    .arg(totalAzError * 60.0 / testCases.size(), 0, 'f', 1));
        debugLog(QString("Average Dec error: %1° (%2 arcmin)")
                    .arg(totalAltError / testCases.size(), 0, 'f', 4)
                    .arg(totalAltError * 60.0 / testCases.size(), 0, 'f', 1));
        debugLog(QString("Maximum RA error: %1° (%2 arcmin)")
                    .arg(maxAzError, 0, 'f', 4)
                    .arg(maxAzError * 60.0, 0, 'f', 1));
        debugLog(QString("Maximum Dec error: %1° (%2 arcmin)")
                    .arg(maxAltError, 0, 'f', 4)
                    .arg(maxAltError * 60.0, 0, 'f', 1));
    }
    
    if (passCount < testCases.size()) {
        debugLog("");
        debugLog("DIAGNOSTIC RECOMMENDATIONS:");
        if (totalAzError / testCases.size() > 5.0) {
            debugLog("• Large RA errors suggest azimuth convention problems");
            debugLog("• Check if azimuth is measured from North vs South");
            debugLog("• Verify East vs West azimuth direction");
        }
        if (totalAltError / testCases.size() > 5.0) {
            debugLog("• Large Dec errors suggest altitude calculation issues");
            debugLog("• Check trigonometric formula implementation");
            debugLog("• Verify coordinate system assumptions");
        }
        if (failCount > passCount) {
            debugLog("• High failure rate suggests fundamental algorithm problem");
            debugLog("• Review spherical astronomy coordinate transformation");
            debugLog("• Check local sidereal time calculation");
        }
    } else {
        debugLog("");
        debugLog("🎉 ALL TESTS PASSED! Your coordinate conversion is working correctly.");
    }
    
    debugLog("=======================================");
    debugLog("");
}

void StellinaProcessor::runRandomTestSubset(int numTests) {
    QRandomGenerator random;
    debugLog(QString("=== RANDOM TEST SUBSET (%1 tests) ===").arg(numTests));
    
    QList<CoordinateTestCase> allTests = getCoordinateTestCases();
        
    for (const CoordinateTestCase &testCase : allTests) {
        if (random.generateDouble() < 0.1) testSingleCoordinate(testCase);
    }
    
    debugLog("Random subset test complete.");
    debugLog("");
}

void StellinaProcessor::testSingleCoordinate(const CoordinateTestCase &testCase) {
    // Save current settings
    QString savedLocation = m_observerLocation;
    m_observerLocation = QString("%1,%2").arg(testCase.lat, 0, 'f', 6).arg(testCase.lon, 0, 'f', 6);
    
    debugLog(QString("Testing: %1").arg(testCase.description));
    
    double calculatedRA, calculatedDec;
    bool success = convertAltAzToRaDec(testCase.expectedAlt, testCase.expectedAz, testCase.date, calculatedRA, calculatedDec);
    
    if (success) {
        double raError = qAbs(calculatedRA - testCase.raOfDate);
        if (raError > 180.0) raError = 360.0 - raError;
        double decError = qAbs(calculatedDec - testCase.decOfDate);
        
        QString status = (raError < 1.0 && decError < 1.0) ? "✓ PASS" : "✗ FAIL";
        debugLog(QString("%1: RA error=%2°, Dec error=%3°")
                    .arg(status)
                    .arg(raError, 0, 'f', 3)
                    .arg(decError, 0, 'f', 3));
    } else {
        debugLog("✗ FAIL: Conversion failed");
    }
    
    // Restore settings
    m_observerLocation = savedLocation;
}

void StellinaProcessor::runAccuracyAnalysis() {
    debugLog("=== ACCURACY ANALYSIS ===");
    debugLog("Testing coordinate conversion accuracy across different sky regions...");
    debugLog("");
    
    QList<CoordinateTestCase> testCases = getCoordinateTestCases();
    
    // Group by different criteria
    QMap<QString, QList<double>> errorsByRegion;
    QMap<QString, QList<double>> errorsByEra;
    QMap<QString, QList<double>> errorsByAltitude;
    
    for (const CoordinateTestCase &testCase : testCases) {
        QString savedLocation = m_observerLocation;
        m_observerLocation = QString("%1,%2").arg(testCase.lat, 0, 'f', 6).arg(testCase.lon, 0, 'f', 6);
        
        double calculatedRA, calculatedDec;
        bool success = convertAltAzToRaDec(testCase.expectedAlt, testCase.expectedAz, testCase.date, calculatedRA, calculatedDec);
        
        if (success) {
            double raError = qAbs(calculatedRA - testCase.raOfDate);
            if (raError > 180.0) raError = 360.0 - raError;
            double decError = qAbs(calculatedDec - testCase.decOfDate);
            double totalError = sqrt(raError*raError + decError*decError);
            
            // Categorize by latitude region
            QString region;
            if (qAbs(testCase.lat) > 66.5) region = "Polar";
            else if (qAbs(testCase.lat) > 45.0) region = "High Latitude";
            else if (qAbs(testCase.lat) > 23.5) region = "Mid Latitude";
            else region = "Tropical";
            
            errorsByRegion[region].append(totalError);
            
            // Categorize by era
            QString era;
            if (testCase.julianDay < 2000000) era = "Ancient";
            else if (testCase.julianDay < 2400000) era = "Medieval";
            else if (testCase.julianDay < 2500000) era = "Modern";
            else era = "Future";
            
            errorsByEra[era].append(totalError);
            
            // Categorize by altitude
            QString altCategory;
            if (qAbs(testCase.expectedAlt) > 60) altCategory = "High Altitude";
            else if (qAbs(testCase.expectedAlt) > 30) altCategory = "Mid Altitude";
            else altCategory = "Low Altitude";
            
            errorsByAltitude[altCategory].append(totalError);
        }
        
        m_observerLocation = savedLocation;
    }
    
    // Report results by category
    auto reportCategory = [this](const QString &categoryName, const QMap<QString, QList<double>> &errorMap) {
        debugLog(QString("--- %1 Analysis ---").arg(categoryName));
        for (auto it = errorMap.begin(); it != errorMap.end(); ++it) {
            if (!it.value().isEmpty()) {
                double avgError = 0.0;
                for (double error : it.value()) avgError += error;
                avgError /= it.value().size();
                
                debugLog(QString("%1: avg error = %2° (%3 tests)")
                            .arg(it.key(), -15)
                            .arg(avgError, 0, 'f', 3)
                            .arg(it.value().size()));
            }
        }
        debugLog("");
    };
    
    reportCategory("Latitude Region", errorsByRegion);
    reportCategory("Historical Era", errorsByEra);
    reportCategory("Altitude Range", errorsByAltitude);
    
    debugLog("Analysis complete.");
    debugLog("");
}
