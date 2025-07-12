#include "StellinaProcessor.h"
#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonArray>
#include <QDateTime>
#include <QSettings>
#include <QSplitter>
#include <QMenuBar>
#include <QStatusBar>
#include <QHeaderView>
#include <QElapsedTimer>
#include <QProcess>
#include <QStandardItemModel>
#include <QTableWidgetItem>
#include <cmath>

StellinaProcessor::StellinaProcessor(QWidget *parent)
    : QMainWindow(parent)
    , m_sirilClient(new SirilClient(this))
    , m_processingTimer(new QTimer(this))
    , m_processing(false)
    , m_processingMode(MODE_BASIC_PLATESOLVE)
    , m_currentImageIndex(0)
    , m_processedCount(0)
    , m_errorCount(0)
    , m_skippedCount(0)
    , m_darkCalibratedCount(0)
    , m_registeredCount(0)
    , m_processingStartTime(0)
    , m_currentStage(STAGE_DARK_CALIBRATION)
    , m_qualityFilter(true)
    , m_debugMode(false)
    , m_focalLength(400.0)
    , m_pixelSize(2.40)
    , m_observerLocation("London")
    , m_autoMatchDarks(true)
    , m_temperatureTolerance(5)
    , m_exposureTolerance(10)
    , m_createMasterDark(true)
    , m_sequenceName("stellina_sequence")
{
    setWindowTitle("Enhanced Stellina Processor for Siril - v2.0");
    setMinimumSize(1000, 800);
    
    // Initialize stacking parameters
    m_stackingParams.method = "median";
    m_stackingParams.rejection = "sigma";
    m_stackingParams.rejectionLow = 3.0;
    m_stackingParams.rejectionHigh = 3.0;
    m_stackingParams.normalizeImages = true;
    m_stackingParams.applyDrizzle = false;
    m_stackingParams.drizzleScale = 1.5;
    m_stackingParams.outputFormat = "fits";
    
    // Setup timer
    m_processingTimer->setSingleShot(false);
    m_processingTimer->setInterval(2000); // Process one image every 2 seconds
    
    setupUI();
    setupMenu();
    connectSignals();
    updateUI();
    
    // Load settings
    QSettings settings;
    m_sourceDirectory = settings.value("sourceDirectory").toString();
    m_outputDirectory = settings.value("outputDirectory").toString();
    m_darkDirectory = settings.value("darkDirectory").toString();
    m_qualityFilter = settings.value("qualityFilter", true).toBool();
    m_focalLength = settings.value("focalLength", 400.0).toDouble();
    m_pixelSize = settings.value("pixelSize", 2.40).toDouble();
    m_observerLocation = settings.value("observerLocation", "London").toString();
    m_processingMode = static_cast<ProcessingMode>(settings.value("processingMode", 0).toInt());
    
    // Update UI with loaded settings
    m_sourceDirectoryEdit->setText(m_sourceDirectory);
    m_outputDirectoryEdit->setText(m_outputDirectory);
    m_darkDirectoryEdit->setText(m_darkDirectory);
    m_qualityFilterCheck->setChecked(m_qualityFilter);
    m_focalLengthSpin->setValue(m_focalLength);
    m_pixelSizeSpin->setValue(m_pixelSize);
    m_observerLocationEdit->setText(m_observerLocation);
    m_processingModeCombo->setCurrentIndex(m_processingMode);
    
    logMessage("Enhanced Stellina Processor started. Connect to Siril and select processing mode.", "blue");
    
    // Scan for dark frames if directory is set
    if (!m_darkDirectory.isEmpty()) {
        scanDarkFrames();
    }
}

StellinaProcessor::~StellinaProcessor() {
    // Save settings
    QSettings settings;
    settings.setValue("sourceDirectory", m_sourceDirectory);
    settings.setValue("outputDirectory", m_outputDirectory);
    settings.setValue("darkDirectory", m_darkDirectory);
    settings.setValue("qualityFilter", m_qualityFilter);
    settings.setValue("focalLength", m_focalLength);
    settings.setValue("pixelSize", m_pixelSize);
    settings.setValue("observerLocation", m_observerLocation);
    settings.setValue("processingMode", static_cast<int>(m_processingMode));
}

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
    
    // Connection group
    m_connectionGroup = new QGroupBox("Siril Connection");
    QHBoxLayout *connectionLayout = new QHBoxLayout(m_connectionGroup);
    
    m_testConnectionButton = new QPushButton("Test Connection");
    m_connectionStatus = new QLabel("Not connected");
    m_connectionStatus->setStyleSheet("color: red; font-weight: bold;");
    
    connectionLayout->addWidget(m_testConnectionButton);
    connectionLayout->addWidget(m_connectionStatus);
    connectionLayout->addStretch();
    
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
    
    // Input group
    m_inputGroup = new QGroupBox("Input/Output Directories");
    QGridLayout *inputLayout = new QGridLayout(m_inputGroup);
    
    inputLayout->addWidget(new QLabel("Source Directory:"), 0, 0);
    m_sourceDirectoryEdit = new QLineEdit;
    m_selectSourceButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_sourceDirectoryEdit, 0, 1);
    inputLayout->addWidget(m_selectSourceButton, 0, 2);
    
    inputLayout->addWidget(new QLabel("Output Directory:"), 1, 0);
    m_outputDirectoryEdit = new QLineEdit;
    m_selectOutputButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_outputDirectoryEdit, 1, 1);
    inputLayout->addWidget(m_selectOutputButton, 1, 2);
    
    inputLayout->addWidget(new QLabel("Dark Frames Directory:"), 2, 0);
    m_darkDirectoryEdit = new QLineEdit;
    m_selectDarkButton = new QPushButton("Browse...");
    inputLayout->addWidget(m_darkDirectoryEdit, 2, 1);
    inputLayout->addWidget(m_selectDarkButton, 2, 2);
    
    QHBoxLayout *darkInfoLayout = new QHBoxLayout;
    m_darkFramesCount = new QLabel("No dark frames loaded");
    m_refreshDarkButton = new QPushButton("Refresh");
    darkInfoLayout->addWidget(m_darkFramesCount);
    darkInfoLayout->addWidget(m_refreshDarkButton);
    darkInfoLayout->addStretch();
    inputLayout->addLayout(darkInfoLayout, 3, 1, 1, 2);
    
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
    layout->addWidget(m_connectionGroup);
    layout->addWidget(m_modeGroup);
    layout->addWidget(m_inputGroup);
    layout->addWidget(m_basicOptionsGroup);
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
    
    m_createMasterDarkCheck = new QCheckBox("Create master dark frames from multiple darks");
    m_createMasterDarkCheck->setChecked(true);
    darkOptionsLayout->addWidget(m_createMasterDarkCheck, 3, 0, 1, 2);
    
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
    toolsMenu->addAction("Test &Connection", this, &StellinaProcessor::onTestConnection);
    toolsMenu->addAction("&Refresh Dark Frames", this, &StellinaProcessor::onRefreshDarkFrames);
    toolsMenu->addSeparator();
    toolsMenu->addAction("&Clear Log", this, &StellinaProcessor::onClearLog);
    
    // Help menu
    QMenu *helpMenu = menuBar->addMenu("&Help");
    helpMenu->addAction("&About", [this]() {
        QMessageBox::about(this, "About Enhanced Stellina Processor",
                          "Enhanced Stellina Processor for Siril v2.0\n\n"
                          "A Qt application for advanced processing of Stellina telescope images\n"
                          "using Siril's capabilities.\n\n"
                          "New Features:\n"
                          "• Dark frame calibration with automatic matching\n"
                          "• Astrometric registration and stacking\n"
                          "• Full processing pipeline\n"
                          "• Advanced rejection algorithms\n"
                          "• Drizzle enhancement support\n"
                          "• Comprehensive processing reports");
    });
}

void StellinaProcessor::connectSignals() {
    // UI signals
    connect(m_selectSourceButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectSourceDirectory);
    connect(m_selectOutputButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectOutputDirectory);
    connect(m_selectDarkButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectDarkDirectory);
    connect(m_startButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStartProcessing);
    connect(m_stopButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStopProcessing);
    connect(m_testConnectionButton, &QPushButton::clicked,
            this, &StellinaProcessor::onTestConnection);
    connect(m_clearLogButton, &QPushButton::clicked,
            this, &StellinaProcessor::onClearLog);
    connect(m_refreshDarkButton, &QPushButton::clicked,
            this, &StellinaProcessor::onRefreshDarkFrames);
    connect(m_previewStackButton, &QPushButton::clicked,
            this, &StellinaProcessor::onPreviewStacking);
    
    // Processing mode combo
    connect(m_processingModeCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &StellinaProcessor::onProcessingModeChanged);
    
    // Siril client signals
    connect(m_sirilClient, &SirilClient::connected,
            this, &StellinaProcessor::onSirilConnected);
    connect(m_sirilClient, &SirilClient::disconnected,
            this, &StellinaProcessor::onSirilDisconnected);
    connect(m_sirilClient, &SirilClient::commandExecuted,
            this, &StellinaProcessor::onSirilCommandExecuted);
    connect(m_sirilClient, &SirilClient::errorOccurred,
            this, &StellinaProcessor::onSirilError);
    
    // Processing timer
    connect(m_processingTimer, &QTimer::timeout,
            this, &StellinaProcessor::onProcessingTimer);
    
    // Settings updates
    connect(m_qualityFilterCheck, &QCheckBox::toggled, [this](bool checked) {
        m_qualityFilter = checked;
    });
    connect(m_debugModeCheck, &QCheckBox::toggled, [this](bool checked) {
        m_debugMode = checked;
    });
    connect(m_focalLengthSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_focalLength = value;
    });
    connect(m_pixelSizeSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_pixelSize = value;
    });
    connect(m_observerLocationEdit, &QLineEdit::textChanged, [this](const QString &text) {
        m_observerLocation = text;
    });
    
    // Dark calibration settings
    connect(m_autoMatchDarksCheck, &QCheckBox::toggled, [this](bool checked) {
        m_autoMatchDarks = checked;
    });
    connect(m_temperatureToleranceSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [this](int value) {
        m_temperatureTolerance = value;
    });
    connect(m_exposureToleranceSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [this](int value) {
        m_exposureTolerance = value;
    });
    connect(m_createMasterDarkCheck, &QCheckBox::toggled, [this](bool checked) {
        m_createMasterDark = checked;
    });
    
    // Stacking settings
    connect(m_stackingMethodCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.method = text;
        onStackingParametersChanged();
    });
    connect(m_rejectionMethodCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.rejection = text;
        onStackingParametersChanged();
    });
    connect(m_rejectionLowSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.rejectionLow = value;
        onStackingParametersChanged();
    });
    connect(m_rejectionHighSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.rejectionHigh = value;
        onStackingParametersChanged();
    });
    connect(m_normalizeCheck, &QCheckBox::toggled, [this](bool checked) {
        m_stackingParams.normalizeImages = checked;
        onStackingParametersChanged();
    });
    connect(m_drizzleCheck, &QCheckBox::toggled, [this](bool checked) {
        m_stackingParams.applyDrizzle = checked;
        m_drizzleScaleSpin->setEnabled(checked);
        onStackingParametersChanged();
    });
    connect(m_drizzleScaleSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.drizzleScale = value;
        onStackingParametersChanged();
    });
    connect(m_outputFormatCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.outputFormat = text;
        onStackingParametersChanged();
    });
}

void StellinaProcessor::onProcessingModeChanged() {
    m_processingMode = static_cast<ProcessingMode>(m_processingModeCombo->currentData().toInt());
    
    QString description;
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        description = "Basic plate solving with coordinate annotation - original functionality";
        break;
    case MODE_DARK_CALIBRATION:
        description = "Dark frame calibration only - removes thermal noise from light frames";
        break;
    case MODE_ASTROMETRIC_STACKING:
        description = "Astrometric registration and stacking - aligns and combines images";
        break;
    case MODE_FULL_PIPELINE:
        description = "Complete processing pipeline - dark calibration, plate solving, registration, and stacking";
        break;
    }
    
    m_modeDescription->setText(description);
    updateUI();
}

void StellinaProcessor::onSelectDarkDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Dark Frames Directory", m_darkDirectory);
    if (!dir.isEmpty()) {
        m_darkDirectory = dir;
        m_darkDirectoryEdit->setText(dir);
        scanDarkFrames();
        updateUI();
    }
}

void StellinaProcessor::onRefreshDarkFrames() {
    if (!m_darkDirectory.isEmpty()) {
        scanDarkFrames();
    }
}

void StellinaProcessor::onPreviewStacking() {
    if (m_imagesToProcess.isEmpty()) {
        if (!findStellinaImages()) {
            QMessageBox::warning(this, "No Images", "No valid Stellina images found for preview.");
            return;
        }
    }
    
    if (m_imagesToProcess.size() < 3) {
        QMessageBox::warning(this, "Insufficient Images", "Need at least 3 images for preview stacking.");
        return;
    }
    
    if (!m_sirilClient->isConnected()) {
        QMessageBox::warning(this, "Connection Error", "Please connect to Siril first.");
        return;
    }
    
    // Take first 5 images for preview
    QStringList previewImages = m_imagesToProcess.mid(0, qMin(5, m_imagesToProcess.size()));
    
    logMessage(QString("Creating preview stack with %1 images...").arg(previewImages.size()), "blue");
    
    // Create a temporary sequence for preview
    QString previewSequence = "preview_stack";
    if (createSequence(previewImages, previewSequence)) {
        if (performGlobalRegistration(previewSequence)) {
            if (performStacking(previewSequence, m_stackingParams)) {
                logMessage("Preview stack created successfully! Check output directory.", "green");
            } else {
                logMessage("Preview stacking failed.", "red");
            }
        } else {
            logMessage("Preview registration failed.", "red");
        }
    } else {
        logMessage("Failed to create preview sequence.", "red");
    }
}

void StellinaProcessor::onStackingParametersChanged() {
    // This function is called when stacking parameters change
    // Update the internal stacking parameters structure
    m_stackingParams.method = m_stackingMethodCombo->currentText();
    m_stackingParams.rejection = m_rejectionMethodCombo->currentText();
    m_stackingParams.rejectionLow = m_rejectionLowSpin->value();
    m_stackingParams.rejectionHigh = m_rejectionHighSpin->value();
    m_stackingParams.normalizeImages = m_normalizeCheck->isChecked();
    m_stackingParams.applyDrizzle = m_drizzleCheck->isChecked();
    m_stackingParams.drizzleScale = m_drizzleScaleSpin->value();
    m_stackingParams.outputFormat = m_outputFormatCombo->currentText();
    
    if (m_debugMode) {
        logMessage(QString("Stacking parameters updated: method=%1, rejection=%2")
                      .arg(m_stackingParams.method)
                      .arg(m_stackingParams.rejection), "gray");
    }
}

void StellinaProcessor::scanDarkFrames() {
    m_darkFrames.clear();
    
    if (m_darkDirectory.isEmpty()) {
        m_darkFramesCount->setText("No dark frames directory selected");
        return;
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        m_darkFramesCount->setText("Dark frames directory does not exist");
        return;
    }
    
    QStringList darkFiles = darkDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    
    logMessage(QString("Scanning %1 potential dark frames...").arg(darkFiles.size()), "blue");
    
    // Group dark frames by exposure, temperature, and binning
    QMap<QString, QStringList> darkGroups;
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperature = extractTemperature(fullPath);
        QString binning = extractBinning(fullPath);
        
        if (exposure > 0) { // Valid dark frame
            QString key = QString("%1_%2_%3").arg(exposure).arg(temperature).arg(binning);
            darkGroups[key].append(fullPath);
        }
    }
    
    // Create DarkFrame entries
    for (auto it = darkGroups.begin(); it != darkGroups.end(); ++it) {
        QStringList parts = it.key().split('_');
        if (parts.size() >= 3) {
            DarkFrame dark;
            dark.filepath = it.value().first(); // Use first file as representative
            dark.exposure = parts[0].toInt();
            dark.temperature = parts[1].toInt();
            dark.binning = parts[2];
            
            m_darkFrames.append(dark);
        }
    }
    
    // Update UI
    m_darkFramesCount->setText(QString("%1 dark frame groups found").arg(m_darkFrames.size()));
    
    // Update dark frames table
    m_darkFramesTable->setRowCount(m_darkFrames.size());
    for (int i = 0; i < m_darkFrames.size(); ++i) {
        const DarkFrame &dark = m_darkFrames[i];
        
        m_darkFramesTable->setItem(i, 0, new QTableWidgetItem(QFileInfo(dark.filepath).fileName()));
        m_darkFramesTable->setItem(i, 1, new QTableWidgetItem(QString::number(dark.exposure)));
        m_darkFramesTable->setItem(i, 2, new QTableWidgetItem(QString::number(dark.temperature)));
        m_darkFramesTable->setItem(i, 3, new QTableWidgetItem(dark.binning));
        
        // Count how many dark frames in this group
        QString key = QString("%1_%2_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning);
        int count = darkGroups[key].size();
        m_darkFramesTable->setItem(i, 4, new QTableWidgetItem(QString::number(count)));
    }
    
    logMessage(QString("Dark frame scan complete: %1 groups found").arg(m_darkFrames.size()), "green");
}

// Dark calibration functions
bool StellinaProcessor::findMatchingDarkFrame(const QString &lightFrame, DarkFrame &darkFrame) {
    if (!m_autoMatchDarks || m_darkFrames.isEmpty()) {
        return false;
    }
    
    int lightExposure = extractExposureTime(lightFrame);
    int lightTemperature = extractTemperature(lightFrame);
    QString lightBinning = extractBinning(lightFrame);
    
    if (lightExposure <= 0) {
        return false;
    }
    
    // Find best matching dark frame
    for (const DarkFrame &dark : m_darkFrames) {
        bool exposureMatch = qAbs(dark.exposure - lightExposure) <= (lightExposure * m_exposureTolerance / 100);
        bool temperatureMatch = qAbs(dark.temperature - lightTemperature) <= m_temperatureTolerance;
        bool binningMatch = (dark.binning == lightBinning);
        
        if (exposureMatch && temperatureMatch && binningMatch) {
            darkFrame = dark;
            return true;
        }
    }
    
    return false;
}

bool StellinaProcessor::applyDarkCalibration(const QString &lightFrame, const QString &darkFrame, const QString &outputFrame) {
    // Load light frame
    if (!m_sirilClient->loadImage(lightFrame)) {
        return false;
    }
    
    // Apply dark subtraction
    QString command = QString("calibrate %1").arg(QFileInfo(darkFrame).baseName());
    if (!m_sirilClient->sendSirilCommand(command)) {
        return false;
    }
    
    // Save calibrated frame
    if (!m_sirilClient->saveImage(outputFrame)) {
        return false;
    }
    
    m_sirilClient->closeImage();
    return true;
}

// Astrometric stacking functions
bool StellinaProcessor::createSequence(const QStringList &imageList, const QString &sequenceName) {
    if (imageList.isEmpty()) {
        return false;
    }
    
    // Change to output directory
    if (!m_sirilClient->changeDirectory(m_outputDirectory)) {
        logMessage("Warning: Could not change to output directory", "orange");
    }
    
    // Copy images to a temporary location with sequential naming
    QString seqDir = QDir(m_outputDirectory).absoluteFilePath(sequenceName);
    QDir().mkpath(seqDir);
    
    logMessage(QString("Creating sequence with %1 images...").arg(imageList.size()), "blue");
    
    for (int i = 0; i < imageList.size(); ++i) {
        QString srcFile = imageList[i];
        QString dstFile = QDir(seqDir).absoluteFilePath(QString("%1_%2.fits").arg(sequenceName).arg(i + 1, 4, 10, QChar('0')));
        
        if (!QFile::copy(srcFile, dstFile)) {
            logMessage(QString("Failed to copy %1 to sequence directory").arg(QFileInfo(srcFile).fileName()), "red");
            return false;
        }
    }
    
    // Create Siril sequence
    QString command = QString("cd \"%1\"").arg(seqDir);
    if (!m_sirilClient->sendSirilCommand(command)) {
        return false;
    }
    
    command = QString("convert %1 -out=%1").arg(sequenceName);
    return m_sirilClient->sendSirilCommand(command);
}

bool StellinaProcessor::performGlobalRegistration(const QString &sequenceName) {
    logMessage("Performing global registration...", "blue");
    m_registrationStatusLabel->setText("Registration: In Progress");
    
    // Load sequence
    QString command = QString("load %1").arg(sequenceName);
    if (!m_sirilClient->sendSirilCommand(command)) {
        m_registrationStatusLabel->setText("Registration: Failed");
        return false;
    }
    
    // Perform global registration
    command = "register pp";
    if (!m_sirilClient->sendSirilCommand(command)) {
        m_registrationStatusLabel->setText("Registration: Failed");
        return false;
    }
    
    m_registrationStatusLabel->setText("Registration: Complete");
    logMessage("Global registration completed", "green");
    return true;
}

bool StellinaProcessor::performStacking(const QString &sequenceName, const StackingParams &params) {
    logMessage(QString("Performing stacking with method: %1").arg(params.method), "blue");
    m_stackingStatusLabel->setText("Stacking: In Progress");
    
    // Build stacking command
    QString command = QString("stack pp_%1 rej %2 %3 %4 -norm=%5")
                         .arg(sequenceName)
                         .arg(params.rejection)
                         .arg(params.rejectionLow)
                         .arg(params.rejectionHigh)
                         .arg(params.normalizeImages ? "additive" : "none");
    
    if (params.method != "sum") {
        command += QString(" -type=%1").arg(params.method);
    }
    
    if (params.applyDrizzle) {
        command += QString(" -drizzle -scale=%1").arg(params.drizzleScale);
    }
    
    if (!m_sirilClient->sendSirilCommand(command)) {
        m_stackingStatusLabel->setText("Stacking: Failed");
        return false;
    }
    
    // Save final stack
    QString outputName = QString("stacked_%1.%2").arg(sequenceName).arg(params.outputFormat);
    if (!m_sirilClient->saveImage(outputName)) {
        m_stackingStatusLabel->setText("Stacking: Failed");
        return false;
    }
    
    m_finalStackedImage = QDir(m_outputDirectory).absoluteFilePath(outputName);
    m_stackingStatusLabel->setText("Stacking: Complete");
    logMessage(QString("Stacking completed: %1").arg(outputName), "green");
    return true;
}

// Enhanced processing pipeline
void StellinaProcessor::startStellinaProcessing() {
    if (!validateProcessingInputs()) {
        return;
    }
    
    logMessage(QString("Starting %1 processing pipeline...")
                  .arg(m_processingModeCombo->currentText()), "blue");
    
    // Find images to process
    if (!findStellinaImages()) {
        logMessage("No valid Stellina images found in source directory.", "red");
        return;
    }
    
    // Initialize processing state
    m_processing = true;
    m_currentImageIndex = 0;
    m_processedCount = 0;
    m_errorCount = 0;
    m_skippedCount = 0;
    m_darkCalibratedCount = 0;
    m_registeredCount = 0;
    m_processingStartTime = QDateTime::currentMSecsSinceEpoch();
    
    // Clear previous results
    m_darkCalibratedFiles.clear();
    m_plateSolvedFiles.clear();
    m_registeredFiles.clear();
    m_finalStackedImage.clear();
    
    // Set up progress tracking
    m_progressBar->setMaximum(m_imagesToProcess.length());
    m_progressBar->setValue(0);
    
    // Determine processing stages based on mode
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        m_currentStage = STAGE_PLATE_SOLVING;
        break;
    case MODE_DARK_CALIBRATION:
        m_currentStage = STAGE_DARK_CALIBRATION;
        break;
    case MODE_ASTROMETRIC_STACKING:
        m_currentStage = STAGE_REGISTRATION;
        break;
    case MODE_FULL_PIPELINE:
        m_currentStage = STAGE_DARK_CALIBRATION;
        break;
    }
    
    updateProcessingStatus();
    updateUI();
    
    logMessage(QString("Found %1 images to process.").arg(m_imagesToProcess.length()), "green");
    
    // Start processing timer
    m_processingTimer->start();
}

void StellinaProcessor::processNextImage() {
    if (m_currentImageIndex >= m_imagesToProcess.length()) {
        // Current stage complete, move to next stage or finish
        switch (m_currentStage) {
        case STAGE_DARK_CALIBRATION:
            if (m_processingMode == MODE_DARK_CALIBRATION) {
                finishProcessing();
                return;
            } else {
                // Move to plate solving stage
                m_currentStage = STAGE_PLATE_SOLVING;
                m_currentImageIndex = 0;
                // Use dark-calibrated files if available
                if (!m_darkCalibratedFiles.isEmpty()) {
                    m_imagesToProcess = m_darkCalibratedFiles;
                }
                updateProcessingStatus();
                return;
            }
            break;
            
        case STAGE_PLATE_SOLVING:
            if (m_processingMode == MODE_BASIC_PLATESOLVE) {
                finishProcessing();
                return;
            } else {
                // Move to registration stage
                m_currentStage = STAGE_REGISTRATION;
                if (performAstrometricStacking()) {
                    m_currentStage = STAGE_COMPLETE;
                    finishProcessing();
                } else {
                    logMessage("Astrometric stacking failed", "red");
                    finishProcessing();
                }
                return;
            }
            break;
            
        case STAGE_REGISTRATION:
            if (performAstrometricStacking()) {
                m_currentStage = STAGE_COMPLETE;
                finishProcessing();
            } else {
                logMessage("Astrometric stacking failed", "red");
                finishProcessing();
            }
            return;
            
        default:
            finishProcessing();
            return;
        }
    }
    
    QString currentFile = m_imagesToProcess[m_currentImageIndex];
    QString baseName = QFileInfo(currentFile).baseName();
    
    logMessage(QString("Processing image %1 of %2: %3 (Stage: %4)")
                  .arg(m_currentImageIndex + 1)
                  .arg(m_imagesToProcess.length())
                  .arg(QFileInfo(currentFile).fileName())
                  .arg(getStageDescription()), "blue");
    
    bool success = false;
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        success = processImageDarkCalibration(currentFile);
        break;
    case STAGE_PLATE_SOLVING:
        success = processImagePlatesolving(currentFile);
        break;
    default:
        success = processImagePlatesolving(currentFile); // fallback
        break;
    }
    
    // Update counters
    if (success) {
        m_processedCount++;
    } else {
        m_errorCount++;
    }
    
    // Update progress
    m_currentImageIndex++;
    m_progressBar->setValue(m_currentImageIndex);
    updateProcessingStatus();
    
    // Update time estimate
    qint64 elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    double avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
    qint64 remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
    
    m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                    .arg(formatProcessingTime(remaining)));
}

bool StellinaProcessor::processImageDarkCalibration(const QString &lightFrame) {
    m_currentTaskLabel->setText("Dark calibration...");
    
    DarkFrame matchingDark;
    if (!findMatchingDarkFrame(lightFrame, matchingDark)) {
        logMessage("No matching dark frame found, skipping dark calibration", "orange");
        m_skippedCount++;
        return true; // Continue processing without dark calibration
    }
    
    QString outputName = QString("dark_calibrated_%1.fits")
                            .arg(QFileInfo(lightFrame).baseName());
    QString outputPath = QDir(m_outputDirectory).absoluteFilePath(outputName);
    
    if (applyDarkCalibration(lightFrame, matchingDark.filepath, outputPath)) {
        m_darkCalibratedFiles.append(outputPath);
        m_darkCalibratedCount++;
        logMessage(QString("Dark calibration successful: %1").arg(outputName), "green");
        return true;
    } else {
        logMessage("Dark calibration failed", "red");
        return false;
    }
}

bool StellinaProcessor::processImagePlatesolving(const QString &fitsPath) {
    m_currentTaskLabel->setText("Plate solving...");
    
    // This is essentially the same as the original plate solving code
    // but adapted for the new pipeline structure
    
    QString baseName = QFileInfo(fitsPath).baseName();
    
    // Find corresponding JSON file
    QDir sourceDir(QFileInfo(fitsPath).dir());
    QString jsonPath;
    
    QStringList jsonCandidates = {
        baseName + ".json",
        baseName + ".JSON",
        baseName + "-stacking.json",
        baseName + "-stacking.JSON",
        QFileInfo(fitsPath).completeBaseName() + ".json",
        QFileInfo(fitsPath).completeBaseName() + ".JSON",
        QFileInfo(fitsPath).completeBaseName() + "-stacking.json",
        QFileInfo(fitsPath).completeBaseName() + "-stacking.JSON"
    };
    
    bool jsonFound = false;
    for (const QString &candidate : jsonCandidates) {
        QString candidatePath = sourceDir.absoluteFilePath(candidate);
        if (QFile::exists(candidatePath)) {
            jsonPath = candidatePath;
            jsonFound = true;
            break;
        }
    }
    
    if (!jsonFound) {
        logMessage(QString("No JSON file found for %1").arg(QFileInfo(fitsPath).fileName()), "red");
        return false;
    }
    
    // Load JSON metadata
    QJsonObject json = loadStellinaJson(jsonPath);
    if (json.isEmpty()) {
        logMessage(QString("Failed to load JSON metadata"), "red");
        return false;
    }
    
    // Extract coordinates
    double alt, az;
    if (!extractCoordinates(json, alt, az)) {
        logMessage("Failed to extract coordinates from JSON", "red");
        return false;
    }
    
    // Convert Alt/Az to RA/Dec
    double ra, dec;
    QString dateObs = json["date_obs"].toString();
    if (!convertAltAzToRaDec(alt, az, dateObs, ra, dec)) {
        logMessage("Failed to convert coordinates", "red");
        return false;
    }
    
    // Load image in Siril
    if (!m_sirilClient->loadImage(fitsPath)) {
        logMessage("Failed to load image in Siril", "red");
        return false;
    }
    
    // Perform plate solving
    bool platesolveSuccess = m_sirilClient->platesolve(ra, dec, m_focalLength, m_pixelSize, true);
    
    if (!platesolveSuccess) {
        logMessage("Plate solving failed (continuing anyway)", "orange");
    }
    
    // Save processed image
    QString outputName = QString("processed_%1").arg(baseName);
    QString outputPath = QDir(m_outputDirectory).absoluteFilePath(outputName);
    
    if (!m_sirilClient->saveImage(outputPath)) {
        logMessage("Failed to save processed image", "red");
        m_sirilClient->closeImage();
        return false;
    }
    
    // Close image
    m_sirilClient->closeImage();
    
    // Add to plate solved files list for potential stacking
    m_plateSolvedFiles.append(outputPath);
    
    return true;
}

bool StellinaProcessor::performAstrometricStacking() {
    if (m_plateSolvedFiles.isEmpty()) {
        logMessage("No plate-solved images available for stacking", "red");
        return false;
    }
    
    if (m_plateSolvedFiles.size() < 3) {
        logMessage("Need at least 3 images for stacking", "red");
        return false;
    }
    
    logMessage(QString("Starting astrometric stacking with %1 images").arg(m_plateSolvedFiles.size()), "blue");
    
    // Create sequence
    m_sequenceName = QString("stellina_stack_%1").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    
    m_subTaskProgressBar->setVisible(true);
    m_subTaskProgressBar->setMaximum(4);
    m_subTaskProgressBar->setValue(0);
    
    // Step 1: Create sequence
    m_currentTaskLabel->setText("Creating image sequence...");
    m_subTaskProgressBar->setValue(1);
    
    if (!createSequence(m_plateSolvedFiles, m_sequenceName)) {
        logMessage("Failed to create image sequence", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 2: Perform registration
    m_currentTaskLabel->setText("Performing astrometric registration...");
    m_subTaskProgressBar->setValue(2);
    
    if (!performGlobalRegistration(m_sequenceName)) {
        logMessage("Failed to perform registration", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 3: Perform stacking
    m_currentTaskLabel->setText("Stacking registered images...");
    m_subTaskProgressBar->setValue(3);
    
    if (!performStacking(m_sequenceName, m_stackingParams)) {
        logMessage("Failed to perform stacking", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 4: Complete
    m_currentTaskLabel->setText("Stacking complete");
    m_subTaskProgressBar->setValue(4);
    m_subTaskProgressBar->setVisible(false);
    
    return true;
}

bool StellinaProcessor::validateProcessingInputs() {
    if (!m_sirilClient->isConnected()) {
        QMessageBox::warning(this, "Connection Error", "Please connect to Siril first.");
        return false;
    }
    
    if (m_sourceDirectory.isEmpty() || m_outputDirectory.isEmpty()) {
        QMessageBox::warning(this, "Directory Error", "Please select both source and output directories.");
        return false;
    }
    
    // Check if dark frames are needed but not available
    if ((m_processingMode == MODE_DARK_CALIBRATION || m_processingMode == MODE_FULL_PIPELINE) 
        && m_autoMatchDarks && m_darkFrames.isEmpty()) {
        int ret = QMessageBox::question(this, "No Dark Frames", 
                                       "Dark calibration is enabled but no dark frames were found. Continue without dark calibration?",
                                       QMessageBox::Yes | QMessageBox::No);
        if (ret == QMessageBox::No) {
            return false;
        }
    }
    
    return true;
}

void StellinaProcessor::updateProcessingStatus() {
    QString stageText = getStageDescription();
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1").arg(
            m_processing ? "In Progress" : "Ready"));
        break;
    case STAGE_PLATE_SOLVING:
        // Keep existing registration/stacking status as ready for now
        break;
    case STAGE_REGISTRATION:
        m_registrationStatusLabel->setText("Registration: In Progress");
        break;
    case STAGE_STACKING:
        m_stackingStatusLabel->setText("Stacking: In Progress");
        break;
    case STAGE_COMPLETE:
        m_darkCalibrationStatusLabel->setText("Dark Calibration: Complete");
        m_registrationStatusLabel->setText("Registration: Complete");
        m_stackingStatusLabel->setText("Stacking: Complete");
        break;
    }
    
    if (m_processing) {
        m_progressLabel->setText(QString("Processing %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                    .arg(m_currentImageIndex + 1)
                                    .arg(m_imagesToProcess.length())
                                    .arg(stageText)
                                    .arg(m_processedCount)
                                    .arg(m_errorCount));
        
        if (m_darkCalibratedCount > 0) {
            m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1 completed").arg(m_darkCalibratedCount));
        }
    }
}

QString StellinaProcessor::getStageDescription() const {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION: return "Dark Calibration";
    case STAGE_PLATE_SOLVING: return "Plate Solving";
    case STAGE_REGISTRATION: return "Registration";
    case STAGE_STACKING: return "Stacking";
    case STAGE_COMPLETE: return "Complete";
    default: return "Unknown";
    }
}

void StellinaProcessor::finishProcessing() {
    m_processing = false;
    m_processingTimer->stop();
    
    qint64 totalTime = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    
    QString completionMessage = QString("Processing complete! Total time: %1\n\n")
                                   .arg(formatProcessingTime(totalTime));
    
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        completionMessage += QString("Plate solved: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_DARK_CALIBRATION:
        completionMessage += QString("Dark calibrated: %1, Skipped: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_skippedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_ASTROMETRIC_STACKING:
        completionMessage += QString("Images processed: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
        
    case MODE_FULL_PIPELINE:
        completionMessage += QString("Dark calibrated: %1, Plate solved: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
    }
    
    logMessage(completionMessage, "green");
    
    // Save processing report
    saveProcessingReport();
    
    updateUI();
    
    // Show completion dialog
    QMessageBox::information(this, "Processing Complete", completionMessage);
}

void StellinaProcessor::saveProcessingReport() {
    QString reportPath = QDir(m_outputDirectory).absoluteFilePath(
        QString("processing_report_%1.txt").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")));
    
    QFile reportFile(reportPath);
    if (!reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return;
    }
    
    QTextStream out(&reportFile);
    
    out << "Enhanced Stellina Processor - Processing Report\n";
    out << "==============================================\n\n";
    out << "Date: " << QDateTime::currentDateTime().toString() << "\n";
    out << "Processing Mode: " << m_processingModeCombo->currentText() << "\n";
    out << "Source Directory: " << m_sourceDirectory << "\n";
    out << "Output Directory: " << m_outputDirectory << "\n";
    out << "Dark Directory: " << m_darkDirectory << "\n\n";
    
    out << "Processing Statistics:\n";
    out << "- Total images: " << m_imagesToProcess.length() << "\n";
    out << "- Successfully processed: " << m_processedCount << "\n";
    out << "- Dark calibrated: " << m_darkCalibratedCount << "\n";
    out << "- Errors: " << m_errorCount << "\n";
    out << "- Skipped: " << m_skippedCount << "\n\n";
    
    if (m_processingMode == MODE_ASTROMETRIC_STACKING || m_processingMode == MODE_FULL_PIPELINE) {
        out << "Stacking Parameters:\n";
        out << "- Method: " << m_stackingParams.method << "\n";
        out << "- Rejection: " << m_stackingParams.rejection << "\n";
        out << "- Rejection thresholds: " << m_stackingParams.rejectionLow << " / " << m_stackingParams.rejectionHigh << "\n";
        out << "- Normalization: " << (m_stackingParams.normalizeImages ? "Yes" : "No") << "\n";
        out << "- Drizzle: " << (m_stackingParams.applyDrizzle ? QString("Yes (scale: %1)").arg(m_stackingParams.drizzleScale) : "No") << "\n";
        out << "- Output format: " << m_stackingParams.outputFormat << "\n\n";
    }
    
    if (!m_finalStackedImage.isEmpty()) {
        out << "Final stacked image: " << m_finalStackedImage << "\n\n";
    }
    
    qint64 totalTime = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    out << "Total processing time: " << formatProcessingTime(totalTime) << "\n";
    
    reportFile.close();
    logMessage(QString("Processing report saved: %1").arg(QFileInfo(reportPath).fileName()), "blue");
}

QString StellinaProcessor::formatProcessingTime(qint64 milliseconds) {
    qint64 seconds = milliseconds / 1000;
    qint64 minutes = seconds / 60;
    qint64 hours = minutes / 60;
    
    seconds %= 60;
    minutes %= 60;
    
    if (hours > 0) {
        return QString("%1h %2m %3s").arg(hours).arg(minutes).arg(seconds);
    } else if (minutes > 0) {
        return QString("%1m %2s").arg(minutes).arg(seconds);
    } else {
        return QString("%1s").arg(seconds);
    }
}

// Utility functions for FITS metadata extraction
int StellinaProcessor::extractExposureTime(const QString &fitsFile) {
    // This would normally use a FITS library like cfitsio
    // For now, return a placeholder value
    // In a real implementation, you'd read the EXPTIME or EXPOSURE header
    Q_UNUSED(fitsFile)
    return 300; // Placeholder: 5 minutes
}

int StellinaProcessor::extractTemperature(const QString &fitsFile) {
    // This would normally read the CCD-TEMP or similar header
    Q_UNUSED(fitsFile)
    return -10; // Placeholder: -10°C
}

QString StellinaProcessor::extractBinning(const QString &fitsFile) {
    // This would normally read the XBINNING/YBINNING headers
    Q_UNUSED(fitsFile)
    return "1x1"; // Placeholder
}

// Rest of the implementation continues with existing functions...
// (onSelectSourceDirectory, onSelectOutputDirectory, etc. remain the same)

void StellinaProcessor::onSelectSourceDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Stellina Images Directory", m_sourceDirectory);
    if (!dir.isEmpty()) {
        m_sourceDirectory = dir;
        m_sourceDirectoryEdit->setText(dir);
        
        QDir testDir(dir);
        logMessage(QString("Selected directory: %1").arg(testDir.absolutePath()), "blue");
        logMessage(QString("Directory exists: %1").arg(testDir.exists() ? "Yes" : "No"), "blue");
        
        QStringList allFiles = testDir.entryList(QDir::Files);
        logMessage(QString("Total files in directory: %1").arg(allFiles.count()), "blue");
        
        updateUI();
    }
}

void StellinaProcessor::onSelectOutputDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Output Directory", m_outputDirectory);
    if (!dir.isEmpty()) {
        m_outputDirectory = dir;
        m_outputDirectoryEdit->setText(dir);
        updateUI();
    }
}

void StellinaProcessor::onStartProcessing() {
    startStellinaProcessing();
}

void StellinaProcessor::onStopProcessing() {
    if (m_processing) {
        m_processingTimer->stop();
        m_processing = false;
        logMessage("Processing stopped by user.", "orange");
        finishProcessing();
    }
}

void StellinaProcessor::onTestConnection() {
    if (m_sirilClient->isConnected()) {
        m_sirilClient->disconnectFromSiril();
    } else {
        logMessage("Attempting to connect to Siril...", "blue");
        if (m_sirilClient->connectToSiril()) {
            logMessage("Successfully connected to Siril!", "green");
        } else {
            logMessage(QString("Failed to connect: %1").arg(m_sirilClient->lastError()), "red");
        }
    }
    updateUI();
}

void StellinaProcessor::onClearLog() {
    m_logTextEdit->clear();
    logMessage("Log cleared.", "gray");
}

void StellinaProcessor::onSirilConnected() {
    logMessage("Connected to Siril successfully!", "green");
    
    QString workDir = m_sirilClient->getWorkingDirectory();
    if (!workDir.isEmpty()) {
        logMessage(QString("Siril working directory: %1").arg(workDir), "blue");
    }
    
    updateUI();
}

void StellinaProcessor::onSirilDisconnected() {
    logMessage("Disconnected from Siril.", "orange");
    updateUI();
}

void StellinaProcessor::onSirilCommandExecuted(const QString &command, bool success) {
    if (m_debugMode) {
        QString status = success ? "SUCCESS" : "FAILED";
        logMessage(QString("Command [%1]: %2").arg(status).arg(command), success ? "gray" : "red");
    }
}

void StellinaProcessor::onSirilError(const QString &error) {
    logMessage(QString("Siril error: %1").arg(error), "red");
}

void StellinaProcessor::onProcessingTimer() {
    processNextImage();
}

void StellinaProcessor::updateUI() {
    bool connected = m_sirilClient->isConnected();
    bool canProcess = connected && !m_sourceDirectory.isEmpty() && !m_outputDirectory.isEmpty() && !m_processing;
    
    m_startButton->setEnabled(canProcess);
    m_stopButton->setEnabled(m_processing);
    m_previewStackButton->setEnabled(canProcess && (m_processingMode == MODE_ASTROMETRIC_STACKING || m_processingMode == MODE_FULL_PIPELINE));
    
    updateConnectionStatus();
}

void StellinaProcessor::updateConnectionStatus() {
    if (m_sirilClient->isConnected()) {
        m_connectionStatus->setText("Connected to Siril");
        m_connectionStatus->setStyleSheet("color: green; font-weight: bold;");
        m_testConnectionButton->setText("Disconnect");
    } else {
        m_connectionStatus->setText("Not connected");
        m_connectionStatus->setStyleSheet("color: red; font-weight: bold;");
        m_testConnectionButton->setText("Test Connection");
    }
}

void StellinaProcessor::logMessage(const QString &message, const QString &color) {
    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    QString formattedMessage = QString("<span style='color: %1'>[%2] %3</span>")
                                  .arg(color)
                                  .arg(timestamp)
                                  .arg(message);
    
    m_logTextEdit->append(formattedMessage);
    
    QTextCursor cursor = m_logTextEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    m_logTextEdit->setTextCursor(cursor);
    
    if (color == "red" || color == "green" || color == "blue") {
        m_statusLabel->setText(message);
    }
}

// Existing functions that remain mostly unchanged...
bool StellinaProcessor::findStellinaImages() {
    m_imagesToProcess.clear();
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: '%1'").arg(m_sourceDirectory), "red");
        return false;
    }
    
    logMessage(QString("Scanning directory: %1").arg(sourceDir.absolutePath()), "blue");
    
    QStringList allFiles = sourceDir.entryList(QDir::Files);
    logMessage(QString("Found %1 total files in directory").arg(allFiles.count()), "gray");
    
    QStringList jsonFiles = sourceDir.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
    QStringList fitFiles = sourceDir.entryList(QStringList() << "*.fit", QDir::Files);
    QStringList allFitsFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
    
    logMessage(QString("File breakdown: %1 .json, %2 .fit, %3 .fits files")
                  .arg(jsonFiles.count()).arg(fitFiles.count()).arg(allFitsFiles.count()), "blue");
    
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.count()), "blue");
    
    int validPairs = 0;
    int jsonMissing = 0;
    int qualityRejected = 0;
    
    for (const QString &fitsFile : fitsFiles) {
        QString baseName = QFileInfo(fitsFile).baseName();
        
        QStringList jsonCandidates = {
            baseName + ".json",
            baseName + ".JSON",
            baseName + "-stacking.json",
            baseName + "-stacking.JSON",
            QFileInfo(fitsFile).completeBaseName() + ".json",
            QFileInfo(fitsFile).completeBaseName() + ".JSON",
            QFileInfo(fitsFile).completeBaseName() + "-stacking.json",
            QFileInfo(fitsFile).completeBaseName() + "-stacking.JSON"
        };
        
        QString jsonFile;
        bool jsonFound = false;
        
        for (const QString &candidate : jsonCandidates) {
            if (sourceDir.exists(candidate)) {
                jsonFile = candidate;
                jsonFound = true;
                break;
            }
        }
        
        if (jsonFound) {
            QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
            QString jsonPath = sourceDir.absoluteFilePath(jsonFile);
            
            if (m_qualityFilter) {
                QJsonObject json = loadStellinaJson(jsonPath);
                if (!json.isEmpty() && !checkStellinaQuality(json)) {
                    qualityRejected++;
                    if (m_debugMode) {
                        logMessage(QString("Rejected %1: failed quality check").arg(fitsFile), "gray");
                    }
                    continue;
                }
            }
            
            m_imagesToProcess.append(fitsPath);
            validPairs++;
        } else {
            jsonMissing++;
            if (m_debugMode && jsonMissing <= 10) {
                logMessage(QString("No JSON file found for %1").arg(fitsFile), "orange");
            }
        }
    }
    
    logMessage(QString("File pairing results: %1 valid pairs, %2 missing JSON, %3 quality rejected")
                  .arg(validPairs).arg(jsonMissing).arg(qualityRejected), "blue");
    
    return !m_imagesToProcess.isEmpty();
}

QJsonObject StellinaProcessor::loadStellinaJson(const QString &jsonPath) {
    QFile file(jsonPath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (m_debugMode) {
            logMessage(QString("Failed to open JSON file: %1").arg(jsonPath), "red");
        }
        return QJsonObject();
    }
    
    QJsonParseError error;
    QJsonDocument doc = QJsonDocument::fromJson(file.readAll(), &error);
    
    if (error.error != QJsonParseError::NoError) {
        if (m_debugMode) {
            logMessage(QString("JSON parse error in %1: %2").arg(jsonPath).arg(error.errorString()), "red");
        }
        return QJsonObject();
    }
    
    return doc.object();
}

bool StellinaProcessor::extractCoordinates(const QJsonObject &json, double &alt, double &az) {
    if (json.contains("motors")) {
        QJsonObject motors = json["motors"].toObject();
        if (motors.contains("ALT") && motors.contains("AZ")) {
            alt = motors["ALT"].toDouble();
            az = motors["AZ"].toDouble();
            return true;
        }
    }
    
    if (json.contains("altitude") && json.contains("azimuth")) {
        alt = json["altitude"].toDouble();
        az = json["azimuth"].toDouble();
        return true;
    }
    
    if (json.contains("alt") && json.contains("az")) {
        alt = json["alt"].toDouble();
        az = json["az"].toDouble();
        return true;
    }
    
    return false;
}

bool StellinaProcessor::convertAltAzToRaDec(double alt, double az, const QString &dateObs, double &ra, double &dec) {
    // Placeholder implementation
    ra = az;
    dec = alt;
    Q_UNUSED(dateObs)
    return true;
}

bool StellinaProcessor::checkStellinaQuality(const QJsonObject &json) {
    if (json.contains("stackingData")) {
        QJsonObject stackingData = json["stackingData"].toObject();
        
        if (stackingData.contains("error") && !stackingData["error"].isNull()) {
            return false;
        }
        
        if (stackingData.contains("liveRegistrationResult")) {
            QJsonObject regResult = stackingData["liveRegistrationResult"].toObject();
            
            if (regResult.contains("status")) {
                int status = regResult["status"].toInt();
                if (status != 0) {
                    return false;
                }
            }
            
            if (regResult.contains("statusMessage")) {
                QString statusMsg = regResult["statusMessage"].toString();
                if (statusMsg != "StackingOk") {
                    return false;
                }
            }
            
            if (regResult.contains("starsUsed")) {
                int starsUsed = regResult["starsUsed"].toInt();
                if (starsUsed < 10) {
                    return false;
                }
            }
        }
    }
    
    if (json.contains("quality")) {
        return json["quality"].toBool();
    }
    
    if (json.contains("used_for_stacking")) {
        return json["used_for_stacking"].toBool();
    }
    
    return true;
}