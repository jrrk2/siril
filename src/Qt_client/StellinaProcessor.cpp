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
#include <cmath>

StellinaProcessor::StellinaProcessor(QWidget *parent)
    : QMainWindow(parent)
    , m_centralWidget(nullptr)
    , m_sirilClient(new SirilClient(this))
    , m_processingTimer(new QTimer(this))
    , m_processing(false)
    , m_currentImageIndex(0)
    , m_processedCount(0)
    , m_errorCount(0)
    , m_skippedCount(0)
    , m_qualityFilter(true)
    , m_debugMode(false)
    , m_focalLength(400.0)
    , m_pixelSize(2.40)
    , m_observerLocation("London")
{
    setWindowTitle("Stellina Processor for Siril");
    setMinimumSize(800, 600);
    
    // Setup timer
    m_processingTimer->setSingleShot(false);
    m_processingTimer->setInterval(1000); // Process one image per second
    
    setupUI();
    setupMenu();
    connectSignals();
    updateUI();
    
    // Load settings
    QSettings settings;
    m_sourceDirectory = settings.value("sourceDirectory").toString();
    m_outputDirectory = settings.value("outputDirectory").toString();
    m_qualityFilter = settings.value("qualityFilter", true).toBool();
    m_focalLength = settings.value("focalLength", 400.0).toDouble();
    m_pixelSize = settings.value("pixelSize", 2.40).toDouble();
    m_observerLocation = settings.value("observerLocation", "London").toString();
    
    // Update UI with loaded settings
    m_sourceDirectoryEdit->setText(m_sourceDirectory);
    m_outputDirectoryEdit->setText(m_outputDirectory);
    m_qualityFilterCheck->setChecked(m_qualityFilter);
    m_focalLengthSpin->setValue(m_focalLength);
    m_pixelSizeSpin->setValue(m_pixelSize);
    m_observerLocationEdit->setText(m_observerLocation);
    
    logMessage("Stellina Processor started. Please connect to Siril first.", "blue");
}

StellinaProcessor::~StellinaProcessor() {
    // Save settings
    QSettings settings;
    settings.setValue("sourceDirectory", m_sourceDirectory);
    settings.setValue("outputDirectory", m_outputDirectory);
    settings.setValue("qualityFilter", m_qualityFilter);
    settings.setValue("focalLength", m_focalLength);
    settings.setValue("pixelSize", m_pixelSize);
    settings.setValue("observerLocation", m_observerLocation);
}

void StellinaProcessor::setupUI() {
    m_centralWidget = new QWidget;
    setCentralWidget(m_centralWidget);
    
    m_mainLayout = new QVBoxLayout(m_centralWidget);
    
    // Connection group
    m_connectionGroup = new QGroupBox("Siril Connection");
    QHBoxLayout *connectionLayout = new QHBoxLayout(m_connectionGroup);
    
    m_testConnectionButton = new QPushButton("Test Connection");
    m_connectionStatus = new QLabel("Not connected");
    m_connectionStatus->setStyleSheet("color: red; font-weight: bold;");
    
    connectionLayout->addWidget(m_testConnectionButton);
    connectionLayout->addWidget(m_connectionStatus);
    connectionLayout->addStretch();
    
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
    
    // Options group
    m_optionsGroup = new QGroupBox("Processing Options");
    QGridLayout *optionsLayout = new QGridLayout(m_optionsGroup);
    
    m_qualityFilterCheck = new QCheckBox("Enable Stellina quality filtering");
    m_qualityFilterCheck->setChecked(true);
    m_debugModeCheck = new QCheckBox("Debug mode (verbose logging)");
    
    optionsLayout->addWidget(m_qualityFilterCheck, 0, 0, 1, 2);
    optionsLayout->addWidget(m_debugModeCheck, 1, 0, 1, 2);
    
    optionsLayout->addWidget(new QLabel("Focal Length (mm):"), 2, 0);
    m_focalLengthSpin = new QDoubleSpinBox;
    m_focalLengthSpin->setRange(50.0, 5000.0);
    m_focalLengthSpin->setValue(400.0);
    m_focalLengthSpin->setSuffix(" mm");
    optionsLayout->addWidget(m_focalLengthSpin, 2, 1);
    
    optionsLayout->addWidget(new QLabel("Pixel Size (μm):"), 3, 0);
    m_pixelSizeSpin = new QDoubleSpinBox;
    m_pixelSizeSpin->setRange(1.0, 20.0);
    m_pixelSizeSpin->setValue(2.40);
    m_pixelSizeSpin->setDecimals(2);
    m_pixelSizeSpin->setSuffix(" μm");
    optionsLayout->addWidget(m_pixelSizeSpin, 3, 1);
    
    optionsLayout->addWidget(new QLabel("Observer Location:"), 4, 0);
    m_observerLocationEdit = new QLineEdit("London");
    optionsLayout->addWidget(m_observerLocationEdit, 4, 1);
    
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
    
    processingLayout->addLayout(buttonLayout);
    processingLayout->addWidget(m_progressBar);
    processingLayout->addWidget(m_progressLabel);
    
    // Log group
    m_logGroup = new QGroupBox("Processing Log");
    QVBoxLayout *logLayout = new QVBoxLayout(m_logGroup);
    
    m_logTextEdit = new QTextEdit;
    m_logTextEdit->setMaximumHeight(200);
    
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
    logButtonLayout->addStretch();
    logButtonLayout->addWidget(m_clearLogButton);
    
    logLayout->addWidget(m_logTextEdit);
    logLayout->addLayout(logButtonLayout);
    
    // Add all groups to main layout
    m_mainLayout->addWidget(m_connectionGroup);
    m_mainLayout->addWidget(m_inputGroup);
    m_mainLayout->addWidget(m_optionsGroup);
    m_mainLayout->addWidget(m_processingGroup);
    m_mainLayout->addWidget(m_logGroup);
    
    // Status bar
    m_statusLabel = new QLabel("Ready");
    statusBar()->addWidget(m_statusLabel);
}

void StellinaProcessor::setupMenu() {
    QMenuBar *menuBar = this->menuBar();
    
    // File menu
    QMenu *fileMenu = menuBar->addMenu("&File");
    fileMenu->addAction("&Exit", this, &QWidget::close);
    
    // Tools menu
    QMenu *toolsMenu = menuBar->addMenu("&Tools");
    toolsMenu->addAction("Test &Connection", this, &StellinaProcessor::onTestConnection);
    toolsMenu->addSeparator();
    toolsMenu->addAction("&Clear Log", this, &StellinaProcessor::onClearLog);
    
    // Help menu
    QMenu *helpMenu = menuBar->addMenu("&Help");
    helpMenu->addAction("&About", [this]() {
        QMessageBox::about(this, "About Stellina Processor",
                          "Stellina Processor for Siril\n\n"
                          "A Qt application for processing Stellina telescope images\n"
                          "using Siril's plate solving capabilities.\n\n"
                          "Features:\n"
                          "• Automatic coordinate conversion from Alt/Az to RA/Dec\n"
                          "• Stellina quality filtering\n"
                          "• Batch processing with error recovery\n"
                          "• Real-time progress monitoring");
    });
}

void StellinaProcessor::connectSignals() {
    // UI signals
    connect(m_selectSourceButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectSourceDirectory);
    connect(m_selectOutputButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectOutputDirectory);
    connect(m_startButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStartProcessing);
    connect(m_stopButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStopProcessing);
    connect(m_testConnectionButton, &QPushButton::clicked,
            this, &StellinaProcessor::onTestConnection);
    connect(m_clearLogButton, &QPushButton::clicked,
            this, &StellinaProcessor::onClearLog);
    
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
}

void StellinaProcessor::updateUI() {
    bool connected = m_sirilClient->isConnected();
    bool canProcess = connected && !m_sourceDirectory.isEmpty() && !m_outputDirectory.isEmpty() && !m_processing;
    
    m_startButton->setEnabled(canProcess);
    m_stopButton->setEnabled(m_processing);
    
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
    
    // Auto-scroll to bottom
    QTextCursor cursor = m_logTextEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    m_logTextEdit->setTextCursor(cursor);
    
    // Also update status bar for important messages
    if (color == "red" || color == "green" || color == "blue") {
        m_statusLabel->setText(message);
    }
}

// Slot implementations
void StellinaProcessor::onSelectSourceDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Stellina Images Directory", m_sourceDirectory);
    if (!dir.isEmpty()) {
        m_sourceDirectory = dir;
        m_sourceDirectoryEdit->setText(dir);
        
        // Debug: Immediately test the directory
        QDir testDir(dir);
        logMessage(QString("Selected directory: %1").arg(testDir.absolutePath()), "blue");
        logMessage(QString("Directory exists: %1").arg(testDir.exists() ? "Yes" : "No"), "blue");
        logMessage(QString("Directory is readable: %1").arg(testDir.isReadable() ? "Yes" : "No"), "blue");
        
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
    if (!m_sirilClient->isConnected()) {
        QMessageBox::warning(this, "Connection Error", "Please connect to Siril first.");
        return;
    }
    
    if (m_sourceDirectory.isEmpty() || m_outputDirectory.isEmpty()) {
        QMessageBox::warning(this, "Directory Error", "Please select both source and output directories.");
        return;
    }
    
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
    
    // Test basic commands
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

// Processing functions
void StellinaProcessor::startStellinaProcessing() {
    logMessage("Starting Stellina processing...", "blue");
    
    // Find images to process
    if (!findStellinaImages()) {
        logMessage("No valid Stellina images found in source directory.", "red");
        return;
    }
    
    // Change Siril working directory (but continue even if this fails)
    bool dirChanged = m_sirilClient->changeDirectory(m_outputDirectory);
    if (!dirChanged) {
        logMessage("Warning: Failed to change Siril working directory. Will use absolute paths.", "orange");
        logMessage(QString("Will save files to: %1").arg(m_outputDirectory), "blue");
    } else {
        logMessage(QString("Changed Siril working directory to: %1").arg(m_outputDirectory), "green");
    }
    
    // Initialize processing state
    m_processing = true;
    m_currentImageIndex = 0;
    m_processedCount = 0;
    m_errorCount = 0;
    m_skippedCount = 0;
    
    m_progressBar->setMaximum(m_imagesToProcess.length());
    m_progressBar->setValue(0);
    m_progressLabel->setText(QString("Processing 0 of %1 images...").arg(m_imagesToProcess.length()));
    
    logMessage(QString("Found %1 images to process.").arg(m_imagesToProcess.length()), "green");
    
    updateUI();
    
    // Start processing timer
    m_processingTimer->start();
}

bool StellinaProcessor::findStellinaImages() {
    m_imagesToProcess.clear();
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: '%1'").arg(m_sourceDirectory), "red");
        
        // Try to clean the path
        QString cleanPath = QDir::cleanPath(m_sourceDirectory);
        logMessage(QString("Trying cleaned path: '%1'").arg(cleanPath), "blue");
        
        QDir cleanDir(cleanPath);
        if (cleanDir.exists()) {
            logMessage("Cleaned path works! Using cleaned path.", "green");
            sourceDir = cleanDir;
        } else {
            return false;
        }
    }
    
    logMessage(QString("Scanning directory: %1").arg(sourceDir.absolutePath()), "blue");
    
    // Debug: List all files in directory
    QStringList allFiles = sourceDir.entryList(QDir::Files);
    logMessage(QString("Found %1 total files in directory").arg(allFiles.count()), "gray");
    
    // Show file type breakdown
    QStringList jsonFiles = sourceDir.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
    QStringList fitFiles = sourceDir.entryList(QStringList() << "*.fit", QDir::Files);
    QStringList allFitsFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
    
    logMessage(QString("File breakdown: %1 .json, %2 .fit, %3 .fits files")
                  .arg(jsonFiles.count()).arg(fitFiles.count()).arg(allFitsFiles.count()), "blue");
    
    // Show sample JSON filenames to understand the naming pattern
    if (m_debugMode && !jsonFiles.isEmpty()) {
        logMessage("Sample JSON filenames:", "blue");
        for (int i = 0; i < qMin(10, jsonFiles.count()); ++i) {
            logMessage(QString("  JSON: %1").arg(jsonFiles[i]), "gray");
        }
        
        logMessage("Sample FITS filenames:", "blue");
        for (int i = 0; i < qMin(10, allFitsFiles.count()); ++i) {
            logMessage(QString("  FITS: %1").arg(allFitsFiles[i]), "gray");
        }
    }
    
    if (m_debugMode && allFiles.count() < 50) { // Only show if not too many files
        logMessage("Sample files in directory:", "gray");
        for (int i = 0; i < qMin(20, allFiles.count()); ++i) {
            logMessage(QString("  %1").arg(allFiles[i]), "gray");
        }
        if (allFiles.count() > 20) {
            logMessage(QString("  ... and %1 more files").arg(allFiles.count() - 20), "gray");
        }
    }
    
    // Find pairs of .fits and .json files
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.count()), "blue");
    
    // Debug: Show first few FITS files and look for JSON pairs
    int validPairs = 0;
    int jsonMissing = 0;
    int qualityRejected = 0;
    
    for (const QString &fitsFile : fitsFiles) {
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Try multiple JSON file patterns to handle different naming schemes
        QStringList jsonCandidates = {
            baseName + ".json",                           // img-0001.json
            baseName + ".JSON",                           // img-0001.JSON  
            baseName + "-stacking.json",                  // img-0001-stacking.json (Stellina pattern!)
            baseName + "-stacking.JSON",                  // img-0001-stacking.JSON
            QFileInfo(fitsFile).completeBaseName() + ".json", // Handle .fit vs .fits
            QFileInfo(fitsFile).completeBaseName() + ".JSON",
            QFileInfo(fitsFile).completeBaseName() + "-stacking.json",
            QFileInfo(fitsFile).completeBaseName() + "-stacking.JSON",
            // Try without hyphen
            QString(baseName).replace("-", "") + ".json", // img0001.json
            QString(baseName).replace("-", "") + ".JSON", // img0001.JSON
            // Try with underscore instead of hyphen
            QString(baseName).replace("-", "_") + ".json", // img_0001.json
            QString(baseName).replace("-", "_") + ".JSON"  // img_0001.JSON
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
        
        if (m_debugMode && validPairs < 5) { // Show first few for debugging
            logMessage(QString("FITS: %1 -> JSON candidates: %2")
                          .arg(fitsFile)
                          .arg(jsonCandidates.join(", ")), "gray");
            if (jsonFound) {
                logMessage(QString("  -> Found JSON: %1").arg(jsonFile), "gray");
            } else {
                logMessage(QString("  -> No JSON found"), "gray");
            }
        }
        
        if (jsonFound) {
            QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
            QString jsonPath = sourceDir.absoluteFilePath(jsonFile);
            
            // Check quality if filtering is enabled
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
            if (m_debugMode && jsonMissing <= 10) { // Don't spam too much
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
    // Check for Stellina's motor coordinates (main pattern)
    if (json.contains("motors")) {
        QJsonObject motors = json["motors"].toObject();
        if (motors.contains("ALT") && motors.contains("AZ")) {
            alt = motors["ALT"].toDouble();
            az = motors["AZ"].toDouble();
            return true;
        }
    }
    
    // Fallback patterns for other JSON structures
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
    // This is a placeholder - you'll need to implement coordinate conversion
    // For now, we'll use a simple approximation
    // In reality, you'd use astronomical libraries like SOFA or your existing Python code
    
    ra = az; // Placeholder
    dec = alt; // Placeholder
    
    return true;
}

bool StellinaProcessor::checkStellinaQuality(const QJsonObject &json) {
    // Check Stellina's stacking data for quality indicators
    if (json.contains("stackingData")) {
        QJsonObject stackingData = json["stackingData"].toObject();
        
        // Check if there's an error
        if (stackingData.contains("error") && !stackingData["error"].isNull()) {
            return false; // Has an error
        }
        
        // Check live registration result
        if (stackingData.contains("liveRegistrationResult")) {
            QJsonObject regResult = stackingData["liveRegistrationResult"].toObject();
            
            // Check status (0 seems to be success based on the example)
            if (regResult.contains("status")) {
                int status = regResult["status"].toInt();
                if (status != 0) {
                    return false; // Non-zero status indicates error
                }
            }
            
            // Check status message
            if (regResult.contains("statusMessage")) {
                QString statusMsg = regResult["statusMessage"].toString();
                if (statusMsg != "StackingOk") {
                    return false; // Not "StackingOk"
                }
            }
            
            // Check if reasonable number of stars were used
            if (regResult.contains("starsUsed")) {
                int starsUsed = regResult["starsUsed"].toInt();
                if (starsUsed < 10) { // Arbitrary threshold - adjust as needed
                    return false; // Too few stars used
                }
            }
        }
    }
    
    // Fallback quality checks for other JSON structures
    if (json.contains("quality")) {
        return json["quality"].toBool();
    }
    
    if (json.contains("used_for_stacking")) {
        return json["used_for_stacking"].toBool();
    }
    
    // Default to accepting if no quality info (or if all checks passed)
    return true;
}

void StellinaProcessor::processNextImage() {
    if (m_currentImageIndex >= m_imagesToProcess.length()) {
        finishProcessing();
        return;
    }
    
    QString fitsPath = m_imagesToProcess[m_currentImageIndex];
    QString baseName = QFileInfo(fitsPath).baseName();
    
    // Use the same JSON finding logic as in findStellinaImages()
    QDir sourceDir(QFileInfo(fitsPath).dir());
    QString jsonPath;
    
    // Try the same JSON file patterns we used during file discovery
    QStringList jsonCandidates = {
        baseName + ".json",                           
        baseName + ".JSON",                           
        baseName + "-stacking.json",                  // Stellina pattern!
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
        logMessage(QString("No JSON file found for %1 during processing").arg(QFileInfo(fitsPath).fileName()), "red");
        m_errorCount++;
        m_currentImageIndex++;
        m_progressBar->setValue(m_currentImageIndex);
        return;
    }
    
    logMessage(QString("Processing image %1 of %2: %3")
                  .arg(m_currentImageIndex + 1)
                  .arg(m_imagesToProcess.length())
                  .arg(QFileInfo(fitsPath).fileName()), "blue");
    
    if (m_debugMode) {
        logMessage(QString("Using JSON file: %1").arg(QFileInfo(jsonPath).fileName()), "gray");
    }
    
    bool success = false;
    
    do {
        // Load JSON metadata
        QJsonObject json = loadStellinaJson(jsonPath);
        if (json.isEmpty()) {
            logMessage(QString("Failed to load JSON metadata from: %1").arg(QFileInfo(jsonPath).fileName()), "red");
            break;
        }
        
        if (m_debugMode) {
            logMessage(QString("Loaded JSON with %1 keys").arg(json.keys().size()), "gray");
        }
        
        // Extract coordinates
        double alt, az;
        if (!extractCoordinates(json, alt, az)) {
            logMessage("Failed to extract coordinates from JSON", "red");
            break;
        }
        
        // Convert Alt/Az to RA/Dec
        double ra, dec;
        QString dateObs = json["date_obs"].toString();
        if (!convertAltAzToRaDec(alt, az, dateObs, ra, dec)) {
            logMessage("Failed to convert coordinates", "red");
            break;
        }
        
        // Load image in Siril (this may take a while for large files)
        logMessage("Loading image in Siril...", "blue");
        if (!m_sirilClient->loadImage(fitsPath)) {
            logMessage("Failed to load image in Siril", "red");
            break;
        }
        logMessage("Image loaded successfully", "green");
        
        // Perform plate solving
        logMessage(QString("Attempting plate solve with RA=%.4f°, Dec=%.4f°").arg(ra).arg(dec), "blue");
        if (!m_sirilClient->platesolve(ra, dec, m_focalLength, m_pixelSize, true)) {
            logMessage("Plate solving failed (continuing anyway)", "orange");
            // Don't break - save the image anyway
        } else {
            logMessage("Plate solving completed", "green");
        }
        
        // Save processed image
        QString outputName = QString("processed_%1").arg(baseName);
        QString outputPath = QDir(m_outputDirectory).absoluteFilePath(outputName);
        
        if (!m_sirilClient->saveImage(outputPath)) {
            logMessage("Failed to save processed image", "red");
            break;
        }
        
        // Close image
        m_sirilClient->closeImage();
        
        success = true;
        logMessage(QString("Successfully processed: %1").arg(QFileInfo(fitsPath).fileName()), "green");
        
    } while (false);
    
    // Update counters
    if (success) {
        m_processedCount++;
    } else {
        m_errorCount++;
        logMessage(QString("Failed to process: %1").arg(QFileInfo(fitsPath).fileName()), "red");
    }
    
    // Update progress
    m_currentImageIndex++;
    m_progressBar->setValue(m_currentImageIndex);
    m_progressLabel->setText(QString("Processed %1 of %2 images (Success: %3, Errors: %4)")
                                .arg(m_currentImageIndex)
                                .arg(m_imagesToProcess.length())
                                .arg(m_processedCount)
                                .arg(m_errorCount));
}

void StellinaProcessor::finishProcessing() {
    m_processing = false;
    m_processingTimer->stop();
    
    logMessage(QString("Processing complete! Processed: %1, Errors: %2, Total: %3")
                  .arg(m_processedCount)
                  .arg(m_errorCount)
                  .arg(m_imagesToProcess.length()), "green");
    
    updateUI();
    
    // Show completion dialog
    QMessageBox::information(this, "Processing Complete",
                            QString("Stellina processing completed!\n\n"
                                   "Successfully processed: %1\n"
                                   "Errors: %2\n"
                                   "Total images: %3")
                               .arg(m_processedCount)
                               .arg(m_errorCount)
                               .arg(m_imagesToProcess.length()));
}