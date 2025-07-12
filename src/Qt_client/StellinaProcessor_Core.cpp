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
    m_processingTimer->setInterval(2000);
    
    setupUI();
    setupMenu();
    connectSignals();
    updateUI();
    loadSettings();
    
    logMessage("Enhanced Stellina Processor started. Connect to Siril and select processing mode.", "blue");
    
    // Scan for dark frames if directory is set
    if (!m_darkDirectory.isEmpty()) {
        scanDarkFrames();
    }
}

StellinaProcessor::~StellinaProcessor() {
    saveSettings();
}

void StellinaProcessor::loadSettings() {
    QSettings settings;
    m_sourceDirectory = settings.value("sourceDirectory").toString();
    m_darkDirectory = settings.value("darkDirectory").toString();
    m_calibratedDirectory = settings.value("calibratedDirectory").toString();
    m_plateSolvedDirectory = settings.value("plateSolvedDirectory").toString();
    m_stackedDirectory = settings.value("stackedDirectory").toString();
    m_qualityFilter = settings.value("qualityFilter", true).toBool();
    m_focalLength = settings.value("focalLength", 400.0).toDouble();
    m_pixelSize = settings.value("pixelSize", 2.40).toDouble();
    m_observerLocation = settings.value("observerLocation", "London").toString();
    m_processingMode = static_cast<ProcessingMode>(settings.value("processingMode", 0).toInt());
    
    // Update UI with loaded settings
    m_sourceDirectoryEdit->setText(m_sourceDirectory);
    m_darkDirectoryEdit->setText(m_darkDirectory);
    m_calibratedDirectoryEdit->setText(m_calibratedDirectory);
    m_plateSolvedDirectoryEdit->setText(m_plateSolvedDirectory);
    m_stackedDirectoryEdit->setText(m_stackedDirectory);
    m_qualityFilterCheck->setChecked(m_qualityFilter);
    m_focalLengthSpin->setValue(m_focalLength);
    m_pixelSizeSpin->setValue(m_pixelSize);
    m_observerLocationEdit->setText(m_observerLocation);
    m_processingModeCombo->setCurrentIndex(m_processingMode);
}

void StellinaProcessor::saveSettings() {
    QSettings settings;
    settings.setValue("sourceDirectory", m_sourceDirectory);
    settings.setValue("darkDirectory", m_darkDirectory);
    settings.setValue("calibratedDirectory", m_calibratedDirectory);
    settings.setValue("plateSolvedDirectory", m_plateSolvedDirectory);
    settings.setValue("stackedDirectory", m_stackedDirectory);
    settings.setValue("qualityFilter", m_qualityFilter);
    settings.setValue("focalLength", m_focalLength);
    settings.setValue("pixelSize", m_pixelSize);
    settings.setValue("observerLocation", m_observerLocation);
    settings.setValue("processingMode", static_cast<int>(m_processingMode));
}

QString StellinaProcessor::getOutputDirectoryForCurrentStage() const {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        return m_calibratedDirectory.isEmpty() ? m_sourceDirectory : m_calibratedDirectory;
    case STAGE_PLATE_SOLVING:
        return m_plateSolvedDirectory.isEmpty() ? m_sourceDirectory : m_plateSolvedDirectory;
    case STAGE_REGISTRATION:
    case STAGE_STACKING:
        return m_stackedDirectory.isEmpty() ? m_plateSolvedDirectory : m_stackedDirectory;
    case STAGE_COMPLETE:
        return m_stackedDirectory.isEmpty() ? m_plateSolvedDirectory : m_stackedDirectory;
    default:
        return m_sourceDirectory;
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
    // This method is now deprecated - we always use findAllMatchingDarkFrames instead
    // Keep it for compatibility but it's not used in the new workflow
    Q_UNUSED(lightFrame)
    Q_UNUSED(darkFrame)
    return false;
}

QStringList StellinaProcessor::findAllMatchingDarkFrames(int targetExposure, int targetTemperature, const QString &targetBinning) {
    QStringList matchingDarks;
    
    if (m_darkDirectory.isEmpty()) {
        return matchingDarks;
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        return matchingDarks;
    }
    
    QStringList darkFiles = darkDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperature = extractTemperature(fullPath);
        QString binning = extractBinning(fullPath);
        
        // Check if this dark frame matches the light frame characteristics
        bool exposureMatch = qAbs(exposure - targetExposure) <= (targetExposure * m_exposureTolerance / 100);
        bool temperatureMatch = qAbs(temperature - targetTemperature) <= m_temperatureTolerance;
        bool binningMatch = (binning == targetBinning);
        
        if (exposureMatch && temperatureMatch && binningMatch) {
            matchingDarks.append(fullPath);
        }
    }
    
    return matchingDarks;
}

bool StellinaProcessor::createMasterDark(const QStringList &darkFrames, const QString &outputPath) {
    if (darkFrames.isEmpty()) {
        logMessage("No dark frames provided for master dark creation", "red");
        return false;
    }
    
    logMessage(QString("Creating master dark from %1 frames...").arg(darkFrames.size()), "blue");
    
    if (darkFrames.size() == 1) {
        // Only one dark frame, just copy it
        if (QFile::copy(darkFrames.first(), outputPath)) {
            logMessage("Master dark created (single frame copy)", "green");
            return true;
        } else {
            logMessage("Failed to copy single dark frame", "red");
            return false;
        }
    }
    
    // Multiple dark frames - stack them using Siril
    QString tempSeqName = QString("temp_dark_seq_%1").arg(QDateTime::currentMSecsSinceEpoch());
    QString tempDir = QDir(m_calibratedDirectory).absoluteFilePath(tempSeqName);
    QDir().mkpath(tempDir);
    
    // Copy dark frames to temporary directory with sequential naming
    for (int i = 0; i < darkFrames.size(); ++i) {
        QString srcFile = darkFrames[i];
        QString dstFile = QDir(tempDir).absoluteFilePath(QString("%1_%2.fits")
                                                            .arg(tempSeqName)
                                                            .arg(i + 1, 4, 10, QChar('0')));
        
        if (!QFile::copy(srcFile, dstFile)) {
            logMessage(QString("Failed to copy dark frame: %1").arg(QFileInfo(srcFile).fileName()), "red");
            return false;
        }
    }
    
    // Change to temp directory
    QString oldWorkDir = m_sirilClient->getWorkingDirectory();
    if (!m_sirilClient->changeDirectory(tempDir)) {
        logMessage("Failed to change to temp directory", "red");
        return false;
    }
    
    // Create sequence
    QString convertCmd = QString("convert %1 -out=%1").arg(tempSeqName);
    if (!m_sirilClient->sendSirilCommand(convertCmd)) {
        logMessage("Failed to create dark sequence", "red");
        m_sirilClient->changeDirectory(oldWorkDir);
        return false;
    }
    
    // Load the sequence
    QString loadCmd = QString("load %1").arg(tempSeqName);
    if (!m_sirilClient->sendSirilCommand(loadCmd)) {
        logMessage("Failed to load dark sequence", "red");
        m_sirilClient->changeDirectory(oldWorkDir);
        return false;
    }
    
    // Stack using median (best for dark frames)
    QString stackCmd = QString("stack %1 rej none -type=median -norm=none").arg(tempSeqName);
    if (!m_sirilClient->sendSirilCommand(stackCmd)) {
        logMessage("Failed to stack dark frames", "red");
        m_sirilClient->changeDirectory(oldWorkDir);
        return false;
    }
    
    // Save the master dark
    QString tempMasterPath = QDir(tempDir).absoluteFilePath(QString("stacked_%1.fits").arg(tempSeqName));
    if (!m_sirilClient->saveImage(tempMasterPath)) {
        logMessage("Failed to save master dark", "red");
        m_sirilClient->changeDirectory(oldWorkDir);
        return false;
    }
    
    m_sirilClient->closeImage();
    
    // Restore working directory
    m_sirilClient->changeDirectory(oldWorkDir);
    
    // Move master dark to final location
    if (QFile::exists(outputPath)) {
        QFile::remove(outputPath);
    }
    
    if (QFile::rename(tempMasterPath, outputPath)) {
        // Clean up temp directory
        QDir tempDirObj(tempDir);
        tempDirObj.removeRecursively();
        
        logMessage(QString("Master dark created successfully from %1 frames").arg(darkFrames.size()), "green");
        return true;
    } else {
        logMessage("Failed to move master dark to final location", "red");
        return false;
    }
}

bool StellinaProcessor::applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame) {
    logMessage(QString("Applying master dark to %1").arg(QFileInfo(lightFrame).fileName()), "blue");
    
    // Load light frame
    if (!m_sirilClient->loadImage(lightFrame)) {
        logMessage("Failed to load light frame", "red");
        return false;
    }
    
    // Try Siril's calibrate command first
    QString calibrateCmd = QString("calibrate_single %1").arg(masterDark);
    bool calibrateSuccess = m_sirilClient->sendSirilCommand(calibrateCmd);
    
    if (!calibrateSuccess) {
        logMessage("calibrate_single failed, using manual subtraction", "orange");
        
        // Manual subtraction method
        QString loadDarkCmd = QString("load %1 $masterdark").arg(masterDark);
        if (!m_sirilClient->sendSirilCommand(loadDarkCmd)) {
            logMessage("Failed to load master dark", "red");
            m_sirilClient->closeImage();
            return false;
        }
        
        // Subtract master dark from light frame
        if (!m_sirilClient->sendSirilCommand("sub $masterdark")) {
            logMessage("Failed to subtract master dark", "red");
            m_sirilClient->closeImage();
            return false;
        }
    }
    
    // Save calibrated frame
    if (!m_sirilClient->saveImage(outputFrame)) {
        logMessage("Failed to save calibrated frame", "red");
        m_sirilClient->closeImage();
        return false;
    }
    
    m_sirilClient->closeImage();
    logMessage("Master dark applied successfully", "green");
    return true;
}

bool StellinaProcessor::processImageDarkCalibration(const QString &lightFrame) {
    m_currentTaskLabel->setText("Dark calibration...");
    
    // Extract light frame characteristics
    int lightExposure = extractExposureTime(lightFrame);
    int lightTemperature = extractTemperature(lightFrame);
    QString lightBinning = extractBinning(lightFrame);
    
    if (lightExposure <= 0) {
        logMessage("Could not determine light frame exposure time", "red");
        return false;
    }
    
    // Find all matching dark frames
    QStringList matchingDarks = findAllMatchingDarkFrames(lightExposure, lightTemperature, lightBinning);
    
    if (matchingDarks.isEmpty()) {
        logMessage(QString("No matching dark frames found for exposure=%1s, temp=%2°C, binning=%3")
                      .arg(lightExposure).arg(lightTemperature).arg(lightBinning), "orange");
        m_skippedCount++;
        return true; // Continue processing without dark calibration
    }
    
    logMessage(QString("Found %1 matching dark frames").arg(matchingDarks.size()), "blue");
    
    // Create master dark for this combination
    QString masterDarkName = QString("master_dark_%1s_%2C_%3.fits")
                                .arg(lightExposure)
                                .arg(lightTemperature)
                                .arg(lightBinning);
    QString masterDarkPath = QDir(m_calibratedDirectory).absoluteFilePath(masterDarkName);
    
    // Check if master dark already exists
    if (!QFile::exists(masterDarkPath)) {
        if (!createMasterDark(matchingDarks, masterDarkPath)) {
            logMessage("Failed to create master dark", "red");
            return false;
        }
        logMessage(QString("Created master dark: %1").arg(masterDarkName), "green");
    } else {
        logMessage(QString("Using existing master dark: %1").arg(masterDarkName), "blue");
    }
    
    // Apply master dark to light frame
    QString outputName = QString("dark_calibrated_%1.fits")
                            .arg(QFileInfo(lightFrame).baseName());
    QString outputPath = QDir(m_calibratedDirectory).absoluteFilePath(outputName);
    
    if (applyMasterDark(lightFrame, masterDarkPath, outputPath)) {
        m_darkCalibratedFiles.append(outputPath);
        m_darkCalibratedCount++;
        logMessage(QString("Dark calibration successful: %1").arg(outputName), "green");
        return true;
    } else {
        logMessage("Dark calibration failed", "red");
        return false;
    }
}

// Astrometric stacking functions
bool StellinaProcessor::createSequence(const QStringList &imageList, const QString &sequenceName) {
    if (imageList.isEmpty()) {
        return false;
    }
    
    QString outputDir = getOutputDirectoryForCurrentStage();
    
    // Change to output directory
    if (!m_sirilClient->changeDirectory(outputDir)) {
        logMessage("Warning: Could not change to output directory", "orange");
    }
    
    // Copy images to a temporary location with sequential naming
    QString seqDir = QDir(outputDir).absoluteFilePath(sequenceName);
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
    
    QString outputDir = getOutputDirectoryForCurrentStage();
    m_finalStackedImage = QDir(outputDir).absoluteFilePath(outputName);
    m_stackingStatusLabel->setText("Stacking: Complete");
    logMessage(QString("Stacking completed: %1").arg(outputName), "green");
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

bool StellinaProcessor::registerImages(const QStringList &imageList, const QString &referenceImage) {
    Q_UNUSED(imageList)
    Q_UNUSED(referenceImage)
    // Implementation would use Siril's registration commands
    return true;
}

bool StellinaProcessor::stackRegisteredImages(const QStringList &registeredImages, const QString &outputStack) {
    Q_UNUSED(registeredImages)
    Q_UNUSED(outputStack)
    // Implementation would use Siril's stacking commands
    return true;
}

// Stellina image processing functions
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
    // Placeholder implementation - would need proper coordinate conversion
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

void StellinaProcessor::saveProcessingReport() {
    QString outputDir = getOutputDirectoryForCurrentStage();
    QString reportPath = QDir(outputDir).absoluteFilePath(
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