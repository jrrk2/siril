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
#include <QThread>
#include <cmath>
#include <vector>

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

// FITS metadata extraction functions - using cfitsio directly
int StellinaProcessor::extractExposureTime(const QString &fitsFile) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsFile.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return 10; // Fallback: 10 seconds
    }
    
    // For Stellina, exposure is in milliseconds under "EXPOSURE" keyword
    int exposure_ms = 10000; // Default: 10 seconds (10000ms)
    
    if (fits_read_key(fptr, TINT, "EXPOSURE", &exposure_ms, nullptr, &status) != 0) {
        // If EXPOSURE doesn't work, try other common keywords
        double exptime_s = 10.0;
        if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exptime_s, nullptr, &status) == 0) {
            exposure_ms = static_cast<int>(exptime_s * 1000); // Convert seconds to ms
        }
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    
    // Convert milliseconds to seconds for consistent usage
    return exposure_ms / 1000;
}

int StellinaProcessor::extractTemperature(const QString &fitsFile) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsFile.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return 284; // Fallback: 11°C = 284K
    }
    
    // For Stellina, temperature is under "TEMP" keyword in Celsius
    double temperature_celsius = 11.2; // Default from your example
    
    if (fits_read_key(fptr, TDOUBLE, "TEMP", &temperature_celsius, nullptr, &status) != 0) {
        // Try other common temperature keywords if TEMP fails
        QStringList tempKeys = {"CCD-TEMP", "CCD_TEMP", "TEMPERAT", "SET-TEMP"};
        for (const QString &key : tempKeys) {
            QByteArray keyBytes = key.toLocal8Bit();
            if (fits_read_key(fptr, TDOUBLE, keyBytes.data(), &temperature_celsius, nullptr, &status) == 0) {
                break;
            }
            status = 0;
        }
    }
    
    fits_close_file(fptr, &status);
    
    // Convert Celsius to Kelvin and round to nearest integer
    // K = °C + 273.15
    double temperature_kelvin = temperature_celsius + 273.15;
    return static_cast<int>(qRound(temperature_kelvin));
}

QString StellinaProcessor::extractBinning(const QString &fitsFile) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsFile.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return "1x1"; // Fallback
    }
    
    // Stellina doesn't seem to have explicit binning keywords in your example
    // Check the actual image dimensions vs sensor size to infer binning
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, &status) == 0) {
        // Your example shows 3072x2080 which appears to be full resolution
        // Stellina sensor is typically 3072x2080 at full res
        if (naxes[0] == 3072 && naxes[1] == 2080) {
            fits_close_file(fptr, &status);
            return "1x1"; // Full resolution
        } else if (naxes[0] == 1536 && naxes[1] == 1040) {
            fits_close_file(fptr, &status);
            return "2x2"; // Half resolution (2x2 binning)
        } else if (naxes[0] == 1024 && naxes[1] == 693) {
            fits_close_file(fptr, &status);
            return "3x3"; // Third resolution (3x3 binning)
        }
    }
    
    // Try explicit binning keywords anyway
    int xbin = 1, ybin = 1;
    QStringList xbinKeys = {"XBINNING", "BINX", "BIN_X"};
    QStringList ybinKeys = {"YBINNING", "BINY", "BIN_Y"};
    
    for (const QString &key : xbinKeys) {
        QByteArray keyBytes = key.toLocal8Bit();
        if (fits_read_key(fptr, TINT, keyBytes.data(), &xbin, nullptr, &status) == 0) {
            break;
        }
        status = 0;
    }
    
    for (const QString &key : ybinKeys) {
        QByteArray keyBytes = key.toLocal8Bit();
        if (fits_read_key(fptr, TINT, keyBytes.data(), &ybin, nullptr, &status) == 0) {
            break;
        }
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    return QString("%1x%2").arg(xbin).arg(ybin);
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
    QMap<QString, DarkFrame> darkInfo; // Store the representative info for each group
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperatureK = extractTemperature(fullPath); // Now returns Kelvin
        QString binning = extractBinning(fullPath);
        
        if (exposure > 0) { // Valid dark frame
            QString key = QString("%1_%2_%3").arg(exposure).arg(temperatureK).arg(binning);
            darkGroups[key].append(fullPath);
            
            // Store representative info for this group (first file sets the pattern)
            if (!darkInfo.contains(key)) {
                DarkFrame dark;
                dark.filepath = fullPath; // Representative file
                dark.exposure = exposure;
                dark.temperature = temperatureK; // Store as Kelvin
                dark.binning = binning;
                darkInfo[key] = dark;
                
                if (m_debugMode) {
                    int temperatureC = temperatureK - 273; // Convert back to Celsius for display
                    logMessage(QString("Dark group %1: %2s, %3K (%4°C), %5 - first file: %6")
                                  .arg(darkInfo.size())
                                  .arg(exposure)
                                  .arg(temperatureK)
                                  .arg(temperatureC)
                                  .arg(binning)
                                  .arg(QFileInfo(darkFile).fileName()), "gray");
                }
            }
        } else {
            if (m_debugMode) {
                logMessage(QString("Skipped invalid dark frame: %1").arg(QFileInfo(darkFile).fileName()), "orange");
            }
        }
    }
    
    // Create DarkFrame entries from the groups
    m_darkFrames.clear();
    for (auto it = darkGroups.begin(); it != darkGroups.end(); ++it) {
        QString key = it.key();
        if (darkInfo.contains(key)) {
            DarkFrame dark = darkInfo[key];
            // Update filepath to point to the first file in the group
            dark.filepath = it.value().first();
            m_darkFrames.append(dark);
            
            int temperatureC = dark.temperature - 273; // Convert for display
            logMessage(QString("Dark group: %1s exposure, %2K (%3°C), %4 binning → %5 frames")
                          .arg(dark.exposure)
                          .arg(dark.temperature)
                          .arg(temperatureC)
                          .arg(dark.binning)
                          .arg(it.value().size()), "blue");
        }
    }
    
    // Update UI
    m_darkFramesCount->setText(QString("%1 dark frame groups found").arg(m_darkFrames.size()));
    
    // Update dark frames table
    m_darkFramesTable->setRowCount(m_darkFrames.size());
    for (int i = 0; i < m_darkFrames.size(); ++i) {
        const DarkFrame &dark = m_darkFrames[i];
        QString key = QString("%1_%2_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning);
        
        int temperatureC = dark.temperature - 273;
        m_darkFramesTable->setItem(i, 0, new QTableWidgetItem(QString("%1s_%2K_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning)));
        m_darkFramesTable->setItem(i, 1, new QTableWidgetItem(QString("%1s").arg(dark.exposure)));
        m_darkFramesTable->setItem(i, 2, new QTableWidgetItem(QString("%1K (%2°C)").arg(dark.temperature).arg(temperatureC)));
        m_darkFramesTable->setItem(i, 3, new QTableWidgetItem(dark.binning));
        
        // Count how many dark frames in this group
        int count = darkGroups[key].size();
        m_darkFramesTable->setItem(i, 4, new QTableWidgetItem(QString::number(count)));
    }
    
    logMessage(QString("Dark frame scan complete: %1 groups found").arg(m_darkFrames.size()), "green");
    
    // Show summary of what was found
    if (!m_darkFrames.isEmpty()) {
        QStringList summary;
        for (const DarkFrame &dark : m_darkFrames) {
            QString key = QString("%1_%2_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning);
            int count = darkGroups[key].size();
            int temperatureC = dark.temperature - 273;
            summary.append(QString("%1×%2s@%3K").arg(count).arg(dark.exposure).arg(dark.temperature));
        }
        logMessage(QString("Found: %1").arg(summary.join(", ")), "blue");
    }
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
    
    // Use direct FITS manipulation for better performance
    logMessage("Using direct FITS processing for master dark creation", "blue");
    return createMasterDarkDirect(darkFrames, outputPath);
}

bool StellinaProcessor::createMasterDarkDirect(const QStringList &darkFrames, const QString &outputPath) {
    if (darkFrames.isEmpty()) {
        return false;
    }
    
    if (darkFrames.size() == 1) {
        return QFile::copy(darkFrames.first(), outputPath);
    }
    
    logMessage(QString("Creating master dark from %1 frames using direct FITS manipulation...").arg(darkFrames.size()), "blue");
    
    fitsfile *firstFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open first dark frame to get dimensions and header info
    QByteArray firstPath = darkFrames.first().toLocal8Bit();
    if (fits_open_file(&firstFits, firstPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open first dark frame: %1 (FITS error: %2)").arg(darkFrames.first()).arg(status), "red");
        return false;
    }
    
    // Get image dimensions
    int naxis;
    long naxes[2];
    if (fits_get_img_dim(firstFits, &naxis, &status) || 
        fits_get_img_size(firstFits, 2, naxes, &status)) {
        logMessage(QString("Failed to get image dimensions (FITS error: %1)").arg(status), "red");
        fits_close_file(firstFits, &status);
        return false;
    }
    
    long totalPixels = naxes[0] * naxes[1];
    logMessage(QString("Image dimensions: %1 x %2 (%3 pixels)").arg(naxes[0]).arg(naxes[1]).arg(totalPixels), "blue");
    
    // Create output FITS file
    QByteArray outputPathBytes = QString("!%1").arg(outputPath).toLocal8Bit(); // ! forces overwrite
    if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
        logMessage(QString("Failed to create output FITS file: %1 (FITS error: %2)").arg(outputPath).arg(status), "red");
        fits_close_file(firstFits, &status);
        return false;
    }
    
    // Copy header from first frame
    if (fits_copy_header(firstFits, outputFits, &status)) {
        logMessage(QString("Failed to copy FITS header (FITS error: %1)").arg(status), "red");
        fits_close_file(firstFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    fits_close_file(firstFits, &status);
    
    // Allocate memory for accumulation (using double for precision)
    std::vector<double> accumulator(totalPixels, 0.0);
    std::vector<float> pixelBuffer(totalPixels);
    
    int successfulFrames = 0;
    
    // Process each dark frame
    for (int i = 0; i < darkFrames.size(); ++i) {
        const QString &darkPath = darkFrames[i];
        fitsfile *darkFits = nullptr;
        
        QByteArray darkPathBytes = darkPath.toLocal8Bit();
        if (fits_open_file(&darkFits, darkPathBytes.data(), READONLY, &status)) {
            logMessage(QString("Failed to open dark frame %1: %2 (FITS error: %3)")
                          .arg(i+1).arg(QFileInfo(darkPath).fileName()).arg(status), "orange");
            status = 0; // Reset status to continue
            continue;
        }
        
        // Read pixel data
        long firstPixel = 1;
        if (fits_read_img(darkFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                         pixelBuffer.data(), nullptr, &status)) {
            logMessage(QString("Failed to read pixel data from frame %1: %2 (FITS error: %3)")
                          .arg(i+1).arg(QFileInfo(darkPath).fileName()).arg(status), "orange");
            fits_close_file(darkFits, &status);
            status = 0;
            continue;
        }
        
        // Add to accumulator
        for (long j = 0; j < totalPixels; ++j) {
            accumulator[j] += pixelBuffer[j];
        }
        
        successfulFrames++;
        fits_close_file(darkFits, &status);
        
        // Update progress
        if (i % 5 == 0 || i == darkFrames.size() - 1) {
            logMessage(QString("Processed dark frame %1 of %2").arg(i+1).arg(darkFrames.size()), "gray");
        }
    }
    
    if (successfulFrames == 0) {
        logMessage("No dark frames could be processed", "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputPath);
        return false;
    }
    
    logMessage(QString("Successfully read %1 of %2 dark frames, computing average...").arg(successfulFrames).arg(darkFrames.size()), "blue");
    
    // Average the accumulated values
    for (long i = 0; i < totalPixels; ++i) {
        pixelBuffer[i] = static_cast<float>(accumulator[i] / successfulFrames);
    }
    
    // Write averaged data to output
    long firstPixel = 1;
    if (fits_write_img(outputFits, TFLOAT, firstPixel, totalPixels, 
                      pixelBuffer.data(), &status)) {
        logMessage(QString("Failed to write master dark data (FITS error: %1)").arg(status), "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputPath);
        return false;
    }
    
    // Update header with processing info
    QString historyComment = QString("Master dark from %1 frames").arg(successfulFrames);
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(outputFits, historyBytes.data(), &status)) {
        // Non-critical error, just log it
        logMessage("Warning: Could not write processing history to header", "orange");
        status = 0;
    }
    
    // Add custom keyword for frame count
    if (fits_write_key(outputFits, TINT, "NFRAMES", &successfulFrames, "Number of frames averaged", &status)) {
        logMessage("Warning: Could not write frame count to header", "orange");
        status = 0;
    }
    
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error occurred during master dark creation (error: %1)").arg(status), "red");
        QFile::remove(outputPath);
        return false;
    }
    
    logMessage(QString("Master dark created successfully from %1 frames (averaged %2 frames)")
                  .arg(darkFrames.size()).arg(successfulFrames), "green");
    return true;
}

bool StellinaProcessor::applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame) {
    logMessage(QString("Applying master dark to %1").arg(QFileInfo(lightFrame).fileName()), "blue");
    
    // Add a small delay to prevent overwhelming Siril
    QThread::msleep(100);
    
    // Try loading the light frame with retry logic
    int maxRetries = 3;
    bool loadSuccess = false;
    
    for (int attempt = 1; attempt <= maxRetries; ++attempt) {
        if (m_sirilClient->loadImage(lightFrame)) {
            loadSuccess = true;
            break;
        } else {
            if (attempt < maxRetries) {
                logMessage(QString("Load attempt %1 failed, retrying...").arg(attempt), "orange");
                QThread::msleep(500); // Wait before retry
            }
        }
    }
    
    if (!loadSuccess) {
        logMessage("Failed to load light frame after retries", "red");
        return false;
    }
    
    // Try Siril's calibrate command first
    QString calibrateCmd = QString("calibrate_single \"%1\"").arg(masterDark);
    bool calibrateSuccess = m_sirilClient->sendSirilCommand(calibrateCmd);
    
    if (!calibrateSuccess) {
        logMessage("calibrate_single failed, using manual subtraction", "orange");
        
        // Manual subtraction method with retry
        QString loadDarkCmd = QString("load \"%1\" $masterdark").arg(masterDark);
        bool darkLoadSuccess = false;
        
        for (int attempt = 1; attempt <= 2; ++attempt) {
            if (m_sirilClient->sendSirilCommand(loadDarkCmd)) {
                darkLoadSuccess = true;
                break;
            } else if (attempt < 2) {
                logMessage("Failed to load master dark, retrying...", "orange");
                QThread::msleep(300);
            }
        }
        
        if (!darkLoadSuccess) {
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
    
    // Save calibrated frame with retry
    bool saveSuccess = false;
    for (int attempt = 1; attempt <= 2; ++attempt) {
        if (m_sirilClient->saveImage(outputFrame)) {
            saveSuccess = true;
            break;
        } else if (attempt < 2) {
            logMessage("Save attempt failed, retrying...", "orange");
            QThread::msleep(200);
        }
    }
    
    if (!saveSuccess) {
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
    int lightTemperatureK = extractTemperature(lightFrame); // Now in Kelvin
    QString lightBinning = extractBinning(lightFrame);
    
    if (lightExposure <= 0) {
        logMessage("Could not determine light frame exposure time", "red");
        return false;
    }
    
    // Find all matching dark frames
    QStringList matchingDarks = findAllMatchingDarkFrames(lightExposure, lightTemperatureK, lightBinning);
    
    if (matchingDarks.isEmpty()) {
        int temperatureC = lightTemperatureK - 273;
        logMessage(QString("No matching dark frames found for exposure=%1s, temp=%2K (%3°C), binning=%4")
                      .arg(lightExposure).arg(lightTemperatureK).arg(temperatureC).arg(lightBinning), "orange");
        m_skippedCount++;
        return true; // Continue processing without dark calibration
    }
    
    logMessage(QString("Found %1 matching dark frames").arg(matchingDarks.size()), "blue");
    
    // Create master dark for this combination using Kelvin in filename
    QString masterDarkName = QString("master_dark_%1s_%2K_%3.fits")
                                .arg(lightExposure)
                                .arg(lightTemperatureK)
                                .arg(lightBinning);
    QString masterDarkPath = QDir(m_calibratedDirectory).absoluteFilePath(masterDarkName);
    
    // Check if master dark already exists
    if (!QFile::exists(masterDarkPath)) {
        if (!createMasterDark(matchingDarks, masterDarkPath)) {
            logMessage("Failed to create master dark", "red");
            return false;
        }
        int temperatureC = lightTemperatureK - 273;
        logMessage(QString("Created master dark: %1 (from %2K/%3°C data)")
                      .arg(masterDarkName).arg(lightTemperatureK).arg(temperatureC), "green");
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
