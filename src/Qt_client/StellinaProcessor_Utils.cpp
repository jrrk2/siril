// =============================================================================
// StellinaProcessor_ReversedImageUtils.cpp
// Utility functions for handling reversed Stellina images
// =============================================================================

#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QDebug>
#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>
#include <cmath>

// =============================================================================
// Bayer Pattern Detection and Handling
// =============================================================================

QString StellinaProcessor::detectBayerPattern(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        logMessage(QString("Warning: Could not open FITS file to detect bayer pattern: %1").arg(fitsPath), "orange");
        return "RGGB";  // Default assumption
    }
    
    // Try multiple possible bayer pattern keywords
    QStringList bayerKeywords = {"BAYERPAT", "COLORTYP", "CFATYPE", "FILTER", "BAYEROFF"};
    
    for (const QString &keyword : bayerKeywords) {
        char bayerValue[FLEN_VALUE];
        int keyStatus = 0;
        
        QByteArray keyBytes = keyword.toLocal8Bit();
        if (fits_read_key(fptr, TSTRING, keyBytes.data(), bayerValue, nullptr, &keyStatus) == 0) {
            QString pattern = QString::fromLatin1(bayerValue).trimmed().remove('\'').remove('"').toUpper();
            
            if (m_debugMode) {
                logMessage(QString("Found bayer keyword %1 = %2").arg(keyword).arg(pattern), "gray");
            }
            
            // Parse various formats
            if (pattern.contains("RGGB") || pattern == "RGGB") {
                fits_close_file(fptr, &status);
                return "RGGB";
            } else if (pattern.contains("GRBG") || pattern == "GRBG") {
                fits_close_file(fptr, &status);
                return "GRBG";
            } else if (pattern.contains("GBRG") || pattern == "GBRG") {
                fits_close_file(fptr, &status);
                return "GBRG";
            } else if (pattern.contains("BGGR") || pattern == "BGGR") {
                fits_close_file(fptr, &status);
                return "BGGR";
            }
        }
    }
    
    fits_close_file(fptr, &status);
    
    // If no bayer pattern found in header, use filename heuristics
    QString filename = QFileInfo(fitsPath).fileName();
    
    if (isReversedStellinaImage(filename)) {
        logMessage(QString("Detected reversed stellina image, assuming rotated bayer pattern"), "blue");
        return "BGGR";  // Reversed stellina images typically have rotated bayer pattern
    }
    
    // Default to RGGB for normal Stellina images
    return "RGGB";
}

bool StellinaProcessor::isReversedStellinaImage(const QString &filename) {
    QString baseName = QFileInfo(filename).baseName();
    
    // Check for patterns like "img-0001r", "image-0001r", etc.
    // Case insensitive matching
    QRegularExpression reversedPattern(R"(^(img|image)-\d+r$)", QRegularExpression::CaseInsensitiveOption);
    QRegularExpressionMatch match = reversedPattern.match(baseName);
    
    if (match.hasMatch()) {
        if (m_debugMode) {
            logMessage(QString("Identified reversed image pattern: %1").arg(filename), "gray");
        }
        return true;
    }
    
    // Also check for other possible reversed patterns
    QRegularExpression altPattern(R"(.*[_-]r(ev|eversed)?$)", QRegularExpression::CaseInsensitiveOption);
    return altPattern.match(baseName).hasMatch();
}

QString StellinaProcessor::getBaseName(const QString &filename) {
    QString baseName = QFileInfo(filename).baseName();
    
    if (isReversedStellinaImage(filename)) {
        // Remove the 'r' suffix: "img-0001r" -> "img-0001"
        if (baseName.endsWith('r', Qt::CaseInsensitive)) {
            return baseName.left(baseName.length() - 1);
        }
        // Handle other reversed patterns
        QRegularExpression reversedSuffix(R"([_-]r(ev|eversed)?$)", QRegularExpression::CaseInsensitiveOption);
        return baseName.replace(reversedSuffix, "");
    }
    
    return baseName;
}

QPoint StellinaProcessor::getBayerPatternOffset(const QString &pattern) {
    // Returns the offset of the Red pixel in the bayer pattern
    if (pattern == "RGGB") return QPoint(0, 0);  // R G
                                                  // G B
    if (pattern == "GRBG") return QPoint(1, 0);  // G R
                                                  // B G
    if (pattern == "GBRG") return QPoint(0, 1);  // G B
                                                  // R G
    if (pattern == "BGGR") return QPoint(1, 1);  // B G
                                                  // G R
    return QPoint(0, 0);  // Default to RGGB
}

bool StellinaProcessor::needsDarkRotation(const QString &lightBayerPattern, const QString &darkBayerPattern) {
    if (lightBayerPattern == darkBayerPattern) {
        return false;
    }
    
    // Log the mismatch
    if (m_debugMode) {
        logMessage(QString("Bayer pattern mismatch: Light=%1, Dark=%2").arg(lightBayerPattern).arg(darkBayerPattern), "gray");
    }
    
    return true;
}

// =============================================================================
// Dark Frame Rotation Implementation
// =============================================================================

bool StellinaProcessor::rotateDarkFrame(const QString &inputDarkPath, const QString &outputDarkPath, 
                                       const QString &fromPattern, const QString &toPattern) {
    
    logMessage(QString("Rotating dark frame from %1 to %2 pattern").arg(fromPattern).arg(toPattern), "blue");
    
    fitsfile *inputFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open input dark frame
    QByteArray inputPath = inputDarkPath.toLocal8Bit();
    if (fits_open_file(&inputFits, inputPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open input dark frame: %1 (FITS error: %2)").arg(inputDarkPath).arg(status), "red");
        return false;
    }
    
    // Get image dimensions
    long naxes[2];
    if (fits_get_img_size(inputFits, 2, naxes, &status)) {
        logMessage(QString("Failed to get dark frame dimensions (FITS error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    long width = naxes[0];
    long height = naxes[1];
    long totalPixels = width * height;
    
    if (m_debugMode) {
        logMessage(QString("Dark frame dimensions: %1 x %2").arg(width).arg(height), "gray");
    }
    
    // Read input image data
    std::vector<float> inputData(totalPixels);
    if (fits_read_img(inputFits, TFLOAT, 1, totalPixels, nullptr, inputData.data(), nullptr, &status)) {
        logMessage(QString("Failed to read input dark frame data (FITS error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    // Create output FITS file
    QByteArray outputPath = QString("!%1").arg(outputDarkPath).toLocal8Bit();
    if (fits_create_file(&outputFits, outputPath.data(), &status)) {
        logMessage(QString("Failed to create rotated dark frame: %1 (FITS error: %2)").arg(outputDarkPath).arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    // Create output image with same dimensions
    if (fits_create_img(outputFits, FLOAT_IMG, 2, naxes, &status)) {
        logMessage(QString("Failed to create output image structure (FITS error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    // Perform the rotation based on bayer pattern transformation
    std::vector<float> outputData(totalPixels);
    
    if (performBayerPatternRotation(inputData, outputData, width, height, fromPattern, toPattern)) {
        logMessage(QString("Successfully rotated bayer pattern %1 -> %2").arg(fromPattern).arg(toPattern), "green");
    } else {
        logMessage(QString("Warning: Unsupported bayer pattern transformation %1 -> %2, copying data").arg(fromPattern).arg(toPattern), "orange");
        outputData = inputData;
    }
    
    // Write rotated data
    if (fits_write_img(outputFits, TFLOAT, 1, totalPixels, outputData.data(), &status)) {
        logMessage(QString("Failed to write rotated dark frame data (FITS error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    // Copy header from input to output
    if (fits_copy_header(inputFits, outputFits, &status)) {
        logMessage("Warning: Failed to copy header to rotated dark frame", "orange");
        status = 0;
    }
    
    // Update bayer pattern in output header
    QByteArray patternBytes = toPattern.toLocal8Bit();
    char* patternPtr = patternBytes.data();
    if (fits_update_key(outputFits, TSTRING, "BAYERPAT", &patternPtr, "Rotated bayer pattern", &status)) {
        status = 0;
        fits_write_key(outputFits, TSTRING, "BAYERPAT", &patternPtr, "Rotated bayer pattern", &status);
    }
    
    // Add processing history
    QString historyComment = QString("HISTORY Rotated dark frame from %1 to %2 pattern for reversed stellina image")
                                .arg(fromPattern).arg(toPattern);
    QByteArray historyBytes = historyComment.toLocal8Bit();
    fits_write_history(outputFits, historyBytes.data(), &status);
    
    // Add rotation method comment
    QString methodComment = QString("HISTORY Rotation method: %1")
                               .arg(getRotationMethodDescription(fromPattern, toPattern));
    QByteArray methodBytes = methodComment.toLocal8Bit();
    fits_write_history(outputFits, methodBytes.data(), &status);
    
    fits_close_file(inputFits, &status);
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error during dark frame rotation: %1").arg(status), "red");
        QFile::remove(outputDarkPath);
        return false;
    }
    
    logMessage(QString("Successfully created rotated dark frame: %1").arg(QFileInfo(outputDarkPath).fileName()), "green");
    return true;
}

// =============================================================================
// Bayer Pattern Rotation Logic
// =============================================================================

bool StellinaProcessor::performBayerPatternRotation(const std::vector<float> &inputData, 
                                                   std::vector<float> &outputData,
                                                   long width, long height,
                                                   const QString &fromPattern, 
                                                   const QString &toPattern) {
    
    if (fromPattern == toPattern) {
        outputData = inputData;
        return true;
    }
    
    // Determine the rotation needed
    int rotationType = getBayerRotationType(fromPattern, toPattern);
    
    switch (rotationType) {
        case 0: // No rotation
            outputData = inputData;
            break;
            
        case 1: // 90° clockwise rotation
            return rotateBayerImage90CW(inputData, outputData, width, height);
            
        case 2: // 180° rotation
            return rotateBayerImage180(inputData, outputData, width, height);
            
        case 3: // 270° clockwise (90° counter-clockwise) rotation
            return rotateBayerImage270CW(inputData, outputData, width, height);
            
        case 4: // Horizontal flip
            return flipBayerImageHorizontal(inputData, outputData, width, height);
            
        case 5: // Vertical flip
            return flipBayerImageVertical(inputData, outputData, width, height);
            
        default:
            logMessage(QString("Unsupported bayer pattern transformation: %1 -> %2").arg(fromPattern).arg(toPattern), "orange");
            outputData = inputData;
            return false;
    }
    
    return true;
}

int StellinaProcessor::getBayerRotationType(const QString &fromPattern, const QString &toPattern) {
    // Define rotation mappings for common bayer pattern transformations
    // This is based on how the bayer pattern changes with image rotation
    
    struct PatternTransform {
        QString from;
        QString to;
        int rotation;
        QString description;
    };
    
    static const QVector<PatternTransform> transforms = {
        {"RGGB", "BGGR", 2, "180° rotation"},
        {"BGGR", "RGGB", 2, "180° rotation"},
        {"GRBG", "GBRG", 2, "180° rotation"},
        {"GBRG", "GRBG", 2, "180° rotation"},
        
        {"RGGB", "GBRG", 1, "90° clockwise"},
        {"GBRG", "BGGR", 1, "90° clockwise"},
        {"BGGR", "GRBG", 1, "90° clockwise"},
        {"GRBG", "RGGB", 1, "90° clockwise"},
        
        {"RGGB", "GRBG", 3, "270° clockwise"},
        {"GRBG", "BGGR", 3, "270° clockwise"},
        {"BGGR", "GBRG", 3, "270° clockwise"},
        {"GBRG", "RGGB", 3, "270° clockwise"},
        
        {"RGGB", "GRGB", 4, "Horizontal flip"},
        {"GRBG", "RGGB", 4, "Horizontal flip"},
        {"GBRG", "BGGR", 4, "Horizontal flip"},
        {"BGGR", "GBRG", 4, "Horizontal flip"},
        
        {"RGGB", "GBRG", 5, "Vertical flip"},
        {"GRBG", "GBRG", 5, "Vertical flip"},
        {"GBRG", "RGGB", 5, "Vertical flip"},
        {"BGGR", "GRBG", 5, "Vertical flip"}
    };
    
    for (const auto &transform : transforms) {
        if (transform.from == fromPattern && transform.to == toPattern) {
            if (m_debugMode) {
                logMessage(QString("Bayer transformation: %1").arg(transform.description), "gray");
            }
            return transform.rotation;
        }
    }
    
    return 0; // No rotation
}

QString StellinaProcessor::getRotationMethodDescription(const QString &fromPattern, const QString &toPattern) {
    int rotationType = getBayerRotationType(fromPattern, toPattern);
    
    switch (rotationType) {
        case 0: return "No rotation";
        case 1: return "90° clockwise rotation";
        case 2: return "180° rotation";
        case 3: return "270° clockwise rotation";
        case 4: return "Horizontal flip";
        case 5: return "Vertical flip";
        default: return "Unknown transformation";
    }
}

// =============================================================================
// Image Rotation Implementations
// =============================================================================

bool StellinaProcessor::rotateBayerImage180(const std::vector<float> &inputData, 
                                           std::vector<float> &outputData,
                                           long width, long height) {
    for (long y = 0; y < height; ++y) {
        for (long x = 0; x < width; ++x) {
            long srcIdx = y * width + x;
            long dstIdx = (height - 1 - y) * width + (width - 1 - x);
            outputData[dstIdx] = inputData[srcIdx];
        }
    }
    return true;
}

bool StellinaProcessor::rotateBayerImage90CW(const std::vector<float> &inputData, 
                                            std::vector<float> &outputData,
                                            long width, long height) {
    // For 90° clockwise rotation, output dimensions are swapped
    for (long y = 0; y < height; ++y) {
        for (long x = 0; x < width; ++x) {
            long srcIdx = y * width + x;
            long dstIdx = x * height + (height - 1 - y);
            outputData[dstIdx] = inputData[srcIdx];
        }
    }
    return true;
}

bool StellinaProcessor::rotateBayerImage270CW(const std::vector<float> &inputData, 
                                             std::vector<float> &outputData,
                                             long width, long height) {
    // For 270° clockwise rotation, output dimensions are swapped
    for (long y = 0; y < height; ++y) {
        for (long x = 0; x < width; ++x) {
            long srcIdx = y * width + x;
            long dstIdx = (width - 1 - x) * height + y;
            outputData[dstIdx] = inputData[srcIdx];
        }
    }
    return true;
}

bool StellinaProcessor::flipBayerImageHorizontal(const std::vector<float> &inputData, 
                                                std::vector<float> &outputData,
                                                long width, long height) {
    for (long y = 0; y < height; ++y) {
        for (long x = 0; x < width; ++x) {
            long srcIdx = y * width + x;
            long dstIdx = y * width + (width - 1 - x);
            outputData[dstIdx] = inputData[srcIdx];
        }
    }
    return true;
}

bool StellinaProcessor::flipBayerImageVertical(const std::vector<float> &inputData, 
                                              std::vector<float> &outputData,
                                              long width, long height) {
    for (long y = 0; y < height; ++y) {
        for (long x = 0; x < width; ++x) {
            long srcIdx = y * width + x;
            long dstIdx = (height - 1 - y) * width + x;
            outputData[dstIdx] = inputData[srcIdx];
        }
    }
    return true;
}

// =============================================================================
// Enhanced Dark Frame Matching
// =============================================================================

QString StellinaProcessor::getMatchingDarkKey(const StellinaImageData &imageData) {
    return QString("%1_%2_%3_%4")
        .arg(imageData.exposureSeconds)
        .arg(imageData.temperatureKelvin)
        .arg(imageData.binning)
        .arg(imageData.bayerPattern);
}

// Updated findAllMatchingDarkFrames to handle bayer pattern matching
QStringList StellinaProcessor::findAllMatchingDarkFrames(int targetExposure, int targetTemperature, 
                                                        const QString &targetBinning, 
                                                        const QString &targetBayerPattern) {
    QStringList matchingDarks;
    
    if (m_darkDirectory.isEmpty()) {
        if (m_debugMode) {
            logMessage("No dark directory specified", "gray");
        }
        return matchingDarks;
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        if (m_debugMode) {
            logMessage(QString("Dark directory does not exist: %1").arg(m_darkDirectory), "gray");
        }
        return matchingDarks;
    }
    
    QStringList darkFiles = darkDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    
    if (m_debugMode) {
        logMessage(QString("Searching %1 dark frames for match: %2s, %3K, %4, %5")
                      .arg(darkFiles.size())
                      .arg(targetExposure)
                      .arg(targetTemperature)
                      .arg(targetBinning)
                      .arg(targetBayerPattern), "gray");
    }
    
    int exactMatches = 0;
    int rotatedMatches = 0;
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperature = extractTemperature(fullPath);
        QString binning = extractBinning(fullPath);
        QString bayerPattern = detectBayerPattern(fullPath);
        
        // Check if this dark frame matches the light frame characteristics
        bool exposureMatch = qAbs(exposure - targetExposure) <= (targetExposure * m_exposureTolerance / 100);
        bool temperatureMatch = qAbs(temperature - targetTemperature) <= m_temperatureTolerance;
        bool binningMatch = (binning == targetBinning);
        
        if (exposureMatch && temperatureMatch && binningMatch) {
            if (bayerPattern == targetBayerPattern) {
                // Exact match - add directly
                matchingDarks.append(fullPath);
                exactMatches++;
                
                if (m_debugMode) {
                    logMessage(QString("Exact match: %1 (%2)").arg(QFileInfo(darkFile).fileName()).arg(bayerPattern), "gray");
                }
                
            } else if (needsDarkRotation(bayerPattern, targetBayerPattern)) {
                // Need rotation - create rotated version
                QString rotatedDarkName = QString("rotated_%1_%2_to_%3.%4")
                    .arg(QFileInfo(darkFile).baseName())
                    .arg(bayerPattern)
                    .arg(targetBayerPattern)
                    .arg(QFileInfo(darkFile).suffix());
                
                QString rotatedDarkPath = QDir(m_calibratedDirectory).absoluteFilePath(rotatedDarkName);
                
                if (!QFile::exists(rotatedDarkPath)) {
                    if (rotateDarkFrame(fullPath, rotatedDarkPath, bayerPattern, targetBayerPattern)) {
                        matchingDarks.append(rotatedDarkPath);
                        rotatedMatches++;
                        logMessage(QString("Created rotated dark frame: %1").arg(rotatedDarkName), "green");
                    } else {
                        logMessage(QString("Failed to rotate dark frame for pattern %1->%2")
                                      .arg(bayerPattern).arg(targetBayerPattern), "orange");
                    }
                } else {
                    matchingDarks.append(rotatedDarkPath);
                    rotatedMatches++;
                    
                    if (m_debugMode) {
                        logMessage(QString("Using existing rotated dark frame: %1").arg(rotatedDarkName), "gray");
                    }
                }
            }
        }
    }
    
    if (m_debugMode) {
        logMessage(QString("Found %1 matching darks (%2 exact, %3 rotated)")
                      .arg(matchingDarks.size()).arg(exactMatches).arg(rotatedMatches), "blue");
    }
    
    return matchingDarks;
}

// =============================================================================
// Enhanced Dark Frame Processing
// =============================================================================

void StellinaProcessor::loadDarkFrames() {
    m_darkFrames.clear();
    m_darkFramesTable->setRowCount(0);
    
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
    
    // Group dark frames by exposure, temperature, binning, and bayer pattern
    QMap<QString, QStringList> darkGroups;
    QMap<QString, DarkFrame> darkInfo;
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperatureK = extractTemperature(fullPath);
        QString binning = extractBinning(fullPath);
        QString bayerPattern = detectBayerPattern(fullPath);
        
        if (exposure > 0) {
            QString key = QString("%1_%2_%3_%4").arg(exposure).arg(temperatureK).arg(binning).arg(bayerPattern);
            darkGroups[key].append(fullPath);
            
            if (!darkInfo.contains(key)) {
                DarkFrame dark;
                dark.filepath = fullPath;
                dark.exposure = exposure;
                dark.temperature = temperatureK;
                dark.binning = binning;
                dark.bayerPattern = bayerPattern;
                darkInfo[key] = dark;
                
                if (m_debugMode) {
                    int temperatureC = temperatureK - 273;
                    logMessage(QString("Dark group %1: %2s, %3K (%4°C), %5, %6 - first file: %7")
                                  .arg(darkInfo.size())
                                  .arg(exposure)
                                  .arg(temperatureK)
                                  .arg(temperatureC)
                                  .arg(binning)
                                  .arg(bayerPattern)
                                  .arg(QFileInfo(darkFile).fileName()), "gray");
                }
            }
        } else {
            if (m_debugMode) {
                logMessage(QString("Skipped invalid dark frame: %1").arg(QFileInfo(darkFile).fileName()), "orange");
            }
        }
    }
    
    // Convert to list and populate table
    m_darkFrames.clear();
    for (auto it = darkInfo.begin(); it != darkInfo.end(); ++it) {
        m_darkFrames.append(it.value());
    }
    
    // Update table
    m_darkFramesTable->setRowCount(m_darkFrames.size());
    
    for (int i = 0; i < m_darkFrames.size(); ++i) {
        const DarkFrame &dark = m_darkFrames[i];
        QString key = QString("%1_%2_%3_%4").arg(dark.exposure).arg(dark.temperature).arg(dark.binning).arg(dark.bayerPattern);
        int count = darkGroups[key].size();
        int temperatureC = dark.temperature - 273;
        
        m_darkFramesTable->setItem(i, 0, new QTableWidgetItem(QFileInfo(dark.filepath).fileName()));
        m_darkFramesTable->setItem(i, 1, new QTableWidgetItem(QString::number(dark.exposure)));
        m_darkFramesTable->setItem(i, 2, new QTableWidgetItem(QString::number(temperatureC)));
        m_darkFramesTable->setItem(i, 3, new QTableWidgetItem(dark.binning));
        m_darkFramesTable->setItem(i, 4, new QTableWidgetItem(dark.bayerPattern));
        m_darkFramesTable->setItem(i, 5, new QTableWidgetItem(QString::number(count)));
    }
    
    // Update status
    if (m_darkFrames.isEmpty()) {
        m_darkFramesCount->setText("No valid dark frames found");
        logMessage("No valid dark frames found", "orange");
    } else {
        m_darkFramesCount->setText(QString("%1 dark frame groups found").arg(m_darkFrames.size()));
        
        // Create summary
        QStringList summary;
        QMap<QString, QStringList> darkGroups;
        for (const DarkFrame &dark : m_darkFrames) {
            QString key = QString("%1_%2_%3_%4").arg(dark.exposure).arg(dark.temperature).arg(dark.binning).arg(dark.bayerPattern);
            int count = darkGroups[key].size();
            int temperatureC = dark.temperature - 273;
            summary.append(QString("%1×%2s@%3K(%4)")
                              .arg(count)
                              .arg(dark.exposure)
                              .arg(dark.temperature)
                              .arg(dark.bayerPattern));
        }
        logMessage(QString("Found: %1").arg(summary.join(", ")), "blue");
        
        // Log bayer pattern distribution
        QMap<QString, int> patternCounts;
        for (const DarkFrame &dark : m_darkFrames) {
            patternCounts[dark.bayerPattern]++;
        }
        
        QStringList patternSummary;
        for (auto it = patternCounts.begin(); it != patternCounts.end(); ++it) {
            patternSummary.append(QString("%1: %2").arg(it.key()).arg(it.value()));
        }
        logMessage(QString("Bayer patterns: %1").arg(patternSummary.join(", ")), "blue");
    }
}


// Static member definitions
QDateTime StellinaProcessor::s_sessionReferenceTime;
qint64 StellinaProcessor::s_sessionReferenceAcqTime = 0;
bool StellinaProcessor::s_sessionTimingInitialized = false;

// Method implementations
void StellinaProcessor::resetSessionTiming() {
    s_sessionReferenceTime = QDateTime();
    s_sessionReferenceAcqTime = 0;
    s_sessionTimingInitialized = false;
}

bool StellinaProcessor::initializeSessionTiming(const QString &sourceDirectory) {
    QDir sourceDir(sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage("Source directory does not exist", "red");
        return false;
    }
    
    // Find first image pair (FITS + JSON)
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
    fitsFiles.sort();
    
    if (fitsFiles.isEmpty()) {
        logMessage("No FITS files found in source directory", "red");
        return false;
    }
    
    // Try each FITS file until we find one with both FITS and JSON
    for (const QString &fitsFile : fitsFiles) {
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Look for corresponding JSON file
        QStringList jsonCandidates = {
            baseName + "-stacking.json",
            baseName + ".json",
            baseName + "-stacking.JSON",
            baseName + ".JSON"
        };
        
        QString jsonPath;
        for (const QString &candidate : jsonCandidates) {
            QString candidatePath = sourceDir.absoluteFilePath(candidate);
            if (QFile::exists(candidatePath)) {
                jsonPath = candidatePath;
                break;
            }
        }
        
        if (jsonPath.isEmpty()) {
            continue; // Try next FITS file
        }
        
        // Extract DATE-OBS from FITS
        QString dateObs = extractDateObs(fitsPath);
        if (dateObs.isEmpty()) {
            continue;
        }
        
        // Extract acqTime from JSON
        QFile jsonFile(jsonPath);
        if (!jsonFile.open(QIODevice::ReadOnly)) {
            continue;
        }
        
        QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll());
        QJsonObject root = doc.object();
        
        if (!root.contains("acqTime")) {
            continue;
        }
        
        qint64 acqTime = root["acqTime"].toVariant().toLongLong();
        
        // Parse DATE-OBS to QDateTime
        QDateTime dateObsTime = QDateTime::fromString(dateObs, "yyyy-MM-ddThh:mm:ss");
        if (!dateObsTime.isValid()) {
            // Try alternative formats
            QStringList formats = {
                "yyyy-MM-ddThh:mm:ss.zzz",
                "yyyy-MM-ddThh:mm:ss.zzzZ",
                "yyyy-MM-dd hh:mm:ss"
            };
            
            for (const QString &format : formats) {
                dateObsTime = QDateTime::fromString(dateObs, format);
                if (dateObsTime.isValid()) break;
            }
        }
        
        if (!dateObsTime.isValid()) {
            continue;
        }
        
        dateObsTime.setTimeSpec(Qt::UTC);
        
        // Successfully found reference pair - initialize session timing
        s_sessionReferenceTime = dateObsTime;
        s_sessionReferenceAcqTime = acqTime;
        s_sessionTimingInitialized = true;
        
        logMessage(QString("Session timing initialized from %1:").arg(baseName), "blue");
        logMessage(QString("  DATE-OBS: %1").arg(dateObsTime.toString("yyyy-MM-ddThh:mm:ss")), "blue");
        logMessage(QString("  acqTime:  %1").arg(acqTime), "blue");
        
        return true;
    }
    
    logMessage("Could not find valid FITS+JSON pair for timing initialization", "red");
    return false;
}

QDateTime StellinaProcessor::convertAcqTimeToUTC(qint64 acqTime) {
    if (!s_sessionTimingInitialized) {
        logMessage("Session timing not initialized - call initializeSessionTiming() first", "red");
        return QDateTime();
    }
    
    // Calculate time difference from reference
    qint64 timeDifferenceMs = acqTime - s_sessionReferenceAcqTime;
    
    // Add difference to reference time
    QDateTime result = s_sessionReferenceTime.addMSecs(timeDifferenceMs);
    
    return result;
}

bool StellinaProcessor::convertAltAzToRaDecWithPreciseTiming(double alt, double az, const QString &jsonPath, double &ra, double &dec) {
    // Initialize session timing if not already done
    if (!s_sessionTimingInitialized) {
        QString sourceDir = QFileInfo(jsonPath).dir().absolutePath();
        if (!initializeSessionTiming(sourceDir)) {
            logMessage("Failed to initialize session timing", "red");
            return false;
        }
    }
    
    // Parse JSON to get acqTime
    QFile file(jsonPath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (m_debugMode) {
            logMessage(QString("Failed to open JSON file: %1").arg(jsonPath), "red");
        }
        return false;
    }
    
    QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
    QJsonObject root = doc.object();
    
    if (!root.contains("acqTime")) {
        if (m_debugMode) {
            logMessage(QString("No acqTime in JSON: %1").arg(jsonPath), "red");
        }
        return false;
    }
    
    // Get precise timing from acqTime using calculated offset
    qint64 acqTime = root["acqTime"].toVariant().toLongLong();
    QDateTime preciseUTC = convertAcqTimeToUTC(acqTime);
    
    if (!preciseUTC.isValid()) {
        return false;
    }
    
    // Apply mount corrections
    double correctedAlt, correctedAz;
    applyMountTiltCorrection(correctedAlt, correctedAz, alt, az);
    
    // Get observer location
    QStringList locationParts = m_observerLocation.split(',');
    double observer_lat = 51.5074;  // Default to London
    double observer_lon = -0.1278;
    
    if (locationParts.size() >= 2) {
        bool ok1, ok2;
        double lat = locationParts[0].trimmed().toDouble(&ok1);
        double lon = locationParts[1].trimmed().toDouble(&ok2);
        if (ok1 && ok2) {
            observer_lat = lat;
            observer_lon = lon;
        }
    }
    
    // Use precise acqTime-derived timestamp for coordinate conversion
    int year = preciseUTC.date().year();
    int month = preciseUTC.date().month();
    int day = preciseUTC.date().day();
    int hour = preciseUTC.time().hour();
    int minute = preciseUTC.time().minute();
    int second = preciseUTC.time().second();
    
    if (m_debugMode) {
        logMessage(QString("Precise timing: acqTime=%1 -> %2 UTC")
                  .arg(acqTime)
                  .arg(preciseUTC.toString("yyyy-MM-ddThh:mm:ss")), "blue");
    }
    
    // Perform coordinate conversion with precise timing
    auto [jd, ra2000, dec2000, raNow, decNow, lst, ha] = 
        CoordinateUtils::calculateRaDec(year, month, day, hour, minute, second,
                                       correctedAlt, correctedAz, observer_lat, observer_lon);
    
    // Apply systematic corrections if enabled
    if (m_mountTilt.enableDriftCorrection) {
        double minutesElapsed = s_sessionReferenceTime.msecsTo(preciseUTC) / 60000.0;
        
        double raOffset = m_mountTilt.driftRA * minutesElapsed / 3600.0;
        double decOffset = m_mountTilt.driftDec * minutesElapsed / 3600.0;
        
        ra2000 += raOffset;
        dec2000 += decOffset;
        
        if (m_debugMode) {
            logMessage(QString("Applied time-based corrections: RA+=%1°, Dec+=%2° at %.1f min")
                      .arg(raOffset, 0, 'f', 6)
                      .arg(decOffset, 0, 'f', 6)
                      .arg(minutesElapsed), "gray");
        }
    }
    
    ra = ra2000;
    dec = dec2000;
    
    if (m_debugMode) {
        logMessage(QString("Precise conversion: Alt=%.4f°,Az=%.4f° -> RA=%.6f°,Dec=%.6f°")
                  .arg(correctedAlt, 0, 'f', 4)
                  .arg(correctedAz, 0, 'f', 4)
                  .arg(ra, 0, 'f', 6)
                  .arg(dec, 0, 'f', 6), "green");
    }
    
    return true;
}

void StellinaProcessor::validateTimingOffset() {
    logMessage("=== VALIDATING TIMING OFFSET CALCULATION ===", "blue");
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory first", "red");
        return;
    }
    
    // Reset and recalculate timing
    resetSessionTiming();
    
    if (!initializeSessionTiming(m_sourceDirectory)) {
        logMessage("Failed to initialize session timing", "red");
        return;
    }
    
    QDir sourceDir(m_sourceDirectory);
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
    fitsFiles.sort();
    
    logMessage("Validating acqTime conversion against DATE-OBS:", "green");
    logMessage("Image       | acqTime    | Converted UTC     | FITS DATE-OBS     | Diff (ms)", "green");
    logMessage("------------|------------|-------------------|-------------------|---------", "gray");
    
    for (int i = 0; i < qMin(fitsFiles.size(), 10); i++) {
        QString fitsFile = fitsFiles[i];
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Find corresponding JSON
        QString jsonPath = sourceDir.absoluteFilePath(baseName + "-stacking.json");
        if (!QFile::exists(jsonPath)) {
            jsonPath = sourceDir.absoluteFilePath(baseName + ".json");
        }
        if (!QFile::exists(jsonPath)) continue;
        
        // Get acqTime from JSON
        QFile jsonFile(jsonPath);
        if (!jsonFile.open(QIODevice::ReadOnly)) continue;
        
        QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll());
        QJsonObject root = doc.object();
        
        if (!root.contains("acqTime")) continue;
        
        qint64 acqTime = root["acqTime"].toVariant().toLongLong();
        
        // Convert acqTime to UTC
        QDateTime convertedUTC = convertAcqTimeToUTC(acqTime);
        
        // Get DATE-OBS from FITS
        QString dateObs = extractDateObs(fitsPath);
        QDateTime fitsUTC = QDateTime::fromString(dateObs, "yyyy-MM-ddThh:mm:ss");
        fitsUTC.setTimeSpec(Qt::UTC);
        
        // Calculate difference
        qint64 diffMs = fitsUTC.msecsTo(convertedUTC);
        
        logMessage(QString("%1 | %2 | %3 | %4 | %5")
                  .arg(baseName, -11)
                  .arg(acqTime, -10)
                  .arg(convertedUTC.toString("yyyy-MM-ddThh:mm:ss"), -17)
                  .arg(fitsUTC.toString("yyyy-MM-ddThh:mm:ss"), -17)
                  .arg(diffMs, 8), 
                  qAbs(diffMs) < 1000 ? "green" : (qAbs(diffMs) < 5000 ? "orange" : "red"));
    }
    
    logMessage("", "gray");
    logMessage("VALIDATION CRITERIA:", "blue");
    logMessage("• Differences should be minimal (< 1000ms) for most images", "blue");
    logMessage("• Large differences may indicate timing issues in FITS headers", "blue");
    logMessage("• Consistent small differences are acceptable", "blue");
}

void StellinaProcessor::compareDriftWithDynamicOffset() {
    logMessage("=== DRIFT COMPARISON WITH DYNAMIC OFFSET ===", "blue");
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory first", "red");
        return;
    }
    
    // Reset and initialize timing
    resetSessionTiming();
    
    if (!initializeSessionTiming(m_sourceDirectory)) {
        logMessage("Failed to initialize session timing", "red");
        return;
    }
    
    QDir sourceDir(m_sourceDirectory);
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
    fitsFiles.sort();
    
    logMessage("Comparing coordinate drift between timing methods:", "green");
    logMessage("Image      | Method    | Minutes | RA (degrees) | RA Drift | Dec (degrees) | Dec Drift", "green");
    logMessage("-----------|-----------|---------|--------------|----------|---------------|----------", "gray");
    
    double firstRA_DateObs = 0, firstDec_DateObs = 0;
    double firstRA_AcqTime = 0, firstDec_AcqTime = 0;
    QDateTime startTime = s_sessionReferenceTime;
    
    for (int i = 0; i < qMin(fitsFiles.size(), 10); i++) {
        QString fitsFile = fitsFiles[i];
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Find corresponding JSON
        QString jsonPath = sourceDir.absoluteFilePath(baseName + "-stacking.json");
        if (!QFile::exists(jsonPath)) {
            jsonPath = sourceDir.absoluteFilePath(baseName + ".json");
        }
        if (!QFile::exists(jsonPath)) continue;
        
        // Get coordinates from JSON
        QJsonObject metadata = loadStellinaJson(jsonPath);
        double alt, az;
        if (!extractCoordinates(metadata, alt, az)) continue;
        
        // METHOD 1: Using DATE-OBS (current method with drift)
        QString dateObs = extractDateObs(fitsPath);
        double ra_dateobs, dec_dateobs;
        bool dateObsSuccess = convertAltAzToRaDec(alt, az, dateObs, ra_dateobs, dec_dateobs);
        
        // METHOD 2: Using acqTime with dynamic offset
        double ra_acqtime, dec_acqtime;
        bool acqTimeSuccess = convertAltAzToRaDecWithPreciseTiming(alt, az, jsonPath, ra_acqtime, dec_acqtime);
        
        if (!dateObsSuccess || !acqTimeSuccess) continue;
        
        // Calculate elapsed time
        QFile jsonFile(jsonPath);
        jsonFile.open(QIODevice::ReadOnly);
        QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll());
        qint64 acqTime = doc.object()["acqTime"].toVariant().toLongLong();
        QDateTime currentTime = convertAcqTimeToUTC(acqTime);
        
        double minutesElapsed = startTime.msecsTo(currentTime) / 60000.0;
        
        // Calculate drift from first image
        if (i == 0) {
            firstRA_DateObs = ra_dateobs;
            firstDec_DateObs = dec_dateobs;
            firstRA_AcqTime = ra_acqtime;
            firstDec_AcqTime = dec_acqtime;
        }
        
        double raDrift_DateObs = ra_dateobs - firstRA_DateObs;
        double decDrift_DateObs = dec_dateobs - firstDec_DateObs;
        double raDrift_AcqTime = ra_acqtime - firstRA_AcqTime;
        double decDrift_AcqTime = dec_acqtime - firstDec_AcqTime;
        
        // Display results
        logMessage(QString("%1 | DATE-OBS  | %2 | %3 | %4 | %5 | %6")
                  .arg(baseName, -10)
                  .arg(minutesElapsed, 7, 'f', 2)
                  .arg(ra_dateobs, 12, 'f', 6)
                  .arg(raDrift_DateObs, 8, 'f', 4)
                  .arg(dec_dateobs, 13, 'f', 6)
                  .arg(decDrift_DateObs, 8, 'f', 4), "orange");
        
        logMessage(QString("%1 | acqTime   | %2 | %3 | %4 | %5 | %6")
                  .arg(baseName, -10)
                  .arg(minutesElapsed, 7, 'f', 2)
                  .arg(ra_acqtime, 12, 'f', 6)
                  .arg(raDrift_AcqTime, 8, 'f', 4)
                  .arg(dec_acqtime, 13, 'f', 6)
                  .arg(decDrift_AcqTime, 8, 'f', 4), "green");
        
        logMessage("", "gray");
    }
    
    logMessage("=== DRIFT COMPARISON COMPLETE ===", "blue");
    logMessage("The acqTime method should show significantly less drift.", "green");
}

void StellinaProcessor::updateProcessingToDynamicOffset() {
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory first", "red");
        return;
    }
    
    // Initialize session timing
    resetSessionTiming();
    if (!initializeSessionTiming(m_sourceDirectory)) {
        logMessage("Failed to initialize session timing for processing", "red");
        return;
    }
    
    logMessage("=== PROCESSING WITH DYNAMIC OFFSET acqTime ===", "blue");
    
    // Process all images with precise timing
    for (auto &imageData : m_stellinaImageData) {
        if (imageData.hasValidCoordinates && !imageData.originalJsonPath.isEmpty()) {
            
            // Use precise acqTime-based conversion with dynamic offset
            double preciseRA, preciseDec;
            bool success = convertAltAzToRaDecWithPreciseTiming(
                imageData.altitude, 
                imageData.azimuth, 
                imageData.originalJsonPath, 
                preciseRA, 
                preciseDec
            );
            
            if (success) {
                imageData.calculatedRA = preciseRA;
                imageData.calculatedDec = preciseDec;
                imageData.hasCalculatedCoords = true;
                
                if (m_debugMode) {
                    logMessage(QString("Image %1: Precise RA=%.6f°, Dec=%.6f°")
                              .arg(QFileInfo(imageData.originalFitsPath).baseName())
                              .arg(preciseRA, 0, 'f', 6)
                              .arg(preciseDec, 0, 'f', 6), "green");
                }
            }
        }
    }
}
