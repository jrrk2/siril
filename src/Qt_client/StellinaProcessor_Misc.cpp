#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
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
#include <QPainter>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>  // for std::fmod if needed

// Add this validation function to verify acqTime preservation
void StellinaProcessor::validateAcqTimePreservation() {
    logMessage("=== VALIDATING ACQTIME PRESERVATION IN DARK CALIBRATED FILES ===", "blue");
    
    if (m_calibratedDirectory.isEmpty()) {
        logMessage("Please set calibrated directory first", "red");
        return;
    }
    
    QDir calibratedDir(m_calibratedDirectory);
    if (!calibratedDir.exists()) {
        logMessage("Calibrated directory does not exist", "red");
        return;
    }
    
    QStringList calibratedFiles = calibratedDir.entryList(
        QStringList() << "dark_calibrated_*.fits" << "*calibrated*.fits", 
        QDir::Files);
    
    if (calibratedFiles.isEmpty()) {
        logMessage("No calibrated files found for validation", "red");
        return;
    }
    
    logMessage("Validating acqTime preservation and boot time calibration:", "green");
    logMessage("File                     | acqTime    | RefAcqTime | RefUTC            | Reconstructed UTC     | Status", "green");
    logMessage("-------------------------|------------|------------|-------------------|-----------------------|-------", "gray");
    
    int validCount = 0;
    int totalCount = 0;
    
    for (const QString &file : calibratedFiles) {
        QString filePath = calibratedDir.absoluteFilePath(file);
        totalCount++;
        
        qint64 acqTime, refAcqTime;
        QDateTime refUTCTime;
        
        if (readAcqTimeFromFits(filePath, acqTime, refAcqTime, refUTCTime)) {
            // Reconstruct the precise UTC time
            QDateTime reconstructedUTC = reconstructUTCFromAcqTime(acqTime, refAcqTime, refUTCTime);
            
            // Check if the reconstruction is valid
            bool isValid = reconstructedUTC.isValid();
            if (isValid) {
                validCount++;
            }
            
            logMessage(QString("%1 | %2 | %3 | %4 | %5 | %6")
                      .arg(QFileInfo(file).fileName(), -24)
                      .arg(acqTime, -10)
                      .arg(refAcqTime, -10)
                      .arg(refUTCTime.toString("yyyy-MM-ddThh:mm:ss"), -17)
                      .arg(reconstructedUTC.toString("yyyy-MM-ddThh:mm:ss"), -21)
                      .arg(isValid ? "OK" : "ERROR"), 
                      isValid ? "green" : "red");
        } else {
            logMessage(QString("%1 | %2 | %3 | %4 | %5 | %6")
                      .arg(QFileInfo(file).fileName(), -24)
                      .arg("N/A", -10)
                      .arg("N/A", -10)
                      .arg("N/A", -17)
                      .arg("N/A", -21)
                      .arg("MISSING"), "orange");
        }
    }
    
    logMessage("", "gray");
    logMessage(QString("SUMMARY: %1/%2 files have valid acqTime preservation (%3%)")
              .arg(validCount)
              .arg(totalCount)
              .arg(validCount * 100.0 / totalCount, 0, 'f', 1), "blue");
    
    if (validCount == totalCount) {
        logMessage("✓ All calibrated files have properly preserved acqTime data", "green");
    } else if (validCount > 0) {
        logMessage("⚠ Some files are missing acqTime data - may have less precise timing", "orange");
    } else {
        logMessage("✗ No files have preserved acqTime data - timing will use DATE-OBS fallback", "red");
    }
    
    logMessage("", "gray");
    logMessage("VALIDATION CRITERIA:", "blue");
    logMessage("• acqTime should be preserved from original JSON files", "blue");
    logMessage("• RefAcqTime and RefUTC should match the session boot time calibration", "blue");
    logMessage("• Reconstructed UTC should be valid and consistent with original timing", "blue");
    logMessage("• Files with 'OK' status will have precise msec timing for plate solving", "blue");
}

// Additional utility function to compare timing accuracy
void StellinaProcessor::compareTimingAccuracy() {
    logMessage("=== COMPARING TIMING ACCURACY: acqTime vs DATE-OBS ===", "blue");
    
    if (m_calibratedDirectory.isEmpty()) {
        logMessage("Please set calibrated directory first", "red");
        return;
    }
    
    QDir calibratedDir(m_calibratedDirectory);
    if (!calibratedDir.exists()) {
        logMessage("Calibrated directory does not exist", "red");
        return;
    }
    
    QStringList calibratedFiles = calibratedDir.entryList(
        QStringList() << "dark_calibrated_*.fits" << "*calibrated*.fits", 
        QDir::Files);
    
    if (calibratedFiles.isEmpty()) {
        logMessage("No calibrated files found for comparison", "red");
        return;
    }
    
    logMessage("Comparing coordinate calculations using different timing methods:", "green");
    logMessage("File                     | acqTime RA  | acqTime Dec | DATE-OBS RA | DATE-OBS Dec | RA Diff | Dec Diff", "green");
    logMessage("-------------------------|-------------|-------------|-------------|--------------|---------|----------", "gray");
    
    int comparisonCount = 0;
    double totalRADiff = 0.0, totalDecDiff = 0.0;
    double maxRADiff = 0.0, maxDecDiff = 0.0;
    
    for (const QString &file : calibratedFiles) {
        QString filePath = calibratedDir.absoluteFilePath(file);
        
        // Get coordinates using acqTime (precise timing)
        double ra_acqTime, dec_acqTime;
        bool hasAcqTime = convertAltAzToRaDecFromCalibratedFits(filePath, ra_acqTime, dec_acqTime);
        
        // Get coordinates using DATE-OBS (fallback timing)
        double ra_dateobs, dec_dateobs;
        bool hasDateObs = false;
        
        fitsfile *fptr = nullptr;
        int status = 0;
        QByteArray pathBytes = filePath.toLocal8Bit();
        
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status) == 0) {
            double alt, az;
            char dateobs[FLEN_VALUE];
            
            if (fits_read_key(fptr, TDOUBLE, "STELLALT", &alt, nullptr, &status) == 0 &&
                fits_read_key(fptr, TDOUBLE, "STELLAZ", &az, nullptr, &status) == 0 &&
                fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) == 0) {
                
                QString dateObsStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
                hasDateObs = convertAltAzToRaDec(alt, az, dateObsStr, ra_dateobs, dec_dateobs);
            }
            fits_close_file(fptr, &status);
        }
        
        if (hasAcqTime && hasDateObs) {
            double raDiff = ra_acqTime - ra_dateobs;
            double decDiff = dec_acqTime - dec_dateobs;
            
            // Handle RA wraparound at 0/360 degrees
            if (raDiff > 180.0) raDiff -= 360.0;
            if (raDiff < -180.0) raDiff += 360.0;
            
            // Convert to arcseconds for display
            double raArcsec = raDiff * 3600.0;
            double decArcsec = decDiff * 3600.0;
            
            totalRADiff += fabs(raArcsec);
            totalDecDiff += fabs(decArcsec);
            maxRADiff = qMax(maxRADiff, fabs(raArcsec));
            maxDecDiff = qMax(maxDecDiff, fabs(decArcsec));
            comparisonCount++;
            
            QString status = (fabs(raArcsec) < 1.0 && fabs(decArcsec) < 1.0) ? "good" : 
                            (fabs(raArcsec) < 5.0 && fabs(decArcsec) < 5.0) ? "fair" : "poor";
            
            logMessage(QString("%1 | %2 | %3 | %4 | %5 | %6\" | %7\"")
                      .arg(QFileInfo(file).fileName(), -24)
                      .arg(ra_acqTime, 11, 'f', 6)
                      .arg(dec_acqTime, 11, 'f', 6)
                      .arg(ra_dateobs, 11, 'f', 6)
                      .arg(dec_dateobs, 12, 'f', 6)
                      .arg(raArcsec, 7, 'f', 2)
                      .arg(decArcsec, 8, 'f', 2),
                      status == "good" ? "green" : (status == "fair" ? "orange" : "red"));
        } else {
            logMessage(QString("%1 | %2 | %3 | %4 | %5 | %6 | %7")
                      .arg(QFileInfo(file).fileName(), -24)
                      .arg(hasAcqTime ? "OK" : "FAIL", -11)
                      .arg(hasAcqTime ? "OK" : "FAIL", -11)
                      .arg(hasDateObs ? "OK" : "FAIL", -11)
                      .arg(hasDateObs ? "OK" : "FAIL", -12)
                      .arg("N/A", -7)
                      .arg("N/A", -8), "orange");
        }
    }
    
    if (comparisonCount > 0) {
        double avgRADiff = totalRADiff / comparisonCount;
        double avgDecDiff = totalDecDiff / comparisonCount;
        
        logMessage("", "gray");
        logMessage(QString("TIMING ACCURACY SUMMARY (%1 files):").arg(comparisonCount), "blue");
        logMessage(QString("  Average RA difference:  %1 arcsec").arg(avgRADiff, 0, 'f', 2), "blue");
        logMessage(QString("  Average Dec difference: %1 arcsec").arg(avgDecDiff, 0, 'f', 2), "blue");
        logMessage(QString("  Maximum RA difference:  %1 arcsec").arg(maxRADiff, 0, 'f', 2), "blue");
        logMessage(QString("  Maximum Dec difference: %1 arcsec").arg(maxDecDiff, 0, 'f', 2), "blue");
        
        if (avgRADiff < 1.0 && avgDecDiff < 1.0) {
            logMessage("✓ Excellent timing accuracy - acqTime provides sub-arcsecond precision", "green");
        } else if (avgRADiff < 5.0 && avgDecDiff < 5.0) {
            logMessage("⚠ Good timing accuracy - small differences may be due to mount drift", "orange");
        } else {
            logMessage("✗ Poor timing accuracy - may indicate timing calibration issues", "red");
        }
    } else {
        logMessage("No files available for timing comparison", "red");
    }
}
// Add this function to your StellinaProcessor class
bool StellinaProcessor::readAcqTimeFromFits(const QString &fitsPath, qint64 &acqTime, qint64 &refAcqTime, QDateTime &refUTCTime) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        logMessage(QString("Failed to open FITS file for acqTime reading: %1 (status: %2)").arg(fitsPath).arg(status), "red");
        return false;
    }
    
    // Read the acqTime value
    if (fits_read_key(fptr, TLONGLONG, "STL-ACQ", &acqTime, nullptr, &status) != 0) {
        logMessage(QString("No acqTime found in FITS header: %1").arg(QFileInfo(fitsPath).fileName()), "orange");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read the reference acqTime for boot time calibration
    status = 0;
    if (fits_read_key(fptr, TLONGLONG, "STL-REF", &refAcqTime, nullptr, &status) != 0) {
        logMessage(QString("No reference acqTime found in FITS header: %1").arg(QFileInfo(fitsPath).fileName()), "orange");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read the reference UTC time string
    char refTimeStr[FLEN_VALUE];
    status = 0;
    if (fits_read_key(fptr, TSTRING, "STL-REFT", refTimeStr, nullptr, &status) != 0) {
        logMessage(QString("No reference UTC time found in FITS header: %1").arg(QFileInfo(fitsPath).fileName()), "orange");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Parse the reference UTC time
    QString refTimeString = QString::fromLatin1(refTimeStr).trimmed().remove('\'').remove('"');
    refUTCTime = QDateTime::fromString(refTimeString, Qt::ISODate);
    if (!refUTCTime.isValid()) {
        logMessage(QString("Invalid reference UTC time format: %1").arg(refTimeString), "orange");
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    
    if (m_debugMode) {
        logMessage(QString("Read timing data from %1: acqTime=%2, refAcqTime=%3, refUTC=%4")
                  .arg(QFileInfo(fitsPath).fileName())
                  .arg(acqTime)
                  .arg(refAcqTime)
                  .arg(refUTCTime.toString(Qt::ISODate)), "gray");
    }
    
    return true;
}

// Function to reconstruct precise UTC time from stored acqTime values
QDateTime StellinaProcessor::reconstructUTCFromAcqTime(qint64 acqTime, qint64 refAcqTime, const QDateTime &refUTCTime) {
    // Calculate time difference from reference
    qint64 timeDifferenceMs = acqTime - refAcqTime;
    
    // Add difference to reference time
    QDateTime result = refUTCTime.addMSecs(timeDifferenceMs);
    
    return result;
}

// Enhanced function to use acqTime for precise coordinate conversion in plate solving
bool StellinaProcessor::convertAltAzToRaDecFromCalibratedFits(const QString &calibratedFitsPath, double &ra, double &dec) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = calibratedFitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        logMessage(QString("Failed to open calibrated FITS file: %1").arg(calibratedFitsPath), "red");
        return false;
    }
    
    // Read Stellina Alt/Az coordinates
    double alt, az;
    if (fits_read_key(fptr, TDOUBLE, "STELLALT", &alt, nullptr, &status) != 0 ||
        fits_read_key(fptr, TDOUBLE, "STELLAZ", &az, nullptr, &status) != 0) {
        logMessage(QString("No Stellina Alt/Az coordinates in calibrated file: %1").arg(QFileInfo(calibratedFitsPath).fileName()), "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Try to read stored acqTime values for precise timing
    qint64 acqTime, refAcqTime;
    QDateTime refUTCTime;
    status = 0;
    
    if (readAcqTimeFromFits(calibratedFitsPath, acqTime, refAcqTime, refUTCTime)) {
        // Use precise acqTime-based timing for coordinate conversion
        QDateTime preciseUTC = reconstructUTCFromAcqTime(acqTime, refAcqTime, refUTCTime);
        
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
        
        // Use precise timing for coordinate conversion
        int year = preciseUTC.date().year();
        int month = preciseUTC.date().month();
        int day = preciseUTC.date().day();
        int hour = preciseUTC.time().hour();
        int minute = preciseUTC.time().minute();
        int second = preciseUTC.time().second();
        
        if (m_debugMode) {
            logMessage(QString("Using precise acqTime timing: %1 ms -> %2 UTC")
                      .arg(acqTime)
                      .arg(preciseUTC.toString("yyyy-MM-ddThh:mm:ss")), "blue");
        }
        
        // Perform coordinate conversion with precise timing
        auto [jd, ra2000, dec2000, raNow, decNow, lst, ha] = 
            CoordinateUtils::calculateRaDec(year, month, day, hour, minute, second,
                                           correctedAlt, correctedAz, observer_lat, observer_lon);
        
        // Apply systematic corrections if enabled
        if (m_mountTilt.enableDriftCorrection) {
            double minutesElapsed = refUTCTime.msecsTo(preciseUTC) / 60000.0;
            
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
            logMessage(QString("Precise coordinate conversion: Alt=%.4f°,Az=%.4f° -> RA=%.6f°,Dec=%.6f°")
                      .arg(correctedAlt, 0, 'f', 4)
                      .arg(correctedAz, 0, 'f', 4)
                      .arg(ra, 0, 'f', 6)
                      .arg(dec, 0, 'f', 6), "green");
        }
        
        return true;
    }
    
    fits_close_file(fptr, &status);
    
    // Fallback to DATE-OBS if acqTime is not available
    logMessage("acqTime not available, falling back to DATE-OBS timing", "orange");
    
    // Read DATE-OBS for fallback timing
    char dateobs[FLEN_VALUE];
    status = 0;
    if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) != 0) {
        logMessage(QString("No DATE-OBS found in calibrated file: %1").arg(QFileInfo(calibratedFitsPath).fileName()), "red");
        return false;
    }
    
    QString dateObsStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
    
    // Use the existing convertAltAzToRaDec function as fallback
    return convertAltAzToRaDec(alt, az, dateObsStr, ra, dec);
}

bool StellinaProcessor::writeStellinaMetadataWithCoordinates(const QString &fitsPath, const StellinaImageData &imageData) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        logMessage(QString("Failed to open FITS file for metadata writing: %1 (status: %2)").arg(fitsPath).arg(status), "red");
        return false;
    }
    
    // Clean existing Stellina keywords first
    fits_close_file(fptr, &status);
    cleanExistingStellinaKeywords(fitsPath);
    
    // Reopen for writing
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        logMessage(QString("Failed to reopen FITS file for metadata writing: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Write original Stellina Alt/Az coordinates
    double alt = imageData.altitude;
    double az = imageData.azimuth;
    
    if (fits_write_key(fptr, TDOUBLE, "STELLALT", &alt, "Stellina altitude (degrees)", &status) ||
        fits_write_key(fptr, TDOUBLE, "STELLAZ", &az, "Stellina azimuth (degrees)", &status)) {
        logMessage("Warning: Could not write Stellina Alt/Az coordinates to FITS header", "orange");
        status = 0; // Continue anyway
    }
    
    // NEW: Write acqTime from JSON - Extract and write the msec since boot timing
    qint64 acqTime = 0;
    if (!imageData.originalJsonPath.isEmpty() && QFile::exists(imageData.originalJsonPath)) {
        QFile jsonFile(imageData.originalJsonPath);
        if (jsonFile.open(QIODevice::ReadOnly)) {
            QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll());
            QJsonObject root = doc.object();
            
            if (root.contains("acqTime")) {
                acqTime = root["acqTime"].toVariant().toLongLong();
                
                // Write acqTime as a long long integer (milliseconds since telescope boot)
                if (fits_write_key(fptr, TLONGLONG, "STL-ACQ", &acqTime, "Stellina acqTime (ms since boot)", &status)) {
                    logMessage("Warning: Could not write acqTime to FITS header", "orange");
                    status = 0;
                } else {
                    if (m_debugMode) {
                        logMessage(QString("Wrote acqTime: %1 ms").arg(acqTime), "gray");
                    }
                }
                
                // Also write the boot time calibration reference values for reconstruction
                if (s_sessionTimingInitialized) {
                    qint64 refAcqTime = s_sessionReferenceAcqTime;
                    if (fits_write_key(fptr, TLONGLONG, "STL-REF", &refAcqTime, "Session reference acqTime (ms)", &status)) {
                        logMessage("Warning: Could not write reference acqTime", "orange");
                        status = 0;
                    }
                    
                    // Write reference UTC time as ISO string for reconstruction
                    QString refUTCStr = s_sessionReferenceTime.toString(Qt::ISODate);
                    QByteArray refUTCBytes = refUTCStr.toLocal8Bit();
                    char* refUTCPtr = refUTCBytes.data();
                    if (fits_write_key(fptr, TSTRING, "STL-REFT", &refUTCPtr, "Session reference UTC time", &status)) {
                        logMessage("Warning: Could not write reference UTC time", "orange");
                        status = 0;
                    }
                    
                    if (m_debugMode) {
                        logMessage(QString("Wrote timing calibration: ref_acqTime=%1, ref_UTC=%2")
                                  .arg(refAcqTime).arg(refUTCStr), "gray");
                    }
                }
            }
        }
    }
    
    // Convert Alt/Az to RA/DEC and write to header
    double ra, dec;
    if (imageData.hasValidCoordinates && !imageData.dateObs.isEmpty()) {
      double obslat, obslong, jd, lst, ha;
      if (convertAltAzToRaDecExt(imageData.altitude, imageData.azimuth, imageData.dateObs, ra, dec, obslat, obslong, jd, lst, ha)) {
	    if (fits_write_key(fptr, TDOUBLE, "OBSLAT", &obslat, "Observer latitude for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write observer latitude", "orange");
                    status = 0;
                }
	    if (fits_write_key(fptr, TDOUBLE, "OBSLONG", &obslong, "Observer longitude for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write observer longitude", "orange");
                    status = 0;
                }
	    if (fits_write_key(fptr, TDOUBLE, "JULIAN", &jd, "Julian date for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write Julian date", "orange");
                    status = 0;
                }
	    if (fits_write_key(fptr, TDOUBLE, "LST", &lst, "Local solar time for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write local solar time", "orange");
                    status = 0;
                }
	    if (fits_write_key(fptr, TDOUBLE, "HOURANG", &ha, "Hour angle for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write hour angle", "orange");
                    status = 0;
                }
            // Write calculated RA/DEC coordinates
            if (fits_write_key(fptr, TDOUBLE, "STELLRA", &ra, "Calculated RA from Alt/Az (degrees)", &status) ||
                fits_write_key(fptr, TDOUBLE, "STELLDEC", &dec, "Calculated Dec from Alt/Az (degrees)", &status)) {
                logMessage("Warning: Could not write calculated RA/DEC to FITS header", "orange");
                status = 0;
            } else {
                if (m_debugMode) {
                    logMessage(QString("Wrote calculated coordinates: RA=%1°, Dec=%2°")
                              .arg(ra, 0, 'f', 6).arg(dec, 0, 'f', 6), "gray");
                }
                
                // Write coordinate conversion metadata
                QString obsLocation = m_observerLocation;
                QByteArray locationBytes = obsLocation.toLocal8Bit();
                char* locationPtr = locationBytes.data();
                if (fits_write_key(fptr, TSTRING, "OBSLOC", &locationPtr, "Observer location for coordinate conversion", &status)) {
                    logMessage("Warning: Could not write observer location", "orange");
                    status = 0;
                }
                
                QString conversion_method = "Alt/Az to RA/Dec from Stellina mount position";
                QByteArray methodBytes = conversion_method.toLocal8Bit();
                char* methodPtr = methodBytes.data();
                if (fits_write_key(fptr, TSTRING, "COORDMET", &methodPtr, "Mount coordinate conversion method", &status)) {
                    logMessage("Warning: Could not write conversion method", "orange");
                    status = 0;
                }
            }
        } else {
            logMessage("Warning: Failed to convert Alt/Az to RA/DEC during metadata writing", "orange");
        }
    } else {
        logMessage("Warning: Cannot convert coordinates - missing valid Alt/Az or observation time", "orange");
    }
    
    // Write original file paths (relative names only for portability)
    QString origFitsRelative = QFileInfo(imageData.originalFitsPath).fileName();
    QString origJsonRelative = QFileInfo(imageData.originalJsonPath).fileName();
    
    QByteArray origFitsBytes = origFitsRelative.toLocal8Bit();
    QByteArray origJsonBytes = origJsonRelative.toLocal8Bit();
    char* origFitsPtr = origFitsBytes.data();
    char* origJsonPtr = origJsonBytes.data();
    
    if (fits_write_key(fptr, TSTRING, "STELLORIG", &origFitsPtr, "Original Stellina FITS file", &status) ||
        fits_write_key(fptr, TSTRING, "STELLJSON", &origJsonPtr, "Original Stellina JSON file", &status)) {
        logMessage("Warning: Could not write original file references", "orange");
        status = 0;
    }
    
    // Write exposure and temperature for easier matching
    int exposure = imageData.exposureSeconds;
    int temperature = imageData.temperatureKelvin;
    
    if (fits_write_key(fptr, TINT, "STELLEXP", &exposure, "Stellina exposure (seconds)", &status) ||
        fits_write_key(fptr, TINT, "STELLTEMP", &temperature, "Stellina temperature (Kelvin)", &status)) {
        logMessage("Warning: Could not write exposure/temperature metadata", "orange");
        status = 0;
    }
    
    // Write processing stage
    QString processingStage = "COORDINATES_CALCULATED";
    QByteArray stageBytes = processingStage.toLocal8Bit();
    char* stagePtr = stageBytes.data();
    
    if (fits_write_key(fptr, TSTRING, "STELLSTG", &stagePtr, "Stellina processing stage", &status)) {
        logMessage("Warning: Could not write processing stage", "orange");
        status = 0;
    }
    
    // Write processing timestamp
    QString timestamp = QDateTime::currentDateTime().toString(Qt::ISODate);
    QByteArray timestampBytes = timestamp.toLocal8Bit();
    char* timestampPtr = timestampBytes.data();
    
    if (fits_write_key(fptr, TSTRING, "STELLTS", &timestampPtr, "Stellina processing timestamp", &status)) {
        logMessage("Warning: Could not write timestamp", "orange");
        status = 0;
    }
    
    // Add comprehensive processing history including acqTime
    QString historyComment = QString("Stellina: Alt=%.2f°, Az=%.2f°, RA=%.6f°, Dec=%.6f°, Exp=%1s, acqTime=%2ms")
                                .arg(exposure)
                                .arg(imageData.altitude)
                                .arg(imageData.azimuth)
                                .arg(ra)
                                .arg(dec)
                                .arg(acqTime);
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(fptr, historyBytes.data(), &status)) {
        // Non-critical error
        logMessage("Warning: Could not write processing history", "orange");
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error during metadata writing (error: %1)").arg(status), "red");
        return false;
    }
    
    if (m_debugMode) {
        logMessage(QString("Wrote complete Stellina metadata with acqTime to: %1").arg(QFileInfo(fitsPath).fileName()), "gray");
    }
    
    return true;
}
