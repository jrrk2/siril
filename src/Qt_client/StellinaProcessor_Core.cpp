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
#include <limits>
#include <algorithm>  // for std::fmod if needed

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
    testLibnovaConversion();
 
    logMessage("Enhanced Stellina Processor started. Connect to Siril and select processing mode.", "blue");
    
    // Scan for dark frames if directory is set
    if (!m_darkDirectory.isEmpty()) {
        scanDarkFrames();
    }
}

StellinaProcessor::~StellinaProcessor() {
    saveSettings();
}

// Add these helper functions to StellinaProcessor_Core.cpp
// (based on your working lstest.cpp code)

// Constants
const double PI = 3.14159265358979323846;
const double DEG_TO_RAD = PI / 180.0;
const double RAD_TO_DEG = 180.0 / PI;

// Function to calculate Julian Date (JD) from a given Gregorian date
double StellinaProcessor::calculateJD(int year, int month, int day, int hour, int minute, int second) {
    int A = (14 - month) / 12;
    int Y = year + 4800 - A;
    int M = month + 12 * A - 3;

    double JD = day + ((153 * M + 2) / 5) + 365 * Y + Y / 4 - Y / 100 + Y / 400 - 32045;
    JD += hour / 24.0 + minute / 1440.0 + second / 86400.0;
    return JD;
}

// Function to calculate Local Sidereal Time (LST) from Julian Date (JD) and observer's longitude
double StellinaProcessor::calculateLST(double JD, double longitude) {
    // Calculate Julian Century (T)
    double T = (JD - 2451545.0) / 36525.0;

    // Calculate Greenwich Sidereal Time (GST)
    double GST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + T * T * (0.000387933 - T / 38710000.0);
    GST = fmod(GST, 360.0);  // Ensure GST is between 0 and 360 degrees

    // Convert to Local Sidereal Time (LST) by adding observer's longitude
    double LST = GST + longitude;  // Longitude in degrees
    LST = fmod(LST, 360.0);  // Ensure LST is between 0 and 360 degrees

    // Convert LST to hours
    LST /= 15.0;  // 15 degrees = 1 hour

    // Ensure LST is between 0 and 24 hours
    if (LST < 0) LST += 24.0;
    return LST;
}

// Function to convert ALT/AZ to RA/DEC (from your working code)
void StellinaProcessor::altAzToRaDec(double alt, double az, double lat, double lst, double &ra, double &dec) {
    // Convert Alt/Az to radians
    alt *= DEG_TO_RAD;
    az *= DEG_TO_RAD;
    lat *= DEG_TO_RAD;

    // Calculate Declination
    dec = asin(sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az));

    // Calculate Hour Angle
    double H = atan2(sin(az), cos(az) * sin(lat) - tan(dec) * cos(lat));

    // Convert back to degrees
    ra = lst * 15.0 - H * RAD_TO_DEG;
    dec *= RAD_TO_DEG;

    // Ensure RA is between 0 and 360 degrees
    if (ra < 0) ra += 360;
    if (ra >= 360) ra -= 360;
}

// Replace the convertAltAzToRaDec function with this proven implementation
bool StellinaProcessor::convertAltAzToRaDec(double alt, double az, const QString &dateObs, double &ra, double &dec) {
    // Parse observer location from settings (default to London if not set)
    QStringList locationParts = m_observerLocation.split(',');
    double observer_lat = 51.5074;  // London latitude (degrees)
    double observer_lon = -0.1278;  // London longitude (degrees)
    
    // Try to parse observer location if it's in "lat,lon" format
    if (locationParts.size() >= 2) {
        bool ok1, ok2;
        double lat = locationParts[0].trimmed().toDouble(&ok1);
        double lon = locationParts[1].trimmed().toDouble(&ok2);
        if (ok1 && ok2) {
            observer_lat = lat;
            observer_lon = lon;
        }
    }
    
    if (m_debugMode) {
        logMessage(QString("Using observer location: lat=%1°, lon=%2°")
                      .arg(observer_lat, 0, 'f', 4)
                      .arg(observer_lon, 0, 'f', 4), "gray");
    }
    
    // Parse observation time from FITS DATE-OBS header
    if (m_debugMode) {
        logMessage(QString("Raw dateObs parameter: '%1' (length: %2)").arg(dateObs).arg(dateObs.length()), "gray");
    }
    
    QDateTime obsTime;
    if (!dateObs.isEmpty()) {
        // Try different date formats that might be in FITS headers
        QStringList formats = {
            "yyyy-MM-ddThh:mm:ss.zzz",
            "yyyy-MM-ddThh:mm:ss.zzzZ",
            "yyyy-MM-ddThh:mm:ss",
            "yyyy-MM-ddThh:mm:ssZ",
            "yyyy-MM-dd hh:mm:ss.zzz",
            "yyyy-MM-dd hh:mm:ss",
            "yyyy-MM-dd"
        };
        
        for (const QString &format : formats) {
            obsTime = QDateTime::fromString(dateObs, format);
            if (obsTime.isValid()) {
                // Ensure we're working in UTC
                obsTime.setTimeSpec(Qt::UTC);
                break;
            }
        }
    }
    
    if (!obsTime.isValid()) {
        logMessage(QString("ERROR: Could not parse observation time '%1' - coordinate conversion requires accurate time").arg(dateObs), "red");
        logMessage("Supported time formats: yyyy-MM-ddThh:mm:ss.zzz, yyyy-MM-ddThh:mm:ss, yyyy-MM-dd hh:mm:ss", "red");
        return false;
    }
    
    if (m_debugMode) {
        logMessage(QString("Observation time: %1").arg(obsTime.toString(Qt::ISODate)), "gray");
    }
    
    // Calculate Julian Date using proven method
    double jd = calculateJD(obsTime.date().year(),
                           obsTime.date().month(),
                           obsTime.date().day(),
                           obsTime.time().hour(),
                           obsTime.time().minute(),
                           obsTime.time().second());
    
    if (m_debugMode) {
        logMessage(QString("Julian Day: %1").arg(jd, 0, 'f', 6), "gray");
    }
    
    // Calculate Local Sidereal Time using proven method
    double lst = calculateLST(jd, observer_lon);
    
    if (m_debugMode) {
        logMessage(QString("Local Sidereal Time: %1 hours (%2°)")
                      .arg(lst, 0, 'f', 4)
                      .arg(lst * 15.0, 0, 'f', 2), "gray");
    }
    
    // Convert Alt/Az to RA/Dec using proven method
    altAzToRaDec(alt, az, observer_lat, lst, ra, dec);
    
    if (m_debugMode) {
        logMessage(QString("Coordinate conversion using proven algorithm:"), "blue");
        logMessage(QString("  Input Alt/Az: %1°, %2°").arg(alt, 0, 'f', 4).arg(az, 0, 'f', 4), "blue");
        logMessage(QString("  Observer: %1°N, %2°E").arg(observer_lat, 0, 'f', 4).arg(observer_lon, 0, 'f', 4), "blue");
        logMessage(QString("  Time: %1 (JD %2)").arg(obsTime.toString(Qt::ISODate)).arg(jd, 0, 'f', 6), "blue");
        logMessage(QString("  LST: %1 hours").arg(lst, 0, 'f', 4), "blue");
        logMessage(QString("  Result RA/Dec: %1°, %2°").arg(ra, 0, 'f', 6).arg(dec, 0, 'f', 6), "blue");
        
        // Also show in hours:minutes:seconds format for RA
        double ra_hours = ra / 15.0;  // Convert degrees to hours
        int h = static_cast<int>(ra_hours);
        int m = static_cast<int>((ra_hours - h) * 60);
        double s = ((ra_hours - h) * 60 - m) * 60;
        logMessage(QString("  RA in HMS: %1h %2m %3s").arg(h).arg(m).arg(s, 0, 'f', 2), "blue");
        
        // Show declination in degrees:arcminutes:arcseconds
        int d = static_cast<int>(dec);
        int am = static_cast<int>(qAbs(dec - d) * 60);
        double as = (qAbs(dec - d) * 60 - am) * 60;
        logMessage(QString("  Dec in DMS: %1° %2' %3\"").arg(d).arg(am).arg(as, 0, 'f', 1), "blue");
    }
    
    // Sanity check the results
    if (ra < 0 || ra >= 360.0) {
        logMessage(QString("Warning: RA out of range: %1°").arg(ra), "orange");
    }
    
    if (dec < -90.0 || dec > 90.0) {
        logMessage(QString("Warning: Dec out of range: %1°").arg(dec), "orange");
    }
    
    return true;
}

// Add these function declarations to StellinaProcessor.h in the private section:
/*
private:
    double calculateJD(int year, int month, int day, int hour, int minute, int second);
    double calculateLST(double JD, double longitude);
    void altAzToRaDec(double alt, double az, double lat, double lst, double &ra, double &dec);
*/

QString StellinaProcessor::extractDateObs(const QString &fitsFile) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsFile.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        if (m_debugMode) {
            logMessage(QString("Failed to open FITS file for DATE-OBS extraction: %1 (status: %2)").arg(fitsFile).arg(status), "red");
        }
        return QString();
    }
    
    // Try to read DATE-OBS keyword
    char dateobs[FLEN_VALUE];  // FLEN_VALUE is typically 71 characters
    char comment[FLEN_COMMENT];
    
    if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, comment, &status) == 0) {
        fits_close_file(fptr, &status);
        QString result = QString::fromLatin1(dateobs).trimmed();
        
        // Remove any surrounding quotes that FITS might include
        if (result.startsWith('\'') && result.endsWith('\'')) {
            result = result.mid(1, result.length() - 2);
        }
        if (result.startsWith('"') && result.endsWith('"')) {
            result = result.mid(1, result.length() - 2);
        }
        
        if (m_debugMode) {
            logMessage(QString("Extracted DATE-OBS: '%1' from %2").arg(result).arg(QFileInfo(fitsFile).fileName()), "gray");
        }
        
        return result;
    }
    
    // If DATE-OBS doesn't exist, try other common date keywords
    QStringList dateKeywords = {"DATE", "DATE_OBS", "DATEOBS", "OBS-DATE"};
    
    for (const QString &keyword : dateKeywords) {
        status = 0; // Reset status
        QByteArray keyBytes = keyword.toLatin1();
        
        if (fits_read_key(fptr, TSTRING, keyBytes.data(), dateobs, comment, &status) == 0) {
            fits_close_file(fptr, &status);
            QString result = QString::fromLatin1(dateobs).trimmed();
            
            // Remove quotes
            if (result.startsWith('\'') && result.endsWith('\'')) {
                result = result.mid(1, result.length() - 2);
            }
            if (result.startsWith('"') && result.endsWith('"')) {
                result = result.mid(1, result.length() - 2);
            }
            
            if (m_debugMode) {
                logMessage(QString("Found date in keyword '%1': '%2' from %3").arg(keyword).arg(result).arg(QFileInfo(fitsFile).fileName()), "gray");
            }
            
            return result;
        }
    }
    
    fits_close_file(fptr, &status);
    
    if (m_debugMode) {
        logMessage(QString("No date keyword found in FITS header: %1").arg(QFileInfo(fitsFile).fileName()), "orange");
    }
    
    return QString();
}

// Also add this helper function to parse observer location more robustly
bool StellinaProcessor::parseObserverLocation(const QString &location, double &lat, double &lon, double &elevation) {
    // Default values (London)
    lat = 51.5074;
    lon = -0.1278;
    elevation = 0.0;
    
    if (location.isEmpty()) {
        return false;
    }
    
    // Try to parse different formats:
    // "London" -> lookup coordinates (not implemented here)
    // "51.5074,-0.1278" -> lat,lon
    // "51.5074,-0.1278,25" -> lat,lon,elevation
    // "51°30'27\"N,0°7'40\"W" -> DMS format (not implemented here)
    
    QStringList parts = location.split(',');
    if (parts.size() >= 2) {
        bool ok1, ok2;
        double parsedLat = parts[0].trimmed().toDouble(&ok1);
        double parsedLon = parts[1].trimmed().toDouble(&ok2);
        
        if (ok1 && ok2) {
            lat = parsedLat;
            lon = parsedLon;
            
            if (parts.size() >= 3) {
                bool ok3;
                double parsedElev = parts[2].trimmed().toDouble(&ok3);
                if (ok3) {
                    elevation = parsedElev;
                }
            }
            return true;
        }
    }
    
    // Could add city name lookup here
    // For now, just return false for unknown formats
    return false;
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
// Replace the applyMasterDark function in StellinaProcessor_Core.cpp with this implementation

bool StellinaProcessor::applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame) {
    logMessage(QString("Applying master dark to %1 using direct FITS processing").arg(QFileInfo(lightFrame).fileName()), "blue");
    
    fitsfile *lightFits = nullptr;
    fitsfile *darkFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open light frame
    QByteArray lightPath = lightFrame.toLocal8Bit();
    if (fits_open_file(&lightFits, lightPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open light frame: %1 (FITS error: %2)").arg(lightFrame).arg(status), "red");
        return false;
    }
    
    // Open master dark
    QByteArray darkPath = masterDark.toLocal8Bit();
    if (fits_open_file(&darkFits, darkPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open master dark: %1 (FITS error: %2)").arg(masterDark).arg(status), "red");
        fits_close_file(lightFits, &status);
        return false;
    }
    
    // Get dimensions and verify they match
    long lightNaxes[2], darkNaxes[2];
    if (fits_get_img_size(lightFits, 2, lightNaxes, &status) ||
        fits_get_img_size(darkFits, 2, darkNaxes, &status)) {
        logMessage(QString("Failed to get image dimensions (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    if (lightNaxes[0] != darkNaxes[0] || lightNaxes[1] != darkNaxes[1]) {
        logMessage(QString("Image dimensions mismatch - Light: %1x%2, Dark: %3x%4")
                      .arg(lightNaxes[0]).arg(lightNaxes[1])
                      .arg(darkNaxes[0]).arg(darkNaxes[1]), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    long totalPixels = lightNaxes[0] * lightNaxes[1];
    logMessage(QString("Processing %1 x %2 image (%3 pixels)").arg(lightNaxes[0]).arg(lightNaxes[1]).arg(totalPixels), "blue");
    
    // Create output FITS file
    QByteArray outputPathBytes = QString("!%1").arg(outputFrame).toLocal8Bit(); // ! forces overwrite
    if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
        logMessage(QString("Failed to create output FITS file: %1 (FITS error: %2)").arg(outputFrame).arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    // Copy header from light frame
    if (fits_copy_header(lightFits, outputFits, &status)) {
        logMessage(QString("Failed to copy FITS header (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    // Allocate memory for pixel data
    std::vector<float> lightPixels(totalPixels);
    std::vector<float> darkPixels(totalPixels);
    std::vector<float> resultPixels(totalPixels);
    
    // Read light frame pixel data
    long firstPixel = 1;
    if (fits_read_img(lightFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                     lightPixels.data(), nullptr, &status)) {
        logMessage(QString("Failed to read light frame pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Read master dark pixel data
    if (fits_read_img(darkFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                     darkPixels.data(), nullptr, &status)) {
        logMessage(QString("Failed to read master dark pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Close input files as we no longer need them
    fits_close_file(lightFits, &status);
    fits_close_file(darkFits, &status);
    
    // Perform dark subtraction with negative truncation
    long negativePixels = 0;
    float minResult = std::numeric_limits<float>::max();
    float maxResult = std::numeric_limits<float>::lowest();
    
    for (long i = 0; i < totalPixels; ++i) {
        float result = lightPixels[i] - darkPixels[i];
        
        // Truncate negative values to 0
        if (result < 0.0f) {
            result = 0.0f;
            negativePixels++;
        }
        
        resultPixels[i] = result;
        
        // Track min/max for statistics
        if (result < minResult) minResult = result;
        if (result > maxResult) maxResult = result;
    }
    
    // Log statistics
    double negativePercent = (static_cast<double>(negativePixels) / totalPixels) * 100.0;
    logMessage(QString("Dark subtraction complete - Negative pixels truncated: %1 (%2%)")
                  .arg(negativePixels).arg(negativePercent, 0, 'f', 2), "blue");
    logMessage(QString("Result range: %1 to %2").arg(minResult, 0, 'f', 2).arg(maxResult, 0, 'f', 2), "gray");
    
    // Write result to output file
    if (fits_write_img(outputFits, TFLOAT, firstPixel, totalPixels, 
                      resultPixels.data(), &status)) {
        logMessage(QString("Failed to write calibrated pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Update header with processing information
    QString historyComment = QString("Dark calibrated using master dark: %1").arg(QFileInfo(masterDark).fileName());
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(outputFits, historyBytes.data(), &status)) {
        // Non-critical error
        logMessage("Warning: Could not write processing history to header", "orange");
        status = 0;
    }
    
    // Add custom keywords for processing info
    int truncatedCount = static_cast<int>(negativePixels);
    if (fits_write_key(outputFits, TINT, "NEGTRUNC", &truncatedCount, "Number of negative pixels truncated", &status)) {
        logMessage("Warning: Could not write truncation count to header", "orange");
        status = 0;
    }
    
    QString processingDate = QDateTime::currentDateTime().toString(Qt::ISODate);
    QByteArray dateBytes = processingDate.toLocal8Bit();
    char* datePtr = dateBytes.data();
    if (fits_write_key(outputFits, TSTRING, "DARKCAL", &datePtr, "Dark calibration date", &status)) {
        logMessage("Warning: Could not write processing date to header", "orange");
        status = 0;
    }
    
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error occurred during dark calibration (error: %1)").arg(status), "red");
        QFile::remove(outputFrame);
        return false;
    }
    
    logMessage(QString("Dark calibration successful: %1").arg(QFileInfo(outputFrame).fileName()), "green");
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
// Add these test functions to StellinaProcessor_Core.cpp
// Call testLibnovaConversion() from your constructor or a menu action for verification

void StellinaProcessor::testLibnovaConversion() {
    logMessage("=== Testing libnova coordinate conversion ===", "blue");
    
    // Test case 1: First image from your log
    // Blind solve result: RA=10.9127702516°, Dec=41.2122118775°
    // Your current result: RA=132.720487°, Dec=22.962432°
    testSingleConversion(
        "Test 1 - First image",
        42.0410,           // Alt from log
        286.8526,          // Az from log  
        "2024-01-09T22:13:29",  // DATE-OBS from log
        10.9127702516,     // Expected RA from blind solve
        41.2122118775,     // Expected Dec from blind solve
        132.720487,        // Current RA result
        22.962432          // Current Dec result
    );
    
    // Test case 2: Known good coordinate conversion
    // Use a well-known object at a specific time for validation
    // M31 (Andromeda Galaxy) coordinates: RA=10.684°, Dec=41.269°
    testSingleConversion(
        "Test 2 - M31 reference",
        45.0,              // Approximate alt for M31 from London
        290.0,             // Approximate az for M31 from London
        "2024-01-09T22:00:00",  // Round time
        10.684,            // M31 RA (known)
        41.269,            // M31 Dec (known)
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    // Test case 3: Different time to check time dependency
    testSingleConversion(
        "Test 3 - Time dependency check",
        42.0410,           // Same Alt as test 1
        286.8526,          // Same Az as test 1
        "2024-01-09T23:13:29",  // 1 hour later
        0.0,               // Expected will be different due to time
        0.0,               // Expected will be different due to time
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    // Test case 4: Different observer location
    testSingleConversion(
        "Test 4 - Different location (Paris)",
        42.0410,           // Same Alt
        286.8526,          // Same Az
        "2024-01-09T22:13:29",  // Same time
        0.0,               // Expected will be different due to location
        0.0,               // Expected will be different due to location
        0.0,               // Will be calculated
        0.0,               // Will be calculated
        48.8566,           // Paris latitude
        2.3522             // Paris longitude
    );
    
    // Test case 5: Zenith pointing (RA should equal LST, Dec should equal latitude)
    // For 2024-01-09T22:13:29 at London, LST ≈ 17.5 hours ≈ 262.5°
    testSingleConversion(
        "Test 5 - Zenith pointing",
        90.0,              // Pointing straight up
        0.0,               // Azimuth doesn't matter at zenith
        "2024-01-09T22:13:29",
        262.5,             // Expected RA ≈ LST (will verify in test)
        51.5074,           // Dec should equal observer latitude
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    logMessage("=== libnova coordinate conversion tests complete ===", "blue");
}

void StellinaProcessor::testSingleConversion(const QString &testName,
                                           double alt, double az, 
                                           const QString &dateObs,
                                           double expectedRA, double expectedDec,
                                           double currentRA, double currentDec,
                                           double testLat, double testLon) {
    
    logMessage(QString("--- %1 ---").arg(testName), "green");
    logMessage(QString("Input: Alt=%1°, Az=%2°, Time=%3")
                  .arg(alt, 0, 'f', 4)
                  .arg(az, 0, 'f', 4)
                  .arg(dateObs), "gray");
    
    // Save current observer location
    QString savedLocation = m_observerLocation;
    
    // Set test location if provided
    if (testLat != 0.0 || testLon != 0.0) {
        m_observerLocation = QString("%1,%2").arg(testLat, 0, 'f', 4).arg(testLon, 0, 'f', 4);
        logMessage(QString("Using test location: %1°N, %2°E")
                      .arg(testLat, 0, 'f', 4).arg(testLon, 0, 'f', 4), "gray");
    }
    
    // Perform conversion
    double calculatedRA, calculatedDec;
    bool success = convertAltAzToRaDec(alt, az, dateObs, calculatedRA, calculatedDec);
    
    if (success) {
        logMessage(QString("Calculated: RA=%1°, Dec=%2°")
                      .arg(calculatedRA, 0, 'f', 6)
                      .arg(calculatedDec, 0, 'f', 6), "blue");
        
        // Show in HMS/DMS format too
        double ra_hours = calculatedRA / 15.0;
        int h = static_cast<int>(ra_hours);
        int m = static_cast<int>((ra_hours - h) * 60);
        double s = ((ra_hours - h) * 60 - m) * 60;
        
        int d = static_cast<int>(calculatedDec);
        int am = static_cast<int>(qAbs(calculatedDec - d) * 60);
        double as = (qAbs(calculatedDec - d) * 60 - am) * 60;
        
        logMessage(QString("          RA=%1h%2m%3s, Dec=%4°%5'%6\"")
                      .arg(h).arg(m, 2, 10, QChar('0')).arg(s, 0, 'f', 1)
                      .arg(d).arg(am, 2, 10, QChar('0')).arg(as, 0, 'f', 1), "blue");
        
        // Compare with expected results if provided
        if (expectedRA != 0.0 || expectedDec != 0.0) {
            double raError = qAbs(calculatedRA - expectedRA);
            double decError = qAbs(calculatedDec - expectedDec);
            
            logMessage(QString("Expected: RA=%1°, Dec=%2°")
                          .arg(expectedRA, 0, 'f', 6)
                          .arg(expectedDec, 0, 'f', 6), "orange");
            logMessage(QString("Error:    RA=%1° (%2 arcmin), Dec=%3° (%4 arcmin)")
                          .arg(raError, 0, 'f', 6).arg(raError * 60, 0, 'f', 1)
                          .arg(decError, 0, 'f', 6).arg(decError * 60, 0, 'f', 1), 
                      (raError < 1.0 && decError < 1.0) ? "green" : "red");
            
            // Check if we're close to expected (within 1 degree)
            if (raError < 1.0 && decError < 1.0) {
                logMessage("✓ PASS: Within 1° tolerance", "green");
            } else {
                logMessage("✗ FAIL: Outside 1° tolerance", "red");
                
                // Try to diagnose the issue
                if (raError > 10.0) {
                    logMessage("Large RA error suggests coordinate system issue", "red");
                }
                if (qAbs(calculatedRA - currentRA) < 1.0) {
                    logMessage("Current result matches - may be systematic error", "orange");
                }
            }
        }
        
        // Compare with current result if provided
        if (currentRA != 0.0 || currentDec != 0.0) {
            logMessage(QString("Previous: RA=%1°, Dec=%2°")
                          .arg(currentRA, 0, 'f', 6)
                          .arg(currentDec, 0, 'f', 6), "gray");
        }
        
    } else {
        logMessage("✗ CONVERSION FAILED", "red");
    }
    
    // Restore original observer location
    m_observerLocation = savedLocation;
    
    logMessage("", "gray");  // Blank line for separation
}
