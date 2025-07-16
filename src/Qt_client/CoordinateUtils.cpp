#include "CoordinateUtils.h"
#include <QDateTime>
#include <cmath>

// Constants based on Stellarium's implementation
const double J2000_EPOCH = 2451545.0;  // Julian day for J2000.0 epoch
const double OBLIQUITY_J2000 = 23.4392911 * DEG_TO_RAD;  // J2000 obliquity in radians
const double JULIAN_CENTURY = 36525.0;
const double ARCSEC_TO_RAD = DEG_TO_RAD / 3600.0;
const double HOUR_TO_RAD = M_PI / 12.0;  // Convert hours to radians
const double RAD_TO_HOUR = 12.0 / M_PI;  // Convert radians to hours

// Transformation matrices (updated during each computation)
static double matJ2000ToEquinoxEqu[9];
static double matEquinoxEquToJ2000[9];
static double matEquinoxEquToAltAz[9];
static double matAltAzToEquinoxEqu[9];

// Helper functions for matrix operations
static void matrixMultiply3x3(const double a[9], const double b[9], double result[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i * 3 + j] = 0;
            for (int k = 0; k < 3; k++) {
                result[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
}

static void matrixTranspose3x3(const double matrix[9], double result[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i * 3 + j] = matrix[j * 3 + i];
        }
    }
}

static void transformVector3(const double matrix[9], double x, double y, double z, 
                           double& outX, double& outY, double& outZ) {
    outX = matrix[0] * x + matrix[1] * y + matrix[2] * z;
    outY = matrix[3] * x + matrix[4] * y + matrix[5] * z;
    outZ = matrix[6] * x + matrix[7] * y + matrix[8] * z;
}

static void setRotationMatrixZ(double angle, double matrix[9]) {
    double c = cos(angle);
    double s = sin(angle);
    
    matrix[0] = c;   matrix[1] = -s;  matrix[2] = 0;
    matrix[3] = s;   matrix[4] = c;   matrix[5] = 0;
    matrix[6] = 0;   matrix[7] = 0;   matrix[8] = 1;
}

static void setRotationMatrixY(double angle, double matrix[9]) {
    double c = cos(angle);
    double s = sin(angle);
    
    matrix[0] = c;   matrix[1] = 0;   matrix[2] = s;
    matrix[3] = 0;   matrix[4] = 1;   matrix[5] = 0;
    matrix[6] = -s;  matrix[7] = 0;   matrix[8] = c;
}

static void setRotationMatrixX(double angle, double matrix[9]) {
    double c = cos(angle);
    double s = sin(angle);
    
    matrix[0] = 1;   matrix[1] = 0;   matrix[2] = 0;
    matrix[3] = 0;   matrix[4] = c;   matrix[5] = -s;
    matrix[6] = 0;   matrix[7] = s;   matrix[8] = c;
}

static void setIdentityMatrix(double matrix[9]) {
    for (int i = 0; i < 9; i++) matrix[i] = 0.0;
    matrix[0] = matrix[4] = matrix[8] = 1.0;
}

double CoordinateUtils::computeJulianDay(int year, int month, int day, 
                                           int hour, int minute, int second)
{
    // Stellarium uses the standard Julian Day formula
    // Based on Jean Meeus "Astronomical Algorithms"
    
    if (month < 3) {
        year -= 1;
        month += 12;
    }
    
    int A = year / 100;
    int B = 2 - A + (A / 4);
    
    if (year < 1582 || (year == 1582 && month < 10) || 
        (year == 1582 && month == 10 && day < 15)) {
        B = 0;  // Julian calendar
    }
    
    double E = floor(365.25 * (year + 4716));
    double F = floor(30.6001 * (month + 1));
    
    double julianDay = B + day + E + F - 1524.5;
    julianDay += (hour + minute / 60.0 + second / 3600.0) / 24.0;
    
    return julianDay;
}

// Stellarium-style precession calculation
static void computePrecessionMatrix(double jd, double matrix[9]) {
    // Based on IAU 2000 precession model (Capitaine et al. 2003)
    // This is what Stellarium uses
    
    double T = (jd - J2000_EPOCH) / JULIAN_CENTURY;
    
    // Precession angles in arcseconds
    double zeta_A = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * T;
    double z_A = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * T + 
                 (0.30188 - 0.000344 * T) * T * T;
    double theta_A = (2004.3109 - 0.85330 * T - 0.000217 * T * T) * T;
    
    // Convert to radians
    zeta_A *= ARCSEC_TO_RAD;
    z_A *= ARCSEC_TO_RAD;
    theta_A *= ARCSEC_TO_RAD;
    
    // Build precession matrix
    double cos_zeta = cos(zeta_A);
    double sin_zeta = sin(zeta_A);
    double cos_z = cos(z_A);
    double sin_z = sin(z_A);
    double cos_theta = cos(theta_A);
    double sin_theta = sin(theta_A);
    
    matrix[0] = cos_zeta * cos_z * cos_theta - sin_zeta * sin_z;
    matrix[1] = -sin_zeta * cos_z * cos_theta - cos_zeta * sin_z;
    matrix[2] = -sin_theta * cos_z;
    matrix[3] = cos_zeta * sin_z * cos_theta + sin_zeta * cos_z;
    matrix[4] = -sin_zeta * sin_z * cos_theta + cos_zeta * cos_z;
    matrix[5] = -sin_theta * sin_z;
    matrix[6] = cos_zeta * sin_theta;
    matrix[7] = -sin_zeta * sin_theta;
    matrix[8] = cos_theta;
}

static void computeLocalTransforms(double latitude, double longitude, double lst) {
    // Create transformation from equatorial to horizontal coordinates
    double latRad = latitude * DEG_TO_RAD;
    double lstRad = lst * HOUR_TO_RAD;
    
    // Rotation matrices
    double rotLat[9], rotLst[9];
    setRotationMatrixY(M_PI_2 - latRad, rotLat);
    setRotationMatrixZ(lstRad, rotLst);
    
    // Combined transformation
    matrixMultiply3x3(rotLat, rotLst, matEquinoxEquToAltAz);
    matrixTranspose3x3(matEquinoxEquToAltAz, matAltAzToEquinoxEqu);
}

double CoordinateUtils::localSiderealTime(double longitude, double jd)
{
    // Stellarium's GMST calculation based on IAU 2000
    double T = (jd - J2000_EPOCH) / JULIAN_CENTURY;
    
    // Calculate Greenwich Mean Sidereal Time (GMST)
    // Formula from Meeus "Astronomical Algorithms" Chapter 12
    double gmst = 280.46061837 + 360.98564736629 * (jd - J2000_EPOCH) +
                  0.000387933 * T * T - T * T * T / 38710000.0;
    
    // Normalize to 0-360 degrees
    gmst = fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;
    
    // Convert to hours and add longitude
    double lst = gmst / 15.0 + longitude / 15.0;
    
    return normalizeHours(lst);
}

std::tuple<double, double, double> CoordinateUtils::raDecToAltAz(
    double ra, double dec, double latitude, double longitude, double lst)
{
    // Convert to Cartesian coordinates
    double raRad = ra * DEG_TO_RAD;
    double decRad = dec * DEG_TO_RAD;
    
    double x = cos(decRad) * cos(raRad);
    double y = cos(decRad) * sin(raRad);
    double z = sin(decRad);
    
    // Set up local transformation matrices
    computeLocalTransforms(latitude, longitude, lst);
    
    // Transform to horizontal coordinates
    double altX, altY, altZ;
    transformVector3(matEquinoxEquToAltAz, x, y, z, altX, altY, altZ);
    
    // Convert back to spherical coordinates
    double alt = asin(altZ) * RAD_TO_DEG;
    double az = atan2(altY, altX) * RAD_TO_DEG;
    
    // Normalize azimuth to 0-360
    az = normalizeDegrees(az);
    
    // Calculate hour angle
    double ha = lst - ra / 15.0;
    ha = normalizeHours(ha + 24.0);
    
    return std::make_tuple(alt, az, ha);
}

std::tuple<double, double, double> CoordinateUtils::altAzToRaDec(
    double alt, double az, double latitude, double longitude, double lst)
{
    // Apply atmospheric refraction correction
    double trueAlt = correctRefraction(alt);
    
    // Convert to Cartesian coordinates
    double altRad = trueAlt * DEG_TO_RAD;
    double azRad = az * DEG_TO_RAD;
    
    double x = cos(altRad) * cos(azRad);
    double y = cos(altRad) * sin(azRad);
    double z = sin(altRad);
    
    // Set up local transformation matrices
    computeLocalTransforms(latitude, longitude, lst);
    
    // Transform to equatorial coordinates
    double eqX, eqY, eqZ;
    transformVector3(matAltAzToEquinoxEqu, x, y, z, eqX, eqY, eqZ);
    
    // Convert back to spherical coordinates
    double ra = atan2(eqY, eqX) * RAD_TO_DEG;
    double dec = asin(eqZ) * RAD_TO_DEG;
    
    // Normalize RA to 0-360
    ra = normalizeDegrees(ra);
    
    // Calculate hour angle
    double haHours = lst - ra / 15.0;
    if (haHours < -12.0) haHours += 24.0;
    else if (haHours > 12.0) haHours -= 24.0;
    
    return std::make_tuple(ra, dec, haHours);
}

std::tuple<double, double> CoordinateUtils::j2000ToJNow(double ra2000, double dec2000)
{
    // Get current time for precession calculation
    QDateTime now = QDateTime::currentDateTimeUtc();
    double currentJD = computeJulianDay(now.date().year(), now.date().month(), now.date().day(),
                                       now.time().hour(), now.time().minute(), now.time().second());
    
    // Calculate precession matrix from J2000 to current epoch
    computePrecessionMatrix(currentJD, matJ2000ToEquinoxEqu);
    
    // Convert to Cartesian coordinates
    double raRad = ra2000 * DEG_TO_RAD;
    double decRad = dec2000 * DEG_TO_RAD;
    
    double x = cos(decRad) * cos(raRad);
    double y = cos(decRad) * sin(raRad);
    double z = sin(decRad);
    
    // Apply precession
    double newX, newY, newZ;
    transformVector3(matJ2000ToEquinoxEqu, x, y, z, newX, newY, newZ);
    
    // Convert back to spherical coordinates
    double raNow = atan2(newY, newX) * RAD_TO_DEG;
    double decNow = asin(newZ) * RAD_TO_DEG;
    
    // Normalize RA
    raNow = normalizeDegrees(raNow);
    
    return std::make_tuple(raNow, decNow);
}

std::tuple<double, double> CoordinateUtils::jNowToJ2000(double raNow, double decNow)
{
    // Get current time for precession calculation
    QDateTime now = QDateTime::currentDateTimeUtc();
    double currentJD = computeJulianDay(now.date().year(), now.date().month(), now.date().day(),
                                       now.time().hour(), now.time().minute(), now.time().second());
    
    // Calculate precession matrix from J2000 to current epoch
    computePrecessionMatrix(currentJD, matJ2000ToEquinoxEqu);
    // Create inverse matrix (transpose for rotation matrix)
    matrixTranspose3x3(matJ2000ToEquinoxEqu, matEquinoxEquToJ2000);
    
    // Convert to Cartesian coordinates
    double raRad = raNow * DEG_TO_RAD;
    double decRad = decNow * DEG_TO_RAD;
    
    double x = cos(decRad) * cos(raRad);
    double y = cos(decRad) * sin(raRad);
    double z = sin(decRad);
    
    // Apply inverse precession
    double newX, newY, newZ;
    transformVector3(matEquinoxEquToJ2000, x, y, z, newX, newY, newZ);
    
    // Convert back to spherical coordinates
    double ra2000 = atan2(newY, newX) * RAD_TO_DEG;
    double dec2000 = asin(newZ) * RAD_TO_DEG;
    
    // Normalize RA
    ra2000 = normalizeDegrees(ra2000);
    
    return std::make_tuple(ra2000, dec2000);
}

std::tuple<double, double, double, double, double, double, double> 
CoordinateUtils::calculateAltAz(
    int year, int month, int day, int hour, int minute, int second,
    double ra, double dec, double latitude, double longitude)
{
    // Calculate Julian date
    double jd = computeJulianDay(year, month, day, hour, minute, second);
    
    // Calculate local sidereal time
    double lst = localSiderealTime(longitude, jd);
    
    // Convert J2000 coordinates to current epoch
    auto [raNow, decNow] = j2000ToJNow(ra, dec);
    
    // Calculate altitude, azimuth, and hour angle
    auto [alt, az, ha] = raDecToAltAz(raNow, decNow, latitude, longitude, lst);
    
    return std::make_tuple(jd, raNow, decNow, alt, az, lst, ha);
}

std::tuple<double, double, double, double, double, double, double>
CoordinateUtils::calculateRaDec(
    int year, int month, int day, int hour, int minute, int second,
    double alt, double az, double latitude, double longitude)
{
    // Calculate Julian date
    double jd = computeJulianDay(year, month, day, hour, minute, second);
    
    // Calculate local sidereal time
    double lst = localSiderealTime(longitude, jd);
    
    // Convert horizontal to equatorial
    auto [raNow, decNow, ha] = altAzToRaDec(alt, az, latitude, longitude, lst);
    
    // Convert current epoch to J2000
    auto [ra2000, dec2000] = jNowToJ2000(raNow, decNow);
    
    return std::make_tuple(jd, ra2000, dec2000, raNow, decNow, lst, ha);
}

double CoordinateUtils::correctRefraction(double apparentAlt)
{
    // Stellarium's refraction model based on Bennett (1982)
    // "The Calculation of Astronomical Refraction in Marine Navigation"
    
    if (apparentAlt < -0.575) {
        // Below -34.5 arcminutes, refraction is approximately constant
        return apparentAlt - 34.5 / 60.0;
    }
    
    double altRad = apparentAlt * DEG_TO_RAD;
    double tanAlt = tan(altRad);
    
    // Bennett's formula for refraction in arcminutes
    double refraction = 1.02 / tanAlt - 0.0019279 + 0.0000001475 * tanAlt;
    
    // For low altitudes, use a more accurate formula
    if (apparentAlt < 15.0) {
        double cotAlt = 1.0 / tanAlt;
        refraction = cotAlt - 0.0033 * cotAlt * cotAlt * cotAlt;
    }
    
    // Convert to degrees and apply correction
    return apparentAlt - refraction / 60.0;
}

QString CoordinateUtils::formatRaAsHMS(double ra)
{
    // Convert RA from degrees to hours
    double raHours = ra / 15.0;
    
    // Extract hours, minutes, and seconds
    int hours = static_cast<int>(floor(raHours));
    double remainingMinutes = (raHours - hours) * 60.0;
    int minutes = static_cast<int>(floor(remainingMinutes));
    double seconds = (remainingMinutes - minutes) * 60.0;
    
    // Format as string
    return QString("%1h %2m %3s")
        .arg(hours, 2, 10, QChar('0'))
        .arg(minutes, 2, 10, QChar('0'))
        .arg(seconds, 4, 'f', 1, QChar('0'));
}

QString CoordinateUtils::formatDecAsDMS(double dec)
{
    // Extract sign
    QString sign = (dec < 0) ? "-" : "+";
    dec = fabs(dec);
    
    // Extract degrees, minutes, and seconds
    int degrees = static_cast<int>(floor(dec));
    double remainingMinutes = (dec - degrees) * 60.0;
    int minutes = static_cast<int>(floor(remainingMinutes));
    double seconds = (remainingMinutes - minutes) * 60.0;
    
    // Format as string
    return QString("%1%2° %3' %4\"")
        .arg(sign)
        .arg(degrees, 2, 10, QChar('0'))
        .arg(minutes, 2, 10, QChar('0'))
        .arg(seconds, 4, 'f', 1, QChar('0'));
}

double CoordinateUtils::parseRaFromHMS(int hours, int minutes, double seconds)
{
    // Convert to decimal hours
    double raHours = hours + minutes / 60.0 + seconds / 3600.0;
    
    // Convert hours to degrees
    return raHours * 15.0;
}

double CoordinateUtils::parseDecFromDMS(int degrees, int minutes, double seconds)
{
    // Convert to decimal degrees
    double decDegrees = abs(degrees) + minutes / 60.0 + seconds / 3600.0;
    
    // Apply sign
    if (degrees < 0) {
        decDegrees = -decDegrees;
    }
    
    return decDegrees;
}

double CoordinateUtils::normalizeHours(double hours)
{
    hours = fmod(hours, 24.0);
    if (hours < 0) hours += 24.0;
    return hours;
}

double CoordinateUtils::normalizeDegrees(double degrees)
{
    degrees = fmod(degrees, 360.0);
    if (degrees < 0) degrees += 360.0;
    return degrees;
}

double CoordinateUtils::centuriesSinceJ2000(double jd)
{
    return (jd - J2000_EPOCH) / JULIAN_CENTURY;
}