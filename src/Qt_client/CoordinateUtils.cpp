#include "CoordinateUtils.h"
#include <QDateTime>
#include <cmath>

double CoordinateUtils::computeJulianDay(int year, int month, int day, 
                                           int hour, int minute, int second)
{
    // Handle month adjustment (January and February are treated as months 13 and 14 of the previous year)
    if (month < 3) {
        year -= 1;
        month += 12;
    }
    
    // Calculate the Julian day based on the Gregorian calendar
    int A = year / 100;
    int B = 2 - A + (A / 4);
    
    double E = std::floor(365.25 * static_cast<double>(year + 4716));
    double F = std::floor(30.6001 * static_cast<double>(month + 1));
    
    // Calculate final Julian day with time component
    double julianDay = static_cast<double>(B + day) + E + F - 1524.5;
    julianDay += (hour + minute / 60.0 + second / 3600.0) / 24.0;
    
    return julianDay;
}

double CoordinateUtils::localSiderealTime(double longitude, double jd)
{
    // Calculate days since J2000.0
    double daysSinceJ2000 = jd - J2000_EPOCH;
    
    // Calculate T (Julian centuries since J2000.0)
    double T = daysSinceJ2000 / 36525.0;
    
    // Calculate Greenwich Mean Sidereal Time
    double midnight = std::floor(jd) + 0.5;
    double daysSinceMidnight = jd - midnight;
    double hoursSinceMidnight = daysSinceMidnight * 24.0;
    double wholeDaysSinceJ2000 = midnight - J2000_EPOCH;
    
    double GMST = 6.697374558 +
                  0.06570982441908 * wholeDaysSinceJ2000 +
                  1.00273790935 * hoursSinceMidnight +
                  0.000026 * T * T;
    
    // Convert to the range 0-24 hours
    GMST = normalizeHours(GMST);
    
    // Add the longitude (in hours) to get local sidereal time
    double LST = GMST + longitude / 15.0;
    
    // Normalize to 0-24 hours
    return normalizeHours(LST);
}

std::tuple<double, double, double> CoordinateUtils::raDecToAltAz(
    double ra, double dec, double latitude, double longitude, double lst)
{
    // Convert dec to radians
    double decRad = dec * DEG_TO_RAD;
    
    // Calculate Hour Angle in hours and convert to radians
    double ha = normalizeHours(lst - ra / 15.0 + 24.0);
    double haRad = ha * HOUR_TO_RAD;
    
    // Convert latitude to radians
    double latRad = latitude * DEG_TO_RAD;
    
    // Calculate intermediate values
    double x = std::cos(haRad) * std::cos(decRad);
    double y = std::sin(haRad) * std::cos(decRad);
    double z = std::sin(decRad);
    
    // Transform to horizontal coordinates
    double xhor = x * std::cos((90.0 - latitude) * DEG_TO_RAD) - z * std::sin((90.0 - latitude) * DEG_TO_RAD);
    double yhor = y;
    double zhor = x * std::sin((90.0 - latitude) * DEG_TO_RAD) + z * std::cos((90.0 - latitude) * DEG_TO_RAD);
    
    // Calculate azimuth and altitude
    double az = std::atan2(yhor, xhor) * RAD_TO_DEG + 180.0;
    double alt = std::asin(zhor) * RAD_TO_DEG;
    
    // Normalize azimuth to 0-360 degrees
    az = normalizeDegrees(az);
    
    return std::make_tuple(alt, az, ha);
}

std::tuple<double, double, double> CoordinateUtils::altAzToRaDec(
    double alt, double az, double latitude, double longitude, double lst)
{
    // Correct for atmospheric refraction
    double trueAlt = correctRefraction(alt);
    
    // Convert to radians
    double altRad = trueAlt * DEG_TO_RAD;
    double azRad = az * DEG_TO_RAD;
    double latRad = latitude * DEG_TO_RAD;
    
    // Calculate declination
    double sinDec = std::sin(altRad) * std::sin(latRad) +
                    std::cos(altRad) * std::cos(latRad) * std::cos(azRad);
    double dec = std::asin(sinDec) * RAD_TO_DEG;
    
    // Calculate hour angle
    double cosDecRad = std::cos(dec * DEG_TO_RAD);
    double cosHA = (std::sin(altRad) - std::sin(latRad) * sinDec) /
                   (std::cos(latRad) * cosDecRad);
    double sinHA = -std::cos(altRad) * std::sin(azRad) / cosDecRad;
    double haRad = std::atan2(sinHA, cosHA);
    double haHours = haRad * RAD_TO_HOUR;
    
    // Ensure hour angle is in -12 to +12 range
    if (haHours < -12.0) haHours += 24.0;
    else if (haHours > 12.0) haHours -= 24.0;
    
    // Calculate RA
    double raHours = lst - haHours;
    
    // Normalize to 0-24 hours
    raHours = normalizeHours(raHours);
    
    // Convert RA to degrees
    double raDeg = raHours * 15.0;
    
    return std::make_tuple(raDeg, dec, haHours);
}

std::tuple<double, double> CoordinateUtils::j2000ToJNow(double ra2000, double dec2000)
{
    // Get current time
    QDateTime now = QDateTime::currentDateTimeUtc();
    QDateTime j2000(QDate(2000, 1, 1), QTime(12, 0, 0), Qt::UTC);
    
    // Calculate T in Julian centuries
    double T = j2000.daysTo(now) / 36525.0;
    
    // Calculate correction factors M and N
    double M = 1.2812323 * T + 0.0003879 * T * T + 0.0000101 * T * T * T;
    double N = 0.5567530 * T - 0.0001185 * T * T + 0.0000116 * T * T * T;
    
    // Calculate delta RA and delta Dec
    double ra2000Rad = ra2000 * DEG_TO_RAD;
    double dec2000Rad = dec2000 * DEG_TO_RAD;
    double deltaRa = M + N * std::sin(ra2000Rad) * std::tan(dec2000Rad);
    double deltaDec = N * std::cos(ra2000Rad);
    
    // Add corrections to get current epoch coordinates
    double raNow = ra2000 + deltaRa;
    double decNow = dec2000 + deltaDec;
    
    return std::make_tuple(raNow, decNow);
}

std::tuple<double, double> CoordinateUtils::jNowToJ2000(double raNow, double decNow)
{
    // Get current time
    QDateTime now = QDateTime::currentDateTimeUtc();
    QDateTime j2000(QDate(2000, 1, 1), QTime(12, 0, 0), Qt::UTC);
    
    // Calculate T in Julian centuries
    double T = j2000.daysTo(now) / 36525.0;
    
    // Calculate correction factors M and N
    double M = 1.2812323 * T + 0.0003879 * T * T + 0.0000101 * T * T * T;
    double N = 0.5567530 * T - 0.0001185 * T * T + 0.0000116 * T * T * T;
    
    // Calculate delta RA and delta Dec
    double raNowRad = raNow * DEG_TO_RAD;
    double decNowRad = decNow * DEG_TO_RAD;
    double deltaRa = M + N * std::sin(raNowRad) * std::tan(decNowRad);
    double deltaDec = N * std::cos(raNowRad);
    
    // Subtract corrections to get J2000 coordinates
    double ra2000 = raNow - deltaRa;
    double dec2000 = decNow - deltaDec;
    
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
    // Convert to radians
    double altRad = apparentAlt * DEG_TO_RAD;
    
    // Calculate refraction correction in arcminutes
    double refractionArcmin;
    
    if (apparentAlt > 15.0) {
        // Above 15° altitude - simpler formula
        refractionArcmin = 1.02 / std::tan(altRad + 10.3 / (apparentAlt + 5.11));
    } else if (apparentAlt >= 0.0) {
        // Low altitude - more complex formula
        double r = 0.1594 + apparentAlt * (0.0196 + 0.00002 * apparentAlt);
        refractionArcmin = r * (1.0 / std::tan(altRad));
    } else {
        // Below horizon - use horizon value
        refractionArcmin = 34.0;
    }
    
    // Convert arcminutes to degrees
    double refractionDegrees = refractionArcmin / 60.0;
    
    // Return true altitude
    return apparentAlt - refractionDegrees;
}

QString CoordinateUtils::formatRaAsHMS(double ra)
{
    // Convert RA from degrees to hours
    double raHours = ra / 15.0;
    
    // Extract hours, minutes, and seconds
    double hours = std::floor(raHours);
    double minutes = std::floor((raHours - hours) * 60.0);
    double seconds = ((raHours - hours) * 60.0 - minutes) * 60.0;
    
    // Format as string
    return QString("%1h %2m %3s")
        .arg(static_cast<int>(hours), 2, 10, QChar('0'))
        .arg(static_cast<int>(minutes), 2, 10, QChar('0'))
        .arg(seconds, 4, 'f', 1, QChar('0'));
}

QString CoordinateUtils::formatDecAsDMS(double dec)
{
    // Extract sign
    QString sign = (dec < 0) ? "-" : "+";
    dec = std::abs(dec);
    
    // Extract degrees, minutes, and seconds
    double degrees = std::floor(dec);
    double minutes = std::floor((dec - degrees) * 60.0);
    double seconds = ((dec - degrees) * 60.0 - minutes) * 60.0;
    
    // Format as string
    return QString("%1%2° %3' %4\"")
        .arg(sign)
        .arg(static_cast<int>(degrees), 2, 10, QChar('0'))
        .arg(static_cast<int>(minutes), 2, 10, QChar('0'))
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
    double decDegrees = std::abs(degrees) + minutes / 60.0 + seconds / 3600.0;
    
    // Apply sign
    if (degrees < 0) {
        decDegrees = -decDegrees;
    }
    
    return decDegrees;
}

double CoordinateUtils::normalizeHours(double hours)
{
    hours = std::fmod(hours, 24.0);
    if (hours < 0) hours += 24.0;
    return hours;
}

double CoordinateUtils::normalizeDegrees(double degrees)
{
    degrees = std::fmod(degrees, 360.0);
    if (degrees < 0) degrees += 360.0;
    return degrees;
}

double CoordinateUtils::centuriesSinceJ2000(double jd)
{
    return (jd - J2000_EPOCH) / 36525.0;
}
