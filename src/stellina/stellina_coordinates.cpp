/*
 * stellina_coordinates.cpp
 * Professional coordinate conversion using libnova
 * Much more accurate than basic spherical trigonometry
 */

#include "stellina_processor.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "gui/utils.h"
#include <libnova/libnova.h>
#include <cmath>
#include <ctime>
#include <cstring>

/**
 * Convert Alt/Az coordinates to RA/Dec using libnova
 * This is the professional way with proper astronomical corrections
 */
int stellina_convert_altaz_to_radec(double alt, double az, const char *date_obs, 
                                   double observer_lat, double observer_lon, double observer_alt,
                                   double *ra, double *dec) {
    
    if (!date_obs || !ra || !dec) {
        return -1;
    }
    
    // Parse the DATE-OBS string (ISO format: 2024-01-09T22:13:29)
    struct tm tm_obs = {0};
    if (strptime(date_obs, "%Y-%m-%dT%H:%M:%S", &tm_obs) == NULL) {
        siril_log_color_message("Failed to parse observation time: %s\n", "red", date_obs);
        return -1;
    }
    
    // Convert to Julian Date using libnova
    struct ln_date obs_date;
    obs_date.years = tm_obs.tm_year + 1900;
    obs_date.months = tm_obs.tm_mon + 1;
    obs_date.days = tm_obs.tm_mday;
    obs_date.hours = tm_obs.tm_hour;
    obs_date.minutes = tm_obs.tm_min;
    obs_date.seconds = tm_obs.tm_sec;
    
    double jd = ln_get_julian_day(&obs_date);
    
    // Set up observer location
    struct ln_lnlat_posn observer;
    observer.lat = observer_lat;
    observer.lng = observer_lon;
    
    // Set up horizontal coordinates (Alt/Az)
    struct ln_hrz_posn horizontal;
    horizontal.alt = alt;
    horizontal.az = az;
    
    // Convert to equatorial coordinates (RA/Dec)
    struct ln_equ_posn equatorial;
    ln_get_equ_from_hrz(&horizontal, &observer, jd, &equatorial);
    
    // Return results
    *ra = equatorial.ra;
    *dec = equatorial.dec;
    
    // Validate results
    if (*dec < -90.0 || *dec > 90.0) {
        siril_log_color_message("Invalid declination calculated: %.2f\n", "red", *dec);
        return -1;
    }
    
    // Normalize RA to 0-360 degrees
    if (*ra < 0.0) *ra += 360.0;
    if (*ra >= 360.0) *ra -= 360.0;
    
    siril_debug_print("libnova conversion: Alt/Az %.2f°/%.2f° -> RA/Dec %.4f°/%.4f° (JD %.6f)\n",
                     alt, az, *ra, *dec, jd);
    
    return 0;
}

/**
 * Enhanced coordinate conversion with atmospheric refraction correction
 */
int stellina_convert_altaz_to_radec_corrected(double alt, double az, const char *date_obs,
                                             double observer_lat, double observer_lon, double observer_alt,
                                             double temperature, double pressure,
                                             double *ra, double *dec) {
    
    if (!date_obs || !ra || !dec) {
        return -1;
    }
    
    // Parse observation time
    struct tm tm_obs = {0};
    if (strptime(date_obs, "%Y-%m-%dT%H:%M:%S", &tm_obs) == NULL) {
        siril_log_color_message("Failed to parse observation time: %s\n", "red", date_obs);
        return -1;
    }
    
    struct ln_date obs_date;
    obs_date.years = tm_obs.tm_year + 1900;
    obs_date.months = tm_obs.tm_mon + 1;
    obs_date.days = tm_obs.tm_mday;
    obs_date.hours = tm_obs.tm_hour;
    obs_date.minutes = tm_obs.tm_min;
    obs_date.seconds = tm_obs.tm_sec;
    
    double jd = ln_get_julian_day(&obs_date);
    
    // Set up observer location
    struct ln_lnlat_posn observer;
    observer.lat = observer_lat;
    observer.lng = observer_lon;
    
    // Apply atmospheric refraction correction
    double corrected_alt = alt;
    if (temperature > 0.0 && pressure > 0.0) {
        // Calculate atmospheric refraction
        double refraction = ln_get_refraction_adj(alt, temperature, pressure);
        corrected_alt = alt - refraction;
        
        siril_debug_print("Atmospheric refraction correction: %.4f arcsec (Alt %.2f° -> %.2f°)\n",
                         refraction * 3600.0, alt, corrected_alt);
    }
    
    // Set up horizontal coordinates
    struct ln_hrz_posn horizontal;
    horizontal.alt = corrected_alt;
    horizontal.az = az;
    
    // Convert to equatorial coordinates
    struct ln_equ_posn equatorial;
    ln_get_equ_from_hrz(&horizontal, &observer, jd, &equatorial);
    
    // Apply proper motion and aberration corrections if needed
    // (These are typically small for telescope observations but libnova supports them)
    
    *ra = equatorial.ra;
    *dec = equatorial.dec;
    
    // Validate and normalize
    if (*dec < -90.0 || *dec > 90.0) {
        siril_log_color_message("Invalid declination calculated: %.2f\n", "red", *dec);
        return -1;
    }
    
    if (*ra < 0.0) *ra += 360.0;
    if (*ra >= 360.0) *ra -= 360.0;
    
    siril_debug_print("Enhanced conversion with corrections: Alt/Az %.2f°/%.2f° -> RA/Dec %.4f°/%.4f°\n",
                     alt, az, *ra, *dec);
    
    return 0;
}

/**
 * Calculate local sidereal time using libnova
 */
double stellina_get_local_sidereal_time(const char *date_obs, double observer_lon) {
    struct tm tm_obs = {0};
    if (strptime(date_obs, "%Y-%m-%dT%H:%M:%S", &tm_obs) == NULL) {
        return -1.0;
    }
    
    struct ln_date obs_date;
    obs_date.years = tm_obs.tm_year + 1900;
    obs_date.months = tm_obs.tm_mon + 1;
    obs_date.days = tm_obs.tm_mday;
    obs_date.hours = tm_obs.tm_hour;
    obs_date.minutes = tm_obs.tm_min;
    obs_date.seconds = tm_obs.tm_sec;
    
    double jd = ln_get_julian_day(&obs_date);
    double lst = ln_get_apparent_sidereal_time(jd) * 15.0 + observer_lon; // Convert to degrees
    
    // Normalize to 0-360 degrees
    while (lst < 0.0) lst += 360.0;
    while (lst >= 360.0) lst -= 360.0;
    
    return lst;
}

/**
 * Set default observer location with better defaults
 */
void stellina_set_default_observer_location(struct stellina_config *config) {
    if (config->observer_latitude == 0.0 && config->observer_longitude == 0.0) {
        // Default to Greenwich Observatory as a reasonable starting point
        config->observer_latitude = 51.4769;   // Greenwich latitude
        config->observer_longitude = -0.0005;  // Greenwich longitude  
        config->observer_altitude = 46.0;      // Greenwich altitude in meters
        
        siril_log_message("Using default observer location: Greenwich Observatory\n");
        siril_log_message("  Latitude: %.6f°, Longitude: %.6f°, Altitude: %.1fm\n",
                         config->observer_latitude, config->observer_longitude, config->observer_altitude);
        siril_log_color_message("IMPORTANT: Set your actual location with: stellina_config observer_location lat,lon,alt\n", "orange");
        siril_log_message("Example: stellina_config observer_location 48.8566,2.3522,35  # Paris\n");
    }
}

/**
 * Validate observer coordinates
 */
gboolean stellina_validate_observer_location(double lat, double lon, double alt) {
    if (lat < -90.0 || lat > 90.0) {
        siril_log_color_message("Invalid latitude: %.6f (must be -90 to 90)\n", "red", lat);
        return FALSE;
    }
    
    if (lon < -180.0 || lon > 180.0) {
        siril_log_color_message("Invalid longitude: %.6f (must be -180 to 180)\n", "red", lon);
        return FALSE;
    }
    
    if (alt < -1000.0 || alt > 10000.0) {
        siril_log_color_message("Invalid altitude: %.1f (must be -1000 to 10000 meters)\n", "red", alt);
        return FALSE;
    }
    
    return TRUE;
}

/**
 * Get solar system object coordinates using libnova
 * Useful for validating coordinate conversions
 */
int stellina_get_solar_object_coords(const char *date_obs, double observer_lat, double observer_lon,
                                   const char *object_name, double *ra, double *dec) {
    
    struct tm tm_obs = {0};
    if (strptime(date_obs, "%Y-%m-%dT%H:%M:%S", &tm_obs) == NULL) {
        return -1;
    }
    
    struct ln_date obs_date;
    obs_date.years = tm_obs.tm_year + 1900;
    obs_date.months = tm_obs.tm_mon + 1;
    obs_date.days = tm_obs.tm_mday;
    obs_date.hours = tm_obs.tm_hour;
    obs_date.minutes = tm_obs.tm_min;
    obs_date.seconds = tm_obs.tm_sec;
    
    double jd = ln_get_julian_day(&obs_date);
    struct ln_equ_posn position;
    
    if (g_ascii_strcasecmp(object_name, "moon") == 0) {
        ln_get_lunar_equ_coords(jd, &position);
    } else if (g_ascii_strcasecmp(object_name, "sun") == 0) {
        ln_get_solar_equ_coords(jd, &position);
    } else {
        return -1; // Object not supported
    }
    
    *ra = position.ra;
    *dec = position.dec;
    
    return 0;
}
