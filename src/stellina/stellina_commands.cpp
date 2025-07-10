/*
 * stellina_commands.cpp
 * Clean final Siril command implementations for Stellina processing
 */

#include "stellina_processor.h"
#include "stellina_commands.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <fitsio.h>

// Define Siril command constants if not available
#ifndef CMD_OK
#define CMD_OK 0
#endif
#ifndef CMD_ARG_ERROR
#define CMD_ARG_ERROR 1
#endif
#ifndef CMD_GENERIC_ERROR
#define CMD_GENERIC_ERROR -1
#endif

// Progress callback function
static void stellina_progress_handler(int current, int total, const char* message) {
    char progress_msg[256];
    snprintf(progress_msg, sizeof(progress_msg), "Stellina: %s (%d/%d)", message, current, total);
    siril_log_message("%s\n", progress_msg);
}

// Global configuration instance
static struct stellina_config g_stellina_config = {
    .quality_filter_enabled = TRUE,
    .platesolve_fallback_enabled = TRUE,
    .add_coordinate_keywords = TRUE,
    .observer_latitude = 0.0,
    .observer_longitude = 0.0,
    .observer_altitude = 0.0,
    .output_prefix = "stellina_"
};

/**
 * Main Stellina processing command
 */
int cmd_process_stellina(int nb) {
    if (nb < 2) {
        siril_log_color_message("Usage: stellina_process input_dir output_dir [quality_filter=true]\n", "red");
        siril_log_message("  input_dir: Directory containing FITS and JSON files\n");
        siril_log_message("  output_dir: Directory for processed files\n");
        siril_log_message("  quality_filter: Enable quality filtering (true/false, default: true)\n");
        siril_log_message("Example: stellina_process ~/stellina_raw ~/stellina_processed\n");
        return CMD_ARG_ERROR;
    }
    
    char *input_dir = word[1];
    char *output_dir = word[2];
    
    if (nb >= 3) {
        if (g_ascii_strcasecmp(word[3], "false") == 0 || g_ascii_strcasecmp(word[3], "0") == 0) {
            g_stellina_config.quality_filter_enabled = FALSE;
        } else {
            g_stellina_config.quality_filter_enabled = TRUE;
        }
    }
    
    siril_log_message("Starting Stellina batch processing...\n");
    siril_log_message("Input directory: %s\n", input_dir);
    siril_log_message("Output directory: %s\n", output_dir);
    siril_log_message("Quality filtering: %s\n", 
                     g_stellina_config.quality_filter_enabled ? "enabled" : "disabled");
    
    stellina_set_progress_callback(stellina_progress_handler);
    
    int result = stellina_process_directory(input_dir, output_dir, &g_stellina_config);
    
    if (result == 0) {
        siril_log_color_message("Stellina processing completed successfully\n", "green");
        return CMD_OK;
    } else {
        siril_log_color_message("Stellina processing failed\n", "red");
        return CMD_GENERIC_ERROR;
    }
}

/**
 * Stellina configuration command
 */
int cmd_stellina_config(int nb) {
    if (nb < 2) {
        siril_log_color_message("Usage: stellina_config parameter value\n", "red");
        siril_log_message("Available parameters:\n");
        siril_log_message("  quality_filter <true|false>     - Enable/disable quality filtering\n");
        siril_log_message("  platesolve_fallback <true|false> - Enable/disable platesolve fallback\n");
        siril_log_message("  add_keywords <true|false>       - Add coordinate keywords to FITS\n");
        siril_log_message("  observer_location <lat,lon,alt>  - Set observer location (degrees,degrees,meters)\n");
        siril_log_message("  output_prefix <prefix>           - Set output filename prefix\n");
        siril_log_message("  show                            - Show current configuration\n");
        siril_log_message("\nExamples:\n");
        siril_log_message("  stellina_config quality_filter false\n");
        siril_log_message("  stellina_config observer_location 48.8566,2.3522,35\n");
        siril_log_message("  stellina_config output_prefix processed_\n");
        return CMD_ARG_ERROR;
    }
    
    char *parameter = word[1];
    
    if (g_ascii_strcasecmp(parameter, "show") == 0) {
        siril_log_color_message("Current Stellina Configuration:\n", "green");
        siril_log_message("  Quality filtering: %s\n", 
                         g_stellina_config.quality_filter_enabled ? "enabled" : "disabled");
        siril_log_message("  Platesolve fallback: %s\n", 
                         g_stellina_config.platesolve_fallback_enabled ? "enabled" : "disabled");
        siril_log_message("  Add coordinate keywords: %s\n", 
                         g_stellina_config.add_coordinate_keywords ? "enabled" : "disabled");
        siril_log_message("  Observer location: %.6f°, %.6f°, %.1fm\n",
                         g_stellina_config.observer_latitude, 
                         g_stellina_config.observer_longitude,
                         g_stellina_config.observer_altitude);
        siril_log_message("  Output prefix: %s\n", g_stellina_config.output_prefix);
        return CMD_OK;
    }
    
    if (nb < 3) {
        siril_log_color_message("Error: Missing value for parameter %s\n", "red", parameter);
        return CMD_ARG_ERROR;
    }
    
    char *value = word[2];
    
    if (g_ascii_strcasecmp(parameter, "quality_filter") == 0) {
        gboolean enabled = (g_ascii_strcasecmp(value, "true") == 0 || g_ascii_strcasecmp(value, "1") == 0);
        g_stellina_config.quality_filter_enabled = enabled;
        siril_log_message("Quality filtering %s\n", enabled ? "enabled" : "disabled");
        
    } else if (g_ascii_strcasecmp(parameter, "platesolve_fallback") == 0) {
        gboolean enabled = (g_ascii_strcasecmp(value, "true") == 0 || g_ascii_strcasecmp(value, "1") == 0);
        g_stellina_config.platesolve_fallback_enabled = enabled;
        siril_log_message("Platesolve fallback %s\n", enabled ? "enabled" : "disabled");
        
    } else if (g_ascii_strcasecmp(parameter, "add_keywords") == 0) {
        gboolean enabled = (g_ascii_strcasecmp(value, "true") == 0 || g_ascii_strcasecmp(value, "1") == 0);
        g_stellina_config.add_coordinate_keywords = enabled;
        siril_log_message("Coordinate keywords %s\n", enabled ? "enabled" : "disabled");
        
    } else if (g_ascii_strcasecmp(parameter, "observer_location") == 0) {
        double lat, lon, alt;
        int parsed = sscanf(value, "%lf,%lf,%lf", &lat, &lon, &alt);
        if (parsed != 3) {
            siril_log_color_message("Error: Invalid observer location format. Use: lat,lon,alt\n", "red");
            siril_log_message("Example: 48.8566,2.3522,35\n");
            return CMD_ARG_ERROR;
        }
        
        if (!stellina_validate_observer_location(lat, lon, alt)) {
            return CMD_ARG_ERROR;
        }
        
        g_stellina_config.observer_latitude = lat;
        g_stellina_config.observer_longitude = lon;
        g_stellina_config.observer_altitude = alt;
        siril_log_color_message("Observer location set to: %.6f°, %.6f°, %.1fm\n", "green", lat, lon, alt);
        
    } else if (g_ascii_strcasecmp(parameter, "output_prefix") == 0) {
        if (strlen(value) >= sizeof(g_stellina_config.output_prefix)) {
            siril_log_color_message("Error: Output prefix too long (max %zu characters)\n", "red", 
                                   sizeof(g_stellina_config.output_prefix) - 1);
            return CMD_ARG_ERROR;
        }
        g_strlcpy(g_stellina_config.output_prefix, value, sizeof(g_stellina_config.output_prefix));
        siril_log_message("Output prefix set to: %s\n", value);
        
    } else {
        siril_log_color_message("Error: Unknown parameter '%s'\n", "red", parameter);
        siril_log_message("Use 'stellina_config' without arguments to see available parameters\n");
        return CMD_ARG_ERROR;
    }
    
    return CMD_OK;
}

/**
 * Simple FITS file loader for stats command
 */
static int load_fits_for_stats(const char* filename, fits* fit) {
    if (!filename || !fit) return -1;
    
    memset(fit, 0, sizeof(fits));
    
    int status = 0;
    fits_open_file(&fit->fptr, filename, READONLY, &status);
    if (status != 0) {
        return -1;
    }
    
    int naxis;
    long naxes[3] = {0, 0, 1};
    fits_get_img_param(fit->fptr, 3, &fit->bitpix, &naxis, naxes, &status);
    if (status != 0) {
        fits_close_file(fit->fptr, &status);
        return -1;
    }
    
    fit->rx = naxes[0];
    fit->ry = naxes[1];
    fit->naxes[2] = (naxis > 2) ? naxes[2] : 1;
    
    if (fit->bitpix == 16) {
        fit->type = DATA_USHORT;
    } else if (fit->bitpix == -32) {
        fit->type = DATA_FLOAT;
    } else {
        fit->type = DATA_USHORT;
    }
    
    return 0;
}

/**
 * Close FITS file for stats
 */
static void close_fits_for_stats(fits* fit) {
    if (fit && fit->fptr) {
        int status = 0;
        fits_close_file(fit->fptr, &status);
        fit->fptr = nullptr;
    }
}

/**
 * Stellina statistics and information command
 */
int cmd_stellina_stats(int nb) {
    if (nb == 1) {
        siril_log_color_message("Stellina Extension Information:\n", "green");
        siril_log_message("  Version: 1.0\n");
        siril_log_message("  Purpose: Process Stellina telescope FITS/JSON pairs\n");
        siril_log_message("  Features:\n");
        siril_log_message("    - Alt/Az to RA/Dec coordinate conversion using libnova\n");
        siril_log_message("    - Quality filtering based on Stellina metadata\n");
        siril_log_message("    - Atmospheric refraction correction\n");
        siril_log_message("    - FITS header annotation\n");
        siril_log_message("    - Batch processing with progress reporting\n");
        siril_log_message("\nCommands:\n");
        siril_log_message("  stellina_process - Process directory of FITS/JSON pairs\n");
        siril_log_message("  stellina_config  - Configure processing parameters\n");
        siril_log_message("  stellina_stats   - Show statistics and test files\n");
        siril_log_message("  stellina_test_coords - Test coordinate conversion\n");
        return CMD_OK;
    }
    
    if (nb == 3) {
        char *fits_file = word[1];
        char *json_file = word[2];
        
        siril_log_color_message("Analyzing Stellina file pair:\n", "green");
        siril_log_message("FITS: %s\n", fits_file);
        siril_log_message("JSON: %s\n", json_file);
        
        if (!g_file_test(fits_file, G_FILE_TEST_EXISTS)) {
            siril_log_color_message("Error: FITS file not found: %s\n", "red", fits_file);
            return CMD_GENERIC_ERROR;
        }
        
        if (!g_file_test(json_file, G_FILE_TEST_EXISTS)) {
            siril_log_color_message("Error: JSON file not found: %s\n", "red", json_file);
            return CMD_GENERIC_ERROR;
        }
        
        struct stellina_metadata *metadata = stellina_parse_json_enhanced(json_file);
        if (!metadata) {
            siril_log_color_message("Failed to parse JSON metadata\n", "red");
            return CMD_GENERIC_ERROR;
        }
        
        siril_log_color_message("\nStellina Metadata:\n", "blue");
        siril_log_message("  Altitude: %.2f°\n", metadata->altitude);
        siril_log_message("  Azimuth: %.2f°\n", metadata->azimuth);
        siril_log_message("  Observation time: %s\n", 
                         strlen(metadata->date_obs) > 0 ? metadata->date_obs : "Not found");
        siril_log_message("  Quality: %s (%s)\n", 
                         metadata->quality_accepted ? "ACCEPTED" : "REJECTED",
                         metadata->quality_reason);
        if (metadata->fwhm > 0.0) {
            siril_log_message("  FWHM: %.2f arcsec\n", metadata->fwhm);
        }
        if (metadata->star_count > 0) {
            siril_log_message("  Star count: %d\n", metadata->star_count);
        }
        siril_log_message("  Used for stacking: %s\n", 
                         metadata->used_for_stacking ? "Yes" : "No");
        
        if (metadata->altitude > 0.0 && metadata->azimuth >= 0.0 && strlen(metadata->date_obs) > 0) {
            stellina_set_default_observer_location(&g_stellina_config);
            
            double ra, dec;
            int result = stellina_convert_altaz_to_radec(
                metadata->altitude, metadata->azimuth, metadata->date_obs,
                g_stellina_config.observer_latitude, 
                g_stellina_config.observer_longitude,
                g_stellina_config.observer_altitude,
                &ra, &dec
            );
            
            if (result == 0) {
                siril_log_color_message("\nCalculated Coordinates:\n", "blue");
                siril_log_message("  RA: %.6f° (%.2fh %.2fm %.2fs)\n", ra, 
                                 ra/15.0, fmod(ra*4.0, 60.0), fmod(ra*240.0, 60.0));
                siril_log_message("  Dec: %.6f° (%+.0f° %.0f' %.1f\")\n", dec,
                                 floor(dec), fmod(fabs(dec)*60.0, 60.0), fmod(fabs(dec)*3600.0, 60.0));
                
                double lst = stellina_get_local_sidereal_time(metadata->date_obs, 
                                                            g_stellina_config.observer_longitude);
                if (lst >= 0.0) {
                    siril_log_message("  Local Sidereal Time: %.2f°\n", lst);
                }
            } else {
                siril_log_color_message("Failed to convert coordinates\n", "red");
            }
        } else {
            siril_log_color_message("Insufficient data for coordinate conversion\n", "orange");
        }
        
        fits fit = { 0 };
        if (load_fits_for_stats(fits_file, &fit) == 0) {
            siril_log_color_message("\nFITS File Information:\n", "blue");
            siril_log_message("  Dimensions: %d x %d\n", fit.rx, fit.ry);
            siril_log_message("  Channels: %d\n", fit.naxes[2]);
            siril_log_message("  Data type: %s\n", fit.type == DATA_USHORT ? "16-bit" : 
                             fit.type == DATA_FLOAT ? "32-bit float" : "Other");
            
            if (fit.fptr) {
                double existing_ra, existing_dec;
                int ra_status = 0, dec_status = 0;
                fits_read_key(fit.fptr, TDOUBLE, "RA", &existing_ra, NULL, &ra_status);
                fits_read_key(fit.fptr, TDOUBLE, "DEC", &existing_dec, NULL, &dec_status);
                
                if (ra_status == 0 && dec_status == 0) {
                    siril_log_message("  Existing RA/Dec: %.6f°, %.6f°\n", existing_ra, existing_dec);
                } else {
                    siril_log_message("  No existing RA/Dec coordinates found\n");
                }
            }
            
            close_fits_for_stats(&fit);
        } else {
            siril_log_color_message("Failed to read FITS file\n", "red");
        }
        
        stellina_metadata_free(metadata);
        return CMD_OK;
    }
    
    siril_log_color_message("Usage: stellina_stats [fits_file json_file]\n", "red");
    siril_log_message("  With no arguments: Show general information\n");
    siril_log_message("  With files: Analyze specific FITS/JSON pair\n");
    return CMD_ARG_ERROR;
}

/**
 * Test coordinate conversion command
 */
int cmd_stellina_test_coords(int nb) {
    if (nb < 4) {
        siril_log_color_message("Usage: stellina_test_coords alt az date_obs [lat lon alt_m]\n", "red");
        siril_log_message("  alt: Altitude in degrees\n");
        siril_log_message("  az: Azimuth in degrees\n");
        siril_log_message("  date_obs: Observation time (ISO format: 2024-01-09T22:13:29)\n");
        siril_log_message("  lat, lon, alt_m: Observer location (optional, uses config if not provided)\n");
        siril_log_message("Example: stellina_test_coords 45.5 180.0 2024-01-09T22:13:29\n");
        return CMD_ARG_ERROR;
    }
    
    double alt = g_ascii_strtod(word[1], NULL);
    double az = g_ascii_strtod(word[2], NULL);
    char *date_obs = word[3];
    
    double obs_lat = g_stellina_config.observer_latitude;
    double obs_lon = g_stellina_config.observer_longitude;
    double obs_alt = g_stellina_config.observer_altitude;
    
    if (nb >= 7) {
        obs_lat = g_ascii_strtod(word[4], NULL);
        obs_lon = g_ascii_strtod(word[5], NULL);
        obs_alt = g_ascii_strtod(word[6], NULL);
    } else {
        stellina_set_default_observer_location(&g_stellina_config);
        obs_lat = g_stellina_config.observer_latitude;
        obs_lon = g_stellina_config.observer_longitude;
        obs_alt = g_stellina_config.observer_altitude;
    }
    
    siril_log_color_message("Testing coordinate conversion:\n", "green");
    siril_log_message("Input:\n");
    siril_log_message("  Altitude: %.2f°\n", alt);
    siril_log_message("  Azimuth: %.2f°\n", az);
    siril_log_message("  Time: %s\n", date_obs);
    siril_log_message("  Observer: %.6f°, %.6f°, %.1fm\n", obs_lat, obs_lon, obs_alt);
    
    double ra, dec;
    int result = stellina_convert_altaz_to_radec(alt, az, date_obs, obs_lat, obs_lon, obs_alt, &ra, &dec);
    
    if (result == 0) {
        siril_log_color_message("Result:\n", "green");
        siril_log_message("  RA: %.6f° (%.2fh %.2fm %.2fs)\n", ra, 
                         ra/15.0, fmod(ra*4.0, 60.0), fmod(ra*240.0, 60.0));
        siril_log_message("  Dec: %.6f° (%+.0f° %.0f' %.1f\")\n", dec,
                         floor(dec), fmod(fabs(dec)*60.0, 60.0), fmod(fabs(dec)*3600.0, 60.0));
        
        double lst = stellina_get_local_sidereal_time(date_obs, obs_lon);
        if (lst >= 0.0) {
            siril_log_message("  Local Sidereal Time: %.2f°\n", lst);
        }
        
        return CMD_OK;
    } else {
        siril_log_color_message("Coordinate conversion failed\n", "red");
        return CMD_GENERIC_ERROR;
    }
}

/**
 * Initialize Stellina extension commands
 */
void stellina_commands_init() {
    siril_log_color_message("Initializing Stellina extension commands...\n", "green");
    siril_log_message("Available Stellina commands:\n");
    siril_log_message("  stellina_process - Process directory of FITS/JSON pairs\n");
    siril_log_message("  stellina_config  - Configure processing parameters\n");
    siril_log_message("  stellina_stats   - Show statistics and analyze files\n");
    siril_log_message("  stellina_test_coords - Test coordinate conversion\n");
    siril_log_message("Use 'help <command>' for detailed usage information\n");
}