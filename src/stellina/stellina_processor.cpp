/*
 * stellina_processor.cpp
 * Main processing logic for Stellina images in Siril
 */

#include "stellina_processor.h"
#include "core/siril.h"
#include "io/single_image.h"
#include "io/fits_keywords.h"
#include "astrometry/platesolve.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include <dirent.h>
#include <fnmatch.h>

static stellina_progress_callback g_progress_callback = NULL;

/**
 * Set progress callback for GUI updates
 */
void stellina_set_progress_callback(stellina_progress_callback callback) {
    g_progress_callback = callback;
}

/**
 * Find matching FITS and JSON files in directory
 */
static GList *stellina_find_matching_files(const char *input_dir) {
    GList *matches = NULL;
    DIR *dir;
    struct dirent *entry;
    
    dir = opendir(input_dir);
    if (!dir) {
        siril_log_color_message("Cannot open directory: %s\n", "red", input_dir);
        return NULL;
    }
    
    // First pass: find all FITS files
    GList *fits_files = NULL;
    rewinddir(dir);
    while ((entry = readdir(dir)) != NULL) {
        if (fnmatch("*.fit*", entry->d_name, FNM_CASEFOLD) == 0) {
            char *full_path = g_build_filename(input_dir, entry->d_name, NULL);
            fits_files = g_list_prepend(fits_files, full_path);
        }
    }
    
    // Second pass: find matching JSON files for each FITS
    for (GList *l = fits_files; l != NULL; l = l->next) {
        char *fits_path = (char *)l->data;
        char *basename = g_path_get_basename(fits_path);
        
        // Remove extension and look for matching JSON
        char *dot = strrchr(basename, '.');
        if (dot) *dot = '\0';
        
        char *json_name = g_strdup_printf("%s.json", basename);
        char *json_path = g_build_filename(input_dir, json_name, NULL);
        
        if (g_file_test(json_path, G_FILE_TEST_EXISTS)) {
            // Create match pair
            char **pair = g_new(char *, 2);
            pair[0] = g_strdup(fits_path);  // FITS path
            pair[1] = g_strdup(json_path);  // JSON path
            matches = g_list_prepend(matches, pair);
        } else {
            siril_debug_print("No matching JSON found for %s\n", fits_path);
        }
        
        g_free(basename);
        g_free(json_name);
        g_free(json_path);
    }
    
    closedir(dir);
    g_list_free_full(fits_files, g_free);
    
    // Sort matches by filename
    matches = g_list_reverse(matches);
    
    siril_log_message("Found %d matching FITS/JSON pairs\n", g_list_length(matches));
    return matches;
}

/**
 * Process a directory of Stellina images
 */
int stellina_process_directory(const char *input_dir, const char *output_dir, struct stellina_config *config) {
    if (!input_dir || !output_dir || !config) {
        return -1;
    }
    
    // Set default observer location if not set
    stellina_set_default_observer_location(config);
    
    // Find matching files
    GList *matches = stellina_find_matching_files(input_dir);
    if (!matches) {
        siril_log_color_message("No matching FITS/JSON pairs found in %s\n", "red", input_dir);
        return -1;
    }
    
    struct stellina_stats stats = {0};
    stats.total_images = g_list_length(matches);
    
    int current = 0;
    for (GList *l = matches; l != NULL; l = l->next) {
        char **pair = (char **)l->data;
        char *fits_path = pair[0];
        char *json_path = pair[1];
        
        current++;
        
        if (g_progress_callback) {
            char *basename = g_path_get_basename(fits_path);
            char *message = g_strdup_printf("Processing %s", basename);
            g_progress_callback(current, stats.total_images, message);
            g_free(basename);
            g_free(message);
        }
        
        // Generate output path
        char *basename = g_path_get_basename(fits_path);
        char *output_name = g_strdup_printf("%s%s", config->output_prefix, basename);
        char *output_path = g_build_filename(output_dir, output_name, NULL);
        
        // Process single image
        int result = stellina_process_single_image(fits_path, json_path, output_path, config);
        
        if (result == 0) {
            stats.processed_images++;
        } else if (result == -2) {
            stats.quality_rejected++;
        } else {
            char *error_msg = g_strdup_printf("Failed to process %s", basename);
            stats.error_messages = g_list_prepend(stats.error_messages, error_msg);
        }
        
        g_free(basename);
        g_free(output_name);
        g_free(output_path);
    }
    
    // Print summary
    siril_log_color_message("Stellina Processing Summary:\n", "green");
    siril_log_message("  Total images: %d\n", stats.total_images);
    siril_log_message("  Successfully processed: %d\n", stats.processed_images);
    siril_log_message("  Quality rejected: %d\n", stats.quality_rejected);
    siril_log_message("  Plate solve succeeded: %d\n", stats.platesolve_succeeded);
    siril_log_message("  Plate solve failed: %d\n", stats.platesolve_failed);
    
    if (stats.error_messages) {
        siril_log_color_message("Errors encountered:\n", "red");
        for (GList *l = stats.error_messages; l != NULL; l = l->next) {
            siril_log_message("  %s\n", (char *)l->data);
        }
    }
    
    // Cleanup
    for (GList *l = matches; l != NULL; l = l->next) {
        char **pair = (char **)l->data;
        g_free(pair[0]);
        g_free(pair[1]);
        g_free(pair);
    }
    g_list_free(matches);
    g_list_free_full(stats.error_messages, g_free);
    
    return (stats.processed_images > 0) ? 0 : -1;
}

/**
 * Process a single Stellina image
 */
int stellina_process_single_image(const char *fits_path, const char *json_path, const char *output_path, struct stellina_config *config) {
    if (!fits_path || !json_path || !output_path || !config) {
        return -1;
    }
    
    // Parse JSON metadata
    struct stellina_metadata *metadata = stellina_parse_json(json_path);
    if (!metadata) {
        siril_log_color_message("Failed to parse JSON: %s\n", "red", json_path);
        return -1;
    }
    
    // Check quality if filtering enabled
    if (config->quality_filter_enabled && !stellina_check_quality(metadata)) {
        siril_debug_print("Skipping %s: quality rejected\n", fits_path);
        stellina_metadata_free(metadata);
        return -2; // Quality rejected
    }
    
    // Load the FITS image
    fits fit = {0};
    if (readfits(fits_path, &fit, NULL, com.pref.force_16bit) != READFITS_OK) {
        siril_log_color_message("Failed to load FITS file: %s\n", "red", fits_path);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    // Convert Alt/Az to RA/Dec
    double ra, dec;
    int coord_result = stellina_convert_altaz_to_radec(
        metadata->altitude, metadata->azimuth, metadata->date_obs,
        config->observer_latitude, config->observer_longitude, config->observer_altitude,
        &ra, &dec);
    
    if (coord_result != 0) {
        siril_log_color_message("Failed coordinate conversion for %s\n", "red", fits_path);
        clearfits(&fit);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    siril_log_message("Processing %s: Alt/Az %.2f°/%.2f° -> RA/Dec %.4f°/%.4f°\n",
                     g_path_get_basename(fits_path), metadata->altitude, metadata->azimuth, ra, dec);
    
    // Attempt plate solving if enabled
    gboolean platesolve_success = FALSE;
    if (config->platesolve_fallback_enabled) {
        // TODO: Integrate with Siril's plate solving
        // This would call Siril's internal plate solving functions
        siril_debug_print("Plate solving not yet implemented in this version\n");
    }
    
    // Add coordinate information to FITS header
    if (config->add_coordinate_keywords) {
        stellina_add_coordinates_to_header(&fit, metadata, ra, dec);
    }
    
    // Save processed image
    if (savefits(output_path, &fit) != 0) {
        siril_log_color_message("Failed to save processed image: %s\n", "red", output_path);
        clearfits(&fit);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    siril_log_message("Saved processed image: %s\n", output_path);
    
    // Cleanup
    clearfits(&fit);
    stellina_metadata_free(metadata);
    
    return 0;
}

/**
 * Add coordinate information to FITS header
 */
int stellina_add_coordinates_to_header(fits *fit, const struct stellina_metadata *metadata, double ra, double dec) {
    if (!fit || !fit->fptr || !metadata) {
        return -1;
    }
    
    int status = 0;
    
    // Add calculated coordinates
    fits_update_key(fit->fptr, TDOUBLE, "RA_CALC", &ra, "Calculated RA (degrees)", &status);
    fits_update_key(fit->fptr, TDOUBLE, "DEC_CALC", &dec, "Calculated Dec (degrees)", &status);
    
    // Add original Alt/Az
    fits_update_key(fit->fptr, TDOUBLE, "ALT_ORIG", &metadata->altitude, "Original altitude (degrees)", &status);
    fits_update_key(fit->fptr, TDOUBLE, "AZ_ORIG", &metadata->azimuth, "Original azimuth (degrees)", &status);
    
    // Add processing information
    char stellina_info[] = "PROCESSED";
    fits_update_key(fit->fptr, TSTRING, "STELLINA", stellina_info, "Processed by Stellina extension", &status);
    
    // Add quality information
    char quality_status[32];
    g_strlcpy(quality_status, metadata->quality_accepted ? "ACCEPTED" : "REJECTED", sizeof(quality_status));
    fits_update_key(fit->fptr, TSTRING, "STLQUAL", quality_status, "Stellina quality assessment", &status);
    
    if (metadata->fwhm > 0.0) {
        fits_update_key(fit->fptr, TDOUBLE, "STL_FWHM", &metadata->fwhm, "Stellina FWHM (pixels)", &status);
    }
    
    if (metadata->star_count > 0) {
        fits_update_key(fit->fptr, TINT, "STL_STARS", &metadata->star_count, "Stellina star count", &status);
    }
    
    if (status != 0) {
        siril_log_color_message("Warning: Some FITS header keywords could not be updated (status=%d)\n", "orange", status);
        return -1;
    }
    
    return 0;
}

/**
 * Free stellina stats structure
 */
void stellina_stats_free(struct stellina_stats *stats) {
    if (stats) {
        g_list_free_full(stats->error_messages, g_free);
        g_list_free_full(stats->warning_messages, g_free);
        stats->error_messages = NULL;
        stats->warning_messages = NULL;
    }
}