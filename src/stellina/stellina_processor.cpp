/*
 * stellina_processor.cpp
 * Final fixed Stellina processing implementation using CFITSIO directly
 */

#include "stellina_processor.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "gui/utils.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <vector>
#include <regex>
#include <cmath>
#include <fitsio.h>

static stellina_progress_callback g_progress_callback = nullptr;

/**
 * Set progress callback
 */
void stellina_set_progress_callback(stellina_progress_callback callback) {
    g_progress_callback = callback;
}

/**
 * C++17 compatible ends_with function
 */
bool string_ends_with(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

/**
 * Simple FITS file loader using CFITSIO directly
 */
int load_fits_file(const char* filename, fits* fit) {
    if (!filename || !fit) return -1;
    
    memset(fit, 0, sizeof(fits));
    
    int status = 0;
    fits_open_file(&fit->fptr, filename, READWRITE, &status);
    if (status != 0) {
        // Try read-only if readwrite fails
        status = 0;
        fits_open_file(&fit->fptr, filename, READONLY, &status);
        if (status != 0) {
            char err_text[FLEN_ERRMSG];
            fits_get_errstatus(status, err_text);
            siril_log_color_message("FITS open error: %s\n", "red", err_text);
            return -1;
        }
    }
    
    // Get image dimensions
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
    fit->naxis = naxis;
    
    // Set data type based on bitpix
    if (fit->bitpix == 16) {
        fit->type = DATA_USHORT;
    } else if (fit->bitpix == -32) {
        fit->type = DATA_FLOAT;
    } else {
        fit->type = DATA_USHORT; // default
    }
    
    return 0;
}

/**
 * Simple FITS file saver using CFITSIO directly
 */
int save_fits_file(const char* output_path, fits* fit) {
    if (!output_path || !fit || !fit->fptr) return -1;
    
    int status = 0;
    fitsfile* outfptr;
    
    // Remove file if it exists
    remove(output_path);
    
    // Create new file
    fits_create_file(&outfptr, output_path, &status);
    if (status != 0) {
        char err_text[FLEN_ERRMSG];
        fits_get_errstatus(status, err_text);
        siril_log_color_message("FITS create error: %s\n", "red", err_text);
        return -1;
    }
    
    // Copy header from input to output
    fits_copy_file(fit->fptr, outfptr, 1, 1, 1, &status);
    if (status != 0) {
        char err_text[FLEN_ERRMSG];
        fits_get_errstatus(status, err_text);
        siril_log_color_message("FITS copy error: %s\n", "red", err_text);
        fits_close_file(outfptr, &status);
        return -1;
    }
    
    // Close output file
    fits_close_file(outfptr, &status);
    return (status == 0) ? 0 : -1;
}

/**
 * Close FITS file
 */
void close_fits_file(fits* fit) {
    if (fit && fit->fptr) {
        int status = 0;
        fits_close_file(fit->fptr, &status);
        fit->fptr = nullptr;
    }
}

/**
 * Find matching FITS and JSON files in directory
 * Modified to handle Stellina's specific naming pattern
 */
std::vector<std::pair<std::string, std::string>> find_matching_files(const std::string& directory) {
    std::vector<std::pair<std::string, std::string>> matches;
    
    try {
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            if (!entry.is_regular_file()) continue;
            
            std::string path = entry.path().string();
            std::string filename = entry.path().filename().string();
            
            // Look for FITS files
            if (string_ends_with(filename, ".fits") || string_ends_with(filename, ".fit")) {
                // Extract base pattern (e.g., "img-0001" from "img-0001.fits")
                std::string base_name = filename.substr(0, filename.find_last_of('.'));
                
                // Generate corresponding JSON filename patterns
                std::string stacking_json_path = directory + "/" + base_name + "-stacking.json";
                std::string output_json_path = directory + "/" + base_name + "-output.json";
                
                // Check for stacking JSON first (preferred)
                if (std::filesystem::exists(stacking_json_path)) {
                    matches.emplace_back(stacking_json_path, path);
                    siril_debug_print("Found matching pair: %s <-> %s\n", 
                                    stacking_json_path.c_str(), path.c_str());
                }
                // Fall back to output JSON if stacking not found
                else if (std::filesystem::exists(output_json_path)) {
                    matches.emplace_back(output_json_path, path);
                    siril_debug_print("Found matching pair: %s <-> %s\n", 
                                    output_json_path.c_str(), path.c_str());
                }
            }
        }
    } catch (const std::exception& e) {
        siril_log_color_message("Error scanning directory %s: %s\n", "red", 
                               directory.c_str(), e.what());
    }
    
    siril_log_message("Found %zu matching FITS/JSON pairs\n", matches.size());
    return matches;
}

/**
 * Extract DATE-OBS from FITS header
 */
bool extract_date_from_fits(const char* fits_path, char* date_obs, size_t date_obs_size) {
    if (!fits_path || !date_obs || date_obs_size == 0) {
        return false;
    }
    
    // Initialize output
    date_obs[0] = '\0';
    
    // Open FITS file
    fitsfile* fptr = nullptr;
    int status = 0;
    
    if (fits_open_file(&fptr, fits_path, READONLY, &status) != 0) {
        char err_text[FLEN_ERRMSG];
        fits_get_errstatus(status, err_text);
        siril_log_color_message("FITS open error: %s\n", "red", err_text);
        return false;
    }
    
    // Read DATE-OBS
    char date_value[256] = {0};
    status = 0;
    
    if (fits_read_key_str(fptr, "DATE-OBS", date_value, nullptr, &status) != 0) {
        siril_log_color_message("DATE-OBS not found in FITS header\n", "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Copy value
    g_strlcpy(date_obs, date_value, date_obs_size);
    
    // Close file
    fits_close_file(fptr, &status);
    
    siril_debug_print("Extracted DATE-OBS from FITS: %s\n", date_obs);
    return true;
}

/**
 * Enhanced JSON parsing with better error handling
 */
struct stellina_metadata *stellina_parse_json_enhanced(const char *json_path) {
    if (!json_path || !g_file_test(json_path, G_FILE_TEST_EXISTS)) {
        siril_log_color_message("JSON file not found: %s\n", "red", json_path ? json_path : "null");
        return nullptr;
    }
    
    struct stellina_metadata *metadata = g_new0(struct stellina_metadata, 1);
    
    // Initialize with defaults
    metadata->quality_accepted = TRUE;
    g_strlcpy(metadata->quality_reason, "No quality info found - assuming accepted", 
              sizeof(metadata->quality_reason));
    metadata->fwhm = 0.0;
    metadata->star_count = 0;
    metadata->used_for_stacking = TRUE;
    
    JsonParser *parser = json_parser_new();
    GError *error = nullptr;
    
    if (!json_parser_load_from_file(parser, json_path, &error)) {
        siril_log_color_message("Failed to parse JSON file %s: %s\n", "red", 
                               json_path, error->message);
        g_error_free(error);
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return nullptr;
    }
    
    JsonNode *root = json_parser_get_root(parser);
    if (!JSON_NODE_HOLDS_OBJECT(root)) {
        siril_log_color_message("Invalid JSON format in %s\n", "red", json_path);
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return nullptr;
    }
    
    JsonObject *root_object = json_node_get_object(root);
    
    // Extract altitude (try multiple field names)
    if (json_object_has_member(root_object, "altitude")) {
        metadata->altitude = json_object_get_double_member(root_object, "altitude");
    } else if (json_object_has_member(root_object, "alt")) {
        metadata->altitude = json_object_get_double_member(root_object, "alt");
    } else if (json_object_has_member(root_object, "Alt")) {
        metadata->altitude = json_object_get_double_member(root_object, "Alt");
    } else if (json_object_has_member(root_object, "motors") && 
               JSON_NODE_HOLDS_OBJECT(json_object_get_member(root_object, "motors"))) {
        JsonObject *motors = json_object_get_object_member(root_object, "motors");
        if (json_object_has_member(motors, "ALT")) {
            metadata->altitude = json_object_get_double_member(motors, "ALT");
        }
    } else if (json_object_has_member(root_object, "telescope") && 
               JSON_NODE_HOLDS_OBJECT(json_object_get_member(root_object, "telescope"))) {
        JsonObject *telescope = json_object_get_object_member(root_object, "telescope");
        if (json_object_has_member(telescope, "motors") && 
            JSON_NODE_HOLDS_OBJECT(json_object_get_member(telescope, "motors"))) {
            JsonObject *motors = json_object_get_object_member(telescope, "motors");
            if (json_object_has_member(motors, "ALT")) {
                metadata->altitude = json_object_get_double_member(motors, "ALT");
            }
        }
    }
    
    // Extract azimuth (try multiple field names)
    if (json_object_has_member(root_object, "azimuth")) {
        metadata->azimuth = json_object_get_double_member(root_object, "azimuth");
    } else if (json_object_has_member(root_object, "az")) {
        metadata->azimuth = json_object_get_double_member(root_object, "az");
    } else if (json_object_has_member(root_object, "Az")) {
        metadata->azimuth = json_object_get_double_member(root_object, "Az");
    } else if (json_object_has_member(root_object, "motors") && 
               JSON_NODE_HOLDS_OBJECT(json_object_get_member(root_object, "motors"))) {
        JsonObject *motors = json_object_get_object_member(root_object, "motors");
        if (json_object_has_member(motors, "AZ")) {
            metadata->azimuth = json_object_get_double_member(motors, "AZ");
        }
    } else if (json_object_has_member(root_object, "telescope") && 
               JSON_NODE_HOLDS_OBJECT(json_object_get_member(root_object, "telescope"))) {
        JsonObject *telescope = json_object_get_object_member(root_object, "telescope");
        if (json_object_has_member(telescope, "motors") && 
            JSON_NODE_HOLDS_OBJECT(json_object_get_member(telescope, "motors"))) {
            JsonObject *motors = json_object_get_object_member(telescope, "motors");
            if (json_object_has_member(motors, "AZ")) {
                metadata->azimuth = json_object_get_double_member(motors, "AZ");
            }
        }
    }
    
    // Extract observation time (try multiple field names)
    const char *date_str = nullptr;
    if (json_object_has_member(root_object, "date_obs")) {
        date_str = json_object_get_string_member(root_object, "date_obs");
    } else if (json_object_has_member(root_object, "time")) {
        date_str = json_object_get_string_member(root_object, "time");
    } else if (json_object_has_member(root_object, "timestamp")) {
        date_str = json_object_get_string_member(root_object, "timestamp");
    } else if (json_object_has_member(root_object, "DATE-OBS")) {
        date_str = json_object_get_string_member(root_object, "DATE-OBS");
    }
    
    if (date_str) {
        g_strlcpy(metadata->date_obs, date_str, sizeof(metadata->date_obs));
    } else {
        // Default to empty string
        metadata->date_obs[0] = '\0';
    }
    
    // Extract quality information
    gboolean quality_found = FALSE;
    
    if (json_object_has_member(root_object, "quality")) {
        metadata->quality_accepted = json_object_get_boolean_member(root_object, "quality");
        g_strlcpy(metadata->quality_reason, 
                  metadata->quality_accepted ? "Accepted by Stellina" : "Rejected by Stellina",
                  sizeof(metadata->quality_reason));
        quality_found = TRUE;
    }
    
    if (json_object_has_member(root_object, "used_for_stacking")) {
        metadata->used_for_stacking = json_object_get_boolean_member(root_object, "used_for_stacking");
        if (!quality_found) {
            metadata->quality_accepted = metadata->used_for_stacking;
            g_strlcpy(metadata->quality_reason,
                      metadata->used_for_stacking ? "Used for stacking" : "Not used for stacking",
                      sizeof(metadata->quality_reason));
            quality_found = TRUE;
        }
    }
    
    // Extract additional quality metrics
    if (json_object_has_member(root_object, "fwhm")) {
        metadata->fwhm = json_object_get_double_member(root_object, "fwhm");
    }
    
    if (json_object_has_member(root_object, "star_count")) {
        metadata->star_count = json_object_get_int_member(root_object, "star_count");
    } else if (json_object_has_member(root_object, "stars")) {
        metadata->star_count = json_object_get_int_member(root_object, "stars");
    }
    
    // Apply additional quality heuristics if no explicit quality flag
    if (!quality_found) {
        metadata->quality_accepted = TRUE;
        g_strlcpy(metadata->quality_reason, "No quality flag - assuming good", 
                  sizeof(metadata->quality_reason));
        
        // Check for obvious problems
        if (metadata->fwhm > 8.0) {
            metadata->quality_accepted = FALSE;
            g_strlcpy(metadata->quality_reason, "Poor FWHM > 8.0", sizeof(metadata->quality_reason));
        } else if (metadata->star_count < 5) {
            metadata->quality_accepted = FALSE;
            g_strlcpy(metadata->quality_reason, "Too few stars < 5", sizeof(metadata->quality_reason));
        } else if (metadata->altitude < 30.0) {
            metadata->quality_accepted = FALSE;
            g_strlcpy(metadata->quality_reason, "Low altitude < 30°", sizeof(metadata->quality_reason));
        }
    }
    
    g_object_unref(parser);
    
    siril_debug_print("Parsed JSON: Alt=%.2f°, Az=%.2f°, Quality=%s (%s)\n",
                     metadata->altitude, metadata->azimuth,
                     metadata->quality_accepted ? "ACCEPTED" : "REJECTED",
                     metadata->quality_reason);
    
    return metadata;
}

/**
 * Add coordinates and metadata to FITS header - fixed for CFITSIO
 */
int stellina_add_coordinates_to_header_enhanced(fits *fit, const struct stellina_metadata *metadata, 
                                              double ra, double dec) {
    if (!fit || !fit->fptr || !metadata) {
        return -1;
    }
    
    int status = 0;
    
    // Add calculated coordinates - cast away const for CFITSIO
    double ra_copy = ra;
    double dec_copy = dec;
    double alt_copy = metadata->altitude;
    double az_copy = metadata->azimuth;
    
    fits_update_key(fit->fptr, TDOUBLE, "RA_CALC", &ra_copy, 
                   "Calculated RA from Alt/Az (degrees)", &status);
    fits_update_key(fit->fptr, TDOUBLE, "DEC_CALC", &dec_copy, 
                   "Calculated Dec from Alt/Az (degrees)", &status);
    
    // Add original Alt/Az
    fits_update_key(fit->fptr, TDOUBLE, "ALT_ORIG", &alt_copy, 
                   "Original altitude from Stellina (degrees)", &status);
    fits_update_key(fit->fptr, TDOUBLE, "AZ_ORIG", &az_copy, 
                   "Original azimuth from Stellina (degrees)", &status);
    
    // Add quality information
    const char *quality_str = metadata->quality_accepted ? "ACCEPTED" : "REJECTED";
    fits_update_key(fit->fptr, TSTRING, "QUALITY", (void*)quality_str, 
                   "Image quality assessment", &status);
    fits_update_key(fit->fptr, TSTRING, "QUAL_RSN", (void*)metadata->quality_reason, 
                   "Quality assessment reason", &status);
    
    // Add FWHM if available
    if (metadata->fwhm > 0.0) {
        double fwhm_copy = metadata->fwhm;
        fits_update_key(fit->fptr, TDOUBLE, "FWHM", &fwhm_copy, 
                       "Full Width Half Maximum (arcsec)", &status);
    }
    
    // Add star count if available
    if (metadata->star_count > 0) {
        int star_count_copy = metadata->star_count;
        fits_update_key(fit->fptr, TINT, "STARCOUNT", &star_count_copy, 
                       "Number of detected stars", &status);
    }
    
    // Add processing marker
    const char *processor = "STELLINA_SIRIL";
    fits_update_key(fit->fptr, TSTRING, "PROCSSED", (void*)processor, 
                   "Processed by Stellina Siril extension", &status);
    
    if (status != 0) {
        char err_text[FLEN_ERRMSG];
        fits_get_errstatus(status, err_text);
        siril_log_color_message("FITS header update error: %s\n", "red", err_text);
        return -1;
    }
    
    siril_debug_print("Added coordinates and metadata to FITS header\n");
    return 0;
}

/**
 * Complete directory processing implementation
 */
int stellina_process_directory(const char *input_dir, const char *output_dir, 
                              struct stellina_config *config) {
    if (!input_dir || !output_dir || !config) {
        return -1;
    }
    
    siril_log_message("Starting Stellina directory processing\n");
    siril_log_message("  Input: %s\n", input_dir);
    siril_log_message("  Output: %s\n", output_dir);
    siril_log_message("  Quality filtering: %s\n", 
                     config->quality_filter_enabled ? "enabled" : "disabled");
    
    // Find matching files
    auto matches = find_matching_files(input_dir);
    if (matches.empty()) {
        siril_log_color_message("No matching FITS/JSON pairs found in %s\n", "orange", input_dir);
        return 0;
    }
    
    // Create output directory
    try {
        std::filesystem::create_directories(output_dir);
    } catch (const std::exception& e) {
        siril_log_color_message("Failed to create output directory %s: %s\n", "red", 
                               output_dir, e.what());
        return -1;
    }
    
    // Set default observer location if needed
    stellina_set_default_observer_location(config);
    
    // Process each pair
    int processed = 0, skipped = 0, errors = 0;
    
    for (size_t i = 0; i < matches.size(); ++i) {
        const std::string& json_path = matches[i].first;
        const std::string& fits_path = matches[i].second;
        
        // Update progress
        if (g_progress_callback) {
            std::filesystem::path fp(fits_path);
            std::string message = "Processing " + fp.filename().string();
            g_progress_callback(i + 1, matches.size(), message.c_str());
        }
        
        // Generate output path
        std::filesystem::path input_fits(fits_path);
        std::filesystem::path output_fits = std::filesystem::path(output_dir) / 
                                          (std::string(config->output_prefix) + input_fits.filename().string());
        
        try {
            int result = stellina_process_single_image(fits_path.c_str(), json_path.c_str(), 
                                                     output_fits.string().c_str(), config);
            if (result == 0) {
                processed++;
            } else if (result == 1) {
                skipped++;
            } else {
                errors++;
            }
        } catch (const std::exception& e) {
            siril_log_color_message("Error processing %s: %s\n", "red", 
                                   fits_path.c_str(), e.what());
            errors++;
        }
    }
    
    // Print summary
    siril_log_color_message("Stellina processing complete:\n", "green");
    siril_log_message("  Total pairs found: %zu\n", matches.size());
    siril_log_message("  Successfully processed: %d\n", processed);
    siril_log_message("  Skipped (quality): %d\n", skipped);
    siril_log_message("  Errors: %d\n", errors);
    
    return processed > 0 ? 0 : -1;
}

/**
 * Complete single image processing implementation
 */
int stellina_process_single_image(const char *fits_path, const char *json_path, 
                                 const char *output_path, struct stellina_config *config) {
    if (!fits_path || !json_path || !output_path || !config) {
        return -1;
    }
    
    siril_debug_print("Processing single image: %s\n", fits_path);
    
    // Parse JSON metadata
    struct stellina_metadata *metadata = stellina_parse_json_enhanced(json_path);
    if (!metadata) {
        siril_log_color_message("Failed to parse JSON metadata: %s\n", "red", json_path);
        return -1;
    }
    
    // Check quality if filtering enabled
    if (config->quality_filter_enabled && !stellina_check_quality(metadata)) {
        siril_log_message("Skipping %s: %s\n", fits_path, metadata->quality_reason);
        stellina_metadata_free(metadata);
        return 1; // Skipped
    }
    
    // If date is not in JSON, try to get it from FITS header
    char date_obs[256];
    if (metadata->date_obs[0] == '\0') {
        siril_log_message("No DATE-OBS in JSON, reading from FITS header\n");
        if (extract_date_from_fits(fits_path, date_obs, sizeof(date_obs))) {
            siril_log_message("Using DATE-OBS from FITS: %s\n", date_obs);
        } else {
            siril_log_color_message("Failed to extract DATE-OBS from FITS\n", "red");
            stellina_metadata_free(metadata);
            return -1;
        }
    } else {
        // Use date from metadata
        g_strlcpy(date_obs, metadata->date_obs, sizeof(date_obs));
    }
    
    // Convert Alt/Az to RA/Dec
    double ra, dec;
    int coord_result = stellina_convert_altaz_to_radec(
        metadata->altitude, metadata->azimuth, date_obs,
        config->observer_latitude, config->observer_longitude, config->observer_altitude,
        &ra, &dec
    );
    
    if (coord_result != 0) {
        siril_log_color_message("Failed to convert coordinates for %s\n", "red", fits_path);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    // Load FITS file using our direct CFITSIO wrapper
    fits fit = { 0 };
    if (load_fits_file(fits_path, &fit) != 0) {
        siril_log_color_message("Failed to load FITS file: %s\n", "red", fits_path);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    // Add coordinates to header
    if (config->add_coordinate_keywords) {
        if (stellina_add_coordinates_to_header_enhanced(&fit, metadata, ra, dec) != 0) {
            siril_log_color_message("Failed to add coordinates to header: %s\n", "red", fits_path);
            close_fits_file(&fit);
            stellina_metadata_free(metadata);
            return -1;
        }
    }
    
    // Save to output path
    if (save_fits_file(output_path, &fit) != 0) {
        siril_log_color_message("Failed to save FITS file: %s\n", "red", output_path);
        close_fits_file(&fit);
        stellina_metadata_free(metadata);
        return -1;
    }
    
    siril_log_message("Successfully processed %s -> %s\n", fits_path, output_path);
    siril_debug_print("  Coordinates: Alt/Az %.2f°/%.2f° -> RA/Dec %.4f°/%.4f°\n",
                     metadata->altitude, metadata->azimuth, ra, dec);
    
    // Cleanup
    close_fits_file(&fit);
    stellina_metadata_free(metadata);
    
    return 0; // Success
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