/*
 * stellina_commands.cpp
 * Siril command implementations for Stellina processing
 */

#include "stellina_processor.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "core/processing.h"
#include "io/path_parse.h"

// Global configuration with defaults
static struct stellina_config g_stellina_config = {
    .quality_filter_enabled = TRUE,
    .platesolve_fallback_enabled = TRUE,
    .add_coordinate_keywords = TRUE,
    .observer_latitude = 0.0,
    .observer_longitude = 0.0,
    .observer_altitude = 0.0,
    .output_prefix = "stellina_"
};

// Progress callback for GUI integration
static void stellina_gui_progress(int current, int total, const char *message) {
    double progress = (double)current / (double)total;
    set_progress_bar_data(message, progress);
    
    // Update GUI
    if (com.headless == FALSE) {
        gtk_main_iteration_do(FALSE);
    }
}

/**
 * process_stellina command implementation
 * Usage: process_stellina input_directory output_directory [options]
 * Options: -no-quality-filter, -no-platesolve-fallback, -observer lat,lon,alt
 */
int cmd_process_stellina(struct command_info *cmd) {
    if (cmd->argc < 3) {
        siril_log_color_message(_("Usage: process_stellina input_directory output_directory [options]\n"), "red");
        siril_log_message(_("Options:\n"));
        siril_log_message(_("  -no-quality-filter    : Process all images regardless of Stellina quality assessment\n"));
        siril_log_message(_("  -no-platesolve        : Skip plate solving entirely\n"));
        siril_log_message(_("  -observer lat,lon,alt : Set observer location (degrees, degrees, meters)\n"));
        return CMD_ARG_ERROR;
    }
    
    const char *input_dir = cmd->argv[1];
    const char *output_dir = cmd->argv[2];
    
    // Parse options
    struct stellina_config config = g_stellina_config; // Copy global config
    
    for (int i = 3; i < cmd->argc; i++) {
        if (strcmp(cmd->argv[i], "-no-quality-filter") == 0) {
            config.quality_filter_enabled = FALSE;
        } else if (strcmp(cmd->argv[i], "-no-platesolve") == 0) {
            config.platesolve_fallback_enabled = FALSE;
        } else if (strcmp(cmd->argv[i], "-observer") == 0) {
            if (i + 1 < cmd->argc) {
                char *observer_str = cmd->argv[++i];
                if (sscanf(observer_str, "%lf,%lf,%lf", 
                          &config.observer_latitude, 
                          &config.observer_longitude, 
                          &config.observer_altitude) != 3) {
                    siril_log_color_message(_("Invalid observer location format. Use: lat,lon,alt\n"), "red");
                    return CMD_ARG_ERROR;
                }
            } else {
                siril_log_color_message(_("-observer option requires lat,lon,alt argument\n"), "red");
                return CMD_ARG_ERROR;
            }
        }
    }
    
    // Validate directories
    if (!g_file_test(input_dir, G_FILE_TEST_IS_DIR)) {
        siril_log_color_message(_("Input directory does not exist: %s\n"), "red", input_dir);
        return CMD_ARG_ERROR;
    }
    
    // Create output directory if it doesn't exist
    if (!g_file_test(output_dir, G_FILE_TEST_IS_DIR)) {
        if (g_mkdir_with_parents(output_dir, 0755) != 0) {
            siril_log_color_message(_("Failed to create output directory: %s\n"), "red", output_dir);
            return CMD_GENERIC_ERROR;
        }
    }
    
    // Set up progress reporting
    stellina_set_progress_callback(stellina_gui_progress);
    
    siril_log_color_message(_("Starting Stellina processing...\n"), "green");
    siril_log_message(_("Input directory: %s\n"), input_dir);
    siril_log_message(_("Output directory: %s\n"), output_dir);
    siril_log_message(_("Quality filtering: %s\n"), config.quality_filter_enabled ? "enabled" : "disabled");
    siril_log_message(_("Plate solve fallback: %s\n"), config.platesolve_fallback_enabled ? "enabled" : "disabled");
    
    if (config.observer_latitude != 0.0 || config.observer_longitude != 0.0) {
        siril_log_message(_("Observer location: %.6f°, %.6f°, %.1fm\n"), 
                         config.observer_latitude, config.observer_longitude, config.observer_altitude);
    }
    
    // Process the directory
    int result = stellina_process_directory(input_dir, output_dir, &config);
    
    if (result == 0) {
        siril_log_color_message(_("Stellina processing completed successfully\n"), "green");
    } else {
        siril_log_color_message(_("Stellina processing failed with errors\n"), "red");
    }
    
    return (result == 0) ? CMD_OK : CMD_GENERIC_ERROR;
}

/**
 * stellina_config command implementation
 * Usage: stellina_config option value
 */
int cmd_stellina_config(struct command_info *cmd) {
    if (cmd->argc < 3) {
        siril_log_color_message(_("Usage: stellina_config option value\n"), "red");
        siril_log_message(_("Options:\n"));
        siril_log_message(_("  quality_filter true|false\n"));
        siril_log_message(_("  platesolve_fallback true|false\n"));
        siril_log_message(_("  add_coordinates true|false\n"));
        siril_log_message(_("  observer_location lat,lon,alt\n"));
        siril_log_message(_("  output_prefix string\n"));
        return CMD_ARG_ERROR;
    }
    
    const char *option = cmd->argv[1];
    const char *value = cmd->argv[2];
    
    if (strcmp(option, "quality_filter") == 0) {
        g_stellina_config.quality_filter_enabled = (strcmp(value, "true") == 0);
        siril_log_message(_("Quality filtering: %s\n"), 
                         g_stellina_config.quality_filter_enabled ? "enabled" : "disabled");
    } else if (strcmp(option, "platesolve_fallback") == 0) {
        g_stellina_config.platesolve_fallback_enabled = (strcmp(value, "true") == 0);
        siril_log_message(_("Plate solve fallback: %s\n"), 
                         g_stellina_config.platesolve_fallback_enabled ? "enabled" : "disabled");
    } else if (strcmp(option, "add_coordinates") == 0) {
        g_stellina_config.add_coordinate_keywords = (strcmp(value, "true") == 0);
        siril_log_message(_("Add coordinate keywords: %s\n"), 
                         g_stellina_config.add_coordinate_keywords ? "enabled" : "disabled");
    } else if (strcmp(option, "observer_location") == 0) {
        if (sscanf(value, "%lf,%lf,%lf", 
                  &g_stellina_config.observer_latitude, 
                  &g_stellina_config.observer_longitude, 
                  &g_stellina_config.observer_altitude) == 3) {
            siril_log_message(_("Observer location set to: %.6f°, %.6f°, %.1fm\n"), 
                             g_stellina_config.observer_latitude, 
                             g_stellina_config.observer_longitude, 
                             g_stellina_config.observer_altitude);
        } else {
            siril_log_color_message(_("Invalid observer location format. Use: lat,lon,alt\n"), "red");
            return CMD_ARG_ERROR;
        }
    } else if (strcmp(option, "output_prefix") == 0) {
        g_strlcpy(g_stellina_config.output_prefix, value, sizeof(g_stellina_config.output_prefix));
        siril_log_message(_("Output prefix set to: %s\n"), g_stellina_config.output_prefix);
    } else {
        siril_log_color_message(_("Unknown option: %s\n"), "red", option);
        return CMD_ARG_ERROR;
    }
    
    return CMD_OK;
}

/**
 * stellina_stats command implementation
 * Shows current processing statistics
 */
int cmd_stellina_stats(struct command_info *cmd) {
    siril_log_message(_("Current Stellina Configuration:\n"));
    siril_log_message(_("  Quality filtering: %s\n"), 
                     g_stellina_config.quality_filter_enabled ? "enabled" : "disabled");
    siril_log_message(_("  Plate solve fallback: %s\n"), 
                     g_stellina_config.platesolve_fallback_enabled ? "enabled" : "disabled");
    siril_log_message(_("  Add coordinate keywords: %s\n"), 
                     g_stellina_config.add_coordinate_keywords ? "enabled" : "disabled");
    siril_log_message(_("  Output prefix: %s\n"), g_stellina_config.output_prefix);
    
    if (g_stellina_config.observer_latitude != 0.0 || g_stellina_config.observer_longitude != 0.0) {
        siril_log_message(_("  Observer location: %.6f°, %.6f°, %.1fm\n"), 
                         g_stellina_config.observer_latitude, 
                         g_stellina_config.observer_longitude, 
                         g_stellina_config.observer_altitude);
    } else {
        siril_log_message(_("  Observer location: not set\n"));
    }
    
    return CMD_OK;
}