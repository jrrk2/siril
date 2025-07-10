/*
 * stellina_processor.cpp
 * Minimal Stellina processing implementation
 */

#include "stellina_processor.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstring>

static stellina_progress_callback g_progress_callback = nullptr;

/**
 * Set progress callback
 */
void stellina_set_progress_callback(stellina_progress_callback callback) {
    g_progress_callback = callback;
}

/**
 * Placeholder for directory processing
 */
int stellina_process_directory(const char *input_dir, const char *output_dir, struct stellina_config *config) {
    if (!input_dir || !output_dir || !config) {
        return -1;
    }
    
    std::cout << "Processing Stellina directory: " << input_dir << " -> " << output_dir << std::endl;
    std::cout << "Quality filtering: " << (config->quality_filter_enabled ? "enabled" : "disabled") << std::endl;
    
    // For now, just return success
    return 0;
}

/**
 * Placeholder for single image processing  
 */
int stellina_process_single_image(const char *fits_path, const char *json_path, const char *output_path, struct stellina_config *config) {
    if (!fits_path || !json_path || !output_path || !config) {
        return -1;
    }
    
    std::cout << "Processing single image: " << fits_path << std::endl;
    
    // For now, just return success
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