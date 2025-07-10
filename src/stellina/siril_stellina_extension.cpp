/*
 * Siril Stellina Extension - C++ Plugin for Robust Stellina Processing
 * 
 * This extension adds native Stellina support to Siril with:
 * - Robust coordinate conversion from Stellina Alt/Az to RA/Dec
 * - Quality filtering based on Stellina metadata
 * - Graceful plate solving with fallbacks
 * - Proper FITS header management
 * - Batch processing with detailed progress reporting
 */

#ifndef SIRIL_STELLINA_EXTENSION_H
#define SIRIL_STELLINA_EXTENSION_H

#include "core/siril.h"
#include "core/proto.h"
#include "io/single_image.h"
#include "registration/registration.h"
#include "gui/callbacks.h"
#include <string>
#include <vector>
#include <memory>

// Forward declarations
struct StellinaMetadata;
struct ProcessingResults;

class StellinaProcessor {
public:
    StellinaProcessor();
    ~StellinaProcessor();
    
    // Main processing functions
    bool processDirectory(const std::string& inputDir, const std::string& outputDir);
    bool processSingleImage(const std::string& fitsPath, const std::string& jsonPath);
    
    // Configuration
    void setQualityFiltering(bool enabled) { qualityFilterEnabled_ = enabled; }
    void setPlatesolveFallback(bool enabled) { platesolveFallback_ = enabled; }
    void setObserverLocation(double lat, double lon, double alt);
    
    // Results and progress
    ProcessingResults getResults() const { return results_; }
    void setProgressCallback(std::function<void(int, int, const std::string&)> callback);

private:
    // Core processing methods
    StellinaMetadata parseJSON(const std::string& jsonPath);
    bool convertAltAzToRADec(double alt, double az, const std::string& dateObs, 
                           double& ra, double& dec);
    bool checkStellinaQuality(const StellinaMetadata& metadata);
    bool loadAndPrepareImage(const std::string& fitsPath);
    bool performPlatesolving(double ra, double dec);
    bool addCoordinatesToHeader(double ra, double dec, double alt, double az, 
                              const std::string& dateObs);
    bool saveProcessedImage(const std::string& outputPath);
    
    // Siril integration
    fits* currentImage_;
    bool qualityFilterEnabled_;
    bool platesolveFallback_;
    ProcessingResults results_;
    std::function<void(int, int, const std::string&)> progressCallback_;
    
    // Observer location for coordinate conversion
    double observerLat_, observerLon_, observerAlt_;
};

// Stellina metadata structure
struct StellinaMetadata {
    double altitude;
    double azimuth;
    std::string dateObs;
    bool qualityAccepted;
    std::string qualityReason;
    double fwhm;
    int starCount;
    bool usedForStacking;
};

// Processing results
struct ProcessingResults {
    int totalImages = 0;
    int processedImages = 0;
    int qualityRejected = 0;
    int platesolveFailed = 0;
    int succeeded = 0;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
};

// Siril command interface functions
extern "C" {
    // Register the extension with Siril
    void siril_stellina_init();
    
    // New Siril commands added by this extension
    int process_stellina_directory(struct command_info* cmd);
    int process_stellina_image(struct command_info* cmd);
    int set_stellina_observer(struct command_info* cmd);
    int stellina_quality_filter(struct command_info* cmd);
}

// Implementation details

StellinaProcessor::StellinaProcessor() 
    : currentImage_(nullptr)
    , qualityFilterEnabled_(true)
    , platesolveFallback_(true)
    , observerLat_(0.0)
    , observerLon_(0.0) 
    , observerAlt_(0.0) {
}

bool StellinaProcessor::processDirectory(const std::string& inputDir, const std::string& outputDir) {
    // Find all FITS/JSON pairs
    std::vector<std::pair<std::string, std::string>> imagePairs;
    
    // Scan directory for matching files (implementation would use filesystem APIs)
    // This would replace your Python find_matching_files function
    
    results_ = ProcessingResults();
    results_.totalImages = imagePairs.size();
    
    for (size_t i = 0; i < imagePairs.size(); ++i) {
        const auto& [fitsPath, jsonPath] = imagePairs[i];
        
        if (progressCallback_) {
            progressCallback_(i + 1, imagePairs.size(), 
                            "Processing " + std::filesystem::path(fitsPath).filename().string());
        }
        
        try {
            if (processSingleImage(fitsPath, jsonPath)) {
                results_.succeeded++;
            }
            results_.processedImages++;
        } catch (const std::exception& e) {
            results_.errors.push_back("Error processing " + fitsPath + ": " + e.what());
            siril_log_color_message("Error processing %s: %s\n", "red", fitsPath.c_str(), e.what());
        }
    }
    
    // Print summary
    siril_log_color_message("Stellina processing complete: %d/%d succeeded\n", 
                           "green", results_.succeeded, results_.totalImages);
    
    return results_.succeeded > 0;
}

bool StellinaProcessor::processSingleImage(const std::string& fitsPath, const std::string& jsonPath) {
    // Parse Stellina metadata
    StellinaMetadata metadata = parseJSON(jsonPath);
    
    // Check quality if filtering enabled
    if (qualityFilterEnabled_ && !checkStellinaQuality(metadata)) {
        results_.qualityRejected++;
        siril_log_message("Skipping %s: %s\n", fitsPath.c_str(), metadata.qualityReason.c_str());
        return false;
    }
    
    // Load the FITS image
    if (!loadAndPrepareImage(fitsPath)) {
        results_.errors.push_back("Failed to load image: " + fitsPath);
        return false;
    }
    
    // Convert Alt/Az to RA/Dec
    double ra, dec;
    if (!convertAltAzToRADec(metadata.altitude, metadata.azimuth, metadata.dateObs, ra, dec)) {
        results_.warnings.push_back("Failed coordinate conversion for: " + fitsPath);
        // Continue processing with original coordinates
    }
    
    // Attempt plate solving
    bool platesolved = performPlatesolving(ra, dec);
    if (!platesolved) {
        results_.platesolveFailed++;
        if (platesolveFallback_) {
            siril_log_message("Plate solving failed for %s, using calculated coordinates\n", 
                            fitsPath.c_str());
        } else {
            results_.errors.push_back("Plate solving failed: " + fitsPath);
            return false;
        }
    }
    
    // Add coordinate information to FITS header
    addCoordinatesToHeader(ra, dec, metadata.altitude, metadata.azimuth, metadata.dateObs);
    
    // Save processed image
    std::string outputPath = generateOutputPath(fitsPath);
    return saveProcessedImage(outputPath);
}

bool StellinaProcessor::performPlatesolving(double ra, double dec) {
    // Use Siril's internal plate solving functions
    // This gives us access to the same algorithms but with better error handling
    
    struct platesolve_params params;
    params.ra = ra;
    params.dec = dec;
    params.force_platesolve = TRUE;
    params.use_hint = TRUE;
    
    // Call Siril's internal plate solving
    int result = platesolve_from_command(&params);
    
    if (result == 0) {
        siril_log_color_message("Plate solving succeeded\n", "green");
        return true;
    } else {
        siril_log_color_message("Plate solving failed, continuing with calculated coordinates\n", "orange");
        return false;
    }
}

bool StellinaProcessor::addCoordinatesToHeader(double ra, double dec, double alt, double az, 
                                             const std::string& dateObs) {
    if (!currentImage_ || !currentImage_->header) {
        return false;
    }
    
    // Add calculated coordinates to FITS header
    // This uses Siril's internal FITS header manipulation functions
    fits_update_key(currentImage_->fptr, TDOUBLE, "RA_CALC", &ra, "Calculated RA (degrees)", nullptr);
    fits_update_key(currentImage_->fptr, TDOUBLE, "DEC_CALC", &dec, "Calculated Dec (degrees)", nullptr);
    fits_update_key(currentImage_->fptr, TDOUBLE, "ALT_ORIG", &alt, "Original altitude (degrees)", nullptr);
    fits_update_key(currentImage_->fptr, TDOUBLE, "AZ_ORIG", &az, "Original azimuth (degrees)", nullptr);
    fits_update_key(currentImage_->fptr, TSTRING, "STELLINA", "PROCESSED", "Processed by Stellina extension", nullptr);
    
    return true;
}

// Siril command implementations
extern "C" {
    void siril_stellina_init() {
        // Register new commands with Siril
        add_command("process_stellina", process_stellina_directory, 
                   "process_stellina input_dir output_dir [quality_filter]");
        add_command("stellina_image", process_stellina_image,
                   "stellina_image fits_file json_file");
        add_command("set_stellina_observer", set_stellina_observer,
                   "set_stellina_observer latitude longitude altitude");
        
        siril_log_color_message("Stellina extension loaded successfully\n", "green");
    }
    
    int process_stellina_directory(struct command_info* cmd) {
        if (cmd->args_count < 2) {
            siril_log_message("Usage: process_stellina input_dir output_dir [quality_filter]\n");
            return CMD_ARG_ERROR;
        }
        
        std::string inputDir = cmd->args[0];
        std::string outputDir = cmd->args[1];
        bool qualityFilter = (cmd->args_count > 2) ? (strcmp(cmd->args[2], "true") == 0) : true;
        
        StellinaProcessor processor;
        processor.setQualityFiltering(qualityFilter);
        
        // Set progress callback to update Siril's progress bar
        processor.setProgressCallback([](int current, int total, const std::string& message) {
            set_progress_bar_data(message.c_str(), (double)current / total);
        });
        
        bool success = processor.processDirectory(inputDir, outputDir);
        
        return success ? CMD_OK : CMD_GENERIC_ERROR;
    }
}

#endif // SIRIL_STELLINA_EXTENSION_H