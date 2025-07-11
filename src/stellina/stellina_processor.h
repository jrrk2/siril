/*
 * stellina_processor.h
 * Stellina processing extension for Siril - Fixed header
 */

#ifndef STELLINA_PROCESSOR_H
#define STELLINA_PROCESSOR_H

// C includes
#include <glib.h>
#include <json-glib/json-glib.h>

// C++ includes (only when compiling as C++)
#ifdef __cplusplus
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <libnova/libnova.h>

// Siril includes (now that we have GTK)
#include "core/siril.h"
#include "core/proto.h"
#include "io/single_image.h"
#include "io/fits_keywords.h"
#include "registration/registration.h"
#endif

// Stellina metadata structure (C-compatible)
struct stellina_metadata {
    double altitude;
    double azimuth;
    char date_obs[256];
    gboolean quality_accepted;
    char quality_reason[256];
    double fwhm;
    int star_count;
    gboolean used_for_stacking;
    char original_filename[256];
};

// Processing statistics
struct stellina_stats {
    int total_images;
    int processed_images;
    int quality_rejected;
    int platesolve_failed;
    int platesolve_succeeded;
    int coordinate_conversion_failed;
    GList *error_messages;
    GList *warning_messages;
};

// Configuration
struct stellina_config {
    gboolean quality_filter_enabled;
    gboolean platesolve_fallback_enabled;
    gboolean add_coordinate_keywords;
    double observer_latitude;
    double observer_longitude;
    double observer_altitude;
    char output_prefix[64];
};

extern "C" {
#include "algos/astrometry_solver.h" // For plate solving
#include "core/siril_world_cs.h" // For coordinate system functions
};

gchar *stellina_platesolve(fits *preffit, SirilWorldCS *target_coords, double forced_focal, double forced_pixsize);

#ifdef __cplusplus
extern "C" {
#endif

// Main processing functions
int stellina_process_directory(const char *input_dir, const char *output_dir, struct stellina_config *config);
int stellina_process_single_image(const char *fits_path, const char *json_path, const char *output_path, struct stellina_config *config);

// Utility functions
struct stellina_metadata *stellina_parse_json(const char *json_path);
struct stellina_metadata *stellina_parse_json_enhanced(const char *json_path);
int stellina_convert_altaz_to_radec(double alt, double az, const char *date_obs, 
                                   double observer_lat, double observer_lon, double observer_alt,
                                   double *ra, double *dec);
gboolean stellina_check_quality(const struct stellina_metadata *metadata);

#ifdef __cplusplus
// Only declare this in C++ mode since 'fits' type needs C++ includes
int stellina_add_coordinates_to_header(fits *fit, const struct stellina_metadata *metadata, 
                                      double ra, double dec);
int stellina_add_coordinates_to_header_enhanced(fits *fit, const struct stellina_metadata *metadata, 
                                              double ra, double dec);
#endif

// Memory management
void stellina_metadata_free(struct stellina_metadata *metadata);
void stellina_stats_free(struct stellina_stats *stats);

// Progress reporting
typedef void (*stellina_progress_callback)(int current, int total, const char *message);
void stellina_set_progress_callback(stellina_progress_callback callback);

// Helper functions
void stellina_set_default_observer_location(struct stellina_config *config);
gboolean stellina_validate_observer_location(double lat, double lon, double alt);

// Coordinate utility functions from stellina_coordinates.cpp
int stellina_convert_altaz_to_radec_corrected(double alt, double az, const char *date_obs,
                                             double observer_lat, double observer_lon, double observer_alt,
                                             double temperature, double pressure,
                                             double *ra, double *dec);
double stellina_get_local_sidereal_time(const char *date_obs, double observer_lon);
int stellina_get_solar_object_coords(const char *date_obs, double observer_lat, double observer_lon,
                                   const char *object_name, double *ra, double *dec);

#ifdef __cplusplus
}
#endif

#endif // STELLINA_PROCESSOR_H
