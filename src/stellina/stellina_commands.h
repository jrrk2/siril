/*
 * stellina_commands.h
 * C-compatible header for Stellina command declarations
 * Updated with correct function signatures for Siril integration
 */

#ifndef STELLINA_COMMANDS_H
#define STELLINA_COMMANDS_H

#ifdef __cplusplus
extern "C" {
#endif

// Forward declare command_info if not already declared
struct command_info;

// Siril command interface functions with correct signatures
// These functions take an int parameter representing the number of arguments
int cmd_process_stellina(int nb);
int cmd_stellina_config(int nb);
int cmd_stellina_stats(int nb);
int cmd_stellina_test_coords(int nb);
int cmd_stellina_platesolve(int nb);
int cmd_stellina_script(int nb);

// Initialization function
void stellina_commands_init(void);

#include "core/siril_world_cs.h" // For coordinate system functions
#include "core/command_line_processor.h"
#include "core/initfile.h"
#include "livestacking/livestacking.h"
#include "gui/registration_preview.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/script_menu.h"
#include "gui/callbacks.h"
#include "gui/utils.h"

#ifdef __cplusplus
}
#endif

// Command descriptions for Siril's help system
#define STR_PROCESS_STELLINA N_("Process Stellina FITS/JSON pairs with coordinate conversion")
#define STR_STELLINA_CONFIG N_("Configure Stellina processing parameters")
#define STR_STELLINA_STATS N_("Show Stellina statistics and analyze files")
#define STR_STELLINA_TEST_COORDS N_("Test coordinate conversion (debugging)")
#define STR_STELLINA_PLATESOLVE N_("Stellina specific platesolve")
#define STR_STELLINA_SCRIPT N_("Stellina specific script")

// Command help strings
#define HELP_PROCESS_STELLINA N_(\
"stellina_process input_dir output_dir [quality_filter]\n\n" \
"Process a directory of Stellina FITS/JSON pairs:\n" \
"- Converts Alt/Az coordinates to RA/Dec using libnova\n" \
"- Applies quality filtering based on Stellina metadata\n" \
"- Adds coordinate information to FITS headers\n" \
"- Organizes processed files in output directory\n\n" \
"Arguments:\n" \
"  input_dir: Directory containing FITS and JSON files\n" \
"  output_dir: Directory for processed files\n" \
"  quality_filter: Enable quality filtering (true/false, default: true)\n\n" \
"Example: stellina_process ~/stellina_raw ~/stellina_processed")

#define HELP_STELLINA_CONFIG N_(\
"stellina_config parameter value\n\n" \
"Configure Stellina processing parameters:\n\n" \
"Parameters:\n" \
"  quality_filter <true|false>     - Enable/disable quality filtering\n" \
"  platesolve_fallback <true|false> - Enable/disable platesolve fallback\n" \
"  add_keywords <true|false>       - Add coordinate keywords to FITS\n" \
"  observer_location <lat,lon,alt>  - Set observer location\n" \
"  output_prefix <prefix>           - Set output filename prefix\n" \
"  show                            - Show current configuration\n\n" \
"Examples:\n" \
"  stellina_config observer_location 48.8566,2.3522,35\n" \
"  stellina_config quality_filter false\n" \
"  stellina_config show")

#define HELP_STELLINA_STATS N_(\
"stellina_stats [fits_file json_file]\n\n" \
"Show Stellina extension information or analyze specific files:\n\n" \
"With no arguments: Show general extension information\n" \
"With file arguments: Analyze specific FITS/JSON pair\n\n" \
"Example: stellina_stats image.fits image.json")

#define HELP_STELLINA_TEST_COORDS N_(\
"stellina_test_coords alt az date_obs [lat lon alt_m]\n\n" \
"Test coordinate conversion from Alt/Az to RA/Dec:\n\n" \
"Arguments:\n" \
"  alt: Altitude in degrees\n" \
"  az: Azimuth in degrees\n" \
"  date_obs: Observation time (ISO format: 2024-01-09T22:13:29)\n" \
"  lat, lon, alt_m: Observer location (optional)\n\n" \
"Example: stellina_test_coords 45.5 180.0 2024-01-09T22:13:29")

#define HELP_STELLINA_PLATESOLVE N_(\
"stellina_platesolve [fits_file]\n\n" \
"Example: stellina_platesolve image.fits")

#define HELP_STELLINA_SCRIPT N_(\
"stellina_script script_file]\n\n" \
"Example: stellina_platesolve myscript.ssf")

#endif // STELLINA_COMMANDS_H
