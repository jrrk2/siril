/*
 * stellina_commands.cpp
 * Clean final Siril command implementations for Stellina processing
 */

#include "core/siril.h"

extern "C" {
#include "io/single_image.h"
};

#include "stellina_processor.h"
#include "stellina_commands.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <fitsio.h>

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

int cmd_stellina_platesolve(int nb) {
    if (nb < 3) {
        siril_log_color_message("Usage: stellina_platesolve ra dec\n", "red");
        siril_log_message("  ra: Guide RA in degrees\n");
        siril_log_message("  dec: Guide DEC in degrees\n");
        siril_log_message("Example: stellina_platesolve 10.776 41.238\n");
        return CMD_ARG_ERROR;
    }

    fits *preffit = &gfit;
    SirilWorldCS *target_coords = siril_world_cs_new_from_objct_ra_dec(word[1], word[2]);
    double forced_focal = 400.0;
    double forced_pixsize = 2.4;
    gchar *err_msg = stellina_platesolve(preffit, target_coords, forced_focal, forced_pixsize);
    puts(err_msg);
    return CMD_OK;
}

static void update_log_icon(gboolean is_running) {
	GtkImage *image = GTK_IMAGE(lookup_widget("image_log"));
	if (is_running)
		gtk_image_set_from_icon_name(image, "gtk-yes", GTK_ICON_SIZE_LARGE_TOOLBAR);
	else
		gtk_image_set_from_icon_name(image, "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);
}

struct log_status_bar_idle_data {
	gchar *myline;
	int line;
};

static gboolean log_status_bar_idle_callback(gpointer p) {
	struct log_status_bar_idle_data *data = (struct log_status_bar_idle_data *) p;

	GtkStatusbar *statusbar_script = GTK_STATUSBAR(lookup_widget("statusbar_script"));
	gchar *status;
	gchar *newline;

	update_log_icon(TRUE);

	newline = g_strdup(data->myline);
	status = g_strdup_printf(_("Processing line %d: %s"), data->line, newline);

	gtk_statusbar_push(statusbar_script, 0, status);
	g_free(newline);
	g_free(status);
	g_free(data->myline);
	free(data);

	return FALSE;	// only run once
}

static void display_command_on_status_bar(int line, char *myline) {
	if (!com.headless) {
		struct log_status_bar_idle_data *data;

		data = (struct log_status_bar_idle_data *)malloc(sizeof(struct log_status_bar_idle_data));
		data->line = line;
		data->myline = myline ? g_strdup(myline) : NULL;
		gdk_threads_add_idle(log_status_bar_idle_callback, data);
	}
}

int check_command_mode() {
	/* until we have a proper implementation of modes, we just forbid other
	 * commands to be run during live stacking */
	int retval = 0;
	if (livestacking_is_started()) {
		retval = g_ascii_strcasecmp(word[0], "livestack") &&
			g_ascii_strcasecmp(word[0], "stop_ls") &&
			g_ascii_strcasecmp(word[0], "exit");
		if (retval)
			siril_log_message(_("This command cannot be run while live stacking is active, ignoring.\n"));

	}
	return retval;
}

static void clear_status_bar() {
	GtkStatusbar *bar = GTK_STATUSBAR(lookup_widget("statusbar_script"));
	gtk_statusbar_remove_all(bar, 0);
	update_log_icon(FALSE);
}

static gboolean end_script(gpointer p) {
	/* GTK+ code is ignored during scripts, this is a good place to redraw everything */
	clear_status_bar();
	gui_function(set_GUI_CWD, NULL);
	gui_function(update_MenuItem, NULL);
	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	update_zoom_label();
	update_display_fwhm();
	display_filename();
	gui_function(new_selection_zone, NULL);
	update_spinCPU(GINT_TO_POINTER(0));
	set_cursor_waiting(FALSE);
	return FALSE;
}

gpointer stellina_execute_script(gpointer p) {
	GInputStream *input_stream = (GInputStream*) p;
	gboolean checked_requires = FALSE;
	gchar *buffer;
	int line = 0, retval = 0;
	int wordnb;
	int startmem, endmem;
	struct timeval t_start, t_end;

	com.script = TRUE;
	com.stop_script = FALSE;
	com.script_thread_exited = FALSE;

	gettimeofday(&t_start, NULL);

	/* Now we want to save the cwd in order to come back after
	 * script execution
	 */
	gchar *saved_cwd = g_strdup(com.wd);
	startmem = get_available_memory() / BYTES_IN_A_MB;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
    int overall_retval = 0;
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
			NULL, NULL))) {
		++line;
		if (com.stop_script) {
			retval = 1;
			g_free (buffer);
			break;
		}
		/* Displays comments */
		if (buffer[0] == '#') {
			siril_log_color_message("%s\n", "blue", buffer);
			g_free (buffer);
			continue;
		}

		/* in Windows case, remove trailing CR */
		remove_trailing_cr(buffer);
		g_strstrip(buffer);

		if (buffer[0] == '\0') {
			g_free (buffer);
			continue;
		}

		display_command_on_status_bar(line, buffer);
		parse_line(buffer, length, &wordnb);
		if (check_requires(&checked_requires, com.pref.script_check_requires)) {
			g_free (buffer);
			break;
		}
		if (check_command_mode()) {
			g_free (buffer);
			continue;
		};

		retval = execute_command(wordnb);
		remove_child_from_children((GPid) -2); // remove the processing thread child
		// reference (speculative - not always necessary, but simplest to
		// call it every time just in case the command ran in the thread.)
		if (retval && retval != CMD_NO_WAIT) {
		  siril_log_message(_("Error in line %d ('%s'): %s.\n"), line, buffer, cmd_err_to_str((cmd_errors)retval));
			siril_log_message(_("Exiting batch processing.\n"));
			g_free (buffer);
			break;
		}
		if (retval != CMD_NO_WAIT && waiting_for_thread()) {
			overall_retval |= 1;
			g_free (buffer);
			break;	// abort script on command failure
		}
		endmem = get_available_memory() / BYTES_IN_A_MB;
		siril_debug_print("End of command %s, memory difference: %d MB\n", word[0], startmem - endmem);
		startmem = endmem;
		memset(word, 0, sizeof word);
		g_free (buffer);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);

	if (!com.headless) {
		com.script = FALSE;
		gui_function(end_script, NULL);
	}

	/* Now we want to restore the saved cwd */
	siril_change_dir(saved_cwd, NULL);
	writeinitfile();
	if (!overall_retval) {
		siril_log_message(_("Script execution finished successfully.\n"));
		gettimeofday(&t_end, NULL);
		show_time_msg(t_start, t_end, _("Total execution time"));
	} else {
		char *msg = siril_log_message(_("Script execution failed.\n"));
		msg[strlen(msg) - 1] = '\0';
		set_progress_bar_data(msg, PROGRESS_DONE);
	}
	g_free(saved_cwd);

	if (com.script_thread) {
		siril_debug_print("Script thread exiting\n");
		com.script_thread_exited = TRUE;
	}
	/* If called from the GUI, re-enable widgets blocked during the script */
	siril_add_idle(script_widgets_idle, NULL);
	return GINT_TO_POINTER(retval);
}

int cmd_stellina_script(int nb) {
  if (nb < 2) {
        siril_log_color_message("Usage: stellina_script file\n", "red");
        return CMD_ARG_ERROR;
  }
  char *script_file = word[1];
  GError *error = NULL;
  GFile *file = g_file_new_for_path(script_file);
  GInputStream *input_stream = (GInputStream*) g_file_read(file, NULL, &error);

  if (input_stream == NULL) {
    if (error != NULL) {
      g_clear_error(&error);
      siril_log_message(_("File [%s] does not exist\n"), script_file);
    }
    g_object_unref(file);
    return CMD_GENERIC_ERROR;
  }
  auto script_thread = g_thread_new("script", stellina_execute_script, input_stream);
  return CMD_OK;
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
