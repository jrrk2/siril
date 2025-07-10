/*
 * stellina_integration.cpp
 * Integration of Stellina commands into Siril's command system
 * 
 * This file shows how to register the Stellina commands with Siril.
 * You would typically add these entries to Siril's main command_list.h
 */

#include "stellina_commands.h"
#include "core/command.h"

// Add these entries to Siril's command_list.h in the appropriate location:

/*
// In the command_word enum (alphabetically ordered):
typedef enum {
    // ... other commands ...
    CMD_STELLINA_CONFIG,
    CMD_STELLINA_PROCESS,
    CMD_STELLINA_STATS,
    CMD_STELLINA_TEST_COORDS,
    // ... other commands ...
} command_word;

// In the commands[] array:
{
    "stellina_config",           // command name
    cmd_stellina_config,         // function pointer
    "stellina_config parameter value", // usage
    STR_STELLINA_CONFIG,         // description
    HELP_STELLINA_CONFIG,        // help text
    FALSE,                       // requires image
    FALSE,                       // requires selection
    CMD_STELLINA_CONFIG          // command enum
},
{
    "stellina_process", 
    cmd_process_stellina,
    "stellina_process input_dir output_dir [quality_filter]",
    STR_PROCESS_STELLINA,
    HELP_PROCESS_STELLINA,
    FALSE,
    FALSE,
    CMD_STELLINA_PROCESS
},
{
    "stellina_stats",
    cmd_stellina_stats,
    "stellina_stats [fits_file json_file]",
    STR_STELLINA_STATS,
    HELP_STELLINA_STATS,
    FALSE,
    FALSE,
    CMD_STELLINA_STATS
},
{
    "stellina_test_coords",
    cmd_stellina_test_coords,
    "stellina_test_coords alt az date_obs [lat lon alt_m]",
    STR_STELLINA_TEST_COORDS,
    HELP_STELLINA_TEST_COORDS,
    FALSE,
    FALSE,
    CMD_STELLINA_TEST_COORDS
},
*/

/**
 * Initialize the Stellina extension
 * Call this from Siril's initialization code
 */
void init_stellina_extension() {
    stellina_commands_init();
    siril_log_color_message("Stellina extension initialized successfully\n", "green");
}

/**
 * Example of how to add menu items for the Stellina commands
 * This would typically go in Siril's GUI initialization
 */
void stellina_add_menu_items() {
    // This is pseudo-code showing how you might add menu items
    // The actual implementation depends on Siril's menu system
    
    /*
    GtkWidget *stellina_menu = gtk_menu_new();
    GtkWidget *stellina_menu_item = gtk_menu_item_new_with_label("Stellina");
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(stellina_menu_item), stellina_menu);
    
    // Add submenu items
    GtkWidget *process_item = gtk_menu_item_new_with_label("Process Directory...");
    g_signal_connect(process_item, "activate", G_CALLBACK(on_stellina_process_activate), NULL);
    gtk_menu_shell_append(GTK_MENU_SHELL(stellina_menu), process_item);
    
    GtkWidget *config_item = gtk_menu_item_new_with_label("Configuration...");
    g_signal_connect(config_item, "activate", G_CALLBACK(on_stellina_config_activate), NULL);
    gtk_menu_shell_append(GTK_MENU_SHELL(stellina_menu), config_item);
    
    // Add to main menu bar
    gtk_menu_shell_append(GTK_MENU_SHELL(menu_bar), stellina_menu_item);
    */
}

/**
 * Build instructions for integrating into Siril
 */
/*
To integrate this Stellina extension into Siril:

1. Copy the stellina/ directory into src/
2. Add the stellina library to the main meson.build:
   - Add stellina_dep to siril_dependencies
   - Make sure the stellina subdirectory is processed

3. In src/core/command_list.h:
   - Add the command enums to command_word
   - Add the command entries to the commands[] array
   - Include "stellina/stellina_commands.h"

4. In src/main.c or appropriate initialization file:
   - Call init_stellina_extension() during startup

5. For menu integration (optional):
   - Add menu items in the GUI initialization code
   - Connect menu items to command execution

6. Build requirements:
   - Ensure libnova is available (already handled in meson.build)
   - json-glib-1.0 dependency (should already be in Siril)
   - Standard CFITSIO and GSL dependencies

Example build on macOS with Homebrew:
$ brew install libnova
$ meson setup build --buildtype=release
$ ninja -C build

The extension will then be available through:
- Command line: stellina_process /path/to/images /path/to/output
- Script: load script.ssf (with stellina commands inside)
- GUI: Through menu items if implemented
*/