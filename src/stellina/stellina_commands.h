/*
 * stellina_commands.h
 * C-only header for Stellina command declarations
 * This can be safely included in C files like command_list.h
 */

#ifndef STELLINA_COMMANDS_H
#define STELLINA_COMMANDS_H

#ifdef __cplusplus
extern "C" {
#endif

// Siril command interface functions - correct signature based on error
// Functions should take int parameter, not struct command_info *
int cmd_process_stellina(int nb);
int cmd_stellina_config(int nb);
int cmd_stellina_stats(int nb);

#ifdef __cplusplus
}
#endif

#define STR_PROCESS N_("Stellina processing")
#define STR_CONFIG N_("Stellina config")
#define STR_STATS N_("Stellina stats")

#endif // STELLINA_COMMANDS_H
