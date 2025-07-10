/*
 * stellina_commands.h
 * C-only header for Stellina command declarations
 * This can be safely included in C files like command_list.h
 */

#ifndef STELLINA_COMMANDS_H
#define STELLINA_COMMANDS_H

#include "core/command.h"

#ifdef __cplusplus
extern "C" {
#endif

// Siril command interface functions - C declarations only
int cmd_process_stellina(struct command_info *cmd);
int cmd_stellina_config(struct command_info *cmd);
int cmd_stellina_stats(struct command_info *cmd);

#ifdef __cplusplus
}
#endif

#endif // STELLINA_COMMANDS_H