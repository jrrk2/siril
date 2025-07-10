/*
 * stellina_commands.cpp
 * Minimal Siril command implementations for Stellina processing
 */

#include "stellina_processor.h"
#include <iostream>
#include <cstring>

// Placeholder command implementations
int cmd_process_stellina(struct command_info *cmd) {
    std::cout << "Stellina: process_stellina command called" << std::endl;
    // For now, just return success
    return 0;
}

int cmd_stellina_config(struct command_info *cmd) {
    std::cout << "Stellina: stellina_config command called" << std::endl;
    return 0;
}

int cmd_stellina_stats(struct command_info *cmd) {
    std::cout << "Stellina: stellina_stats command called" << std::endl;
    return 0;
}