QT += core widgets network

CONFIG += c++17

TARGET = StellinaProcessor
TEMPLATE = app

# Version info
VERSION = 1.0.0
QMAKE_TARGET_PRODUCT = "Stellina Processor for Siril"
QMAKE_TARGET_DESCRIPTION = "Qt application for processing Stellina telescope images using Siril"
QMAKE_TARGET_COPYRIGHT = "Copyright (C) 2025"

# macOS specific settings
macx {
    QMAKE_INFO_PLIST = Info.plist
    ICON = stellina_processor.icns
    
    # macOS deployment target
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.15
    
    # Bundle identifier
    QMAKE_TARGET_BUNDLE_PREFIX = org.stellinatools
}

# Windows specific settings
win32 {
    RC_ICONS = stellina_processor.ico
    VERSION_PE_HEADER = $$VERSION
}

# Linux specific settings
unix:!macx {
    target.path = /usr/local/bin
    INSTALLS += target
}

# Include directories
INCLUDEPATH += .

# Source files
SOURCES += \
    main.cpp \
    StellinaProcessor_Core.cpp \
    stellinaprocessor_plugin_import.cpp \
    StellinaProcessor_Processing.cpp \
    StellinaProcessor_Slots.cpp \
    StellinaProcessor_UI.cpp \

# Header files
HEADERS += \
    StellinaProcessor.h \

# Resources (if you add icons, etc.)
# RESOURCES += resources.qrc

# Compiler flags
QMAKE_CXXFLAGS += -Wall -Wextra

# Debug/Release specific settings
CONFIG(debug, debug|release) {
    DEFINES += DEBUG_BUILD
    TARGET = $${TARGET}_debug
} else {
    DEFINES += RELEASE_BUILD
    QMAKE_CXXFLAGS += -O2
}

# Enable all Qt warnings
CONFIG += warn_on

# Additional libraries (if needed)
LIBS += -lm -lcfitsio

# Deployment settings for macOS
macx {
    # Create app bundle with frameworks
    CONFIG += app_bundle
    
    # For distribution
    # QMAKE_POST_LINK += macdeployqt $$OUT_PWD/$${TARGET}.app
    INCLUDEPATH += /opt/homebrew/include
    LIBS += -L/opt/homebrew/lib -lnova
    }

# Build output directories
CONFIG(debug, debug|release) {
    DESTDIR = build/debug
    OBJECTS_DIR = build/debug/obj
    MOC_DIR = build/debug/moc
    RCC_DIR = build/debug/rcc
    UI_DIR = build/debug/ui
} else {
    DESTDIR = build/release
    OBJECTS_DIR = build/release/obj
    MOC_DIR = build/release/moc
    RCC_DIR = build/release/rcc
    UI_DIR = build/release/ui
}
