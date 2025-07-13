#include <QApplication>
#include <QStyleFactory>
#include <QDir>
#include <QStandardPaths>
#include <QDebug>

#include "StellinaProcessor.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    // Set application properties
    app.setApplicationName("Stellina Processor");
    app.setApplicationDisplayName("Stellina Processor for Siril");
    app.setApplicationVersion("1.0.0");
    app.setOrganizationName("StellinaTools");
    app.setOrganizationDomain("stellinatools.org");
    
    // Set a modern style
    app.setStyle(QStyleFactory::create("Fusion"));
    
    // Apply dark palette (optional)
    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window, QColor(53, 53, 53));
    darkPalette.setColor(QPalette::WindowText, Qt::white);
    darkPalette.setColor(QPalette::Base, QColor(25, 25, 25));
    darkPalette.setColor(QPalette::AlternateBase, QColor(53, 53, 53));
    darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
    darkPalette.setColor(QPalette::ToolTipText, Qt::white);
    darkPalette.setColor(QPalette::Text, Qt::white);
    darkPalette.setColor(QPalette::Button, QColor(53, 53, 53));
    darkPalette.setColor(QPalette::ButtonText, Qt::white);
    darkPalette.setColor(QPalette::BrightText, Qt::red);
    darkPalette.setColor(QPalette::Link, QColor(42, 130, 218));
    darkPalette.setColor(QPalette::Highlight, QColor(42, 130, 218));
    darkPalette.setColor(QPalette::HighlightedText, Qt::black);
    
    // Comment out the next line if you prefer the system theme
    // app.setPalette(darkPalette);
    
    // Create and show main window
    StellinaProcessor window;
    window.show();
    
    qDebug() << "Stellina Processor started";
    qDebug() << "Qt version:" << QT_VERSION_STR;
    qDebug() << "Application data location:" << QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    
    return app.exec();
}
