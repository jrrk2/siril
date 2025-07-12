#ifndef STELLINAPROCESSOR_H
#define STELLINAPROCESSOR_H

#include <QMainWindow>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QTextEdit>
#include <QProgressBar>
#include <QLabel>
#include <QGroupBox>
#include <QCheckBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QTimer>
#include <QStringList>
#include <QStandardPaths>
#include <QJsonObject>
#include <QJsonDocument>

#include "SirilClient.h"

class StellinaProcessor : public QMainWindow {
    Q_OBJECT

public:
    explicit StellinaProcessor(QWidget *parent = nullptr);
    ~StellinaProcessor();

private slots:
    // UI slots
    void onSelectSourceDirectory();
    void onSelectOutputDirectory();
    void onStartProcessing();
    void onStopProcessing();
    void onTestConnection();
    void onClearLog();
    
    // Siril connection slots
    void onSirilConnected();
    void onSirilDisconnected();
    void onSirilCommandExecuted(const QString &command, bool success);
    void onSirilError(const QString &error);
    
    // Processing slots
    void onProcessingTimer();

private:
    void setupUI();
    void setupMenu();
    void connectSignals();
    void updateUI();
    void updateConnectionStatus();
    void logMessage(const QString &message, const QString &color = "black");
    
    // Processing functions
    void startStellinaProcessing();
    void processNextImage();
    void finishProcessing();
    bool findStellinaImages();
    QJsonObject loadStellinaJson(const QString &jsonPath);
    bool extractCoordinates(const QJsonObject &json, double &alt, double &az);
    bool convertAltAzToRaDec(double alt, double az, const QString &dateObs, double &ra, double &dec);
    bool checkStellinaQuality(const QJsonObject &json);
    
    // UI components
    QWidget *m_centralWidget;
    QVBoxLayout *m_mainLayout;
    
    // Connection group
    QGroupBox *m_connectionGroup;
    QPushButton *m_testConnectionButton;
    QLabel *m_connectionStatus;
    
    // Input group
    QGroupBox *m_inputGroup;
    QLineEdit *m_sourceDirectoryEdit;
    QPushButton *m_selectSourceButton;
    QLineEdit *m_outputDirectoryEdit;
    QPushButton *m_selectOutputButton;
    
    // Options group
    QGroupBox *m_optionsGroup;
    QCheckBox *m_qualityFilterCheck;
    QCheckBox *m_debugModeCheck;
    QDoubleSpinBox *m_focalLengthSpin;
    QDoubleSpinBox *m_pixelSizeSpin;
    QLineEdit *m_observerLocationEdit;
    
    // Processing group
    QGroupBox *m_processingGroup;
    QPushButton *m_startButton;
    QPushButton *m_stopButton;
    QProgressBar *m_progressBar;
    QLabel *m_progressLabel;
    
    // Log group
    QGroupBox *m_logGroup;
    QTextEdit *m_logTextEdit;
    QPushButton *m_clearLogButton;
    
    // Status bar
    QLabel *m_statusLabel;
    
    // Core components
    SirilClient *m_sirilClient;
    QTimer *m_processingTimer;
    
    // Processing state
    bool m_processing;
    QStringList m_imagesToProcess;
    int m_currentImageIndex;
    int m_processedCount;
    int m_errorCount;
    int m_skippedCount;
    
    // Settings
    QString m_sourceDirectory;
    QString m_outputDirectory;
    bool m_qualityFilter;
    bool m_debugMode;
    double m_focalLength;
    double m_pixelSize;
    QString m_observerLocation;
};

#endif // STELLINAPROCESSOR_H