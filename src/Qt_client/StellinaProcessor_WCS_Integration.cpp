// StellinaProcessor_WCS_Integration.cpp
// Implementation file for WCS astrometric stacker integration

#include "StellinaProcessor.h"
#include "WcsAstrometricStacker.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QTextEdit>
#include <QDialog>
#include <QSettings>
#include <QDateTime>
#include <QDir>
#include <QMenu>
#include <QMenuBar>
#include <QFileInfo>
#include <QRegularExpression>

void StellinaProcessor::initializeWCSStacker() {
    // Initialize WCS stacker
    m_wcsStacker = new WCSAstrometricStacker(this);
    
    // Set up default WCS stacking parameters
    m_wcsStackingParams.combination = StackingParams::WEIGHTED_MEAN;
    m_wcsStackingParams.rejection = StackingParams::SIGMA_CLIPPING;
    m_wcsStackingParams.sigma_low = 3.0;
    m_wcsStackingParams.sigma_high = 3.0;
    m_wcsStackingParams.normalize_exposure = true;
    m_wcsStackingParams.create_weight_map = true;
    m_wcsStackingParams.output_format = "fits";
    
    // Connect WCS stacker signals
    connect(m_wcsStacker, &WCSAstrometricStacker::progressUpdated,
            this, &StellinaProcessor::onWCSProgressUpdated);
    connect(m_wcsStacker, &WCSAstrometricStacker::statusUpdated,
            this, &StellinaProcessor::onWCSStatusUpdated);
    connect(m_wcsStacker, &WCSAstrometricStacker::stackingComplete,
            this, &StellinaProcessor::onWCSStackingComplete);
}

void StellinaProcessor::setupWCSStackingUI() {
    // WCS-based astrometric stacking group
    m_wcsStackingGroup = new QGroupBox("Advanced WCS Astrometric Stacking");
    QGridLayout *wcsLayout = new QGridLayout(m_wcsStackingGroup);
    
    // Add explanation
    QLabel *wcsInfo = new QLabel("WCS-based stacking uses world coordinate systems for precise "
                                "astrometric alignment. Requires plate-solved images with WCS headers.");
    wcsInfo->setWordWrap(true);
    wcsInfo->setStyleSheet("color: blue; font-style: italic; margin-bottom: 10px;");
    wcsLayout->addWidget(wcsInfo, 0, 0, 1, 4);
    
    // Combination method
    wcsLayout->addWidget(new QLabel("Combination Method:"), 1, 0);
    m_wcsCombinationMethodCombo = new QComboBox;
    m_wcsCombinationMethodCombo->addItem("Weighted Mean", static_cast<int>(StackingParams::WEIGHTED_MEAN));
    m_wcsCombinationMethodCombo->addItem("Sigma Clipped Mean", static_cast<int>(StackingParams::SIGMA_CLIPPED_MEAN));
    m_wcsCombinationMethodCombo->addItem("Median", static_cast<int>(StackingParams::MEDIAN));
    m_wcsCombinationMethodCombo->addItem("Mean", static_cast<int>(StackingParams::MEAN));
    wcsLayout->addWidget(m_wcsCombinationMethodCombo, 1, 1, 1, 3);
    
    // Rejection method
    wcsLayout->addWidget(new QLabel("Rejection Method:"), 2, 0);
    m_wcsRejectionMethodCombo = new QComboBox;
    m_wcsRejectionMethodCombo->addItem("No Rejection", static_cast<int>(StackingParams::NO_REJECTION));
    m_wcsRejectionMethodCombo->addItem("Sigma Clipping", static_cast<int>(StackingParams::SIGMA_CLIPPING));
    m_wcsRejectionMethodCombo->addItem("Percentile Clipping", static_cast<int>(StackingParams::PERCENTILE_CLIPPING));
    m_wcsRejectionMethodCombo->setCurrentIndex(1); // Default to sigma clipping
    wcsLayout->addWidget(m_wcsRejectionMethodCombo, 2, 1, 1, 3);
    
    // Sigma clipping parameters
    wcsLayout->addWidget(new QLabel("Sigma Low:"), 3, 0);
    m_wcsSigmaLowSpin = new QDoubleSpinBox;
    m_wcsSigmaLowSpin->setRange(0.1, 10.0);
    m_wcsSigmaLowSpin->setValue(3.0);
    m_wcsSigmaLowSpin->setDecimals(1);
    wcsLayout->addWidget(m_wcsSigmaLowSpin, 3, 1);
    
    wcsLayout->addWidget(new QLabel("Sigma High:"), 3, 2);
    m_wcsSigmaHighSpin = new QDoubleSpinBox;
    m_wcsSigmaHighSpin->setRange(0.1, 10.0);
    m_wcsSigmaHighSpin->setValue(3.0);
    m_wcsSigmaHighSpin->setDecimals(1);
    wcsLayout->addWidget(m_wcsSigmaHighSpin, 3, 3);
    
    // Advanced options
    m_wcsNormalizeExposureCheck = new QCheckBox("Normalize by exposure time");
    m_wcsNormalizeExposureCheck->setChecked(true);
    wcsLayout->addWidget(m_wcsNormalizeExposureCheck, 4, 0, 1, 2);
    
    m_wcsCreateWeightMapCheck = new QCheckBox("Create weight map");
    m_wcsCreateWeightMapCheck->setChecked(true);
    wcsLayout->addWidget(m_wcsCreateWeightMapCheck, 4, 2, 1, 2);
    
    // Output dimensions
    wcsLayout->addWidget(new QLabel("Output Width:"), 5, 0);
    m_wcsOutputWidthSpin = new QSpinBox;
    m_wcsOutputWidthSpin->setRange(0, 10000);
    m_wcsOutputWidthSpin->setValue(0); // 0 = auto
    m_wcsOutputWidthSpin->setSpecialValueText("Auto");
    wcsLayout->addWidget(m_wcsOutputWidthSpin, 5, 1);
    
    wcsLayout->addWidget(new QLabel("Output Height:"), 5, 2);
    m_wcsOutputHeightSpin = new QSpinBox;
    m_wcsOutputHeightSpin->setRange(0, 10000);
    m_wcsOutputHeightSpin->setValue(0); // 0 = auto
    m_wcsOutputHeightSpin->setSpecialValueText("Auto");
    wcsLayout->addWidget(m_wcsOutputHeightSpin, 5, 3);
    
    // Output pixel scale
    wcsLayout->addWidget(new QLabel("Pixel Scale (\"/px):"), 6, 0);
    m_wcsOutputPixelScaleSpin = new QDoubleSpinBox;
    m_wcsOutputPixelScaleSpin->setRange(0.0, 10.0);
    m_wcsOutputPixelScaleSpin->setValue(0.0); // 0 = auto
    m_wcsOutputPixelScaleSpin->setDecimals(2);
    m_wcsOutputPixelScaleSpin->setSpecialValueText("Auto");
    wcsLayout->addWidget(m_wcsOutputPixelScaleSpin, 6, 1, 1, 2);
    
    // Control buttons
    QHBoxLayout *wcsButtonLayout = new QHBoxLayout;
    m_startWCSStackingButton = new QPushButton("Start WCS Stacking");
    m_startWCSStackingButton->setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; }");
    m_saveWCSResultButton = new QPushButton("Save WCS Result");
    m_saveWCSResultButton->setEnabled(false);
    
    wcsButtonLayout->addWidget(m_startWCSStackingButton);
    wcsButtonLayout->addWidget(m_saveWCSResultButton);
    wcsButtonLayout->addStretch();
    
    wcsLayout->addLayout(wcsButtonLayout, 7, 0, 1, 4);
    
    // Add to stacking tab layout
    QVBoxLayout *stackingTabLayout = qobject_cast<QVBoxLayout*>(m_stackingTab->layout());
    if (stackingTabLayout) {
        stackingTabLayout->addWidget(m_wcsStackingGroup);
    }
    
    // Connect signals
    connect(m_wcsCombinationMethodCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsRejectionMethodCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsSigmaLowSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsSigmaHighSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsNormalizeExposureCheck, &QCheckBox::toggled,
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsCreateWeightMapCheck, &QCheckBox::toggled,
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsOutputWidthSpin, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsOutputHeightSpin, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    connect(m_wcsOutputPixelScaleSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StellinaProcessor::onWCSParametersChanged);
    
    connect(m_startWCSStackingButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStartWCSStacking);
    connect(m_saveWCSResultButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSaveWCSResult);
}

void StellinaProcessor::onWCSParametersChanged() {
    // Update WCS stacking parameters from UI
    m_wcsStackingParams.combination = static_cast<StackingParams::CombinationMethod>(
        m_wcsCombinationMethodCombo->currentData().toInt());
    m_wcsStackingParams.rejection = static_cast<StackingParams::RejectionMethod>(
        m_wcsRejectionMethodCombo->currentData().toInt());
    m_wcsStackingParams.sigma_low = m_wcsSigmaLowSpin->value();
    m_wcsStackingParams.sigma_high = m_wcsSigmaHighSpin->value();
    m_wcsStackingParams.normalize_exposure = m_wcsNormalizeExposureCheck->isChecked();
    m_wcsStackingParams.create_weight_map = m_wcsCreateWeightMapCheck->isChecked();
    m_wcsStackingParams.output_width = m_wcsOutputWidthSpin->value();
    m_wcsStackingParams.output_height = m_wcsOutputHeightSpin->value();
    m_wcsStackingParams.output_pixel_scale = m_wcsOutputPixelScaleSpin->value();
    
    // Update the stacker with new parameters
    m_wcsStacker->setStackingParameters(m_wcsStackingParams);
    
    if (m_debugMode) {
        logMessage(QString("WCS stacking parameters updated: method=%1, rejection=%2")
                      .arg(m_wcsCombinationMethodCombo->currentText())
                      .arg(m_wcsRejectionMethodCombo->currentText()), "gray");
    }
}

void StellinaProcessor::onStartWCSStacking() {
    if (m_plateSolvedDirectory.isEmpty()) {
        QMessageBox::warning(this, "Directory Error", 
                            "Please select the plate-solved images directory first.");
        return;
    }
    
    QDir plateSolvedDir(m_plateSolvedDirectory);
    if (!plateSolvedDir.exists()) {
        QMessageBox::warning(this, "Directory Error", 
                            "Plate-solved directory does not exist.");
        return;
    }
    
    // Find plate-solved FITS files
    QStringList plateSolvedFiles = plateSolvedDir.entryList(
        QStringList() << "plate_solved_*.fits" << "*solved*.fits", 
        QDir::Files);
    
    if (plateSolvedFiles.isEmpty()) {
        QMessageBox::warning(this, "No Images", 
                            "No plate-solved FITS files found in the selected directory.\n"
                            "Please run plate solving first.");
        return;
    }
    
    if (plateSolvedFiles.size() < 3) {
        QMessageBox::warning(this, "Insufficient Images", 
                            "Need at least 3 plate-solved images for effective stacking.");
        return;
    }
    
    logMessage(QString("Starting WCS astrometric stacking with %1 images...")
                  .arg(plateSolvedFiles.size()), "blue");
    
    // Set up progress tracking
    m_wcsStacker->setProgressWidgets(m_progressBar, m_currentTaskLabel);
    
    // Load images into WCS stacker
    int successfullyLoaded = 0;
    for (const QString &fileName : plateSolvedFiles) {
        QString fullPath = plateSolvedDir.absoluteFilePath(fileName);
        
        // Try to find corresponding StellinaImageData
        StellinaImageData *imageData = findImageDataByPath(fullPath);
        if (imageData) {
            if (m_wcsStacker->addImageWithMetadata(fullPath, *imageData)) {
                successfullyLoaded++;
            }
        } else {
            // If no StellinaImageData found, try loading just the FITS file
            if (m_wcsStacker->addImage(fullPath)) {
                successfullyLoaded++;
            }
        }
    }
    
    if (successfullyLoaded < 3) {
        QMessageBox::critical(this, "Loading Error", 
                             QString("Only %1 images could be loaded successfully.\n"
                                   "Need at least 3 images for stacking.").arg(successfullyLoaded));
        return;
    }
    
    logMessage(QString("Successfully loaded %1 images for WCS stacking").arg(successfullyLoaded), "green");
    
    // Update UI state
    m_startWCSStackingButton->setEnabled(false);
    m_startWCSStackingButton->setText("Stacking...");
    m_saveWCSResultButton->setEnabled(false);
    
    // Start the stacking process
    m_wcsStacker->startStacking();
}

void StellinaProcessor::onWCSStackingComplete(bool success) {
    // Restore UI state
    m_startWCSStackingButton->setEnabled(true);
    m_startWCSStackingButton->setText("Start WCS Stacking");
    
    if (success) {
        logMessage("WCS astrometric stacking completed successfully!", "green");
        
        // Enable save button
        m_saveWCSResultButton->setEnabled(true);
        
        // Show completion statistics
        QString stats = QString("Stacking Statistics:\n"
                               "- Images processed: %1\n"
                               "- Total exposure time: %2 seconds\n"
                               "- Average quality: %3\n"
                               "- Output dimensions: %4x%5 pixels")
                           .arg(m_wcsStacker->getImageCount())
                           .arg(m_wcsStacker->getTotalExposureTime(), 0, 'f', 1)
                           .arg(m_wcsStacker->getAverageQuality(), 0, 'f', 3)
                           .arg(m_wcsStacker->getStackedImage().cols)
                           .arg(m_wcsStacker->getStackedImage().rows);
        
        logMessage(stats, "blue");
        
        // Generate and log quality report
        QString qualityReport = m_wcsStacker->getQualityReport();
        logMessage("=== WCS Stacking Quality Report ===", "blue");
        QStringList reportLines = qualityReport.split('\n');
        for (const QString &line : reportLines) {
            if (!line.isEmpty()) {
                logMessage(line, "gray");
            }
        }
        
        QMessageBox::information(this, "Stacking Complete", 
                                "WCS astrometric stacking completed successfully!\n\n" + stats);
    } else {
        logMessage("WCS astrometric stacking failed.", "red");
        m_saveWCSResultButton->setEnabled(false);
        
        QMessageBox::critical(this, "Stacking Failed", 
                             "WCS astrometric stacking failed. Check the log for details.");
    }
    
    // Update processing status
    updateProcessingStatus();
}

void StellinaProcessor::onWCSProgressUpdated(int percentage) {
    // Update progress bar
    if (m_progressBar) {
        m_progressBar->setValue(percentage);
    }
}

void StellinaProcessor::onWCSStatusUpdated(const QString &message) {
    // Update status label and log
    if (m_currentTaskLabel) {
        m_currentTaskLabel->setText(message);
    }
    logMessage(message, "blue");
}

void StellinaProcessor::onSaveWCSResult() {
    if (m_wcsStacker->getStackedImage().empty()) {
        QMessageBox::warning(this, "No Result", "No stacked image available to save.");
        return;
    }
    
    // Get output directory - use stacked directory if set, otherwise plate-solved directory
    QString defaultDir = m_stackedDirectory.isEmpty() ? m_plateSolvedDirectory : m_stackedDirectory;
    QString defaultName = QString("wcs_stacked_%1.fits")
                            .arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    QString defaultPath = QDir(defaultDir).absoluteFilePath(defaultName);
    
    QString outputPath = QFileDialog::getSaveFileName(
        this,
        "Save WCS Stacked Image",
        defaultPath,
        "FITS Images (*.fits *.fit);;All Files (*)"
    );
    
    if (outputPath.isEmpty()) {
        return;
    }
    
    logMessage(QString("Saving WCS stacked result to: %1").arg(outputPath), "blue");
    
    if (m_wcsStacker->saveResult(outputPath)) {
        logMessage("WCS stacked image saved successfully!", "green");
        
        // Also save quality report
        QString reportPath = outputPath;
        reportPath.replace(QRegularExpression("\\.(fits|fit)$"), "_quality_report.txt");
        
        QFile reportFile(reportPath);
        if (reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&reportFile);
            out << m_wcsStacker->getQualityReport();
            reportFile.close();
            
            logMessage(QString("Quality report saved to: %1").arg(reportPath), "blue");
        }
        
        QMessageBox::information(this, "Save Complete", 
                                QString("WCS stacked image saved successfully to:\n%1\n\n"
                                       "Quality report saved to:\n%2").arg(outputPath).arg(reportPath));
    } else {
        logMessage("Failed to save WCS stacked image.", "red");
        QMessageBox::critical(this, "Save Failed", "Failed to save the WCS stacked image.");
    }
}

bool StellinaProcessor::performAstrometricStackingEnhanced() {
    if (m_plateSolvedFiles.isEmpty()) {
        logMessage("No plate-solved images available for stacking", "red");
        return false;
    }
    
    if (m_plateSolvedFiles.size() < 3) {
        logMessage("Need at least 3 images for stacking", "red");
        return false;
    }
    
        logMessage("Using WCS-based astrometric stacking", "blue");
        
        // Set up WCS stacker with current plate-solved files
        m_wcsStacker->setProgressWidgets(m_subTaskProgressBar, m_currentTaskLabel);
        
        int loaded = 0;
        for (const QString &filePath : m_plateSolvedFiles) {
            StellinaImageData *imageData = findImageDataByPath(filePath);
            if (imageData) {
                if (m_wcsStacker->addImageWithMetadata(filePath, *imageData)) {
                    loaded++;
                }
            } else {
                if (m_wcsStacker->addImage(filePath)) {
                    loaded++;
                }
            }
        }
        
        if (loaded < 3) {
            logMessage(QString("Only loaded %1 images for WCS stacking").arg(loaded), "red");
            return false;
        }
        
        // Start WCS stacking
        if (!m_wcsStacker->stackImages()) {
            logMessage("WCS stacking failed", "red");
            return false;
        }
        
        // Save result automatically
        QString outputName = QString("wcs_stacked_%1.fits")
                               .arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
        QString outputPath = QDir(getOutputDirectoryForCurrentStage()).absoluteFilePath(outputName);
        
        if (m_wcsStacker->saveResult(outputPath)) {
            m_finalStackedImage = outputPath;
            logMessage(QString("WCS stacking completed: %1").arg(outputName), "green");
            return true;
        } else {
            logMessage("Failed to save WCS stacked result", "red");
            return false;
        }
        
}

void StellinaProcessor::addWCSMenuItems() {
    // This should be called from your setupMenu() function
    QMenu *wcsMenu = menuBar()->addMenu("WCS &Stacking");
    
    wcsMenu->addAction("&Start WCS Stacking", this, &StellinaProcessor::onStartWCSStacking);
    wcsMenu->addAction("&Save WCS Result", this, &StellinaProcessor::onSaveWCSResult);
    wcsMenu->addSeparator();
    wcsMenu->addAction("&Test WCS Loading", [this]() {
        QString testFile = QFileDialog::getOpenFileName(this, "Select Plate-Solved FITS", 
                                                       m_plateSolvedDirectory, "FITS Files (*.fits *.fit)");
        if (!testFile.isEmpty()) {
            WCSAstrometricStacker testStacker;
            if (testStacker.addImage(testFile)) {
                logMessage(QString("Successfully loaded WCS from: %1").arg(QFileInfo(testFile).fileName()), "green");
                logMessage(QString("Image count: %1").arg(testStacker.getImageCount()), "blue");
            } else {
                logMessage(QString("Failed to load WCS from: %1").arg(QFileInfo(testFile).fileName()), "red");
            }
        }
    });
    
    wcsMenu->addAction("WCS &Quality Analysis", [this]() {
        if (m_wcsStacker->getImageCount() > 0) {
            m_wcsStacker->analyzeImageQuality();
            QString report = m_wcsStacker->getQualityReport();
            
            QDialog *reportDialog = new QDialog(this);
            reportDialog->setWindowTitle("WCS Stacking Quality Report");
            reportDialog->setModal(true);
            reportDialog->resize(800, 600);
            
            QVBoxLayout *layout = new QVBoxLayout(reportDialog);
            QTextEdit *reportText = new QTextEdit;
            reportText->setFont(QFont("monospace", 9));
            reportText->setText(report);
            reportText->setReadOnly(true);
            
            QPushButton *closeButton = new QPushButton("Close");
            connect(closeButton, &QPushButton::clicked, reportDialog, &QDialog::accept);
            
            layout->addWidget(reportText);
            layout->addWidget(closeButton);
            
            reportDialog->exec();
        } else {
            QMessageBox::information(this, "No Data", "No images loaded in WCS stacker for analysis.");
        }
    });
}

void StellinaProcessor::updateWCSUI() {
    bool hasPlatesolved = !m_plateSolvedDirectory.isEmpty() && 
                         QDir(m_plateSolvedDirectory).exists();
    bool hasWCSResult = false; // !m_wcsStacker->getStackedImage().empty();
    
    if (m_startWCSStackingButton) {
        m_startWCSStackingButton->setEnabled(hasPlatesolved && !m_processing);
    }
    if (m_saveWCSResultButton) {
        m_saveWCSResultButton->setEnabled(hasWCSResult && !m_processing);
    }
    
    // Update combination method enable/disable based on rejection method
    bool usingSigmaClipping = (m_wcsStackingParams.rejection == StackingParams::SIGMA_CLIPPING);
    if (m_wcsSigmaLowSpin) {
        m_wcsSigmaLowSpin->setEnabled(usingSigmaClipping);
    }
    if (m_wcsSigmaHighSpin) {
        m_wcsSigmaHighSpin->setEnabled(usingSigmaClipping);
    }
}

void StellinaProcessor::loadWCSSettings() {
    QSettings settings;
    
    int combination = settings.value("wcs/combinationMethod", 
                                    static_cast<int>(StackingParams::WEIGHTED_MEAN)).toInt();
    if (m_wcsCombinationMethodCombo) {
        m_wcsCombinationMethodCombo->setCurrentIndex(
            m_wcsCombinationMethodCombo->findData(combination));
    }
    
    int rejection = settings.value("wcs/rejectionMethod", 
                                  static_cast<int>(StackingParams::SIGMA_CLIPPING)).toInt();
    if (m_wcsRejectionMethodCombo) {
        m_wcsRejectionMethodCombo->setCurrentIndex(
            m_wcsRejectionMethodCombo->findData(rejection));
    }
    
    if (m_wcsSigmaLowSpin) {
        m_wcsSigmaLowSpin->setValue(settings.value("wcs/sigmaLow", 3.0).toDouble());
    }
    if (m_wcsSigmaHighSpin) {
        m_wcsSigmaHighSpin->setValue(settings.value("wcs/sigmaHigh", 3.0).toDouble());
    }
    if (m_wcsNormalizeExposureCheck) {
        m_wcsNormalizeExposureCheck->setChecked(settings.value("wcs/normalizeExposure", true).toBool());
    }
    if (m_wcsCreateWeightMapCheck) {
        m_wcsCreateWeightMapCheck->setChecked(settings.value("wcs/createWeightMap", true).toBool());
    }
    if (m_wcsOutputWidthSpin) {
        m_wcsOutputWidthSpin->setValue(settings.value("wcs/outputWidth", 0).toInt());
    }
    if (m_wcsOutputHeightSpin) {
        m_wcsOutputHeightSpin->setValue(settings.value("wcs/outputHeight", 0).toInt());
    }
    if (m_wcsOutputPixelScaleSpin) {
        m_wcsOutputPixelScaleSpin->setValue(settings.value("wcs/outputPixelScale", 0.0).toDouble());
    }
    
    // Update parameters
    onWCSParametersChanged();
}

void StellinaProcessor::saveWCSSettings() {
    QSettings settings;
    
    if (m_wcsCombinationMethodCombo) {
        settings.setValue("wcs/combinationMethod", m_wcsCombinationMethodCombo->currentData().toInt());
    }
    if (m_wcsRejectionMethodCombo) {
        settings.setValue("wcs/rejectionMethod", m_wcsRejectionMethodCombo->currentData().toInt());
    }
    if (m_wcsSigmaLowSpin) {
        settings.setValue("wcs/sigmaLow", m_wcsSigmaLowSpin->value());
    }
    if (m_wcsSigmaHighSpin) {
        settings.setValue("wcs/sigmaHigh", m_wcsSigmaHighSpin->value());
    }
    if (m_wcsNormalizeExposureCheck) {
        settings.setValue("wcs/normalizeExposure", m_wcsNormalizeExposureCheck->isChecked());
    }
    if (m_wcsCreateWeightMapCheck) {
        settings.setValue("wcs/createWeightMap", m_wcsCreateWeightMapCheck->isChecked());
    }
    if (m_wcsOutputWidthSpin) {
        settings.setValue("wcs/outputWidth", m_wcsOutputWidthSpin->value());
    }
    if (m_wcsOutputHeightSpin) {
        settings.setValue("wcs/outputHeight", m_wcsOutputHeightSpin->value());
    }
    if (m_wcsOutputPixelScaleSpin) {
        settings.setValue("wcs/outputPixelScale", m_wcsOutputPixelScaleSpin->value());
    }
}
