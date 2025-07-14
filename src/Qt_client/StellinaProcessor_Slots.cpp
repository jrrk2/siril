#include "StellinaProcessor.h"

void StellinaProcessor::connectSignals() {
    // UI signals
    connect(m_selectSourceButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectSourceDirectory);
    connect(m_selectDarkButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectDarkDirectory);
    connect(m_selectCalibratedButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectCalibratedDirectory);
    connect(m_selectPlateSolvedButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectPlateSolvedDirectory);
    connect(m_selectStackedButton, &QPushButton::clicked,
            this, &StellinaProcessor::onSelectStackedDirectory);
    connect(m_startButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStartProcessing);
    connect(m_stopButton, &QPushButton::clicked,
            this, &StellinaProcessor::onStopProcessing);
    connect(m_testConnectionButton, &QPushButton::clicked,
            this, &StellinaProcessor::onTestConnection);
    connect(m_clearLogButton, &QPushButton::clicked,
            this, &StellinaProcessor::onClearLog);
    connect(m_refreshDarkButton, &QPushButton::clicked,
            this, &StellinaProcessor::onRefreshDarkFrames);
    connect(m_previewStackButton, &QPushButton::clicked,
            this, &StellinaProcessor::onPreviewStacking);
    
    // Processing mode combo
    connect(m_processingModeCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &StellinaProcessor::onProcessingModeChanged);
    
    // Siril client signals
    connect(m_sirilClient, &SirilClient::connected,
            this, &StellinaProcessor::onSirilConnected);
    connect(m_sirilClient, &SirilClient::disconnected,
            this, &StellinaProcessor::onSirilDisconnected);
    connect(m_sirilClient, &SirilClient::commandExecuted,
            this, &StellinaProcessor::onSirilCommandExecuted);
    connect(m_sirilClient, &SirilClient::errorOccurred,
            this, &StellinaProcessor::onSirilError);
    
    // Processing timer
    connect(m_processingTimer, &QTimer::timeout,
            this, &StellinaProcessor::onProcessingTimer);
    
    // Settings updates
    connect(m_qualityFilterCheck, &QCheckBox::toggled, [this](bool checked) {
        m_qualityFilter = checked;
    });
    connect(m_debugModeCheck, &QCheckBox::toggled, [this](bool checked) {
        m_debugMode = checked;
    });
    connect(m_focalLengthSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_focalLength = value;
    });
    connect(m_pixelSizeSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_pixelSize = value;
    });
    connect(m_observerLocationEdit, &QLineEdit::textChanged, [this](const QString &text) {
        m_observerLocation = text;
    });
    
    // Dark calibration settings
    connect(m_autoMatchDarksCheck, &QCheckBox::toggled, [this](bool checked) {
        m_autoMatchDarks = checked;
    });
    connect(m_temperatureToleranceSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [this](int value) {
        m_temperatureTolerance = value;
    });
    connect(m_exposureToleranceSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [this](int value) {
        m_exposureTolerance = value;
    });
    
    // Stacking settings
    connect(m_stackingMethodCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.method = text;
        onStackingParametersChanged();
    });
    connect(m_rejectionMethodCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.rejection = text;
        onStackingParametersChanged();
    });
    connect(m_rejectionLowSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.rejectionLow = value;
        onStackingParametersChanged();
    });
    connect(m_rejectionHighSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.rejectionHigh = value;
        onStackingParametersChanged();
    });
    connect(m_normalizeCheck, &QCheckBox::toggled, [this](bool checked) {
        m_stackingParams.normalizeImages = checked;
        onStackingParametersChanged();
    });
    connect(m_drizzleCheck, &QCheckBox::toggled, [this](bool checked) {
        m_stackingParams.applyDrizzle = checked;
        m_drizzleScaleSpin->setEnabled(checked);
        onStackingParametersChanged();
    });
    connect(m_drizzleScaleSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [this](double value) {
        m_stackingParams.drizzleScale = value;
        onStackingParametersChanged();
    });
    connect(m_outputFormatCombo, &QComboBox::currentTextChanged, [this](const QString &text) {
        m_stackingParams.outputFormat = text;
        onStackingParametersChanged();
    });

    // Mount tilt correction signals
    connect(m_enableTiltCorrectionCheck, &QCheckBox::toggled, [this](bool checked) {
	m_mountTilt.enableCorrection = checked;
	updateTiltUI();
	saveMountTiltToSettings();
    });

    connect(m_northTiltSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), 
	    [this](double value) {
	m_mountTilt.northTilt = value;
	saveMountTiltToSettings();
    });

    connect(m_eastTiltSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), 
	    [this](double value) {
	m_mountTilt.eastTilt = value;
	saveMountTiltToSettings();
    });

    connect(m_calibrateTiltButton, &QPushButton::clicked, 
	    this, &StellinaProcessor::calibrateMountTilt);

    connect(m_testTiltButton, &QPushButton::clicked, 
	    this, &StellinaProcessor::testMountTiltCorrection);

    connect(m_enableDriftCorrectionCheck, &QCheckBox::toggled, [this](bool checked) {
	m_mountTilt.enableDriftCorrection = checked;
	updateTiltUI();
	saveMountTiltToSettings();
    });

    connect(m_driftRASpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), 
	    [this](double value) {
	m_mountTilt.driftRA = value;
	saveMountTiltToSettings();
    });

    connect(m_driftDecSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), 
	    [this](double value) {
	m_mountTilt.driftDec = value;
	saveMountTiltToSettings();
    });
    
}

// UI Event Handlers
void StellinaProcessor::onProcessingModeChanged() {
    m_processingMode = static_cast<ProcessingMode>(m_processingModeCombo->currentData().toInt());
    
    QString description;
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        description = "Basic plate solving with coordinate annotation - original functionality";
        break;
    case MODE_DARK_CALIBRATION:
        description = "Dark frame calibration only - removes thermal noise from light frames";
        break;
    case MODE_ASTROMETRIC_STACKING:
        description = "Astrometric registration and stacking - aligns and combines images";
        break;
    case MODE_FULL_PIPELINE:
        description = "Complete processing pipeline - dark calibration, plate solving, registration, and stacking";
        break;
    }
    
    m_modeDescription->setText(description);
    updateUI();
}

void StellinaProcessor::onSelectSourceDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Raw Light Frames Directory", m_sourceDirectory);
    if (!dir.isEmpty()) {
        m_sourceDirectory = dir;
        m_sourceDirectoryEdit->setText(dir);
        
        QDir testDir(dir);
        logMessage(QString("Selected raw lights directory: %1").arg(testDir.absolutePath()), "blue");
        logMessage(QString("Directory exists: %1").arg(testDir.exists() ? "Yes" : "No"), "blue");
        
        QStringList allFiles = testDir.entryList(QDir::Files);
        logMessage(QString("Total files in directory: %1").arg(allFiles.count()), "blue");
        
        updateUI();
    }
}

void StellinaProcessor::onSelectDarkDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Dark Frames Directory", m_darkDirectory);
    if (!dir.isEmpty()) {
        m_darkDirectory = dir;
        m_darkDirectoryEdit->setText(dir);
        scanDarkFrames();
        updateUI();
    }
}

void StellinaProcessor::onSelectCalibratedDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Dark-Calibrated Lights Directory", m_calibratedDirectory);
    if (!dir.isEmpty()) {
        m_calibratedDirectory = dir;
        m_calibratedDirectoryEdit->setText(dir);
        
        QDir testDir(dir);
        logMessage(QString("Selected calibrated lights directory: %1").arg(testDir.absolutePath()), "blue");
        
        updateUI();
    }
}

void StellinaProcessor::onSelectPlateSolvedDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Plate-Solved Lights Directory", m_plateSolvedDirectory);
    if (!dir.isEmpty()) {
        m_plateSolvedDirectory = dir;
        m_plateSolvedDirectoryEdit->setText(dir);
        
        QDir testDir(dir);
        logMessage(QString("Selected plate-solved lights directory: %1").arg(testDir.absolutePath()), "blue");
        
        updateUI();
    }
}

void StellinaProcessor::onSelectStackedDirectory() {
    QString dir = QFileDialog::getExistingDirectory(this, "Select Stacked Images Directory", m_stackedDirectory);
    if (!dir.isEmpty()) {
        m_stackedDirectory = dir;
        m_stackedDirectoryEdit->setText(dir);
        
        QDir testDir(dir);
        logMessage(QString("Selected stacked images directory: %1").arg(testDir.absolutePath()), "blue");
        
        updateUI();
    }
}

void StellinaProcessor::onStartProcessing() {
    startStellinaProcessing();
}

void StellinaProcessor::onStopProcessing() {
    if (m_processing) {
        m_processingTimer->stop();
        m_processing = false;
        logMessage("Processing stopped by user.", "orange");
        finishProcessing();
    }
}

void StellinaProcessor::onTestConnection() {
    if (m_sirilClient->isConnected()) {
        m_sirilClient->disconnectFromSiril();
    } else {
        logMessage("Attempting to connect to Siril...", "blue");
        if (m_sirilClient->connectToSiril()) {
            logMessage("Successfully connected to Siril!", "green");
        } else {
            logMessage(QString("Failed to connect: %1").arg(m_sirilClient->lastError()), "red");
        }
    }
    updateUI();
}

void StellinaProcessor::onClearLog() {
    m_logTextEdit->clear();
    logMessage("Log cleared.", "gray");
}

void StellinaProcessor::onRefreshDarkFrames() {
    if (!m_darkDirectory.isEmpty()) {
        scanDarkFrames();
    }
}

void StellinaProcessor::onPreviewStacking() {
    if (m_imagesToProcess.isEmpty()) {
        if (!findStellinaImages()) {
            QMessageBox::warning(this, "No Images", "No valid Stellina images found for preview.");
            return;
        }
    }
    
    if (m_imagesToProcess.size() < 3) {
        QMessageBox::warning(this, "Insufficient Images", "Need at least 3 images for preview stacking.");
        return;
    }
    
    if (!m_sirilClient->isConnected()) {
        QMessageBox::warning(this, "Connection Error", "Please connect to Siril first.");
        return;
    }
    
    // Take first 5 images for preview
    QStringList previewImages = m_imagesToProcess.mid(0, qMin(5, m_imagesToProcess.size()));
    
    logMessage(QString("Creating preview stack with %1 images...").arg(previewImages.size()), "blue");
    
    // Create a temporary sequence for preview
    QString previewSequence = "preview_stack";
    if (createSequence(previewImages, previewSequence)) {
        if (performGlobalRegistration(previewSequence)) {
            if (performStacking(previewSequence, m_stackingParams)) {
                logMessage("Preview stack created successfully! Check output directory.", "green");
            } else {
                logMessage("Preview stacking failed.", "red");
            }
        } else {
            logMessage("Preview registration failed.", "red");
        }
    } else {
        logMessage("Failed to create preview sequence.", "red");
    }
}

void StellinaProcessor::onStackingParametersChanged() {
    // This function is called when stacking parameters change
    // Update the internal stacking parameters structure
    m_stackingParams.method = m_stackingMethodCombo->currentText();
    m_stackingParams.rejection = m_rejectionMethodCombo->currentText();
    m_stackingParams.rejectionLow = m_rejectionLowSpin->value();
    m_stackingParams.rejectionHigh = m_rejectionHighSpin->value();
    m_stackingParams.normalizeImages = m_normalizeCheck->isChecked();
    m_stackingParams.applyDrizzle = m_drizzleCheck->isChecked();
    m_stackingParams.drizzleScale = m_drizzleScaleSpin->value();
    m_stackingParams.outputFormat = m_outputFormatCombo->currentText();
    
    if (m_debugMode) {
        logMessage(QString("Stacking parameters updated: method=%1, rejection=%2")
                      .arg(m_stackingParams.method)
                      .arg(m_stackingParams.rejection), "gray");
    }
}

// Siril Event Handlers
void StellinaProcessor::onSirilConnected() {
    logMessage("Connected to Siril successfully!", "green");
    
    QString workDir = m_sirilClient->getWorkingDirectory();
    if (!workDir.isEmpty()) {
        logMessage(QString("Siril working directory: %1").arg(workDir), "blue");
    }
    
    updateUI();
}

void StellinaProcessor::onSirilDisconnected() {
    logMessage("Disconnected from Siril.", "orange");
    updateUI();
}

void StellinaProcessor::onSirilCommandExecuted(const QString &command, bool success) {
    if (m_debugMode) {
        QString status = success ? "SUCCESS" : "FAILED";
        logMessage(QString("Command [%1]: %2").arg(status).arg(command), success ? "gray" : "red");
    }
}

void StellinaProcessor::onSirilError(const QString &error) {
    logMessage(QString("Siril error: %1").arg(error), "red");
}

void StellinaProcessor::onProcessingTimer() {
    processNextImage();
}

// UI Update Functions
void StellinaProcessor::updateUI() {
    bool connected = m_sirilClient->isConnected();
    bool hasSource = !m_sourceDirectory.isEmpty();
    bool hasRequiredOutputs = false;
    
    // Check if required output directories are set based on processing mode
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        hasRequiredOutputs = !m_plateSolvedDirectory.isEmpty();
        break;
    case MODE_DARK_CALIBRATION:
        hasRequiredOutputs = !m_calibratedDirectory.isEmpty();
        break;
    case MODE_ASTROMETRIC_STACKING:
        hasRequiredOutputs = !m_stackedDirectory.isEmpty();
        break;
    case MODE_FULL_PIPELINE:
        hasRequiredOutputs = !m_calibratedDirectory.isEmpty() && 
                           !m_plateSolvedDirectory.isEmpty() && 
                           !m_stackedDirectory.isEmpty();
        break;
    }
    
    bool canProcess = connected && hasSource && hasRequiredOutputs && !m_processing;
    
    m_startButton->setEnabled(canProcess);
    m_stopButton->setEnabled(m_processing);
    m_previewStackButton->setEnabled(canProcess && (m_processingMode == MODE_ASTROMETRIC_STACKING || m_processingMode == MODE_FULL_PIPELINE));
    
    updateConnectionStatus();
}

// Update UI to show connection type
void StellinaProcessor::updateConnectionStatus() {
    if (m_sirilClient->isConnected()) {
        QString connectionType = m_sirilClient->isUsingCLI() ? "siril-cli" : "Siril GUI";
        m_connectionStatus->setText(QString("Connected to %1").arg(connectionType));
        m_connectionStatus->setStyleSheet("color: green; font-weight: bold;");
        m_testConnectionButton->setText("Disconnect");
    } else {
        m_connectionStatus->setText("Not connected");
        m_connectionStatus->setStyleSheet("color: red; font-weight: bold;");
        m_testConnectionButton->setText("Test Connection");
    }
}

void StellinaProcessor::logMessage(const QString &message, const QString &color) {
    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    QString formattedMessage = QString("<span style='color: %1'>[%2] %3</span>")
                                  .arg(color)
                                  .arg(timestamp)
                                  .arg(message);
    
    m_logTextEdit->append(formattedMessage);
    
    QTextCursor cursor = m_logTextEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    m_logTextEdit->setTextCursor(cursor);
    
    if (color == "red" || color == "green" || color == "blue") {
        m_statusLabel->setText(message);
    }
}

void StellinaProcessor::updateTiltUI() {
    if (!m_enableTiltCorrectionCheck) return;
    
    bool tiltEnabled = m_mountTilt.enableCorrection;
    bool driftEnabled = m_mountTilt.enableDriftCorrection;
    
    // Update tilt UI
    m_northTiltSpin->setEnabled(tiltEnabled);
    m_eastTiltSpin->setEnabled(tiltEnabled);
    m_calibrateTiltButton->setEnabled(tiltEnabled);
    m_testTiltButton->setEnabled(tiltEnabled);
    
    // Update drift UI
    if (m_enableDriftCorrectionCheck) {
        m_enableDriftCorrectionCheck->setEnabled(tiltEnabled);
        m_driftRASpin->setEnabled(tiltEnabled && driftEnabled);
        m_driftDecSpin->setEnabled(tiltEnabled && driftEnabled);
        
        m_enableDriftCorrectionCheck->setChecked(driftEnabled);
        m_driftRASpin->setValue(m_mountTilt.driftRA);
        m_driftDecSpin->setValue(m_mountTilt.driftDec);
    }
    
    // Update status labels
    if (tiltEnabled) {
        m_tiltStatusLabel->setText(QString("Tilt correction enabled: θ_N=%1°, θ_E=%2°")
                                      .arg(m_mountTilt.northTilt, 0, 'f', 4)
                                      .arg(m_mountTilt.eastTilt, 0, 'f', 4));
        m_tiltStatusLabel->setStyleSheet("color: green; font-weight: bold;");
        
        if (m_driftStatusLabel) {
            if (driftEnabled) {
                m_driftStatusLabel->setText(QString("Drift correction: %1°/h RA, %2°/h Dec")
                                              .arg(m_mountTilt.driftRA, 0, 'f', 3)
                                              .arg(m_mountTilt.driftDec, 0, 'f', 3));
                m_driftStatusLabel->setStyleSheet("color: green; font-weight: bold;");
            } else {
                m_driftStatusLabel->setText("Drift correction disabled");
                m_driftStatusLabel->setStyleSheet("color: gray; font-style: italic;");
            }
        }
    } else {
        m_tiltStatusLabel->setText("Tilt correction disabled");
        m_tiltStatusLabel->setStyleSheet("color: gray; font-style: italic;");
        if (m_driftStatusLabel) {
            m_driftStatusLabel->setText("Drift correction disabled");
            m_driftStatusLabel->setStyleSheet("color: gray; font-style: italic;");
        }
    }
}
