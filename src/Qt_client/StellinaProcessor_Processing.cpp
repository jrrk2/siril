#include "StellinaProcessor.h"
#include <QApplication>
#include <QDir>
#include <QProcess>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonDocument>
#include <QDateTime>
#include <QMessageBox>

void StellinaProcessor::startStellinaProcessing() {
    if (!validateProcessingInputs()) {
        return;
    }
    
    logMessage(QString("Starting %1 processing pipeline...")
                  .arg(m_processingModeCombo->currentText()), "blue");
    
    // Find images to process
    if (!findStellinaImages()) {
        logMessage("No valid Stellina images found in source directory.", "red");
        return;
    }
    
    // Change Siril working directory
    QString outputDir = getOutputDirectoryForCurrentStage();
    bool dirChanged = m_sirilClient->changeDirectory(outputDir);
    if (!dirChanged) {
        logMessage("Warning: Failed to change Siril working directory. Will use absolute paths.", "orange");
        logMessage(QString("Will save files to: %1").arg(outputDir), "blue");
    } else {
        logMessage(QString("Changed Siril working directory to: %1").arg(outputDir), "green");
    }
    
    // Initialize processing state
    m_processing = true;
    m_currentImageIndex = 0;
    m_processedCount = 0;
    m_errorCount = 0;
    m_skippedCount = 0;
    m_darkCalibratedCount = 0;
    m_registeredCount = 0;
    m_processingStartTime = QDateTime::currentMSecsSinceEpoch();
    
    // Clear previous results
    m_darkCalibratedFiles.clear();
    m_plateSolvedFiles.clear();
    m_registeredFiles.clear();
    m_finalStackedImage.clear();
    
    // Set up progress tracking
    m_progressBar->setMaximum(m_imagesToProcess.length());
    m_progressBar->setValue(0);
    
    // Determine processing stages based on mode
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        m_currentStage = STAGE_PLATE_SOLVING;
        break;
    case MODE_DARK_CALIBRATION:
        m_currentStage = STAGE_DARK_CALIBRATION;
        break;
    case MODE_ASTROMETRIC_STACKING:
        m_currentStage = STAGE_REGISTRATION;
        break;
    case MODE_FULL_PIPELINE:
        m_currentStage = STAGE_DARK_CALIBRATION;
        break;
    }
    
    updateProcessingStatus();
    updateUI();
    
    logMessage(QString("Found %1 images to process.").arg(m_imagesToProcess.length()), "green");
    
    // Start processing timer
    m_processingTimer->start();
}

void StellinaProcessor::processNextImage() {
    if (m_currentImageIndex >= m_imagesToProcess.length()) {
        // Current stage complete, move to next stage or finish
        switch (m_currentStage) {
        case STAGE_DARK_CALIBRATION:
            if (m_processingMode == MODE_DARK_CALIBRATION) {
                finishProcessing();
                return;
            } else {
                // Move to plate solving stage
                m_currentStage = STAGE_PLATE_SOLVING;
                m_currentImageIndex = 0;
                // Use dark-calibrated files if available
                if (!m_darkCalibratedFiles.isEmpty()) {
                    m_imagesToProcess = m_darkCalibratedFiles;
                }
                updateProcessingStatus();
                return;
            }
            break;
            
        case STAGE_PLATE_SOLVING:
            if (m_processingMode == MODE_BASIC_PLATESOLVE) {
                finishProcessing();
                return;
            } else {
                // Move to registration stage
                m_currentStage = STAGE_REGISTRATION;
                if (performAstrometricStacking()) {
                    m_currentStage = STAGE_COMPLETE;
                    finishProcessing();
                } else {
                    logMessage("Astrometric stacking failed", "red");
                    finishProcessing();
                }
                return;
            }
            break;
            
        case STAGE_REGISTRATION:
            if (performAstrometricStacking()) {
                m_currentStage = STAGE_COMPLETE;
                finishProcessing();
            } else {
                logMessage("Astrometric stacking failed", "red");
                finishProcessing();
            }
            return;
            
        default:
            finishProcessing();
            return;
        }
    }
    
    QString currentFile = m_imagesToProcess[m_currentImageIndex];
    
    logMessage(QString("Processing image %1 of %2: %3 (Stage: %4)")
                  .arg(m_currentImageIndex + 1)
                  .arg(m_imagesToProcess.length())
                  .arg(QFileInfo(currentFile).fileName())
                  .arg(getStageDescription()), "blue");
    
    bool success = false;
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        success = processImageDarkCalibration(currentFile);
        break;
    case STAGE_PLATE_SOLVING:
        success = processImagePlatesolving(currentFile);
        break;
    default:
        success = processImagePlatesolving(currentFile); // fallback
        break;
    }
    
    // Update counters
    if (success) {
        m_processedCount++;
    } else {
        m_errorCount++;
    }
    
    // Update progress
    m_currentImageIndex++;
    m_progressBar->setValue(m_currentImageIndex);
    updateProcessingStatus();
    
    // Update time estimate
    qint64 elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    double avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
    qint64 remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
    
    m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                    .arg(formatProcessingTime(remaining)));
}

bool StellinaProcessor::validateProcessingInputs() {
    if (!m_sirilClient->isConnected()) {
        QMessageBox::warning(this, "Connection Error", "Please connect to Siril first.");
        return false;
    }
    
    if (m_sourceDirectory.isEmpty()) {
        QMessageBox::warning(this, "Directory Error", "Please select the raw light frames directory.");
        return false;
    }
    
    // Check required directories based on processing mode
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        if (m_plateSolvedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the plate-solved output directory.");
            return false;
        }
        break;
        
    case MODE_DARK_CALIBRATION:
        if (m_calibratedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the calibrated lights output directory.");
            return false;
        }
        if (m_autoMatchDarks && m_darkFrames.isEmpty()) {
            int ret = QMessageBox::question(this, "No Dark Frames", 
                                           "Dark calibration is enabled but no dark frames were found. Continue without dark calibration?",
                                           QMessageBox::Yes | QMessageBox::No);
            if (ret == QMessageBox::No) {
                return false;
            }
        }
        break;
        
    case MODE_ASTROMETRIC_STACKING:
        if (m_stackedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the stacked images output directory.");
            return false;
        }
        break;
        
    case MODE_FULL_PIPELINE:
        if (m_calibratedDirectory.isEmpty() || m_plateSolvedDirectory.isEmpty() || m_stackedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", 
                                 "Full pipeline requires all output directories:\n"
                                 "- Calibrated Lights\n"
                                 "- Plate-Solved Lights\n" 
                                 "- Stacked Images");
            return false;
        }
        if (m_autoMatchDarks && m_darkFrames.isEmpty()) {
            int ret = QMessageBox::question(this, "No Dark Frames", 
                                           "Dark calibration is enabled but no dark frames were found. Continue without dark calibration?",
                                           QMessageBox::Yes | QMessageBox::No);
            if (ret == QMessageBox::No) {
                return false;
            }
        }
        break;
    }
    
    return true;
}

void StellinaProcessor::updateProcessingStatus() {
    QString stageText = getStageDescription();
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1").arg(
            m_processing ? "In Progress" : "Ready"));
        break;
    case STAGE_PLATE_SOLVING:
        // Keep existing registration/stacking status as ready for now
        break;
    case STAGE_REGISTRATION:
        m_registrationStatusLabel->setText("Registration: In Progress");
        break;
    case STAGE_STACKING:
        m_stackingStatusLabel->setText("Stacking: In Progress");
        break;
    case STAGE_COMPLETE:
        m_darkCalibrationStatusLabel->setText("Dark Calibration: Complete");
        m_registrationStatusLabel->setText("Registration: Complete");
        m_stackingStatusLabel->setText("Stacking: Complete");
        break;
    }
    
    if (m_processing) {
        m_progressLabel->setText(QString("Processing %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                    .arg(m_currentImageIndex + 1)
                                    .arg(m_imagesToProcess.length())
                                    .arg(stageText)
                                    .arg(m_processedCount)
                                    .arg(m_errorCount));
        
        if (m_darkCalibratedCount > 0) {
            m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1 completed").arg(m_darkCalibratedCount));
        }
    }
}

QString StellinaProcessor::getStageDescription() const {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION: return "Dark Calibration";
    case STAGE_PLATE_SOLVING: return "Plate Solving";
    case STAGE_REGISTRATION: return "Registration";
    case STAGE_STACKING: return "Stacking";
    case STAGE_COMPLETE: return "Complete";
    default: return "Unknown";
    }
}

void StellinaProcessor::finishProcessing() {
    m_processing = false;
    m_processingTimer->stop();
    
    qint64 totalTime = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    
    QString completionMessage = QString("Processing complete! Total time: %1\n\n")
                                   .arg(formatProcessingTime(totalTime));
    
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        completionMessage += QString("Plate solved: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_DARK_CALIBRATION:
        completionMessage += QString("Dark calibrated: %1, Skipped: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_skippedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_ASTROMETRIC_STACKING:
        completionMessage += QString("Images processed: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
        
    case MODE_FULL_PIPELINE:
        completionMessage += QString("Dark calibrated: %1, Plate solved: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
    }
    
    logMessage(completionMessage, "green");
    
    // Save processing report
    saveProcessingReport();
    
    updateUI();
    
    // Show completion dialog
    QMessageBox::information(this, "Processing Complete", completionMessage);
}

bool StellinaProcessor::checkSolveFieldInstalled() {
    QProcess process;
    process.start("/opt/homebrew/bin/solve-field", QStringList() << "--help");
    bool finished = process.waitForFinished(3000);
    return finished && process.exitCode() == 0;
}

QStringList StellinaProcessor::getAstrometryPaths() {
    QStringList astrometryPaths;
    
    // Common directories where astrometry.net tools are installed
    QStringList searchDirs = {
        "/opt/homebrew/bin",              // Homebrew Apple Silicon
        "/opt/homebrew/libexec/astrometry.net/bin",  // Homebrew support tools
        "/usr/local/bin",                 // Homebrew Intel / manual install
        "/usr/local/astrometry/bin",      // Custom install location
        "/usr/bin",                       // System package (Linux)
        "/opt/astrometry.net/bin",        // Alternative install
        "/usr/local/libexec/astrometry.net/bin"  // Support tools location
    };
    
    // Check which directories actually exist and contain astrometry tools
    for (const QString &dir : searchDirs) {
        QDir directory(dir);
        if (directory.exists()) {
            // Check for common astrometry.net tools
            QStringList requiredTools = {"solve-field", "fits2fits", "pnmfile", "image2xy"};
            bool hasTools = false;
            
            for (const QString &tool : requiredTools) {
                if (QFile::exists(directory.absoluteFilePath(tool))) {
                    hasTools = true;
                    break;
                }
            }
            
            if (hasTools) {
                astrometryPaths.append(dir);
                if (m_debugMode) {
                    logMessage(QString("Found astrometry tools in: %1").arg(dir), "gray");
                }
            }
        }
    }
    
    return astrometryPaths;
}

QProcessEnvironment StellinaProcessor::createSolveFieldEnvironment() {
    // Start with current environment
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    
    // Get astrometry.net tool directories
    QStringList astrometryPaths = getAstrometryPaths();
    
    // Get current PATH
    QString currentPath = env.value("PATH");
    QStringList pathList = currentPath.split(QDir::listSeparator(), Qt::SkipEmptyParts);
    
    // Add astrometry paths to the beginning (higher priority)
    for (const QString &astroPath : astrometryPaths) {
        if (!pathList.contains(astroPath)) {
            pathList.prepend(astroPath);
        }
    }
    
    // Set the enhanced PATH
    QString enhancedPath = pathList.join(QDir::listSeparator());
    env.insert("PATH", enhancedPath);
    
    if (m_debugMode) {
        logMessage(QString("Enhanced PATH for solve-field: %1").arg(enhancedPath), "gray");
    }
    
    // Add other useful environment variables for astrometry.net
    
    // Set data directory hints (if not already set)
    if (!env.contains("ASTROMETRY_NET_DATA_DIR")) {
        QStringList dataDirs = {
            "/opt/homebrew/share/astrometry",
            "/usr/local/share/astrometry",
            "/usr/share/astrometry",
            QDir::homePath() + "/.local/share/astrometry"
        };
        
        for (const QString &dataDir : dataDirs) {
            if (QDir(dataDir).exists()) {
                env.insert("ASTROMETRY_NET_DATA_DIR", dataDir);
                if (m_debugMode) {
                    logMessage(QString("Set ASTROMETRY_NET_DATA_DIR: %1").arg(dataDir), "gray");
                }
                break;
            }
        }
    }
    
    // Ensure solve-field can find its configuration
    if (!env.contains("SOLVE_FIELD_CONFIG")) {
        QStringList configPaths = {
            "/opt/homebrew/etc/astrometry.cfg",
            "/usr/local/etc/astrometry.cfg",
            "/etc/astrometry.cfg"
        };
        
        for (const QString &configPath : configPaths) {
            if (QFile::exists(configPath)) {
                env.insert("SOLVE_FIELD_CONFIG", configPath);
                if (m_debugMode) {
                    logMessage(QString("Found astrometry config: %1").arg(configPath), "gray");
                }
                break;
            }
        }
    }
    
    return env;
}

bool StellinaProcessor::runSolveField(const QString &fitsPath, const QString &outputPath, double ra, double dec) {
    logMessage(QString("Running solve-field with coordinate hints: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
    
    // Check if solve-field is available
    if (!checkSolveFieldInstalled()) {
        logMessage("solve-field not found - install astrometry.net", "red");
        return false;
    }
    
    QProcess solveProcess;
    QStringList arguments;
    
    // solve-field arguments with coordinate hints and binning
    arguments << fitsPath;
    arguments << "--overwrite";           // Overwrite existing files
    arguments << "--no-plots";            // Don't create plot files
    arguments << "--new-fits" << outputPath;  // Output solved FITS
    arguments << "--scale-units" << "arcsecperpix";
    arguments << "--scale-low" << "1.0";   // Original pixel scale range
    arguments << "--scale-high" << "1.5";  // solve-field will auto-adjust for downsampling
    arguments << "--ra" << QString::number(ra, 'f', 6);     // Use calculated RA
    arguments << "--dec" << QString::number(dec, 'f', 6);   // Use calculated Dec
    arguments << "--radius" << "1.0";      // Small search radius since we have good hints
    arguments << "--cpulimit" << "60";     // 1 minute timeout (faster with hints)
    arguments << "--no-verify";           // Skip verification step
    arguments << "--crpix-center";        // Set reference pixel at center
    arguments << "-z" << "2";             // Debayer the CFA image (CRITICAL for Stellina)
    arguments << "--downsample" << "2";   // 2x2 binning - solve-field auto-adjusts pixel scale
    
    if (m_debugMode) {
        logMessage(QString("solve-field command: solve-field %1").arg(arguments.join(" ")), "gray");
    }
    
    // Start solve-field process
    QProcessEnvironment enhancedEnv = createSolveFieldEnvironment();
    solveProcess.setProcessEnvironment(enhancedEnv);
    
    // Start solve-field process with absolute path and enhanced environment
    solveProcess.start("/opt/homebrew/bin/solve-field", arguments);
    
    if (!solveProcess.waitForStarted(5000)) {
        logMessage("Failed to start solve-field process", "red");
        return false;
    }
    
    // Monitor progress and keep UI responsive
    QElapsedTimer timer;
    timer.start();
    bool finished = false;
    
    while (!finished && timer.elapsed() < 60000) { // 1 minute timeout with hints
        finished = solveProcess.waitForFinished(1000);
        QApplication::processEvents(); // Keep UI responsive
        
        // Check if user wants to stop
        if (!m_processing) {
            solveProcess.kill();
            return false;
        }
        
        // Show some progress
        if (timer.elapsed() % 10000 == 0) { // Every 10 seconds
            logMessage(QString("solve-field running... (%1s elapsed)").arg(timer.elapsed() / 1000), "gray");
        }
    }
    
    if (!finished) {
        logMessage("solve-field timed out after 1 minute", "orange");
        solveProcess.kill();
        return false;
    }
    
    int exitCode = solveProcess.exitCode();
    QString output = solveProcess.readAllStandardOutput();
    QString errors = solveProcess.readAllStandardError();
    
    if (m_debugMode && !output.isEmpty()) {
        logMessage(QString("solve-field output: %1").arg(output.left(200)), "gray");
    }
    
    if (exitCode == 0) {
        // Check if output file was created
        if (QFile::exists(outputPath)) {
            logMessage("solve-field succeeded - image plate solved!", "green");
            
            // Extract coordinates from solve-field output
            QRegularExpression raRegex(R"(Field center: \(RA,Dec\) = \(([\d\.]+), ([\d\.]+)\))");
            QRegularExpressionMatch match = raRegex.match(output);
            if (match.hasMatch()) {
                double ra = match.captured(1).toDouble();
                double dec = match.captured(2).toDouble();
                logMessage(QString("solve-field found coordinates: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
            }
            
            return true;
        } else {
            logMessage(QString("solve-field completed but output file not found: %1").arg(outputPath), "red");
            return false;
        }
    } else {
        logMessage(QString("solve-field failed with exit code %1").arg(exitCode), "red");
        if (!errors.isEmpty()) {
            logMessage(QString("solve-field errors: %1").arg(errors.left(200)), "red");
        }
        return false;
    }
}
