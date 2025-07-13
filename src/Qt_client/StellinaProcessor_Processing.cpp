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
    process.start("solve-field", QStringList() << "--help");
    bool finished = process.waitForFinished(3000);
    return finished && process.exitCode() == 0;
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
    solveProcess.start("solve-field", arguments);
    
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

// Update processImagePlatesolving to use fallback
bool StellinaProcessor::processImagePlatesolving(const QString &fitsPath) {
    m_currentTaskLabel->setText("Plate solving...");
    
    // IMMEDIATE CONNECTION CHECK
    if (!m_sirilClient->isConnected()) {
        logMessage("ERROR: Not connected to Siril - cannot process image", "red");
        return false;
    }
    
    QString baseName = QFileInfo(fitsPath).baseName();
    
    // Find corresponding JSON file
    QDir sourceDir(QFileInfo(fitsPath).dir());
    QString jsonPath;
    
    QStringList jsonCandidates = {
        baseName + ".json",
        baseName + ".JSON",
        baseName + "-stacking.json",
        baseName + "-stacking.JSON",
        QFileInfo(fitsPath).completeBaseName() + ".json",
        QFileInfo(fitsPath).completeBaseName() + ".JSON",
        QFileInfo(fitsPath).completeBaseName() + "-stacking.json",
        QFileInfo(fitsPath).completeBaseName() + "-stacking.JSON"
    };
    
    bool jsonFound = false;
    for (const QString &candidate : jsonCandidates) {
        QString candidatePath = sourceDir.absoluteFilePath(candidate);
        if (QFile::exists(candidatePath)) {
            jsonPath = candidatePath;
            jsonFound = true;
            break;
        }
    }
    
    if (!jsonFound) {
        logMessage(QString("No JSON file found for %1").arg(QFileInfo(fitsPath).fileName()), "red");
        return false;
    }
    
    // Load JSON metadata
    QJsonObject json = loadStellinaJson(jsonPath);
    if (json.isEmpty()) {
        logMessage(QString("Failed to load JSON metadata"), "red");
        return false;
    }
    
    // Extract coordinates
    double alt, az;
    if (!extractCoordinates(json, alt, az)) {
        logMessage("Failed to extract coordinates from JSON", "red");
        return false;
    }
    
    // Convert Alt/Az to RA/Dec
    double ra, dec;
    QString dateObs = extractDateObs(fitsPath);
    if (!convertAltAzToRaDec(alt, az, dateObs, ra, dec)) {
        logMessage("Failed to convert coordinates", "red");
        return false;
    }
    
    QString outputName = QString("processed_%1").arg(baseName);
    QString outputPath = QDir(m_plateSolvedDirectory).absoluteFilePath(outputName);
    
    // First try: Siril with coordinate hints and binning
    bool platesolveSuccess = false;
    bool usedSiril = false;
    
    if (m_sirilClient->isConnected()) {
        logMessage(QString("Trying Siril plate solve with hints: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
        
        // Load image in Siril
        if (m_sirilClient->loadImage(fitsPath)) {
            // Apply 2x2 CFA-aware binning before plate solving to improve SNR
            logMessage("Applying 2x2 CFA-aware binning to improve plate solving", "blue");
            
            // Try CFA-specific binning commands first
            bool binningSuccess = false;
            
            // Option 1: Try CFA binning if available
            if (m_sirilClient->sendSirilCommand("bin 2 2 -cfa")) {
                logMessage("CFA-aware binning applied successfully", "green");
                binningSuccess = true;
            }
            // Option 2: Try generic binning 
            else if (m_sirilClient->sendSirilCommand("bin 2 2")) {
                logMessage("Generic 2x2 binning applied successfully", "green");
                binningSuccess = true;
            }
            // Option 3: Try extract_green to get single channel before binning
            else if (m_sirilClient->sendSirilCommand("extract_green") && 
                     m_sirilClient->sendSirilCommand("bin 2 2")) {
                logMessage("Extracted green channel and applied binning", "green");
                binningSuccess = true;
            }
            
            if (!binningSuccess) {
                logMessage("Binning failed, continuing with original CFA resolution", "orange");
                logMessage("WARNING: Plate solving CFA images without binning may have lower success rate", "orange");
            }
            
            // Perform plate solving with hints (adjust pixel size if binning succeeded)
            double effectivePixelSize = binningSuccess ? (m_pixelSize * 2.0) : m_pixelSize;
            platesolveSuccess = m_sirilClient->platesolve(ra, dec, m_focalLength, effectivePixelSize, true);
            
            if (platesolveSuccess) {
                logMessage("Siril plate solve succeeded with coordinate hints", "green");
                usedSiril = true;
                
                // Save processed image
                if (!m_sirilClient->saveImage(outputPath)) {
                    logMessage("Failed to save Siril processed image", "red");
                    m_sirilClient->closeImage();
                    return false;
                }
            } else {
                logMessage(QString("Siril plate solve failed: %1").arg(m_sirilClient->lastError()), "orange");
                logMessage("Will try solve-field fallback", "blue");
            }
            
            // Close image in Siril
            m_sirilClient->closeImage();
        } else {
            logMessage("Failed to load image in Siril, will try solve-field", "orange");
        }
    }
    
    // Second try: solve-field with coordinate hints
    if (!platesolveSuccess) {
        m_currentTaskLabel->setText("Plate solving with solve-field...");
        logMessage("Falling back to solve-field with coordinate hints", "blue");
        
        platesolveSuccess = runSolveField(fitsPath, outputPath, ra, dec);
        
        if (platesolveSuccess) {
            logMessage("solve-field succeeded with coordinate hints!", "green");
        } else {
            logMessage("Both Siril and solve-field failed to plate solve image", "red");
            return false;
        }
    }
    
    // Add to plate solved files list for potential stacking
    if (platesolveSuccess) {
        m_plateSolvedFiles.append(outputPath);
        
        QString method = usedSiril ? "Siril (with hints)" : "solve-field (with hints)";
        logMessage(QString("Image successfully plate solved using %1").arg(method), "green");
        return true;
    }
    
    return false;
}
