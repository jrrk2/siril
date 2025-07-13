#include "SirilClient.h"
#include <QDir>
#include <QFileInfo>
#include <QStandardPaths>
#include <QDebug>
#include <QDataStream>
#include <QLocalSocket>
#include <QTimer>
#include <QThread>
#include <QProcess>
#include <QApplication>
#include <QElapsedTimer>

SirilClient::SirilClient(QObject *parent)
    : QObject(parent)
    , m_socket(new QLocalSocket(this))
    , m_connectionTimer(new QTimer(this))
    , m_connected(false)
    , m_useCLI(false)  // Add this
    , m_commandPipe(nullptr)  // Add this
    , m_responsePipe(nullptr)  // Add this
{
    // Set up connection timer
    m_connectionTimer->setSingleShot(true);
    m_connectionTimer->setInterval(30000); // 5 second timeout
    
    // Connect signals
    connect(m_socket, &QLocalSocket::connected, 
            this, &SirilClient::onSocketConnected);
    connect(m_socket, &QLocalSocket::disconnected, 
            this, &SirilClient::onSocketDisconnected);
    connect(m_socket, &QLocalSocket::errorOccurred,
            this, &SirilClient::onSocketError);
    connect(m_connectionTimer, &QTimer::timeout,
            this, &SirilClient::onConnectionTimeout);
}

SirilClient::~SirilClient() {
    disconnectFromSiril();
}

bool SirilClient::connectToCLI() {
    m_commandPipePath = "/tmp/siril_command.in";
    m_responsePipePath = "/tmp/siril_command.out";
    
    // Check if siril-cli is running by checking for pipes
    if (!QFile::exists(m_commandPipePath) || !QFile::exists(m_responsePipePath)) {
        setError("siril-cli not running. Start with: siril-cli");
        return false;
    }
    
    // Open command pipe for writing
    m_commandPipe = new QFile(m_commandPipePath, this);
    if (!m_commandPipe->open(QIODevice::WriteOnly)) {
        setError(QString("Failed to open command pipe: %1").arg(m_commandPipePath));
        return false;
    }
    
    // Open response pipe for reading
    m_responsePipe = new QFile(m_responsePipePath, this);
    if (!m_responsePipe->open(QIODevice::ReadOnly)) {
        setError(QString("Failed to open response pipe: %1").arg(m_responsePipePath));
        m_commandPipe->close();
        return false;
    }
    
    m_connected = true;
    m_useCLI = true;
    
    qDebug() << "Connected to siril-cli via named pipes";
    return true;
}


bool SirilClient::waitForPlatesolvingResult() {
    if (!m_responsePipe) {
        return false;
    }
    
    QElapsedTimer timer;
    timer.start();
    const int PLATESOLVE_TIMEOUT = 30000; // 30 seconds
    
    QTextStream in(m_responsePipe);
    QString line;
    
    // Read response lines until we get a result
    while (timer.elapsed() < PLATESOLVE_TIMEOUT) {
        if (m_responsePipe->bytesAvailable() > 0 || m_responsePipe->waitForReadyRead(100)) {
            while (in.readLineInto(&line)) {
                qDebug() << "siril-cli response:" << line;
                
                // Check for success indicators
                if (line.contains("log: Siril solve succeeded")) {
                    qDebug() << "Plate solving succeeded!";
                    return true;
                }
                
                // Check for failure indicators
                if (line.contains("log: Plate solving failed") || 
                    line.contains("could not be aligned with the reference stars") ||
                    line.contains("Transformation matrix is invalid")) {
                    qDebug() << "Plate solving failed:" << line;
                    setError(QString("Plate solving failed: %1").arg(line));
                    return false;
                }
                
                // Check for other error conditions
                if (line.contains("Error") || line.contains("error")) {
                    qDebug() << "Possible error:" << line;
                    // Don't fail immediately - might be a warning
                }
            }
        }
        
        QApplication::processEvents(); // Keep UI responsive
        
        // Check if user wants to stop
        if (!m_connected) {
            return false;
        }
    }
    
    // Timeout - assume failure
    setError("Timeout waiting for plate solving result");
    return false;
}

bool SirilClient::sendCLICommand(const QString &command) {
    if (!m_connected || !m_commandPipe) {
        setError("Not connected to siril-cli");
        return false;
    }
    
    qDebug() << "Sending CLI command:" << command;
    
    // Write command to pipe
    QTextStream out(m_commandPipe);
    out << command << "\n";
    out.flush();
    
    if (!m_commandPipe->flush()) {
        setError("Failed to flush command to pipe");
        return false;
    }
    
    // For plate solving commands, we need to check if it actually succeeded
    if (command.startsWith("platesolve")) {
        return waitForPlatesolvingResult();
    }
    
    // For other commands, assume success if write worked
    return true;
}


// Update sendSirilCommand to use CLI when available
bool SirilClient::sendSirilCommand(const QString &command) {
    if (m_useCLI) {
        return sendCLICommand(command);
    } else {
        // Use existing socket protocol
        QByteArray payload = command.toUtf8();
        auto result = sendCommand(CMD_SEND_COMMAND, payload);
        bool success = (result.first == STATUS_OK);
        emit commandExecuted(command, success);
        return success;
    }
}

// Update disconnectFromSiril
void SirilClient::disconnectFromSiril() {
    if (m_useCLI) {
        if (m_commandPipe) {
            m_commandPipe->close();
            m_commandPipe = nullptr;
        }
        if (m_responsePipe) {
            m_responsePipe->close();
            m_responsePipe = nullptr;
        }
    } else {
        // Existing socket disconnection
        if (m_socket->state() != QLocalSocket::UnconnectedState) {
            m_socket->disconnectFromServer();
            m_socket->waitForDisconnected(3000);
        }
    }
    m_connected = false;
}
// Update the main connectToSiril function
bool SirilClient::connectToSiril() {
    if (m_connected) {
        return true;
    }
    
    // Try CLI first (more stable)
    if (connectToCLI()) {
        emit connected();
        return true;
    }
    
    // Fall back to socket connection
    QString socketPath = findSirilSocket();
    if (socketPath.isEmpty()) {
        setError("Could not find Siril socket and siril-cli not available.\n"
                "Start Siril GUI or run: siril-cli");
        return false;
    }
    
    // Rest of existing socket connection code...
    m_useCLI = false;
    qDebug() << "Connecting to Siril at:" << socketPath;
    
    // Verify the socket file exists and is accessible
    QFileInfo socketInfo(socketPath);
    if (!socketInfo.exists()) {
        setError(QString("Socket file does not exist: %1").arg(socketPath));
        return false;
    }
    
    if (!socketInfo.isReadable()) {
        setError(QString("Socket file is not readable: %1").arg(socketPath));
        return false;
    }
    
    qDebug() << "Socket file verified, attempting connection...";
    
    m_connectionTimer->start();
    m_socket->connectToServer(socketPath);
    
    return m_socket->waitForConnected(5000);
}

QString SirilClient::findSirilSocket() {
    qDebug() << "Looking for running Siril process...";
    
#ifdef Q_OS_WIN
    // Windows approach using tasklist
    QProcess process;
    process.start("tasklist", QStringList() << "/FO" << "CSV");
    
    if (!process.waitForFinished(3000)) {
        qDebug() << "Failed to execute tasklist command";
        return QString();
    }
    
    QString output = process.readAllStandardOutput();
    QStringList lines = output.split('\n', Qt::SkipEmptyParts);
    
    for (const QString &line : lines) {
        if (line.contains("siril", Qt::CaseInsensitive)) {
            // Parse CSV format: "Image Name","PID","Session Name","Session#","Mem Usage"
            QStringList parts = line.split(',');
            if (parts.size() >= 2) {
                QString pidStr = parts[1].remove('"').trimmed();
                bool ok;
                qint64 pid = pidStr.toLongLong(&ok);
                if (ok) {
                    QString socketPath = QString("\\\\.\\pipe\\siril_%1").arg(pid);
                    qDebug() << "Checking Windows pipe:" << socketPath;
                    return socketPath; // On Windows, assume the pipe exists if process exists
                }
            }
        }
    }
#else
    // Unix/Linux/macOS approach using ps -ef (more reliable than ps -axo)
    QProcess process;
    process.start("ps", QStringList() << "-ef");
    
    if (!process.waitForFinished(3000)) {
        qDebug() << "Failed to execute ps command";
        return QString();
    }
    
    QString output = process.readAllStandardOutput();
    QStringList lines = output.split('\n', Qt::SkipEmptyParts);
    
    QList<qint64> sirilPids;
    
    for (const QString &line : lines) {
        // Skip the header line
        if (line.contains("UID") && line.contains("PID")) {
            continue;
        }
        
        // ps -ef format: UID PID PPID C STIME TTY TIME CMD
        QStringList parts = line.simplified().split(' ');
        if (parts.size() >= 8) {
            QString pid = parts[1];
            QString command = parts[7]; // The CMD field
            
            qDebug() << "Checking process:" << pid << command;
            
            // Check if this is a Siril process (check the actual command, not just the name)
            if (command.contains("siril", Qt::CaseInsensitive) && 
                !command.contains("StellinaProcessor")) { // Exclude our own process
                bool ok;
                qint64 pidNum = pid.toLongLong(&ok);
                if (ok) {
                    sirilPids.append(pidNum);
                    qDebug() << "Found Siril process with PID:" << pidNum << "Command:" << command;
                }
            }
        }
    }
    
    if (sirilPids.isEmpty()) {
        qDebug() << "No running Siril processes found";
        return QString();
    }
    
    // Try each PID to find a valid socket (prefer the most recent)
    std::sort(sirilPids.begin(), sirilPids.end(), std::greater<qint64>());
    
    for (qint64 pid : sirilPids) {
        QString socketPath = QString("/tmp/siril_%1.sock").arg(pid);
        qDebug() << "Checking socket path:" << socketPath;
        
        if (QFile::exists(socketPath)) {
            qDebug() << "Found valid socket:" << socketPath;
            return socketPath;
        }
    }
#endif
    
    qDebug() << "No valid sockets found for running Siril processes";
    return QString();
}

QByteArray SirilClient::createCommandHeader(int commandId, int payloadLength) {
    QByteArray header;
    QDataStream stream(&header, QIODevice::WriteOnly);
    stream.setByteOrder(QDataStream::BigEndian);
    
    stream << static_cast<quint8>(commandId);
    stream << static_cast<quint32>(payloadLength);
    
    return header;
}

QPair<int, QByteArray> SirilClient::sendCommand(int commandId, const QByteArray &payload) {
    if (!m_connected) {
        setError("Not connected to Siril");
        return qMakePair(STATUS_ERROR, QByteArray());
    }
    
    // Create and send command
    QByteArray header = createCommandHeader(commandId, payload.length());
    QByteArray fullCommand = header + payload;
    
    qint64 written = m_socket->write(fullCommand);
    if (written != fullCommand.length()) {
        setError("Failed to write complete command to socket");
        return qMakePair(STATUS_ERROR, QByteArray());
    }
    
    if (!m_socket->waitForBytesWritten(30000)) {
        setError("Timeout writing to socket");
        return qMakePair(STATUS_ERROR, QByteArray());
    }
    
    // Add a small delay with UI processing after sending
    QApplication::processEvents();
    
    // Read response
    return readResponse();
}

QPair<int, QByteArray> SirilClient::readResponse() {
    // Read response header (5 bytes: status + length) with UI polling
    QElapsedTimer timer;
    timer.start();
    const int timeoutMs = 120000; // 2 minutes
    
    while (!m_socket->bytesAvailable() && timer.elapsed() < timeoutMs) {
        if (!m_socket->waitForReadyRead(100)) { // Short 100ms waits
            QApplication::processEvents(); // Keep UI responsive
            
            // Check if socket has any data now
            if (m_socket->bytesAvailable() > 0) {
                break;
            }
        }
    }
    
    if (m_socket->bytesAvailable() < 5 && timer.elapsed() >= timeoutMs) {
        setError("Timeout waiting for response");
        return qMakePair(STATUS_ERROR, QByteArray());
    }
    
    QByteArray headerData = m_socket->read(5);
    if (headerData.length() != 5) {
        setError(QString("Invalid response header length: %1").arg(headerData.length()));
        return qMakePair(STATUS_ERROR, QByteArray());
    }
    
    QDataStream headerStream(headerData);
    headerStream.setByteOrder(QDataStream::BigEndian);
    
    quint8 status;
    quint32 responseLength;
    headerStream >> status >> responseLength;
    
    // Read response data if any
    QByteArray responseData;
    if (responseLength > 0) {
        if (!m_socket->waitForReadyRead(1000)) { // Increased to 2 minutes
            setError("Timeout waiting for response data");
            return qMakePair(STATUS_ERROR, QByteArray());
        }
        
        responseData = m_socket->read(responseLength);
        if (static_cast<quint32>(responseData.length()) != responseLength) {
            setError(QString("Incomplete response data: expected %1, got %2")
                     .arg(responseLength).arg(responseData.length()));
            return qMakePair(STATUS_ERROR, QByteArray());
        }
    }
    
    return qMakePair(static_cast<int>(status), responseData);
}

// High-level Siril operations
QString SirilClient::getWorkingDirectory() {
    auto result = sendCommand(CMD_GET_WORKING_DIRECTORY);
    if (result.first == STATUS_OK) {
        return QString::fromUtf8(result.second);
    }
    return QString();
}

bool SirilClient::isImageLoaded() {
    auto result = sendCommand(CMD_GET_IS_IMAGE_LOADED);
    return result.first == STATUS_OK;
}

bool SirilClient::isSequenceLoaded() {
    auto result = sendCommand(CMD_GET_IS_SEQUENCE_LOADED);
    return result.first == STATUS_OK;
}

QString SirilClient::getUserDataDir() {
    auto result = sendCommand(CMD_GET_USERDATA_DIR);
    if (result.first == STATUS_OK) {
        return QString::fromUtf8(result.second);
    }
    return QString();
}

bool SirilClient::changeDirectory(const QString &path) {
    // Use the sendSirilCommand method which sends via CMD_SEND_COMMAND
    QString command = QString("cd \"%1\"").arg(path); // Quote the path in case of spaces
    qDebug() << "Sending cd command:" << command;
    return sendSirilCommand(command);
}

bool SirilClient::loadImage(const QString &filepath) {
    return sendSirilCommand(QString("load %1").arg(filepath));
}

bool SirilClient::saveImage(const QString &filepath) {
    return sendSirilCommand(QString("save %1").arg(filepath));
}

bool SirilClient::closeImage() {
    return sendSirilCommand("close");
}

bool SirilClient::platesolve(double ra, double dec, double focal, double pixelSize, bool force) {
    QString command = QString("platesolve %1 %2 -focal=%3 -pixelsize=%4")
                          .arg(ra, 0, 'f', 6)
                          .arg(dec, 0, 'f', 6)
                          .arg(focal, 0, 'f', 1)
                          .arg(pixelSize, 0, 'f', 2);
    
    if (force) {
        command += " -force";
    }
    
    return sendSirilCommand(command);
}

// Private slots
void SirilClient::onSocketConnected() {
    m_connectionTimer->stop();
    m_connected = true;
    qDebug() << "Connected to Siril";
    emit connected();
}

void SirilClient::onSocketDisconnected() {
    m_connected = false;
    qDebug() << "Disconnected from Siril";
    emit disconnected();
}

void SirilClient::onSocketError() {
    m_connectionTimer->stop();
    QString error = QString("Socket error: %1").arg(m_socket->errorString());
    qDebug() << error;
    setError(error);
}

void SirilClient::onConnectionTimeout() {
    m_socket->abort();
    setError("Connection timeout");
}

void SirilClient::setError(const QString &error) {
    m_lastError = error;
    emit errorOccurred(error);
}
