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
{
    // Set up connection timer
    m_connectionTimer->setSingleShot(true);
    m_connectionTimer->setInterval(5000); // 5 second timeout
    
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

bool SirilClient::connectToSiril() {
    if (m_connected) {
        return true;
    }
    
    QString socketPath = findSirilSocket();
    if (socketPath.isEmpty()) {
        setError("Could not find Siril socket. Is Siril running?\n"
                "Expected socket path like: /tmp/siril_XXXXX.sock");
        return false;
    }
    
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

void SirilClient::disconnectFromSiril() {
    if (m_socket->state() != QLocalSocket::UnconnectedState) {
        m_socket->disconnectFromServer();
        m_socket->waitForDisconnected(3000);
    }
    m_connected = false;
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
        if (!m_socket->waitForReadyRead(120000)) { // Increased to 2 minutes
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

bool SirilClient::sendSirilCommand(const QString &command) {
    QByteArray payload = command.toUtf8();
    qDebug() << "Sending Siril command:" << command;
    qDebug() << "Payload size:" << payload.size();
    
    auto result = sendCommand(CMD_SEND_COMMAND, payload);
    
    bool success = (result.first == STATUS_OK);
    
    qDebug() << "Command result - Status:" << result.first << "Response size:" << result.second.size();
    
    emit commandExecuted(command, success);
    
    if (!success && !result.second.isEmpty()) {
        setError(QString::fromUtf8(result.second));
    }
    
    return success;
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