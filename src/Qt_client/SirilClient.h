#ifndef SIRILCLIENT_H
#define SIRILCLIENT_H

#include <QObject>
#include <QLocalSocket>
#include <QString>
#include <QByteArray>
#include <QPair>
#include <QTimer>

// Command IDs from Siril
enum SirilCommand {
    CMD_SEND_COMMAND = 1,
    CMD_LOG_MESSAGE = 2,
    CMD_UPDATE_PROGRESS = 3,
    CMD_GET_WORKING_DIRECTORY = 4,
    CMD_GET_FILENAME = 5,
    CMD_GET_DIMENSIONS = 6,
    CMD_GET_PIXELDATA = 7,
    CMD_GET_PIXELDATA_REGION = 8,
    CMD_RELEASE_SHM = 9,
    CMD_SET_PIXELDATA = 10,
    CMD_GET_IMAGE_STATS = 11,
    CMD_GET_KEYWORDS = 12,
    CMD_GET_ICC_PROFILE = 13,
    CMD_GET_FITS_HEADER = 14,
    CMD_GET_FITS_HISTORY = 15,
    CMD_GET_FITS_UNKNOWN_KEYS = 16,
    CMD_GET_IMAGE = 17,
    CMD_GET_PSFSTARS = 18,
    CMD_GET_SEQ_STATS = 19,
    CMD_GET_SEQ_REGDATA = 20,
    CMD_GET_SEQ_IMGDATA = 21,
    CMD_GET_SEQ_PIXELDATA = 22,
    CMD_GET_SEQ_IMAGE = 23,
    CMD_GET_SEQ = 24,
    CMD_GET_CONFIG = 25,
    CMD_GET_USERCONFIG_DIR = 26,
    CMD_GET_IS_IMAGE_LOADED = 27,
    CMD_GET_IS_SEQUENCE_LOADED = 28,
    CMD_GET_SELECTION = 29,
    CMD_SET_SELECTION = 30,
    CMD_GET_ACTIVE_VPORT = 31,
    CMD_GET_STAR_IN_SELECTION = 32,
    CMD_GET_STATS_FOR_SELECTION = 33,
    CMD_PIX2WCS = 34,
    CMD_WCS2PIX = 35,
    CMD_UNDO_SAVE_STATE = 36,
    CMD_GET_BUNDLE_PATH = 37,
    CMD_ERROR_MESSAGEBOX = 38,
    CMD_ERROR_MESSAGEBOX_MODAL = 39,
    CMD_PLOT = 40,
    CMD_CLAIM_THREAD = 41,
    CMD_RELEASE_THREAD = 42,
    CMD_SEQ_FRAME_SET_PIXELDATA = 43,
    CMD_REQUEST_SHM = 44,
    CMD_SET_SEQ_FRAME_INCL = 45,
    CMD_GET_USERDATA_DIR = 46,
    CMD_GET_SYSTEMDATA_DIR = 47
};

// Status codes
enum SirilStatus {
    STATUS_OK = 0,
    STATUS_NONE = 1,
    STATUS_ERROR = 255
};

class SirilClient : public QObject {
    Q_OBJECT

public:
    explicit SirilClient(QObject *parent = nullptr);
    ~SirilClient();
    
    // Core connection management
    bool connectToSiril();
    void disconnectFromSiril();
    bool isConnected() const { return m_connected; }
    
    // High-level Siril operations for Stellina processing
    QString getWorkingDirectory();
    bool isImageLoaded();
    bool isSequenceLoaded();
    QString getUserDataDir();
    bool changeDirectory(const QString &path);
    bool loadImage(const QString &filepath);
    bool saveImage(const QString &filepath);
    bool closeImage();
    bool platesolve(double ra, double dec, double focal = 400.0, double pixelSize = 2.40, bool force = true);
    bool sendSirilCommand(const QString &command);
    
    // Error handling
    QString lastError() const { return m_lastError; }

signals:
    void connected();
    void disconnected();
    void commandExecuted(const QString &command, bool success);
    void errorOccurred(const QString &error);
    void logMessage(const QString &message);

private slots:
    void onSocketConnected();
    void onSocketDisconnected();
    void onSocketError();
    void onConnectionTimeout();

private:
    QLocalSocket *m_socket;
    QTimer *m_connectionTimer;
    bool m_connected;
    QString m_lastError;
    
    // Protocol helpers
    QString findSirilSocket();
    QPair<int, QByteArray> sendCommand(int commandId, const QByteArray &payload = QByteArray());
    QByteArray createCommandHeader(int commandId, int payloadLength);
    QPair<int, QByteArray> readResponse();
    
    void setError(const QString &error);
};

#endif // SIRILCLIENT_H