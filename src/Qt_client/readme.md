# Stellina Processor for Siril

A Qt-based desktop application for processing Stellina telescope images using Siril's external IPC interface.

## Features

- **Real-time Siril Integration**: Communicates directly with Siril via IPC (no script generation needed)
- **Stellina Quality Filtering**: Automatically filters images based on Stellina's quality assessment
- **Coordinate Conversion**: Converts Alt/Az coordinates from Stellina JSON to RA/Dec for plate solving
- **Batch Processing**: Process multiple images with error recovery and progress tracking
- **Modern Qt Interface**: Clean, responsive user interface with real-time logging
- **Cross-platform**: Works on macOS, Linux, and Windows

## Prerequisites

1. **Siril**: Must have Siril 1.4.0+ with external IPC enabled (your modified version)
2. **Qt**: Qt 5.15+ or Qt 6.x with Widgets and Network modules
3. **C++17 Compiler**: GCC, Clang, or MSVC with C++17 support

## Building

### macOS / Linux
```bash
# Clone or extract the source code
cd StellinaProcessor

# Generate Makefile
qmake StellinaProcessor.pro

# Build
make

# Run
./StellinaProcessor
```

### Using Qt Creator
1. Open `StellinaProcessor.pro` in Qt Creator
2. Configure your kit (compiler + Qt version)
3. Build and run

### Command Line with specific Qt version
```bash
# If you have multiple Qt versions installed
/path/to/qt/bin/qmake StellinaProcessor.pro
make
```

## Usage

### 1. Start Siril
First, start your modified version of Siril with external IPC enabled. You should see:
```
log: External IPC enabled on /tmp/siril_XXXX.sock
```

### 2. Launch Stellina Processor
Start the application and click "Test Connection" to verify communication with Siril.

### 3. Configure Processing
- **Source Directory**: Select folder containing Stellina .fits and .json file pairs
- **Output Directory**: Select where processed images should be saved
- **Processing Options**:
  - Enable quality filtering to skip images rejected by Stellina
  - Set focal length (default: 400mm for Stellina)
  - Set pixel size (default: 2.40μm for Stellina)
  - Set observer location for coordinate conversion

### 4. Process Images
Click "Start Processing" to begin batch processing. The application will:
1. Find all .fits/.json pairs in the source directory
2. Filter by quality (if enabled)
3. For each image:
   - Extract Alt/Az coordinates from JSON
   - Convert to RA/Dec coordinates
   - Load image in Siril
   - Perform plate solving with coordinate hints
   - Save processed image
   - Continue to next image even if plate solving fails

## File Structure

```
StellinaProcessor/
├── StellinaProcessor.pro     # qmake project file
├── main.cpp                  # Application entry point
├── SirilClient.h/cpp         # Siril IPC communication
├── StellinaProcessor.h/cpp   # Main application window
└── README.md                 # This file
```

## Siril Command Protocol

The application uses Siril's binary IPC protocol:

### Command Header
```cpp
struct CommandHeader {
    uint8_t command;          // Command ID
    uint32_t length;          // Payload length (big-endian)
} __attribute__((packed));
```

### Key Commands Used
- `CMD_GET_WORKING_DIRECTORY (4)`: Get current working directory
- `CMD_GET_IS_IMAGE_LOADED (27)`: Check if image is loaded
- `CMD_SEND_COMMAND (1)`: Send text commands like "load", "save", "platesolve"

### Response Format
```cpp
struct ResponseHeader {
    uint8_t status;           // 0=OK, 1=NONE, 255=ERROR
    uint32_t length;          // Response data length (big-endian)
} __attribute__((packed));
```

## Troubleshooting

### "Could not find Siril socket"
- Ensure Siril is running with external IPC enabled
- Check that `/tmp/siril_*.sock` exists (macOS/Linux)
- Verify you're using the modified Siril version

### "Connection timeout"
- Siril may be busy or unresponsive
- Try restarting Siril
- Check Siril logs for errors

### "Plate solving failed"
- This is normal for some images - processing continues
- Check coordinate conversion is working correctly
- Verify focal length and pixel size settings

### No images found
- Ensure .fits and .json files are paired (same basename)
- Check file permissions
- Enable debug mode for detailed logging

## Known Limitations

1. **Coordinate Conversion**: Currently uses placeholder coordinate conversion. You'll need to implement proper Alt/Az to RA/Dec conversion using astronomical libraries or your existing Python code.

2. **Quality Assessment**: Basic quality checking implementation. You may need to customize based on your specific Stellina JSON format.

3. **Error Recovery**: While the application continues processing after individual failures, some errors may require manual intervention.

## Development Notes

### Adding Coordinate Conversion
Replace the placeholder in `StellinaProcessor::convertAltAzToRaDec()` with:
- Port your existing Python coordinate conversion code
- Use astronomical libraries like SOFA
- Call external coordinate conversion tools

### Customizing Quality Assessment
Modify `StellinaProcessor::checkStellinaQuality()` based on your JSON structure:
```cpp
bool StellinaProcessor::checkStellinaQuality(const QJsonObject &json) {
    // Add your quality criteria here
    if (json.contains("your_quality_field")) {
        return json["your_quality_field"].toBool();
    }
    return true; // Default to accept
}
```

### Adding New Siril Commands
Add new command IDs to the `SirilCommand` enum and implement wrapper functions in `SirilClient`.

## License

This project is designed to work with your existing Stellina processing pipeline and Siril modifications. Adjust licensing as appropriate for your use case.