run: Release/StellinaProcessor.app/Contents/MacOS/StellinaProcessor
	$<

Release/StellinaProcessor.app/Contents/MacOS/StellinaProcessor: qmake
	xcodebuild

qmake:
	qmake -spec macx-xcode StellinaProcessor.pro

moc:
