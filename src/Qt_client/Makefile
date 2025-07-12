run: Release/OriginSimulator.app/Contents/MacOS/OriginSimulator
	$<

Release/OriginSimulator.app/Contents/MacOS/OriginSimulator: qmake
	xcodebuild

qmake:
	qmake -spec macx-xcode StellinaProcessor.pro

moc:
