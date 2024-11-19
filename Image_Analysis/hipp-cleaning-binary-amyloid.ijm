
function cleaningBinaryAmyloid (thresholdedAmyloid, outputFolder) {
	amyloidImages = getFileList(thresholdedAmyloid);
	print(amyloidImages.length);
	
	for (i = 0; i < amyloidImages.length; i++) {
		currentImage = amyloidImages[i];
		print(currentImage);
		sampleName = currentImage.substring(0, currentImage.length - 22);
		print(sampleName);
		
		open(thresholdedAmyloid + currentImage);
		//run("8-bit");
		original = getTitle();
		
		selectImage(original);

		Image.removeScale();
		setOption("BlackBackground", true);
		run("Dilate");
		run("Dilate");
		run("Erode");
		run("Erode");
		run("Analyze Particles...", "size=543-200000 show=Masks display clear add");
		run("Invert LUT");
		saveAs("tif", outputFolder + sampleName + "-amyloid-manual-threshold-ddee-filtered543-100000");
		print("Dilated, eroded, filtered, and saved");
		roiManager("save", outputFolder + sampleName + "-cleaned-amyloid-ROIs.zip");
		print("ROIs saved");
		
		close("*");
		close("ROI Manager");
		close("Results");
	}
}

binaryAmyloid = "/Users/ksd0849/b1042/alex/an1792-rebuttal-hippocampal/redone-amyloid-manual-threshold/";
output = "/Users/ksd0849/b1042/alex/an1792-rebuttal-hippocampal/redone-amyloid-thresholded-cleaned/";
File.makeDirectory(output);

cleaningBinaryAmyloid(binaryAmyloid, output);