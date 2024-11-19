// Cleaning manually thresholded amyloid IF images

function overlayingIba1BinaryOnBleached(iba1Bin, iba1Orig, outputDir) {
	print("Function running...");
	iba1BinaryFileList = getFileList(iba1Bin);
	iba1OriginalFileList = getFileList(iba1Orig);
	roiCounter = 0;
	
	for (i = 0; i < iba1BinaryFileList.length; i++) {
		open(iba1Bin + iba1BinaryFileList[i]);
		current = getTitle();
		print(current);
		sampleName = current.substring(0, 23);
		print(sampleName);
		
		run("Create Selection");
		open(iba1Orig + iba1OriginalFileList[i]);
		run("Restore Selection");
		//run("Flatten");
	    
	    saveAs("tiff", outputDir + sampleName + "-bleached-iba1-w-binary-selection-overlay");
		print("Done with " + sampleName);
		close("*");
	}
}

// run it!
// NEW CAA CONTROLs
iba1CleanedBinary = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-iba1-cleaned-bleach-corrected-background-subtracted/";
iba1Bleached = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1-bleach-corrected/remontaged/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/bleached-iba1-w-binary-selection-overlay/";
File.makeDirectory(output);
overlayingIba1BinaryOnBleached(iba1CleanedBinary, iba1Bleached, output);

//OLD CAA CONTROLS
//amyloidBinary = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/cleaned-binary-amyloid/";
//amyloidOriginal = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/raw-amyloid/";
//output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/orig-amyloid-w-binary-overlay/";
//File.makeDirectory(output);
//cleaningAmyloidBinary(amyloidBinary, amyloidOriginal, output);

