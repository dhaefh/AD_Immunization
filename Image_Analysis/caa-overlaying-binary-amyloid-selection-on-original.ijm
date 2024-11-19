// Cleaning manually thresholded amyloid IF images

function cleaningAmyloidBinary(amyloidBinary, amyloidOriginal, outputDir) {
	print("Function running...");
	amyloidBinaryFileList = getFileList(amyloidBinary);
	amyloidOriginalFileList = getFileList(amyloidOriginal);
	roiCounter = 0;
	
	for (i = 0; i < amyloidBinaryFileList.length; i++) {
		open(amyloidBinary + amyloidBinaryFileList[i]);
		current = getTitle();
		sampleName = current.substring(0, current.length - 12);
		print(sampleName);
		
		run("Create Selection");
		open(amyloidOriginal + amyloidOriginalFileList[i]);
		run("Restore Selection");
		//run("Flatten");
	    
	    saveAs("tiff", outputDir + sampleName + "original-amyloid-w-binary-selection-overlay");
		print("Done with " + sampleName);
		close("*");
	}
}

// run it!
// NEW CAA CONTROLs
amyloidBinary = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-amyloid-cleaned/final/alex-manual-threshold-plus5/";
amyloidOriginal = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/amyloid/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/orig-amyloid-w-cleaned-binary-overlay/";
File.makeDirectory(output);
cleaningAmyloidBinary(amyloidBinary, amyloidOriginal, output);

//OLD CAA CONTROLS
//amyloidBinary = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/cleaned-binary-amyloid/";
//amyloidOriginal = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/raw-amyloid/";
//output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/old-caa-controls/orig-amyloid-w-binary-overlay/";
//File.makeDirectory(output);
//cleaningAmyloidBinary(amyloidBinary, amyloidOriginal, output);

