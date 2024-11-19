// Cleaning manually thresholded amyloid IF images

function cleaningAmyloidBinary(amyloidBinary, outputDir) {
	print("Function running...");
	amyloidFileList = getFileList(amyloidBinary);
	roiCounter = 0;
	
	for (i = 0; i < amyloidFileList.length; i++) {
		open(amyloidBinary + amyloidFileList[i]);
		current = getTitle();
		sampleName = current.substring(0, current.length - 35);
		
		// dilation erosion
	    run("Dilate");
	    run("Dilate");
	    run("Erode");
	    run("Erode");
	    print("Ran dilation/erosion");
	    
	    run("Analyze Particles...", "size=438-infinity show=Masks clear");
	    	// adjust lower bound based on image resolution
	    run("Invert LUT");
	    print("Filtered");
	    
	    saveAs("tiff", outputDir + sampleName + "_cleaned_binary");
		print("Done with " + sampleName);
		close("*");
	}
}

// run it!
amyloidDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-amyloid/alex-plus5/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-amyloid-cleaned/final/";

cleaningAmyloidBinary(amyloidDirectory, outputDirectory);
print("Script finished");