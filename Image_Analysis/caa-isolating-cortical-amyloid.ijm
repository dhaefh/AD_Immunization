// ImageJ macro for using the cleaned amyloid vessel channels to isolate
// the amyloid plaque signals in the original binary IF amyloid images

function isolatingAmyloidPlaques(cleanedVesselProbMap, amyloidBinary, outputDir) {
	print("Function running...");
	vessel_FL = getFileList(cleanedVesselProbMap);
	allAmyloid_FL = getFileList(amyloidBinary);
	
	for (i = 0; i < vessel_FL.length; i++) {
		open(cleanedVesselProbMap + vessel_FL[i]);
	    current = vessel_FL[i];
	    sampleName = current.substring(0, 23);
	    print("Processing " + sampleName + "...");
	    print("Cleaned vessel image opened");
	    
	    for (j = 0; j < 18; j++) {
	    	run("Dilate");
	    }
	    
	    run("Create Selection");
	    print("Dilated vessels and created selection");
	    
	    open(amyloidBinary + allAmyloid_FL[i]);
	    print("Opened comprehensive binary amyloid image");
	    
	    run("Restore Selection");
	    run("Clear");
	    print("Cleared vessels from binary image to isolate plaques");
		
	    saveAs("tiff", outputDir + sampleName + "_isolated_plaque_signal.tif");
	    print("Saved image");
	    print("Done with " + sampleName + "\n");
	    
	    close("*");
	}
}

// run it!
probabilityMaps = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid/";
binaryAmyloidImages = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-amyloid-cleaned/final/alex-manual-threshold-plus5/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-cortical-amyloid/";

isolatingAmyloidPlaques(probabilityMaps, binaryAmyloidImages, outputDirectory);
print("Script finished");