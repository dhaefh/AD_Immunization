// Auto-thresholding DAPI images and running Watershed on result for CAA controls
// Saves both images created in new separate output folders within the original DAPI folder

function dapiThresholding(currentDir) {
    fileList = getFileList(currentDir);
    roiCounter = 0;
    for (i = 0; i < fileList.length; i++) {
    	open(currentDir + fileList[i]);
        current = getTitle();
		sampleName = current.substring(0, current.length - 4); // change if necessary
		print(sampleName);
		
		// Create a new directory with the sample name
		binaryImages = currentDir + "binary/";
		watershedImages = currentDir + "binary-cleaned/";
		File.makeDirectory(binaryImages);
		File.makeDirectory(watershedImages);
		
		// Rename the original image
		rename(sampleName + ".TIF");
		
		// Process image
		selectImage(sampleName + ".TIF");
		run("Duplicate...", "title=" + sampleName + "-dapi-binary.TIF");
		setAutoThreshold("Li dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("8-bit");
		saveAs("tif", binaryImages + "/" + sampleName + "-dapi-binary.tif");
		print("Saved " + binaryImages + "/" + sampleName + "-dapi-binary.tif");
		
		// Analyze Particles
		//selectImage(sampleName + "-dapi-binary.tif");
		run("Analyze Particles...", "size=0-17500 pixel show=Masks clear overlay");
			// adjust upper pixel bound depending on image resolution
		run("Invert");
		//run("Watershed");
		saveAs("tif", watershedImages + "/" + sampleName + "-dapi-particles-cleaned.tif");
		print("Saved " + watershedImages + "/" + sampleName + "-dapi-particles-cleaned.tif");
		
		close("*"); 
    }
}

// run it!
dapiDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/dapi/"
	// If necessary, split tifs to isolate DAPI images and add them all to the same folder
dapiThresholding(dapiDir);

print("Script finished");