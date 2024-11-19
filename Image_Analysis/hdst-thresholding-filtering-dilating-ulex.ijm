// Alex Edwards, 7-10-2024
// Description: Takes in stitched unprocessed amyloid images as well as stitched, background subtracted ULEX images
// 				Thresholds and filters ULEX images to establish ground truth
//				Creates distance maps for ULEX (only one layer => x µm)
//				Saves coordinates for all ROIs per sample
//				Uses ULEX distance map to distinguish vascular vs cortical amyloid
//				Saves all new data (listed below) in newly created output folders
//					Binarized, filtered ULEX images
//					Undilated ULEX roi coordinates
//					ULEX distance map images
//					ULEX distance map ROI coordinates
//					Cortical amyloid images
//					Vascular amyloid images
// Notes: script assumes all folders have images for each sample that are listed in the same order

function ulexBinarizing (ulexRB) {
	allImages = getFileList(ulexRB);
	
	for (i = allImages.length - 1; i < allImages.length; i++) {
		currentImage = allImages[i];
		print(currentImage);
		open(ulexRB + currentImage);
		// get sample name
		sampleName = currentImage.substring(0, currentImage.length - 21);
		setAutoThreshold("Triangle dark no-reset");
		print("Thresholded");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Analyze Particles...", "size=10-Infinity show=Masks clear");
		print("Filtered");
		run("Invert LUT");
		print("Inverted LUT");
		//getROIcoordinates("undilated"); // save coordinates of all ULEX signals
		//print("Saved undilated ULEX ROI coordinates");
		
		// save binarized ULEX
		saveAs("tif", ulexBinary + sampleName + "-ulex-triangleThreshold-filter10");
		print("Saved binary filtered ULEX image");
		
		rename("thresholded");
	}
}

function ulexDistanceMaps (ulexBin) {
	print("In function");
	allImages = getFileList(ulexBin);
	print(allImages.length);
	
	for (i = 0; i < allImages.length; i++) {
		print("In loop");
		currentImage = allImages[i];
		print(currentImage);
		open(ulexBinary + currentImage);
		// get sample name
		sampleName = currentImage.substring(0, currentImage.length - 36);
		
		for (j = 0; j < 30; j++) {
			run("Dilate");
		}
		print("Made ULEX distance map");
		
		// save ULEX distance map => expanded by x µm
		saveAs("tif", ulexDM + sampleName + "-ulex-triangleThreshold-filter10-distance-map");
		print("Saved ULEX distance map");
		
		rename("distance map");
	}
}

function isolateVascularAmyloid (ulexBin, ulexDistFolder, amyloidOrig, amyloidBin) {
	ulexBinaryImages = getFileList(ulexBin);
	ulexDMImages = getFileList(ulexDistFolder);
	origAmyloidImages = getFileList(amyloidOrig);
	cleanedBinaryAmyloidImages = getFileList(amyloidBin);
	
	for (i = 0; i < ulexBinaryImages.length; i++) {
		open(ulexBin + ulexBinaryImages[i]);
		currentImage = getTitle();
		rename("thresholded");
		sampleName = currentImage.substring(0, currentImage.length - 36);
		//getROIcoordinates("not-expanded"); // save coordinates of all regular binary ULEX signals
		//print("Saved ULEX distance map ROI coordinates");
		
		open(ulexDistFolder + ulexDMImages[i]);
		rename("distance map");
		//getROIcoordinates("expanded"); // save coordinates of all dilated ULEX signals
		//print("Saved ULEX distance map ROI coordinates");
		
		open(amyloidOrig + origAmyloidImages[i]);
		vascular = getTitle();
		run("Duplicate...", "title=cortical-amyloid");
		cortical = getTitle();
		
		open(amyloidBin + cleanedBinaryAmyloidImages[i]);
		vascularBin = getTitle();
		run("Duplicate...", "title=cortical-amyloid-binary");
		corticalBin = getTitle();
		
		//selectImage("distance map");
		//run("Create Selection");
		//selectImage(vascular);
		//run("Restore Selection");
		//run("Clear Outside");
		//run("Select None");
		//print("Isolated vascular amyloid");
		//saveAs("tif", vascAmyloid + sampleName + "-vascular-amyloid");
		
		selectImage("distance map");
		run("Create Selection");
		selectImage(cortical);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		print("Isolated cortical amyloid");
		saveAs("tif", cortAmyloid + sampleName + "-cortical-amyloid");
		
		//selectImage("distance map");
		//run("Create Selection");
		//selectImage(vascularBin);
		//run("Restore Selection");
		//run("Clear Outside");
		//run("Select None");
		//print("Isolated binary vascular amyloid");
		//saveAs("tif", vascBinaryAmyloid + sampleName + "-vascular-amyloid-thresholded-filtered");
		
		selectImage("distance map");
		run("Create Selection");
		selectImage(corticalBin);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		print("Isolated binary cortical amyloid");
		saveAs("tif", cortBinaryAmyloid + sampleName + "-cortical-amyloid-thresholded-filtered");
		
		close("*");
		
	}
}

function getROIcoordinates (type) {
	run("Create Selection");
	print("Created selection");
	roiManager("Add");
	//roiManager("Split");
	roiManager("save", outputFolder + sampleName + "-vasculature-" + type + "-rois.zip");
}

// takes the parameter passed in the bash script (path to the main folder containing stitched images)
//args = getArgument();
args = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/";

// make sorted folder
ulexBackgroundSubtracted = args + "ulex-background-subtracted/";
amyloidStitched = args + "amyloid/";
cleanedBinaryAmyloidStitched = args + "cleaned-binary-amyloid/";

// makes all necessary output folders
ulexBinary = args + "ulex-thresholded-filtered/";
ulexDM = args + "ulex-distance-map/";
vascAmyloid = args + "vascular-amyloid/";
cortAmyloid = args + "cortical-amyloid/";
vascBinaryAmyloid = args + "vascular-amyloid-thresholded/";
cortBinaryAmyloid = args + "cortical-amyloid-thresholded/";

File.makeDirectory(ulexBinary);
File.makeDirectory(ulexDM);
File.makeDirectory(vascAmyloid);
File.makeDirectory(cortAmyloid);
File.makeDirectory(vascBinaryAmyloid);
File.makeDirectory(cortBinaryAmyloid);

print("\nFunction running");
//ulexBinarizing(ulexBackgroundSubtracted);
//ulexDistanceMaps(ulexBinary);
isolateVascularAmyloid(ulexBinary, ulexDM, amyloidStitched, cleanedBinaryAmyloidStitched);