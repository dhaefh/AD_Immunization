// creating distance maps of cleaned binary amyloid images for neighboring spot detection -- cohort 3
function amyloidDistanceChannels(amyloidDir, outputDir) {
	imageFileList = getFileList(amyloidDir);
	
	for (i = 16; i < 17; i++) {
		// open motor neurons
		open(amyloidDir + imageFileList[i]);
		run("8-bit");
		original = getTitle();
		var sampleName = original.substring(0, original.length - 42); // change if necessary
		print("Processing " + sampleName + "...");
		roiCounter = 0;
		
		// create all dilated images
		setAutoThreshold("Default dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		current = getTitle();
		nextDilation(outputDir, original, "dilated_65X", 1);
		//saveAs("tiff", dilated65 + sampleName + "_dilated65");
		d65 = getTitle();
		print("Created 65X dilation");
		
		nextDilation(outputDir, d65, "dilated_130X", 2);
		//saveAs("tiff", dilated130 + sampleName + "_dilated130");
		d130 = getTitle();
		print("Created 130X dilation");
		
		nextDilation(outputDir, d130, "dilated_195X", 3);
		//saveAs("tiff", dilated195 + sampleName + "_dilated195");
		d195 = getTitle();
		print("Created 195X dilation");
		
		nextDilation(outputDir, d195, "dilated_260X", 4);
		//saveAs("tiff", dilated260 + sampleName + "_dilated260");
		d260 = getTitle();
		print("Created 260X dilation");
		
		nextDilation(outputDir, d260, "dilated_325X", 5);
		//saveAs("tiff", dilated325 + sampleName + "_dilated325");
		d325 = getTitle();
		print("Created 325X dilation");
		
		// ensure no overlapping regions
		selectImage(d260);
		run("Create Selection");
		selectImage(d325);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d325 = getTitle();
		
		selectImage(d195);
		run("Create Selection");
		selectImage(d260);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d260 = getTitle();
		
		selectImage(d130);
		run("Create Selection");
		selectImage(d195);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d195 = getTitle();
		
		selectImage(d65);
		run("Create Selection");
		selectImage(d130);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d130 = getTitle();
		
		selectImage(original);
		run("Create Selection");
		selectImage(d65);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d65 = getTitle();
		
		selectImage(d65);
		run("Divide...", "value=1.2");
		selectImage(d130);
		run("Divide...", "value=1.5");
		selectImage(d195);
		run("Divide...", "value=2");
		selectImage(d260);
		run("Divide...", "value=3");
		selectImage(d325);
		run("Divide...", "value=6");
		
		imageCalculator("Add create", current, d65);
		result = "Result of " + current;
		imageCalculator("Add create", result, d130);
		result = "Result of " + result;
		imageCalculator("Add create", result, d195);
		result = "Result of " + result;
		imageCalculator("Add create", result, d260);
		result = "Result of " + result;
		imageCalculator("Add create", result, d325);
		saveAs("tiff", outputDir + sampleName + "_cleaned_map");
		print("Done with " + sampleName);
		close("*");
	}
	
	for (i = 19; i < 20; i++) {
		// open motor neurons
		open(amyloidDir + imageFileList[i]);
		run("8-bit");
		original = getTitle();
		var sampleName = original.substring(0, original.length - 42); // change if necessary
		print("Processing " + sampleName + "...");
		roiCounter = 0;
		
		// create all dilated images
		setAutoThreshold("Default dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		current = getTitle();
		nextDilation(outputDir, original, "dilated_65X", 1);
		//saveAs("tiff", dilated65 + sampleName + "_dilated65");
		d65 = getTitle();
		print("Created 65X dilation");
		
		nextDilation(outputDir, d65, "dilated_130X", 2);
		//saveAs("tiff", dilated130 + sampleName + "_dilated130");
		d130 = getTitle();
		print("Created 130X dilation");
		
		nextDilation(outputDir, d130, "dilated_195X", 3);
		//saveAs("tiff", dilated195 + sampleName + "_dilated195");
		d195 = getTitle();
		print("Created 195X dilation");
		
		nextDilation(outputDir, d195, "dilated_260X", 4);
		//saveAs("tiff", dilated260 + sampleName + "_dilated260");
		d260 = getTitle();
		print("Created 260X dilation");
		
		nextDilation(outputDir, d260, "dilated_325X", 5);
		//saveAs("tiff", dilated325 + sampleName + "_dilated325");
		d325 = getTitle();
		print("Created 325X dilation");
		
		// ensure no overlapping regions
		selectImage(d260);
		run("Create Selection");
		selectImage(d325);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d325 = getTitle();
		
		selectImage(d195);
		run("Create Selection");
		selectImage(d260);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d260 = getTitle();
		
		selectImage(d130);
		run("Create Selection");
		selectImage(d195);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d195 = getTitle();
		
		selectImage(d65);
		run("Create Selection");
		selectImage(d130);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d130 = getTitle();
		
		selectImage(original);
		run("Create Selection");
		selectImage(d65);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d65 = getTitle();
		
		selectImage(d65);
		run("Divide...", "value=1.2");
		selectImage(d130);
		run("Divide...", "value=1.5");
		selectImage(d195);
		run("Divide...", "value=2");
		selectImage(d260);
		run("Divide...", "value=3");
		selectImage(d325);
		run("Divide...", "value=6");
		
		imageCalculator("Add create", current, d65);
		result = "Result of " + current;
		imageCalculator("Add create", result, d130);
		result = "Result of " + result;
		imageCalculator("Add create", result, d195);
		result = "Result of " + result;
		imageCalculator("Add create", result, d260);
		result = "Result of " + result;
		imageCalculator("Add create", result, d325);
		saveAs("tiff", outputDir + sampleName + "_cleaned_map");
		print("Done with " + sampleName);
		close("*");
	}
}

function nextDilation(outputDir, latestImage, newestImage, iteration){
	layer = iteration*65;
	
	if (!(File.exists(outputDir + "dilated" + layer + "_layer/" + sampleName + "_dilated" + layer + ".tif"))) {
		run("Duplicate...", "title=[dilate_" + iteration*65 + "x]");
		duplicate = getTitle();
		print("Duplicated " + latestImage);
		
		//selectImage(latestImage);
		//run("Create Selection");
		//roiManager("Add");
		//print("Created selection");
		//run("Select None");
		
		//selectImage(duplicate);
		for (j=0; j < 65; j++) {
			run("Dilate");
			print(j);
		}
		print("Completed dilations");
		
		rename(newestImage);
		saveAs("tiff", outputDir + "dilated" + layer + "_layer/" + sampleName + "_dilated" + layer);
	}
	else {
		open(outputDir + "dilated" + layer + "_layer/" + sampleName + "_dilated" + layer + ".tif");
	}
}

function clearPreviousAndDivide(counter, n){
	run("Clear");
	run("Select None");
	run("Divide...", "value=" + n);
}

// run it!
binaryVascAmyloidDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/vascular-amyloid-thresholded/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/vascular-amyloid-distance-maps/";
File.makeDirectory(outputDirectory);

dilated65 = outputDirectory + "dilated65_layer/";
dilated130 = outputDirectory + "dilated130_layer/";
dilated195 = outputDirectory + "dilated195_layer/";
dilated260 = outputDirectory + "dilated260_layer/";
dilated325 = outputDirectory + "dilated325_layer/";
File.makeDirectory(dilated65);
File.makeDirectory(dilated130);
File.makeDirectory(dilated195);
File.makeDirectory(dilated260);
File.makeDirectory(dilated325);

var sampleName = " ";
//amyloidDistanceChannels(binaryVascAmyloidDirectory, outputDirectory);

binaryCortAmyloidDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/cortical-amyloid-thresholded/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/cortical-amyloid-distance-maps/";
File.makeDirectory(outputDirectory);

dilated65 = outputDirectory + "dilated65_layer/";
dilated130 = outputDirectory + "dilated130_layer/";
dilated195 = outputDirectory + "dilated195_layer/";
dilated260 = outputDirectory + "dilated260_layer/";
dilated325 = outputDirectory + "dilated325_layer/";
File.makeDirectory(dilated65);
File.makeDirectory(dilated130);
File.makeDirectory(dilated195);
File.makeDirectory(dilated260);
File.makeDirectory(dilated325);

amyloidDistanceChannels(binaryCortAmyloidDirectory, outputDirectory);
print("Script finished");