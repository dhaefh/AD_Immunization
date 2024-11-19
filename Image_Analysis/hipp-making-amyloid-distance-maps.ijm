// creating distance maps of cleaned binary amyloid images for neighboring spot detection -- cohort 3
function amyloidDistanceChannels(amyloidDir, outputDir) {
	imageFileList = getFileList(amyloidDir);
	
	for (i = 0; i < imageFileList.length; i++) {
		// open motor neurons
		open(amyloidDir + imageFileList[i]);
		run("8-bit");
		original = getTitle();
		var sampleName = original.substring(0, original.length - 46); // change if necessary
		print("Processing " + sampleName + "...");
		roiCounter = 0;
		
		// create all dilated images
		setAutoThreshold("Default dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");

		current = getTitle();
		nextDilation(outputDir, original, "dilated_35X", 1);
		//saveAs("tiff", dilated35 + sampleName + "_dilated35");
		d35 = getTitle();
		print("Created 35X dilation");
		
		nextDilation(outputDir, d35, "dilated_70X", 2);
		//saveAs("tiff", dilated70 + sampleName + "_dilated70");
		d70 = getTitle();
		print("Created 70X dilation");
		
		nextDilation(outputDir, d70, "dilated_105X", 3);
		//saveAs("tiff", dilated105 + sampleName + "_dilated105");
		d105 = getTitle();
		print("Created 105X dilation");
		
		nextDilation(outputDir, d105, "dilated_140X", 4);
		//saveAs("tiff", dilated140 + sampleName + "_dilated140");
		d140 = getTitle();
		print("Created 140X dilation");
		
		nextDilation(outputDir, d140, "dilated_175X", 5);
		//saveAs("tiff", dilated175 + sampleName + "_dilated175");
		d175 = getTitle();
		print("Created 175X dilation");
		
		// ensure no overlapping regions
		selectImage(d140);
		run("Create Selection");
		selectImage(d175);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d175 = getTitle();
		
		selectImage(d105);
		run("Create Selection");
		selectImage(d140);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d140 = getTitle();
		
		selectImage(d70);
		run("Create Selection");
		selectImage(d105);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d105 = getTitle();
		
		selectImage(d35);
		run("Create Selection");
		selectImage(d70);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d70 = getTitle();
		
		selectImage(original);
		run("Create Selection");
		selectImage(d35);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d35 = getTitle();
		
		selectImage(d35);
		run("Divide...", "value=1.2");
		selectImage(d70);
		run("Divide...", "value=1.5");
		selectImage(d105);
		run("Divide...", "value=2");
		selectImage(d140);
		run("Divide...", "value=3");
		selectImage(d175);
		run("Divide...", "value=6");
		
		imageCalculator("Add create", current, d35);
		result = "Result of " + current;
		imageCalculator("Add create", result, d70);
		result = "Result of " + result;
		imageCalculator("Add create", result, d105);
		result = "Result of " + result;
		imageCalculator("Add create", result, d140);
		result = "Result of " + result;
		imageCalculator("Add create", result, d175);
		saveAs("tiff", outputDir + sampleName + "_cleaned_map");
		print("Done with " + sampleName);
		close("*");
	}
}

function nextDilation(outputDir, latestImage, newestImage, iteration){
	layer = iteration*35;
	
	if (!(File.exists(outputDir + "dilated" + layer + "_layer/" + sampleName + "_dilated" + layer + ".tif"))) {
		run("Duplicate...", "title=[dilate_" + iteration*35 + "x]");
		duplicate = getTitle();
		print("Duplicated " + latestImage);
		
		//selectImage(latestImage);
		//run("Create Selection");
		//roiManager("Add");
		//print("Created selection");
		//run("Select None");
		
		//selectImage(duplicate);
		for (j=0; j < 35; j++) {
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
binaryCortAmyloidDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/redone-amyloid-thresholded-cleaned/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/redone-amyloid-distance-maps/";
File.makeDirectory(outputDirectory);

dilated35 = outputDirectory + "dilated35_layer/";
dilated70 = outputDirectory + "dilated70_layer/";
dilated105 = outputDirectory + "dilated105_layer/";
dilated140 = outputDirectory + "dilated140_layer/";
dilated175 = outputDirectory + "dilated175_layer/";
File.makeDirectory(dilated35);
File.makeDirectory(dilated70);
File.makeDirectory(dilated105);
File.makeDirectory(dilated140);
File.makeDirectory(dilated175);

var sampleName = " ";
amyloidDistanceChannels(binaryCortAmyloidDirectory, outputDirectory);

print("Script finished");