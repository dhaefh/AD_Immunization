//creating distance maps of cleaned binary amyloid images for neighboring spot detection -- cohort 8
function amyloidDistanceChannels(amyloidDir, outputDir) {
	imageFileList = getFileList(amyloidDir);
	
	for (i = 0; i < imageFileList.length; i++) {
		// open motor neurons
		//run("Bio-Formats Importer", "open=[" + mnDir + neuronFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		open(amyloidDir + imageFileList[i]);
		original = getTitle();
		sampleName = original.substring(0, original.length - 10); // change if necessary
		print("Processing " + sampleName + "...");
		
		// create all dilated images
		setAutoThreshold("Default dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		current = getTitle();
		nextDilation(roiCounter, original, "dilated_35X", 1);
		saveAs("tiff", dilated35 + sampleName + "_dilated35");
		d35 = getTitle();
		print("Created 35X dilation");
		
		nextDilation(roiCounter, d20, "dilated_70X", 2);
		saveAs("tiff", dilated40 + sampleName + "_dilated70");
		d70 = getTitle();
		print("Created 70X dilation");
		
		nextDilation(roiCounter, d40, "dilated_105X", 3);
		saveAs("tiff", dilated60 + sampleName + "_dilated105");
		d105 = getTitle();
		print("Created 105X dilation");
		
		nextDilation(roiCounter, d60, "dilated_140X", 4);
		saveAs("tiff", dilated80 + sampleName + "_dilated140");
		d140 = getTitle();
		print("Created 140X dilation");
		
		nextDilation(roiCounter, d80, "dilated_175X", 5);
		saveAs("tiff", dilated100 + sampleName + "_dilated175");
		d175 = getTitle();
		print("Created 175X dilation");
		
		// ensure no overlapping regions
		selectImage(original);
		run("Create Selection");
		selectImage(d35);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d35 = getTitle();
		
		selectImage(d35);
		run("Create Selection");
		selectImage(d70);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d70 = getTitle();
		
		selectImage(d70);
		run("Create Selection");
		selectImage(d105);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d105 = getTitle();
		
		selectImage(d105);
		run("Create Selection");
		selectImage(d140);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d140 = getTitle();
		
		selectImage(d140);
		run("Create Selection");
		selectImage(d175);
		run("Restore Selection");
		run("Clear");
		run("Select None");
		d175 = getTitle();
		
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
function nextDilation(counter, latestImage, newestImage, iteration){
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
}
function clearPreviousAndDivide(counter, n){
	run("Clear");
	run("Select None");
	run("Divide...", "value=" + n);
}
// run it!
binaryAmyloidDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid/";
outputDirectory = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid-distance-maps/";
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
		
motorNeuronDistanceChannels(binaryAmyloidDirectory, outputDirectory);
print("Script finished");