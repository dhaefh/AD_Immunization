// Performs size correction on images to account for alteration in dimensions 
// of microglia images during bleach correction
// Ensures all multipage components will be the same size and are thus compatible
// to be combined into one file

function sizeCorrection(inputFolder) {
	print("Function running");
	images_FL = getFileList(inputFolder);
	print(images_FL.length);
 	
	for (i = 0; i < images_FL.length; i++) {
		open(inputFolder + images_FL[i]);
		print(images_FL[i]);
		current = getTitle();
		sampleName = current.substring(0, 12);
		print("Processing " + sampleName);
		
		// make output folders
		outputFolder = inputFolder + "un-re-stitched/";
		File.makeDirectory(outputFolder);
		print("Made output folder");

		selectImage(current);
		run("Montage to Stack...", "columns=25 rows=25 border=0");
		print("Made stack");
		
		// re stitch into one image
		run("Make Montage...", "columns=25 rows=25 scale=1");
		
		saveAs("tiff", outputFolder + sampleName + "-un-re-stitched");
		close();
		print("Saved image that has been unstitched and restitched");
		
		print("Done with " + sampleName + "\n");
		
		close("*");
	}
}

// run it!

amyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/amyloid-thresholded-cleaned/";
amyloidDM = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/amyloid-distance-maps/";

//imagesOfInterest = getDirectory("Navigate to the folder holding the images you wish to size correct ");
sizeCorrection(amyloid);
sizeCorrection(amyloidDM);

print("Script finished");