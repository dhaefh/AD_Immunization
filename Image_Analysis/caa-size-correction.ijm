// Performs size correction on images to account for alteration in dimensions 
// of microglia images during bleach correction
// Ensures all multipage components will be the same size and are thus compatible
// to be combined into one file

function sizeCorrection(inputFolder) {
	print("Function running");
	images_FL = getFileList(inputFolder);
	print(images_FL.length);
 	
	for (i = 0; i < 8; i++) {
		open(inputFolder + images_FL[i]);
		print(images_FL[i]);
		current = getTitle();
		sampleName = current.substring(0, 25);
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

corticalAmyloidDM = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid-distance-maps/";
corticalAmyloidOuterLayer = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid-distance-maps/dilated325_layer/";
corticalAmyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid/";

//imagesOfInterest = getDirectory("Navigate to the folder holding the images you wish to size correct ");

sizeCorrection(corticalAmyloidDM);
sizeCorrection(corticalAmyloidOuterLayer);
sizeCorrection(corticalAmyloid);


print("Script finished");