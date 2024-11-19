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

function cleanAndExpandAmyloid () {
	amyloidImages = getFileList(thresholdedAmyloid);
	
	for (i = 14; i < amyloidImages.length; i++) {		
		open(thresholdedAmyloid + amyloidImages[i]);
		currentImage = getTitle();
		sampleName = currentImage.substring(0, currentImage.length - 18);
		print(sampleName);
		
		run("Dilate");
		run("Dilate");
		run("Erode");
		run("Erode");
		run("Analyze Particles...", "size=1000-Infinity show=Masks");
		run("Invert LUT");
		saveAs("tif", cleanedAmyloid + sampleName + "-amyloid-manual-threshold-ddee-filtered1000");
		print("Cleaned binary amyloid and saved\n");
		
		for (j = 0; j < 30; j++) {
			run("Dilate");
		}
		
		// save ULEX distance map => expanded by x µm
		saveAs("tif", expandedAmyloid + sampleName + "-amyloid-manual-threshold-cleaned-expanded");
		print("Saved expanded amyloid signal");
		close("*");
		
		close("*");
		
	}
}

// takes the parameter passed in the bash script (path to the main folder containing stitched images)
thresholdedAmyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-thresholded/manually-thresholded-082024/";
cleanedAmyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-manual-threshold-cleaned/";
expandedAmyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-manual-threshold-cleaned-expanded/";
File.makeDirectory(cleanedAmyloid);
File.makeDirectory(expandedAmyloid);

print("\nFunction running");
cleanAndExpandAmyloid();