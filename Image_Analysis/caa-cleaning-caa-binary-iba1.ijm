// Cleaning binary Iba1 images and isolating included selections on bleached, 
// background subtracted images
// All filtering values updated for cohort 8 (new CAA controls)

function iba1Processing(bleach, raw, binaryPath, output) {
	print("Function running...");
	bleachedImages = getFileList(bleach);
	rawImages = getFileList(raw);
	binaryImages = getFileList(binaryPath);
	
	cleanedIba1BleachedOutput = output + "bleached-iba1-w-cleaned-binary-minus5-overlay/";
	cleanedIba1RawOutput = output + "raw-iba1-w-cleaned-binary-minus5-overlay/";
	cleanedBinaryIba1Output = output + "binary-iba1-cleaned-minus5-bleach-corrected-background-subtracted/"
	
	File.makeDirectory(cleanedIba1BleachedOutput);
	File.makeDirectory(cleanedIba1RawOutput);
	File.makeDirectory(cleanedBinaryIba1Output);
	print("Made output folders");
	
	for (i = 4; i < bleachedImages.length; i++) {
		// open bleached iba1 image
		//run("Bio-Formats Importer", "open=[" + bleach + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		open(bleach + bleachedImages[i]);
		bleachCurrent = getTitle();
		sampleName = bleachCurrent.substring(0, bleachCurrent.length - 32); // change if necessary
		print("Processing " + sampleName + "...");
		
		open(raw + rawImages[i]);
		rawCurrent = getTitle();
		
		// open binary iba1 image
		//run("Bio-Formats Importer", "open=[" + binaryPath + thresholdList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		open(binaryPath + binaryImages[i]);
		//run("Invert");
		orig_binary = getTitle();
		run("Duplicate...", "title=[" + sampleName + "_binary_duplicate]");
		binary_duplicate = getTitle();
		selectImage(orig_binary);
		
		// clean up through a series of dilations, erosions, and filtering
		dilationsErosions();
		de1 = getTitle();
		run("Analyze Particles...", "size=18-Infinity pixel show=Masks clear");
		run("Invert LUT");
		
		dilationsErosions();
		de2 = getTitle();
		run("Analyze Particles...", "size=26-Infinity pixel show=Masks clear");
	    run("Invert LUT");
		
		dilationsErosions();
		de3 = getTitle();
		run("Analyze Particles...", "size=88-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		mask = getTitle();
		print("Completed cleaning of binary image");
		
		// create selection on binary image
		selectImage(mask);
		run("Create Selection");
		close();
		
		//restore selection on duplicated binary
		selectImage(binary_duplicate);
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		run("Analyze Particles...", "size=51-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		saveAs("tiff", cleanedBinaryIba1Output + sampleName + "-binary-cleaned-selection-on-bleached.tif");
		cleanedMask = getTitle();
		print("Created cleaned binary image");
		
		// create selection on saved binary image
		selectImage(cleanedMask);
		run("Create Selection");
		selectImage(sampleName + "-binary-cleaned-selection-on-bleached.tif"); 
		close();
		
		//restore selection on duplicated iba1
		selectImage(bleachCurrent);
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1BleachedOutput + sampleName + "-bleached-cleaned.tif");
		print("Created cleaned bleached image");
		
		//restore selection on duplicated iba1
		selectImage(rawCurrent);
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1RawOutput + sampleName + "-bleached-cleaned.tif");
		print("Created cleaned bleached image");
		
		close("*");
		
		print("Done with " + sampleName + "\n");
		
	}
}

function dilationsErosions() {
	for (j=0; j < 2; j++) {
			run("Dilate");
		}
		for (j=0; j < 2; j++) {
			run("Erode");
		}
}

// run it!
//bleachedDir = getDirectory("Navigate to your bleached images ");
//binaryDir = getDirectory("Navigate to your binary Iba1 images ");
//outputDir = getDirectory("Navigate to your output folder ");

rawDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1/";
bleachedDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1-bleach-corrected/remontaged/";
binaryDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-iba1-cleaned-bleach-corrected-background-subtracted/minus5/";
outputDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/";


iba1Processing(bleachedDir, rawDir, binaryDir, outputDir);
print("Script finished");