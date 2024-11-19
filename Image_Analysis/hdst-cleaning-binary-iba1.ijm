// Cleaning binary Iba1 images and isolating included selections on bleached, 
// background subtracted images
// All filtering values updated for cohort 8 (new CAA controls)

function iba1Processing(bleach, raw, binaryPath, output) {
	print("Function running...");
	bleachedImages = getFileList(bleach);
	rawImages = getFileList(raw);
	binaryImages = getFileList(binaryPath);
	
	cleanedIba1BleachedOutput = output + "bleached-iba1-w-cleaned-binary-overlay-redone-again/";
	cleanedIba1RawOutput = output + "raw-iba1-w-cleaned-binary-overlay-redone-again/";
	cleanedBinaryIba1Output = output + "cleaned-binary-iba1-bleach-corrected-background-subtracted-redone-again/"
	
	File.makeDirectory(cleanedIba1BleachedOutput);
	File.makeDirectory(cleanedIba1RawOutput);
	File.makeDirectory(cleanedBinaryIba1Output);
	print("Made output folders");
	
	for (i = 0; i < bleachedImages.length; i++) {
		// open bleached iba1 image
		open(bleach + bleachedImages[i]);
		run("8-bit");
		bleachCurrent = getTitle();
		sampleName = bleachCurrent.substring(0, bleachCurrent.length - 37); // change if necessary
		print("Processing " + sampleName + "...");
		
		// open binary iba1 image
		open(binaryPath + binaryImages[i]);
		run("8-bit");
		orig_binary = getTitle();
		run("Duplicate...", "title=[" + sampleName + "_binary_duplicate]");
		binary_duplicate = getTitle();
		selectImage(orig_binary);
		
		// clean up through a series of dilations, erosions, and filtering
		dilationsErosions();
		de1 = getTitle();
		run("Analyze Particles...", "size=33-Infinity pixel show=Masks clear");
		run("Invert LUT");
		
		dilationsErosions();
		de2 = getTitle();
		run("Analyze Particles...", "size=48-Infinity pixel show=Masks clear");
	    run("Invert LUT");
	    selectImage(de1);
	    close();
		
		dilationsErosions();
		de3 = getTitle();
		run("Analyze Particles...", "size=164-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		mask = getTitle();
		print("Completed cleaning of binary image");
		selectImage(de2);
		close();
		selectImage(de3);
	    close();
		
		// create selection on binary image
		selectImage(mask);
		run("Create Selection");
		close();
		
		//restore selection on duplicated binary
		selectImage(binary_duplicate);
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		run("Analyze Particles...", "size=94-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		saveAs("tiff", cleanedBinaryIba1Output + sampleName + "-binary-cleaned.tif");
		cleanedMask = getTitle();
		print("Created cleaned binary image");
		
		// create selection on saved binary image
		selectImage(cleanedMask);
		run("Create Selection");
		selectImage(sampleName + "-binary-cleaned.tif"); 
		close();
		
		//restore selection on duplicated iba1
		selectImage(bleachCurrent);
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1BleachedOutput + sampleName + "-bleached-cleaned.tif");
		print("Created cleaned bleached image");
		
		//restore selection on duplicated iba1
		open(raw + rawImages[i]);
		run("8-bit");
		run("Restore Selection");
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1RawOutput + sampleName + "-raw-cleaned.tif");
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

rawDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1/redone-again/";
bleachedDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1-bleach-correction-fixed-restitching/iba1-bleach-corrected-after-stitching/redone-again/";
binaryDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1-bleach-corrected-after-stitching-background-subtracted-manual-threshold/redone-again/";
outputDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/";

iba1Processing(bleachedDir, rawDir, binaryDir, outputDir);
print("Script finished");