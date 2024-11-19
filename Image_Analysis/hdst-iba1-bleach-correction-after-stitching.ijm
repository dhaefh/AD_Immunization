// Standardizing microglia signal by performing bleach correction and background subtraction

function bleachCorrection(origDir, outputDir) {
	print("Function running");
	//orig_FL = getFileList(origDir);
	
	//for (i = 0; i < orig_FL.length; i++) {
	//open(origDir + orig_FL[i]);
	open(origDir);
	original = getTitle();
	sampleName = original.substring(0, original.length - 17);
	print("Processing " + sampleName);
	
	// make output folders
	bleachedImages = outputDir + "iba1-bleach-corrected-after-stitching/";
	backgroundSubImages = outputDir + "iba1-bleach-corrected-after-stitching-background-subtracted/";
	File.makeDirectory(bleachedImages);
	File.makeDirectory(backgroundSubImages);
	print("Made output folders");

	run("Montage to Stack...", "columns=25 rows=25 border=0");
	originalStack = getTitle();
	print("Made stack");
	
	// make middle slides first to ensure the standard is representative of the tissue 
	run("Select All");
	run("Make Substack...", "slices=312-313");
	rename("standardsForBleachCorrection");
	run("Concatenate...", "open image1=standardsForBleachCorrection image2=" + originalStack + " image3=[-- None --]");
	bigStack = getTitle();
	print("Reordered slices");
	
	// standardize signal to noise in all panels
	run("Bleach Correction", "correction=[Histogram Matching]");
	print("Bleaching done");
	
	// re order slices by label
	setSlice(1);
	run("Delete Slice");
	run("Delete Slice");

	print("Stack returned to original order");

	// re stitch into one image
	run("Make Montage...", "columns=25 rows=25 scale=1");
	saveAs("tiff", bleachedImages + sampleName + "-bleach-corrected-after-stitching.tif");
	print("Bleached image saved");
	
	// subtract background
	run("Subtract Background...", "rolling=7");
	saveAs("tiff", backgroundSubImages + sampleName + "-bleach-corrected-after-stitching-background-subtracted.tif");
	
	close("*");
	//}
}

hd1021 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1/HD_102-1-ch2-stitched.tif";
hdA1 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1/NMA22-300-A1-ch2-stitched.tif";
outputFolder = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/iba1/iba1-bleached-after-restitching/";
File.makeDirectory(outputFolder);


bleachCorrection(hd1021, outputFolder);
bleachCorrection(hdA1, outputFolder);

print("Script finished");