
function ptauBackgroundSubtraction(origDir, outputDir) {
	print("Function running");
	orig_FL = getFileList(origDir);
	
	for (i = 0; i < orig_FL.length; i++) {
		open(origDir + orig_FL[i]);
		current = getTitle();
		sampleName = current.substring(0, current.length - 37);
		print("Processing " + sampleName);
		// subtract background
		run("Subtract Background...", "rolling=10");
		saveAs("tiff", outputDir + sampleName + "-bleach-corrected-after-stitching-background-subtracted-r10.tif");
		print("Background subtracted image saved");
		close("*");
	}
}

ptauOrig = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/iba1-3-bleach-corrected-after-stitching/";
ptauOut = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hippocampal/iba1-3-bleached-after-stitching-background-subtraction-r10/";
File.makeDirectory(ptauOut);

ptauBackgroundSubtraction(ptauOrig, ptauOut);

print("Script finished");