// Standardizing microglia signal by performing bleach correction and background subtraction

function bleachCorrection(origDir, outputDir) {
	print("Function running");
	orig_FL = getFileList(origDir);
	
	for (i = 0; i < orig_FL.length; i++) {
		open(origDir + orig_FL[i]);
		original = getTitle();
		sampleName = original.substring(0, original.length - 17);
		print("Processing " + sampleName);
		
		// make output folders
		bleachedImages = outputDir + "iba1-bleach-corrected/";
		backgroundSubImages = outputDir + "iba1-bleach-corrected-background-subtracted/";
		File.makeDirectory(bleachedImages);
		File.makeDirectory(backgroundSubImages);
		print("Made output folders");

		run("Montage to Stack...", "columns=25 rows=25 border=0");
		print("Made stack");
		
		numSlices = nSlices;
		
		// label slices in order
		for (j = 1; j <= numSlices; j++) {
			setSlice(j);
			setMetadata("Label","" + j);
		}
		print("Labeled slices");
		
		// make middle slides first to ensure the standard is representative of the tissue 
		run("Select All");
		run("Make Substack...", "slices=312-313 delete");
		run("Concatenate...", "  title=shuffled keep image1=[Substack (312-313)] image2=Stack");
		bigStack = getTitle();
		print("Reordered slices");
		
		// standardize signal to noise in all panels
		//run("Bleach Correction", "correction=[Histogram Matching]");
		//print("Bleaching done");
		
		// re order slices by label
		run("Make Substack...", "slices=1-2 delete");
		rename("stack312-313");
		
		selectImage(bigStack);
		run("Make Substack...", "slices=1-311 delete");
		rename("stack1-311");
		
		run("Concatenate...", "  title=reordered keep image1=[stack1-311] image2=[stack312-313]");
		run("Concatenate...", "  title=final keep image1=[reordered] image2=" + bigStack);
		print("Stack returned to original order");

		// re stitch into one image
		run("Make Montage...", "columns=25 rows=25 scale=1");
		saveAs("tiff", bleachedImages + sampleName + "-bleach-corrected-redone.tif");
		print("Bleached image saved");
		
		// subtract background
		run("Subtract Background...", "rolling=5");
		saveAs("tiff", backgroundSubImages + sampleName + "-bleach-corrected-background-subtracted-redone.tif");
		
		close("*");
	}
}

iba1Orig = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1/";
outputFolder = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched-fixed-081324/";
File.makeDirectory(outputFolder);

bleachCorrection(iba1Orig, outputFolder);

print("Script finished");