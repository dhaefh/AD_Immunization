// fixing improperly stitched iba1 images after bleaching

function reorderMontage(origDir, outputDir) {
	orig_FL = getFileList(origDir);
	
	bleachedImages = outputDir + "iba1-bleach-corrected/";
	backgroundSubImages = outputDir + "iba1-bleach-corrected-background-subtracted/";
	
	for (i = 0; i < orig_FL.length; i++) {
		open(origDir + orig_FL[i]);
		original = getTitle();
		sampleName = original.substring(0, original.length - 17);
		print("Processing " + sampleName);
		
		run("Montage to Stack...", "columns=25 rows=25 border=0");
		bigStack = getTitle();
		print("Made stack");
		numSlices = nSlices;
		
		run("Stack to Images");

		for (j = 0; j < numSlices; j++) {
			currentName = getTitle();
			print(currentName);
			close();
		}
		
		break
		
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
		saveAs("tiff", bleachedImages + sampleName + "-bleach-corrected-remontaged.tif");
		print("Bleached image saved");
		
		// subtract background
		run("Subtract Background...", "rolling=5");
		saveAs("tiff", backgroundSubImages + sampleName + "-bleach-corrected-background-subtracted-remontaged.tif");
		
		close("*");
	}
}
iba1Orig = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1/";
outputFolder = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/";

reorderMontage(iba1Orig, outputFolder);