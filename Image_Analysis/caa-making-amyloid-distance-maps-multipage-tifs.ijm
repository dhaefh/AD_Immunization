alignment = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/alignment-images/un-re-stitched/";
vascAmyloidDM = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-vascular-amyloid-distance-maps/un-re-stitched/";
corticalAmyloidDM = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-cortical-amyloid-distance-maps/un-re-stitched/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/amyloid-multipages/";
File.makeDirectory(output);

alignmentImages = getFileList(alignment);
vascImages = getFileList(vascAmyloidDM);
cortImages = getFileList(corticalAmyloidDM);

for (i = 0; i < 1; i++) {
	open(alignment + alignmentImages[i]);
	run("Invert");
	run("8-bit");
	alignch1 = getTitle();
	//print(alignch1);
	getDimensions(width, height, channels, slices, frames);
	print(width + " by " + height);
	newWidth = round(width*0.5);
	newHeight = round(height*0.5);
	print(newWidth + " by " + newHeight);
	run("Scale...", "x=0.5 y=0.5 width=" + newWidth + " height=" + newHeight + " interpolation=Bilinear average create title=[alignment-downsampled.tif]");
	//print("Alignment dimensions: " + width + " x " + height);
	run("Duplicate...", "title=[alignment-duplicate-1]");
	run("Duplicate...", "title=[alignment-duplicate-2]");
	
	sampleName = alignch1.substring(0, alignch1.length - 19);
	print(sampleName);
	
	open(vascAmyloidDM + vascImages[i]);
	rename("Vasc-amyloid-DM");
	vascch2 = getTitle();
	print(vascch2);
	getDimensions(width, height, channels, slices, frames);
	print(width + " by " + height);
	newWidth = round(width*0.5);
	newHeight = round(height*0.5);
	print(newWidth + " by " + newHeight);
	run("Scale...", "x=0.5 y=0.5 width=" + newWidth + " height=" + newHeight + " interpolation=Bilinear average create title=[vascular-downsampled.tif]");
	//print("Vascular distance map dimensions: " + width + " x " + height);
	run("Duplicate...", "title=[vasc-DM-duplicate-1]");
	
	open(corticalAmyloidDM + cortImages[i]);
	rename("Cort-amyloid-DM");
	cortch3 = getTitle();
	print(cortch3);
	getDimensions(width, height, channels, slices, frames);
	print(width + " by " + height);
	newWidth = round(width*0.5);
	newHeight = round(height*0.5);
	print(newWidth + " by " + newHeight);
	run("Scale...", "x=0.5 y=0.5 width=" + newWidth + " height=" + newHeight + " interpolation=Bilinear average create title=[cortical-downsampled.tif]");
	//print("Cortical distance map dimensions: " + width + " x " + height);
	run("Duplicate...", "title=[cort-DM-duplicate-1]");
	
	run("Merge Channels...", "c1=[" + alignch1 + "] c2=[" + vascch2 + "] c3=[" + cortch3 + "] create");
	print("Ran merge channels");
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(1);
	Stack.setChannel(2);
	Stack.setChannel(3);
	
	saveAs("tif", output + sampleName + "alignment-and-amyloid-distance-maps-multipage-smaller-test");
	print("Saved alignment + cortical + vascular multipage"); 
	close("*");
	
//	run("Merge Channels...", "c1=[alignment-duplicate-1] c2=[vasc-DM-duplicate-1] create");
//	Stack.setDisplayMode("grayscale");
//	Stack.setChannel(1);
//	Stack.setChannel(2);
	
//	saveAs("tif", output + sampleName + "alignment-and-vascular-distance-map-multipage");
//	print("Saved alignment + vascular multipage"); 
	
//	run("Merge Channels...", "c1=[alignment-duplicate-2] c2=[cort-DM-duplicate-1] create");
//	Stack.setDisplayMode("grayscale");
//	Stack.setChannel(1);
//	Stack.setChannel(2);
	
//	saveAs("tif", output + sampleName + "alignment-and-cortical-distance-map-multipage");
//	print("Saved alignment + cortical multipage"); 
//	close("*");
}
