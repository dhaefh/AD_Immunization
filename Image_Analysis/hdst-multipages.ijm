alignment = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/alignments/first-batch/";
dapi = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/dapi-background-subtracted/first-batch/";
amyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-manual-threshold-cleaned/first-batch/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/multipage-tifs/";

alignList = getFileList(alignment);
dapiList = getFileList(dapi);
amyloidList = getFileList(amyloid);

for (i = 0; i < alignList.length; i++) {
	open(alignment + alignList[i]);
	run("Enhance Contrast", "saturated=0.35");
	run("Invert");
	first = getTitle();
	
	open(dapi + dapiList[i]);
	run("8-bit");
	second = getTitle();
	
	open(amyloid + amyloidList[i]);
	run("8-bit");
	third = getTitle();
	
	run("Merge Channels...", "c1=[" + first + "] c2=[" + second + "] c3=[" + third + "] create");
	
	sampleName = first.substring(10, first.length - 14);
	print(sampleName);
	saveAs("tif", output + sampleName + "-alignment-backgroundSubtractedDAPI-binaryAmyloid");
	close("*");
}