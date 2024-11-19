dapi = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/dapi/original/";
iba1 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1/";
amyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/amyloid/";
alignment = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/alignment-images/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/raw-multipages/";

dapiImages = getFileList(dapi);
iba1Images = getFileList(iba1);
amyloidImages = getFileList(amyloid);
alignmentImages = getFileList(alignment);

for (i = 6; i < 8; i++) {
	open(dapi + dapiImages[i]);
	dapich1 = getTitle();
	run("8-bit");
	
	sampleName = dapich1.substring(0, dapich1.length - 17);
	print(sampleName);
	
	open(alignment + alignmentImages[i]);
	run("Invert");
	run("8-bit");
	alignmentch4 = getTitle();
	
	open(iba1 + iba1Images[i]);
	iba1ch2 = getTitle();
	run("8-bit");
	
	open(amyloid + amyloidImages[i]);
	amyloidch3 = getTitle();
	run("8-bit");
	
	run("Merge Channels...", "c1=[" + alignmentch4 + "] c2=[" + dapich1 + "] create");
	print("Ran merge channels");
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(1);
	Stack.setChannel(2);
	//Stack.setChannel(3);
	
	saveAs("tif", output + sampleName + "-raw-channels-alignment-and-dapi");
	print("Saved multipage"); 
}
