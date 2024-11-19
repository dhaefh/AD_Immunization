dapi = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/dapi/original/";
iba1 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/iba1/";
amyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/amyloid/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/alignment-images/";

dapiImages = getFileList(dapi);
iba1Images = getFileList(iba1);
amyloidImages = getFileList(amyloid);

for (i = 0; i < dapiImages.length; i++) {
	open(dapi + dapiImages[i]);
	dapich1 = getTitle();
	run("8-bit");
	
	sampleName = dapich1.substring(0, dapich1.length - 17);
	print(sampleName);
	
	open(iba1 + iba1Images[i]);
	iba1ch2 = getTitle();
	run("8-bit");
	
	open(amyloid + amyloidImages[i]);
	amyloidch3 = getTitle();
	run("8-bit");
	
	imageCalculator("Add create", dapich1, iba1ch2);
	v1 = getTitle();
	imageCalculator("Add create", v1, amyloidch3);
	final = getTitle();
	
	saveAs("tif", output + sampleName + "-flattened");
	print("Saved alignment image");
	close("*");
}
