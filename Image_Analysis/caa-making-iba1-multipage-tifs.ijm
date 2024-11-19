alignment = "R:/Neurology/Gate_Lab/projects/AN1792-vacc/04 raw-data/07 Image analysis/rebuttal/caa-controls/alignment-images/";
iba1CleanedBinary = "R:/Neurology/Gate_Lab/projects/AN1792-vacc/04 raw-data/07 Image analysis/rebuttal/caa-controls/cleaned-binary-iba1/";
output = "R:/Neurology/Gate_Lab/projects/AN1792-vacc/04 raw-data/07 Image analysis/rebuttal/caa-controls/iba1-multipages/";
File.makeDirectory(output);

alignmentImages = getFileList(alignment);
iba1Images = getFileList(iba1CleanedBinary);

for (i = 0; i < 8; i++) {
	open(alignment + alignmentImages[i]);
	run("Invert");
	run("8-bit");
	rename("Alignment");
	alignch1 = getTitle();
	print(alignch1);
	
	sampleName = alignch1.substring(0, alignch1.length - 18);
	print(sampleName);
	
	open(iba1CleanedBinary + iba1Images[i]);
	rename("Iba1");
	iba1ch2 = getTitle();
	print(iba1ch2);
	
	run("Merge Channels...", "c1=[" + alignch1 + "] c2=[" + iba1ch2 + "] create");
	print("Ran merge channels");
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(1);
	Stack.setChannel(2);
	
	saveAs("tif", output + sampleName + "-alignment-and-iba1-multipage");
	print("Saved multipage"); 
}
