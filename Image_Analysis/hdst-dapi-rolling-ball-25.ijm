path = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/dapi-need-background-subtraction/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/dapi-background-subtraction/";
File.makeDirectory(output);
dapiImages = getFileList(path);

for (i=0; i < dapiImages.length; i++){
	current = dapiImages[i];
	open(path + current);
	run("8-bit");
	run("Subtract Background...", "rolling=25");
	newName = current.substring(0, current.length-4) + "-background-subtracted";
	saveAs("tif", output + newName);
	close("*");
}

