// determining Iba1+ and DAPI+ density in regions overlapping amyloid, 
// neighboring amyloid, and not neighboring amyloid
// new CAA controls

function iba1DensitiesRegardingAmyloid() {
	print("Function running...");
	iba1FileList = getFileList(binaryIba1);
	grayMatterRoiFileList = getFileList(grayMatterROIs);
	amyloidFileList = getFileList(binaryAmyloid);
	outerLayerAmyloidFileList = getFileList(outerLayerAmyloidDistance);
		
	
	for (i = 0; i < iba1FileList.length; i++) {
		// open cleaned binary iba1
		open(binaryIba1 + iba1FileList[i]);
	    currentIba1 = getTitle();
	    sample = iba1FileList[i];
	    sampleName = sample.substring(0, 19); // ADJUST
	    print("Processing " + sampleName + "...");
		
		// open gray matter ROI + add to manager (0)
		open(grayMatterROIs + grayMatterRoiFileList[i]);
		//roiManager("Add");
		run("Select None");
		selectImage(currentIba1);
		rename(sampleName + "-gray-matter");
		//roiManager("Select", 0);
		run("Restore Selection");
		run("Clear Outside"); // isolating gray matter area on iba1 image
		// measure value of total gray matter area
		run("Measure");
		print("Measured whole gray matter area");
		
		// measure total Iba1 coverage in gray matter
		open(binaryIba1 + iba1FileList[i]);
		rename(sampleName + "-iba1-gray-matter");
		selectImage(sampleName + "-gray-matter");
		selectImage(sampleName + "-iba1-gray-matter");
		run("Restore Selection");
		run("Clear Outside");
		run("Select None");
		run("Create Selection");
		run("Measure");
		print("Measured total area of Iba1 in gray matter");
		
		// DIRECTLY OVERLAPPING
		// open cleaned binary amyloid
		open(binaryAmyloid + amyloidFileList[i]);
	    rename(sampleName + "-cortical-binary-amyloid-in-gray-matter");
	    //roiManager("Select", 0);
	    selectImage(sampleName + "-gray-matter");
	    selectImage(sampleName + "-cortical-binary-amyloid-in-gray-matter");
	    run("Restore Selection");
	    run("Clear Outside"); // isolating gray matter area
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid present in gray matter
		//roiManager("Add");
		// measure value for total amyloid area (1)
		run("Measure");
		print("Measured total area of cortical amyloid in gray matter");
		run("Select None");
		
		rename(sampleName + "-iba1-colocalizing-with-cortical-amyloid-in-gray-matter");
		selectImage(sampleName + "-iba1-gray-matter");
		run("Select None");
		run("Create Selection");
		selectImage(sampleName + "-iba1-colocalizing-with-cortical-amyloid-in-gray-matter");
		run("Restore Selection");
		run("Clear Outside");
		run("Select None");
		run("Create Selection");
		run("Measure");
		print("Measured colocalization of iba1 and cortical amyloid in gray matter");
		
		saveAs("Results", outputDir + sampleName + "-results.csv");
		
		close("*");
		close("Results");
		print("\n");
	}
}

// run it!
//binaryIba1 = getDirectory("Navigate to folder containing the cleaned binary Iba1 images");
//grayMatterROIs = getDirectory("Navigate to gray matter traces (saved as ROIs) ");
//binaryAmyloid = getDirectory("Navigate to folder containing the cleaned binary amyloid images");
//outerLayerAmyloidDistance = getDirectory("Navigate to folder containing layer 100 from the distance maps of the cleaned binary amyloid images");
//outputDir = getDirectory("Navigate to the desired output folder ");

binaryIba1 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/final-cleaned-binary-iba1/";
grayMatterROIs = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/gray-matter-traces/";
binaryAmyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-cortical-amyloid/un-re-stitched/";
outerLayerAmyloidDistance = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/binary-cortical-amyloid-distance-maps/dilated175_layer/un-re-stitched/";
outputDir = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-fig4-remake/sorted-072424-FINAL/all-stitched/quantification-output/";
File.makeDirectory(outputDir);

iba1DensitiesRegardingAmyloid();
print("Script Finished");