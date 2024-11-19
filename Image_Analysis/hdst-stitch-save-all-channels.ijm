// Alex Edwards, 7-10-2024
// Description: Takes in individual tiles from EVOS imaging for multiple channels for one sample
//				Duplicates and organizes files by channel
//				Stitches individual channels using Grid/Collection stitching plugin
//					Duplicates and renames each tile with its relative position in the grid
//						(1, 1) is the top left corner
//					If requested for a specific channel, will perform bleach correction on tiles before stitching
//						If bleach correction is performed, metadata from the original tiles are transferred to the bleached tiles
//					If requested for a specific channel, will perform background subtraction on image after stitching
// 						with the desired sigma value		
//				Saves all versions of each stitched image

// Checks if an array contains a certain value
// Returns the index of the value in the array if so
// Returns 500 if not
function contains(array, value) {
    for (arrayNum=0; arrayNum < array.length; arrayNum++) 
        if ( array[arrayNum] == value ) {
        	return arrayNum;
        	print(arrayNum);
        }
    return 500;
}

// Returns the array a (full of strings) as a new array full of integers
function extract_digits(strArray) {
	arr2 = newArray; // return array containing digits
	for (strIndex = 0; strIndex < strArray.length; strIndex++) {
		str = "" + strArray[strIndex];
		digits = "";
		for (strChar = 0; strChar < str.length; strChar++) {
			ch = str.substring(strChar, strChar+1);
			if(!isNaN(parseInt(ch)))
				digits += ch;
		}
		arr2[strIndex] = parseInt(digits);
	}
	return arr2;
}
	
// Creates two arrays, one for positive numbers and one for negatives
// Takes values from array passed as parameter
// Orders both arrays smallest to largest
// Concatenates them to put them in numerical order
function splitPosAndNeg(origArray) {
	neg = newArray();
	pos = newArray();
	
	for (i = 0; i < origArray.length; i++) {
		element = parseInt(origArray[i]);
		if (element < 0) {
			neg = Array.concat(neg, element);
		}
		else {
			pos = Array.concat(pos, element);
		}
	}
	
	// sort positives
	arr_numPos = extract_digits(pos);
	Array.sort(arr_numPos, pos);
	
	// sort negatives and invert
	arr_numNeg = extract_digits(neg);
	Array.sort(arr_numNeg, neg);
	Array.reverse(neg);
	
	finalSorted = Array.concat(neg, pos);
	return finalSorted;
}

// Returns the coordinate value for the designated axis: x or y
// Helper function for establishesTilePositions()
function getCoord(axis) {
	coordStart = indexOf(originalInfo, "Position" + axis + "=");
	coord = originalInfo.substring(coordStart + 11, coordStart + 28);
	
	numCharsBeforeSpace = 0;
	for (c = 0; c < coord.length; c++) {
	    if (fromCharCode(charCodeAt(coord, c)) == ".") {
	        break;
	    }
	    numCharsBeforeSpace++;
	}

	coord = coord.substring(0, numCharsBeforeSpace);
	return coord;
}

// Copies original image's metadata onto bleached image
// Allows for coordinates to be pulled from bleached images for stitching
function transferMetadata(originalFolder, bleachedFolder, channel) {
	origFL = getFileList(originalFolder);
	print(origFL.length);
	firstElement = origFL[0];
	
	numCharsBeforeF = 0;
	for (c = 0; c < firstElement.length; c++) {
	    if (fromCharCode(charCodeAt(firstElement, c)) == "f") {
	        break;
	    }
	    numCharsBeforeF++;
	}
		
	imageNameStart = firstElement.substring(0, numCharsBeforeF + 1);
	
	for (i = 0; i < origFL.length; i++) {
		if (i < 10) {
			num = "0" + i;
		}
		else {
			num = i;
		}
		
		open(bleachedFolder + imageNameStart + num + channel + ".tif");
		bleached = getTitle();
		
		open(originalFolder + imageNameStart + num + channel + ".TIF");
		original = getTitle();
		originalInfo = getMetadata("Info");
		
		selectImage(bleached);
		setMetadata("Info", originalInfo);
		
		selectImage(bleached);
		run("Show Info...");
		
		saveAs("Tif", bleachedFolder + imageNameStart + num + channel + ".tif");
		close("*");
		close("Info for " + bleached);
		close("Info for " + original);
	}
}

// Bleaches images for one channel
// Helper function for stitchChannels()
function bleaching(channelDir, mainDir, channelName) {
	imagesList = getFileList(channelDir);
	print("Started bleaching");
	
	for (i = 0; i < imagesList.length; i++) {
		open(channelDir + imagesList[i]);
	}
	run("Images to Stack", "use");
	run("8-bit");
	
	//File.openSequence(channelDir);
	print("Opened images as stack");
	
	bleachedImages = mainDir + channelName + "-bleached/";
	File.makeDirectory(bleachedImages);
	print("Made " + channelName + " bleached output folder");
	
	numSlices = nSlices;
	
	// standardize signal to noise in all panels
	run("Bleach Correction", "correction=[Histogram Matching]");
	print("Bleaching done");
	
	run("Stack to Images");
	print("Split stack into individual images");
	
	// save all images with same name as original file in bleached output folder
	for (panel = 0; panel < numSlices; panel++) {
		saving = getTitle();
		saveAs("Tiff", bleachedImages + saving);
		close();
	}
	print("Saved all images in bleached output folder with same original file name");
	close("*");
}

// Determines x and y coordinates of tiles
// Helper function for stitchingChannels()
function detemineCoord(fileList, coordinate) {
	values = newArray();
	
	for (i = 0; i < fileList.length; i++) {
		coord = fileList[i];
		
			if (coord.substring(coord.length - 4, coord.length) == ".tif") {
			
			numCharsBeforeComma = 0;
			for (c = 0; c < coord.length; c++) {
			    if (fromCharCode(charCodeAt(coord, c)) == ",") {
			        break;
			    }
			    numCharsBeforeComma++;
			}
			// finding x coord
			if (coordinate == "x") {
				xco = coord.substring(1, numCharsBeforeComma);
				
				if (contains(values, xco) == 500) {
					values = Array.concat(values, xco);
				}
			}
			
			// finding y coord
			if (coordinate == "y") {
				numCharsBeforeEndPar = 0;
				for (c = 0; c < coord.length; c++) {
				    if (fromCharCode(charCodeAt(coord, c)) == ")") {
				        break;
				    }
				    numCharsBeforeEndPar++;
				}
				
				yco = coord.substring(numCharsBeforeComma + 2, numCharsBeforeEndPar);
				
				if (contains(values, yco) == 500) {
					values = Array.concat(values, yco);
				}
			}
		}
	}
	
	Array.getStatistics(values, min, max, mean, stdDev);
	size = max;
	return size;
}

// Renames images for one channel as tile positions
function establishesTilePositions(tilesDir, outputDir, sampleName, channel) {
	tilesFL = getFileList(tilesDir);

	outputFolder = outputDir + sampleName + "-" + channel + "-renamed-w-coords/";
	File.makeDirectory(outputFolder);
	
	xValues = newArray();
	yValues = newArray();
	
	// establishing possible x and y coordinates
	for (i = 0; i < tilesFL.length; i++) {
		current = tilesFL[i];
		open(tilesDir + current);

		
		originalInfo = getMetadata("Info");
		
		xcoord = getCoord("X");
		ycoord = getCoord("Y");
		
		print(current + ": (" + xcoord + ", " + ycoord + ")");
		
		if (contains(xValues, xcoord) == 500) {
			//print(contains(xValues, currentX));
			xValues = Array.concat(xValues, xcoord);
		}
		
		if (contains(yValues, ycoord) == 500) {
			//print(contains(yValues, currentY));
			yValues = Array.concat(yValues, ycoord);
		}
		close();
		close("Info for " + current);
	}		
	// sort x array
	xValues = splitPosAndNeg(xValues);
	//Array.print(xValues);
	
	yValues = splitPosAndNeg(yValues);
	//Array.print(yValues);
		
	// renaming images with tile position
	for (tile = 0; tile < tilesFL.length; tile++) {
		open(tilesDir + tilesFL[tile]);
		current = getTitle();
		print(current);
		
		originalInfo = getMetadata("Info");
		
		xcoord = getCoord("X");
		ycoord = getCoord("Y");
		
		tileX = contains(xValues, xcoord) + 1;
		tileY = contains(yValues, ycoord) + 1;
		
		saveAs("tif", outputFolder + "(" + tileX + ", " + tileY + ")");
		close();
		
		print("(" + tileX + ", " + tileY + ")");
		print("");
	}
}

// Performs bleach correction on iba1 tiles
// Stitches tiles into a full image
// Performs background subtraction on dapi, iba1, and ULEX channels
// Saves images
function bleachingStitchingRollingBall(channel) {
	dvalue = channel - 1;
	fileSuffix = "-stitched";
	
	//if (channel == 2) { // for iba1 images, perform bleach correction
	//	bleaching(sortedSampleFolder + "ch" + channel + "/", sortedSampleFolder, "ch" + channel);
	//	print("Bleaching done");
	//	transferMetadata(sortedSampleFolder + "ch" + channel + "/", sortedSampleFolder + "ch" + channel + "-bleached/", "d" + dvalue);
	//	print("Metadata transferred");
	//	fileSuffix = "-bleachcorrected" + fileSuffix;
	//}
	
	establishesTilePositions(sortedSampleFolder + "ch" + channel + "/", sortedSampleFolder, sampleName, "ch" + channel);
	print("Images named with tile positions");
	
	// get grid size
	tiledImages = sortedSampleFolder + sampleName + "-ch" + channel + "-renamed-w-coords/";
	print(tiledImages);
	coordImagesList = getFileList(tiledImages);
	xSize = detemineCoord(coordImagesList, "x");
	ySize = detemineCoord(coordImagesList, "y");
	
	run("Grid/Collection stitching", "type=[Filename defined position] order=[Defined by filename         ] grid_size_x=" + xSize + " grid_size_y=" + ySize + " tile_overlap=20 first_file_index_x=1 first_file_index_y=1 directory=[" + tiledImages + "] file_names=[({x}, {y}).tif] output_textfile_name=" + sampleName + "_coord.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ignore_z_stage subpixel_accuracy computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
	print("Stitched ch" + channel);
	saveAs("tif", sortedSampleFolder + sampleName + "-ch" + channel + fileSuffix);
	print("Saved " + sampleName + fileSuffix);
	
	// background subtraction for all channels except amyloid
	// values set for HDST samples
	if (channel == 1) { // dapi
		run("Subtract Background...", "rolling=25");
		print("Performed background subtraction on DAPI");
		fileSuffix = fileSuffix + "-rollingball";
	}
	else if (channel == 2) { // iba1
		run("Subtract Background...", "rolling=5");
		print("Performed background subtraction on Iba1");
		fileSuffix = fileSuffix + "-rollingball";
	}
	else if (channel == 4) { // ulex
		run("Subtract Background...", "rolling=20");
		print("Performed background subtraction on ULEX");
		fileSuffix = fileSuffix + "-rollingball";
	}
	
	run("8-bit");
	saveAs("tif", sortedSampleFolder + sampleName + "-ch" + channel + fileSuffix);
	print("Saved " + sampleName + fileSuffix);
	print("");
	close("*");

}

// sorts a folder that contains all exported files for a single image
// creates subfolders for each channel as well as the .stitch files
// keeps original files untouched; instead copies all files into a duplicate sorted folder
function reorganizingFiles(sampleFolder) {
	print("Inside function");
	imageList = getFileList(sampleFolder);
	print(imageList.length);
		
	// iterates through every file
	// copies it to the appropriate subfolder
	for (i = 0; i < imageList.length; i++) {
		//print("Inside for loop");
		currentImage = imageList[i];
		currentPath = sampleFolder + currentImage;
		
		if (indexOf(currentImage, "TR") != -1) {
			File.copy(currentPath, sortedSampleFolder + currentImage);
		}
		else if (indexOf(currentImage, "TM") != -1) {
			File.copy(currentPath, sortedSampleFolder + currentImage);
		}
		else if (currentImage.endsWith(".stitch")) {
			File.copy(currentPath, stitch + currentImage);
		}
		else if (indexOf(currentImage, "d0") != -1) {
			if (indexOf(currentImage, "_R_") != -1) {
				File.copy(currentPath, ch1 + currentImage);
			}
			else {
				File.copy(currentPath, merged + currentImage);
			}
		}
		else if (indexOf(currentImage, "d1") != -1) {
			File.copy(currentPath, ch2 + currentImage);
		}
		else if (indexOf(currentImage, "d2") != -1) {
			File.copy(currentPath, ch3 + currentImage);
		}
		else if (indexOf(currentImage, "d3") != -1) {
			File.copy(currentPath, ch4 + currentImage);
		}
		print("Moved " + currentImage);
	}
}

//// RUNNING THE FUNCTIONS

// takes the parameter passed in the bash script (path to the original sample folder)
args = getArgument();
print(args);
File.isDirectory(args);
File.isFile(args);
File.makeDirectory(args);

items = getFileList(args);
print(items.length);

// make sorted folder
indexOfLastForwardSlash = lastIndexOf(args, "/");
secondLastForwardSlash = lastIndexOf(args.substring(0, indexOfLastForwardSlash), "/");
basePath = args.substring(0, secondLastForwardSlash);
sortedFolder = basePath + "/sorted-082124-FINAL/";
File.makeDirectory(sortedFolder);

// make sample specific sorted folder
sampleName = args.substring(secondLastForwardSlash + 1, args.length - 21); // removes the time and date info from the sample name
print("Currently processing: " + sampleName);

sortedSampleFolder = sortedFolder + sampleName + "/";
//print(sortedSampleFolder);
File.makeDirectory(sortedSampleFolder);

// make all output subfolders
ch1 = sortedSampleFolder + "ch1/";
ch2 = sortedSampleFolder + "ch2/";
ch3 = sortedSampleFolder + "ch3/";
ch4 = sortedSampleFolder + "ch4/";
stitch = sortedSampleFolder + "stitch/";
merged = sortedSampleFolder + "merged/";

File.makeDirectory(ch1);
File.makeDirectory(ch2);
File.makeDirectory(ch3);
File.makeDirectory(ch4);
File.makeDirectory(stitch);
File.makeDirectory(merged);

reorganizingFiles(args);
print("Files copied to new folder and sorted");

bleachingStitchingRollingBall(1);
bleachingStitchingRollingBall(2);
bleachingStitchingRollingBall(3);
bleachingStitchingRollingBall(4);
close("*");