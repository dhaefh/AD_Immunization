// quantification calculations
//- [ ]  Quantify Iba1 coverage in gray matter (excluding expanded ULEX)
            // [2]/[1]*100
//- [ ]  Quantify amount of Iba1 in vessels in gray matter (colocalizing with expanded ULEX)
            // [4]/[1]*100
//- [ ]  Quantify vascular density in gray matter
            // [5]/[1]*100
//- [ ]  Quantify vascular amyloid coverage in gray matter
            // [6]/[1]*100
//- [ ]  Quantify cortical amyloid coverage in gray matter
            // [7]/[1]*100
//- [ ]  Quantify amount of Iba1 overlapping cortical amyloid in gray matter
            // [9]/[2]*100 => percent of microglia that are overlapping cortical amyloid
            // [9]/[6]*100 => percent of cortical amyloid that have microglia
//- [ ]  Quantify amount of Iba1 overlapping vascular amyloid in gray matter
            // [8]/[2]*100 => percent of microglia that are overlapping vascular amyloid
            // [8]/[7]*100  => percent of vascular amyloid that have microglia
//- [ ]  Quantify percent of expanded ULEX covered by amyloid in gray matter
            // [6]/[3]*100 => vascular amyloid is the defined as amyloid colocalized with expanded ULEX in the gray matter
//- [ ]  Quantify microglial recruitment to vascular amyloid in gray matter
            // [11]/[2]*100 => percentage of total microglia being recruited
//- [ ]  Quantify microglial recruitment to cortical amyloid in gray matter
            // [12]/[2]*100 => percentage of total microglia being recruited
            
// values needed
    // 1  gray matter area
    // 2  iba1 within gray /Users/ksd0849/b1169/alex/scripts/AN1792-rebuttal/hdst/hdst-size-correction.shmatter area
    // 3  expanded ulex in gray matter area
    // 4  iba1 colocalizing with expanded ULEX in gray matter area
    // 5  binary cleaned ULEX within gray matter area
    // 6  vascular amyloid within gray matter area
    // 7  cortical amyloid within gray matter area
    // 8  iba1 colocalizing with vascular amyloid in gray matter area
    // 9  iba1 colocalizing with cortical amyloid in gray matter area
    // 10 iba1 colocalizing with outermost layer of vascular amyloid distance map in gray matter area
          // clear actual vascular amyloid signal
    // 11 iba1 colocalizing with outermost layer of cortical amyloid distance map in gray matter area
          // clear actual cortical amyloid signal
          
function measuringValuesNeededForQuantification () {
  // roi manager contents
        // (0) gray matter area
        // (1) iba1 within gray matter area
        // (2) expanded ULEX within gray matter area
        // (3) vascular amyloid within gray matter area
        // (4) cortical amyloid within gray matter area
        
  traces = getFileList(grayMatterTraces);
  iba1Images = getFileList(cleanedBinaryIba1);
  regularULEXImages = getFileList(cleanedBinaryULEX);
  expandedULEXImages = getFileList(expandedULEX);
  amyloidImages = getFileList(amyloid);
  amyloidDMImages = getFileList(amyloidExpanded);
  
  for (i = 0 ; i < 20; i++) {
    // 1  gray matter area
    open(grayMatterTraces + traces[i]);
    currentSample = getTitle();
    sampleName = currentSample.substring(0, 12);
    print("Processing " + sampleName);
    rename(sampleName + " grayMatterTrace");
    
    // 3  expanded ulex in gray matter area
    selectImage(sampleName + " grayMatterTrace");
    open(expandedULEX + expandedULEXImages[i]);
    rename(sampleName + " gray matter area");
    run("Restore Selection");
    run("Measure"); // measures gray matter area
    print("Measured gray matter");
    run("Clear Outside"); // removes all non-gray matter ulex signal
    run("Select None");
    rename(sampleName + " gray matter expanded ulex signal");
    run("Create Selection");
    run("Measure"); // measures gray matter expanded ulex coverage
    print("Measured gray matter expanded ulex coverage");
    
    // 5  binary cleaned ULEX within gray matter area
    //selectImage(sampleName + " grayMatterTrace");
    //open(cleanedBinaryULEX + regularULEXImages[i]);
    //run("Restore Selection");
    //run("Clear Outside"); // removes all non-gray matter ulex signal
    //run("Select None");
    //rename(sampleName + " gray matter binary un-expanded ulex signal");
    //run("Create Selection");
    //run("Measure"); // measures gray matter regular ulex (vessel) coverage
    //print("Measured gray matter binary cleaned ULEX coverage");
    
    // 2  iba1 within gray matter area
    //selectImage(sampleName + " grayMatterTrace");
    //open(cleanedBinaryIba1 + iba1Images[i]);
    //run("Restore Selection");
    //run("Clear Outside"); // removes all non-gray matter non-ulex iba1 signal
    //run("Select None");
    //rename(sampleName + " gray matter iba1 coverage");
    //run("Create Selection");
    //run("Measure"); // measures gray matter iba1 coverage
    //print("Measured gray matter iba1 coverage");
    
    // 4  iba1 colocalizing with expanded ULEX in gray matter area
    //selectImage(sampleName + " gray matter expanded ulex signal");
    //rename(sampleName + " iba1 coverage of expanded ulex signal in gray matter");
    //run("Select None");
    //selectImage(sampleName + " gray matter iba1 coverage");
    //run("Select None");
    //run("Create Selection");
    //selectImage(sampleName + " iba1 coverage of expanded ulex signal in gray matter");
    //run("Restore Selection");
    //run("Clear Outside"); // removes all iba1 not colocalizing with gray matter expanded ulex signal
    //run("Select None");
    //run("Create Selection");
    //run("Measure");
    //print("Measured gray matter iba1 colocalizing with expanded ULEX");
    //selectImage(sampleName + " iba1 coverage of expanded ulex signal in gray matter");
    //close();
    
    //  6  amyloid within gray matter area
    selectImage(sampleName + " grayMatterTrace");
    open(amyloid + amyloidImages[i]);
    run("Restore Selection");
    run("Clear Outside"); // removes all non-gray matter amyloid signal
    run("Select None");
    rename(sampleName + " gray matter amyloid signal");
    run("Create Selection");
    run("Measure"); // measures gray matter amyloid coverage
    print("Measured gray matter amyloid coverage");
        
    // 4  amyloid colocalizing with expanded ULEX in gray matter area
    selectImage(sampleName + " gray matter expanded ulex signal");
    rename(sampleName + " amyloid coverage of expanded ulex signal in gray matter");
    run("Select None");
    selectImage(sampleName + " gray matter amyloid signal");
    run("Select None");
    run("Create Selection");
    selectImage(sampleName + " amyloid coverage of expanded ulex signal in gray matter");
    run("Restore Selection");
    run("Clear Outside"); // removes all amyloid not colocalizing with gray matter expanded ulex signal
    run("Select None");
    run("Create Selection");
    run("Measure");
    print("Measured gray matter amyloid colocalizing with expanded ULEX");
    selectImage(sampleName + " amyloid coverage of expanded ulex signal in gray matter");
    close();
    
    // 8  iba1 colocalizing with amyloid in gray matter area
    //selectImage(sampleName + " gray matter amyloid signal");
    //open(cleanedBinaryIba1 + iba1Images[i]);
    //rename(sampleName + " iba1 coverage of amyloid in gray matter");
    //run("Restore Selection");
    //run("Clear Outside"); // removes all iba1 signal not colocalizing with gray matter amyloid pixels
    //run("Select None");
    //run("Create Selection");
    //run("Measure"); // measures iba1 coverage of amyloid in gray matter
    //print("Measured gray matter iba1 colocalizing with amyloid");
    //close();
        
    // 11 iba1 colocalizing with outermost layer of amyloid distance map in gray matter area
          // clear actual amyloid signal
    //selectImage(sampleName + " grayMatterTrace");
    //open(amyloid + amyloidImages[i]);
    //run("Restore Selection");
    //run("Clear Outside"); // removes all non-gray matter amyloid signal
    //run("Select None");
    //rename(sampleName + " gray matter amyloid signal");
    //run("Create Selection");
    //open(amyloidExpanded + amyloidDMImages[i]);
    //run("Restore Selection");
    //run("Clear");
    //run("Select None");
    //rename(sampleName + " iba1 coverage neighboring gray matter amyloid");
    //selectImage(sampleName + " gray matter iba1 coverage");
    //selectImage(sampleName + " iba1 coverage neighboring gray matter amyloid");
    //run("Restore Selection");
    //run("Clear Outside");
    //run("Select None");
    //run("Create Selection");
    //run("Measure"); // measures iba1 coverage of the expanded amyloid signal (excluding actual amyloid signal) in gray matter
    //print("Measured gray matter iba1 neighboring amyloid");
    
    saveAs("Results", output + sampleName + "-results.csv");
    close("*");
    close("Results");
  }
}

grayMatterTraces = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/gray-matter-with-vasculature/";
cleanedBinaryIba1 = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/final-cleaned-binary-iba1/";
cleanedBinaryULEX = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/ulex-thresholded-filtered/un-re-stitched/";
expandedULEX = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/ulex-distance-map/un-re-stitched/";
amyloid = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-manual-threshold-cleaned/un-re-stitched/";
amyloidExpanded = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/amyloid-manual-threshold-cleaned-expanded/un-re-stitched/";
output = "/projects/b1042/Gate_Lab/alex/an1792-rebuttal-hdst/all-stitched/finals/quantification-results-091124/";
File.makeDirectory(output);

measuringValuesNeededForQuantification();
print("Script finished");