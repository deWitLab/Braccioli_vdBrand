//Macro for Noemi: analysis of intensity; it includes measurements of the overlapping areas (colocalization)
//Date: 2023-09-01
//Modified the previous version to fit the new files: Bf is ch3 and we ignore it.  

//This macro segments the areas postive for the different stainings and measures the intentisity and area of each staining in all the channels. 

//With Noise Correction

#@ File (label="Choose the Source Directory", description="Source directory", style="directory") dir1
#@ File (label="Choose the Destination Directory", description="Destination directory", style="directory") dir2
#@ String (label="Name of the results file:", style="text field") Results_name
#@ Integer (label="Lower Threshold (whole gastruloid)", min=0, max=255, value=25) lower_whole
#@ Integer (label="Upper Threshold (whole gastruloid)", min=0, max=255, value=255) upper_whole
#@ boolean(label = "Ch1") Ch1
#@ String (label="Staining Chanel 1", style="text field") Ch1_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch1_LUT
#@ String (label="ROI_color", choices={"red", "green", "yellow", "magenta", "blue", "cyan", "grays"}, style="listBox") Ch1_ROI_color
#@ Integer (label="Lower Threshold (Ch1)", min=0, max=255, value=25) lower_ch1
#@ Integer (label="Upper Threshold (Ch1)", min=0, max=255, value=255) upper_ch1
#@ boolean(label = "Ch2") Ch2
#@ String (label="Staining Chanel 2", style="text field") Ch2_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch2_LUT
#@ String (label="ROI_color", choices={"red", "green", "yellow", "magenta", "blue", "cyan", "grays"}, style="listBox") Ch2_ROI_color
#@ Integer (label="Lower Threshold (Ch2)", min=0, max=255, value=25) lower_ch2
#@ Integer (label="Upper Threshold (Ch2)", min=0, max=255, value=255) upper_ch2
#@ boolean(label = "Ch3") Ch3
#@ String (label="Staining Chanel 3", style="text field") Ch3_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch3_LUT
#@ String (label="ROI_color", choices={"red", "green", "yellow", "magenta", "blue", "cyan", "grays"}, style="listBox") Ch3_ROI_color
#@ Integer (label="Lower Threshold (Ch3)", min=0, max=255, value=25) lower_ch3
#@ Integer (label="Upper Threshold (Ch3)", min=0, max=255, value=255) upper_ch3
#@ boolean(label = "Ch4") Ch4
#@ String (label="Staining Chanel 4", style="text field") Ch4_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch4_LUT
#@ String (label="ROI_color", choices={"red", "green", "yellow", "magenta", "blue", "cyan", "grays"}, style="listBox") Ch4_ROI_color
#@ Integer (label="Lower Threshold (Ch4)", min=0, max=255, value=25) lower_ch4
#@ Integer (label="Upper Threshold (Ch4)", min=0, max=255, value=255) upper_ch4

function contains( array, value ) {
    for (i=0; i<array.length; i++) 
        if ( array[i] == value ) return true;
    return false;
}

function Masking_zstack(Window_name, lower_th, upper_th) { 
	// Creates a binary mask from a channel with a z-stack 
	selectWindow(Window_name);
	run("Duplicate...", "duplicate");
	setAutoThreshold("Otsu dark");
	setThreshold(lower_th, upper_th);
	run("Convert to Mask", "method=Otsu background=Dark black");
	rename("Mask_" + Window_name);
}

function Merge_Masking_zstack(Window_name,lower_th_m, upper_th_m) { 
	// Creates a binary mask from a channel with a z-stack 
	selectWindow(Window_name);
	run("Duplicate...", "duplicate");
	setAutoThreshold("Otsu dark");
	setThreshold(lower_th_m, upper_th_m);
	run("Convert to Mask", "method=Otsu background=Dark calculate black create");
	rename("Mask_" + Window_name);
}

function ROI_selection(Mask, Staining, ROI, ROI_color, slices, id, ch) { 
	// This function iterates through the slices of the called image, applies a supervised threshold and creates a ROI from the area. 
	// It sets the ROI line widht to 1 and it changes the ROI color to the one called in the function. finally it renames the ROI refering back to the staining and the slice. 
	selectWindow(Mask);
	for (s = 1; s <= slices; s++) {
		Stack.setSlice(s);
		run("Create Selection");
		type = selectionType(); 
		if (type == -1) {
			print("No selection");
			setResult("Label", nResults, id + ":" + Staining + "_area_slice_" + s);
			n=nResults;
			setResult("Area", n-1, 0);
			setResult("Mean", n-1, 0);
			setResult("StdDev", n-1, 0);
			setResult("Min", n-1, 0);
			setResult("Max", n-1, 0);
			setResult("Ch", n-1, ch);
			setResult("Slice", n-1, s);
			updateResults();
		}
		else {
		roiManager("add");
		roiManager("Select", ROI);
		roiManager("Set Line Width", 1);
		roiManager("Set Color", ROI_color);
		roiManager("Rename", Staining + "_area_slice_" + s);
		ROI=ROI+1;
		}
	}
}


function Colocalization(Mask_1, Mask_2, Staining_1, Staining_2, ROI, slices, id, ch) {
	// This function takes two binary masks as input and calculates the common area ("AND") between them. 
	// Then it stores the common areas as a ROI (refered as ovlap_area).
	imageCalculator("AND create stack", Mask_1, Mask_2);
	rename(Staining_1 + "_vs_" + Staining_2 + "_" + id);

	for (s = 1; s <= slices; s++) {
		selectWindow(Staining_1 + "_vs_" + Staining_2 + "_" + id);
		Stack.setSlice(s);
		run("Create Selection");
		type = selectionType();
		if (type == -1) {
			setResult("Label", nResults, id + ":" + Staining_1 + "_" + Staining_2 + "_ovlap_area_" + s);
			n=nResults;
			setResult("Area", n-1, 0);
			setResult("Mean", n-1, 0);
			setResult("StdDev", n-1, 0);
			setResult("Min", n-1, 0);
			setResult("Max", n-1, 0);
			setResult("Ch", n-1, ch);
			setResult("Slice", n-1, s);
			updateResults();
		}
		else {
		roiManager("add");
		roiManager("Select", ROI);
		roiManager("Rename", Staining_1 + "_" + Staining_2 + "_ovlap_area_" + s);
		ROI=ROI+1;
		}
	}
}



//This creates a list with all the files in that folder
//We will use it to iterate through it and open the files that need analysis
list = getFileList(dir1);

//This hides images, you can not follow what the macro is doing till the end but it prevents the images from
//flashing on your screen. Consider deleting it out depending on your purposes. 
setBatchMode(true);

//Select the different parameters that we want to measure; select shape descriptors for gastruloid morphology characterization
run("Set Measurements...", "area mean standard integrated display min max redirect=None decimal=3");

//Close and reset everything so we can start with a white page
run("Clear Results");
close("Empty_areas");
roiManager("reset"); 
run("Close All");
run("Conversions...", " ");
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], ".tif")){
		path = dir1 + "/" +list[i];
		run("Bio-Formats Importer", "open=[" + path  + "]");
		Stack.getDimensions(width, height, channels, slices, frames);
		for (s = 0; s <= slices; s++) {
			Stack.setSlice(s);
			resetMinAndMax();
		}
		rename(i);
		original = getTitle();
		
		//setBatchMode("show"); //This will only show the original image, so we can follow part of what it is happening 
		run("Duplicate...", "title=[" + list[i] + "_dup] duplicate");
		run("32-bit");

		//To smooth the image and correct for the noise
		//Sigma can be adjusted depending on the images
		run("Gaussian Blur...", "sigma=1.5 stack"); 
		run("Median...", "radius=4 stack");
		run("8-bit");
		
		//Each channel corresponds to a respective staining; so first we separate them
		run("Split Channels");

		//Merge all the channels to obtain the shape of the whole gastruloid
		if (Ch3) {
			
			run("Merge Channels...", "c1=[C1-" + list[i] + "_dup] c2=[C2-" + list[i] + "_dup] c3=[C4-" + list[i] + "_dup] c4=[C3-" + list[i] + "_dup] keep");
			run("RGB Color", "slices keep");
			run("8-bit");
			rename("Merge_" + list[i]);
			Stack.getDimensions(width, height, channels, slices, frames);
		}
		else {
			run("Merge Channels...", "c1=[C1-" + list[i] + "_dup] c2=[C2-" + list[i] + "_dup] c3=[C4-" + list[i] + "_dup]  keep");
			run("RGB Color", "slices keep");
			run("8-bit");
			rename("Merge_" + list[i]);
			Stack.getDimensions(width, height, channels, slices, frames);
		}

		// Create masks from all the channels. This masks will be used to calculate the overlapping areas among the channels. 
		Merge_Masking_zstack("Merge_" + list[i], lower_whole, upper_whole);
		Masking_zstack("C1-" + list[i] + "_dup", lower_ch1, upper_ch1);
		Masking_zstack("C2-" + list[i] + "_dup", lower_ch2, upper_ch2);
		Masking_zstack("C4-" + list[i] + "_dup", lower_ch4, upper_ch4);
		
		if (Ch3) {
			Masking_zstack("C3-"+list[i] + "_dup", lower_ch3, upper_ch3);
		}
		
		// Make selections for the postive areas in the merge image and all the channels slice by slice
		activeROIs = roiManager("count");
		ROI_selection("Mask_"  + "Merge_" + list[i], "Total", activeROIs, "white", slices, i, 0);
		activeROIs = roiManager("count");
		ROI_selection("Mask_" + "C1-" + list[i] + "_dup", Ch1_st, activeROIs, Ch1_ROI_color, slices, i, 1);
		activeROIs = roiManager("count");
		ROI_selection("Mask_" + "C2-" + list[i] + "_dup", Ch2_st, activeROIs, Ch2_ROI_color, slices, i, 2);
		activeROIs = roiManager("count");
		ROI_selection("Mask_" + "C4-" + list[i] + "_dup", Ch4_st, activeROIs, Ch4_ROI_color, slices, i, 3);
		activeROIs = roiManager("count");
	
		if (Ch3) {
			ROI_selection("Mask_" + "C3-"+list[i] + "_dup", Ch3_st, activeROIs, Ch3_ROI_color, slices, i, 4);
			activeROIs = roiManager("count");
		}
		
		//Ch1 vs Ch2
		Colocalization("Mask_" + "C1-"+list[i] + "_dup", "Mask_" + "C2-"+list[i] + "_dup", Ch1_st, Ch2_st, activeROIs, slices, i, 0);
		activeROIs = roiManager("count");

		//Ch1 vs Ch4
		Colocalization("Mask_" + "C1-"+list[i] + "_dup","Mask_" + "C4-"+list[i] + "_dup", Ch1_st, Ch4_st, activeROIs, slices, i, 0);
		activeROIs = roiManager("count");
		
		//Ch1 vs Total
		Colocalization("Mask_" + "C1-"+list[i] + "_dup", "Mask_" + "Merge_" + list[i], Ch1_st, "Total", activeROIs, slices, i, 0);
		activeROIs = roiManager("count");
		
		//Ch2 vs Ch4
		Colocalization("Mask_" + "C2-"+list[i] + "_dup","Mask_" + "C4-"+ list[i] + "_dup", Ch2_st, Ch4_st, activeROIs, slices, i, 0);
		activeROIs = roiManager("count");
		
		//Ch2 vs Total
		Colocalization("Mask_" + "C2-"+list[i] + "_dup", "Mask_" + "Merge_" + list[i], Ch2_st, "Total", activeROIs, slices, i, 0);
		activeROIs = roiManager("count");
		
		//Ch4 vs Total
		Colocalization("Mask_" + "C4-"+list[i] + "_dup", "Mask_" + "Merge_" + list[i], Ch4_st, "Total", activeROIs, slices, i, 0);
		activeROIs = roiManager("count");
	

		//If Ch3
		if (Ch3) {
			//Ch3 vs Ch1
			Colocalization("Mask_" + "C3-"+list[i] + "_dup","Mask_" + "C1-"+list[i] + "_dup", Ch3_st, Ch1_st, activeROIs, slices, i, 0);
			activeROIs = roiManager("count");
	
			//Ch3 vs Ch2
			Colocalization("Mask_" + "C3-"+list[i] + "_dup","Mask_" + "C2-"+list[i] + "_dup", Ch3_st, Ch2_st, activeROIs, slices, i, 0);
			activeROIs = roiManager("count");
	
			//Ch3 vs Ch4
			Colocalization("Mask_" + "C4-"+list[i] + "_dup","Mask_" + "C3-"+list[i] + "_dup", Ch4_st, Ch3_st, activeROIs, slices, i, 0);
			activeROIs = roiManager("count");
			
			//Ch3 vs Total
			Colocalization("Mask_" + "C3-"+list[i] + "_dup", "Mask_" + "Merge_" + list[i], Ch3_st, "Total", activeROIs, slices, i, 0);
			activeROIs = roiManager("count");
			}
			
		if (isOpen("Results")) {
			Table.rename("Results", "Empty_areas");
		}

		//Now we measure all the ROIs in every channel
		selectWindow(original);
		roiManager("Deselect");
		roiManager("multi-measure measure_all append");	
		
		//save the Empty_areas table
		if (isOpen("Empty_areas")) {
		selectWindow("Empty_areas");
		Table.save(dir2 + "/" + Results_name +  "_" + i +  "_Empty_areas.csv");
		close("Empty_areas");	
		}
		//save the results table
		selectWindow("Results");
		saveAs("Results",dir2 + "/" + Results_name + "_"+ i + ".csv");
		run("Close");
		
		//Next we set the LUTs
		Stack.setChannel(1);
		run(Ch1_LUT);
		run("Enhance Contrast", "saturated=0.35");
		Stack.setChannel(2);
		run(Ch2_LUT);
		run("Enhance Contrast", "saturated=0.35");
		Stack.setChannel(3);
		run(Ch4_LUT);
		run("Enhance Contrast", "saturated=0.35");
	
		 if (Ch3){
		 	Stack.setChannel(4);
			run(Ch3_LUT);
			run("Enhance Contrast", "saturated=0.35");
		 }

		//Make composite
		Stack.setDisplayMode("composite");
		roiManager("Show All without labels");
		saveAs("Tif", dir2 + "/" + list[i] + "_analyzed");
		
		//Reset RoiManager for next gastruloid
		roiManager("deselect");
		roiManager("Save",dir2 + "/" + list[i] + "_ROIset.zip");
		roiManager("reset");
		
		//Save the merge
		selectWindow("Merge_" + list[i]);
		saveAs("Tif", dir2 + "/" + list[i] + "_merge");
		
		//Close all images windows: this prevents having windows with the same name open (e.g after merging) and prevents runing out of memory. 
		close("*"); 
	}
}


run("Close All");
//setBatchMode("exit and display");

