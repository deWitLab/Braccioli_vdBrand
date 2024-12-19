//Macro for Charis: cell type distribution along the gastruloid, wide field images
//Version: 20240423
//Makes a max projection when you have a 

#@ File (label="Choose the Source Directory", description="Source directory", style="directory") dir1
#@ File (label="Choose the Destination Directory", description="Destination directory", style="directory") dir2
#@ String (label="Staining Channel 1", style="text field") Ch1_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch1_LUT
#@ String (label="Staining Channel 2", style="text field") Ch2_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch2_LUT
#@ String (label="Staining Channel 3", style="text field") Ch3_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch3_LUT
#@ boolean(label = "Ch4") Ch4
#@ String (label="Staining Channel 4", style="text field") Ch4_st
#@ String (label="LUT", choices={"Red", "Green", "Yellow", "Magenta", "Blue", "Cyan", "Grays"}, style="listBox") Ch4_LUT

//setBatchMode(true);
list = getFileList(dir1);

run("Close All");
run("Clear Results");

setTool("polyline");
run("Line Width...", "line=150");

for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	if (endsWith(list[i], ".tif")){
		path = dir1 + "/" +list[i];
		run("Bio-Formats Importer", "open=[" + path  + "]");
		original = getTitle();
		run("Duplicate...", "title=[" + list[i] + "_dup] duplicate");
		duplicate=getTitle();
		Stack.getDimensions(width, height, channels, slices, frames);
		selectWindow(duplicate);
		if (slices>1) {
			run("Z Project...", "projection=[Max Intensity]");
		}

		
		rename(duplicate + "_max");
		
		//Save the individual channels and the merge
		selectWindow(duplicate + "_max");
		run("Split Channels");
		selectWindow("C1-"+list[i] + "_dup_max");
		Ch1_=getTitle();
		run(Ch1_LUT);
		run("Enhance Contrast", "saturated=0.35");
		saveAs("JPEG", dir2 + "/" + list[i] + "_" + Ch1_st);

		selectWindow("C2-"+list[i] + "_dup_max");
		Ch2_=getTitle();
		run(Ch2_LUT);
		run("Enhance Contrast", "saturated=0.35");
		saveAs("JPEG", dir2 + "/" + list[i] + "_" + Ch2_st);
		
		selectWindow("C3-"+list[i] + "_dup_max");
		Ch3_=getTitle();
		run(Ch3_LUT);
		run("Enhance Contrast", "saturated=0.35");
		saveAs("JPEG", dir2 + "/" + list[i] + "_" + Ch3_st);
		
		if (Ch4) {
			selectWindow("C4-"+list[i] + "_dup_max");
			Ch4_=getTitle();
			run(Ch4_LUT);
			run("Enhance Contrast", "saturated=0.35");
			saveAs("JPEG", dir2 + "/" + list[i] + "_" + Ch4_st);
	
			run("Merge Channels...", "c2=[" + Ch4_ + "] c3=[" + Ch3_ + "] c4=[" + Ch1_ + "] c5=[" + Ch2_ + "] create keep");
			saveAs("JPEG", dir2 + "/" + list[i] + "_" + "merge_4_chanels");
			setBatchMode("show"); //This will only show the original image, so we can follow part of what it is happening
			waitForUser("Draw a line ROI from Anterior to Posterior axis of the gastruloid.");
			roiManager("Add");
			
		}
		else {

	
		run("Merge Channels...", "c4=[" + Ch3_ + "] c5=[" + Ch1_ + "] c6=[" + Ch2_ + "] create keep");
		saveAs("JPEG", dir2 + "/" + list[i] + "_" + "merge_2_chanels");
		setBatchMode("show"); //This will only show the original image, so we can follow part of what it is happening
		waitForUser("Draw a line ROI from Anterior to Posterior axis of the gastruloid.");
		roiManager("Add");
		
		}

			//Channel 1
		selectWindow(original);
		Stack.setChannel(1);
		roiManager("Select", 0);
		run("Plots...", "width=530 height=300 font=12 draw draw_ticks list minimum=0 maximum=0 interpolate");
		run("Plot Profile");
		Plot.getValues(x1, y1);
		for (a=0; a<x1.length; a++){
		  	setResult(Ch1_st + "_Distance", a, x1[a]);
		  	setResult(Ch1_st, a, y1[a]);
		  	updateResults();
		}
		
			//Channel 2
		selectWindow(original);
		roiManager("Select", 0);
		Stack.setChannel(2);
		run("Plots...", "width=530 height=300 font=12 draw draw_ticks list minimum=0 maximum=0 interpolate");
		run("Plot Profile");
		Plot.getValues(x2, y2);
		for (b=0; b<x2.length; b++){
		  	setResult(Ch2_st + "_Distance", b, x2[b]);
		  	setResult(Ch2_st, b, y2[b]);
		  	updateResults();
		}
		
		  	//Channel 3
		selectWindow(original);
		roiManager("Select", 0);
		Stack.setChannel(3);
		run("Plots...", "width=530 height=300 font=12 draw draw_ticks list minimum=0 maximum=0 interpolate");
		run("Plot Profile");
		Plot.getValues(x3, y3);
		for (e=0; e<x3.length; e++){
		  	setResult(Ch3_st + "_Distance", e, x3[e]);
		  	setResult(Ch3_st, e, y3[e]);
		  	updateResults();
		}
		
			//Channel 4
		if (Ch4) {
		selectWindow(original);
		roiManager("Select", 0);
		Stack.setChannel(4);
		run("Plots...", "width=530 height=300 font=12 draw draw_ticks list minimum=0 maximum=0 interpolate");
		run("Plot Profile");
		Plot.getValues(x4, y4);
		for (d=0; d<x4.length; d++){
		  	setResult(Ch4_st + "_Distance", d, x4[d]);
		  	setResult(Ch4_st, d, y4[d]);
		  	updateResults();
		}
		}
		
		saveAs("Results", dir2 + "/" + list[i] + "_data.csv");
		run("Clear Results");
		roiManager("reset");
		
		
	}
}

run("Close All");
run("Clear Results");