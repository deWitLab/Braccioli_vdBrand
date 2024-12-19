//Macro for Charis: cell type distribution along the gastruloid, wide field images
//Version: 20240423
//Makes a max projection when you have a 

#@ File (label="Choose the Source Directory", description="Source directory", style="directory") dir1
#@ File (label="Choose the Destination Directory", description="Destination directory", style="directory") dir2
#@ String (label="Name of the results file:", style="text field") Results_name
#@ String (label="Staining Channel 1", style="text field") Ch1_st
#@ String (label="Staining Channel 2", style="text field") Ch2_st
#@ String (label="Staining Channel 3", style="text field") Ch3_st
#@ boolean(label = "Ch4") Ch4
#@ String (label="Staining Channel 4", style="text field") Ch4_st

setBatchMode(true);

list = getFileList(dir1);

run("Close All");
run("Clear Results")

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
		selectWindow(duplicate + "_max");
		run("Median...", "radius=2 stack");
		run("Split Channels");
		
		selectWindow("C1-"+list[i] + "_dup_max");
		rename(Ch1_st + list[i]);
		Ch1_=getTitle();
		
		selectWindow("C2-"+list[i] + "_dup_max");
		rename(Ch2_st + list[i]);
		Ch2_=getTitle();
		
		selectWindow("C3-"+list[i] + "_dup_max");
		rename(Ch3_st + list[i]);
		Ch3_=getTitle();
		
		run("Merge Channels...", "c1=[" + Ch1_ + "] c2=[" + Ch2_ + "] c3=[" + Ch3_ + "] keep");
		run("8-bit");
		run("Duplicate...", "title=[" + list[i] + "_mask] duplicate");
		Mask_1=getTitle();
		Mask_="Mask";
		run("Gaussian Blur...", "sigma=5");
		setBatchMode("show");
		setAutoThreshold("Otsu dark");
		waitForUser("Adjust Threshold. Try to cover all the gastruloid But NO areas outside.\n Do not worry if there are areas iside the gastruloid that the mask is not covering, the macro will fill the holes later.");
		run("Convert to Mask", "method=Otsu background=Dark calculate black create");
		run("Kill Borders");
		setBatchMode("show");
		Mask_2=getTitle();
		run("Erode");
		run("Dilate");
		run("Fill Holes");
		run("Divide...", "value=255.000");
		setBatchMode("hide");
		
		Table.create(list[i]);
		
		run("JACoP ", "imga=[" + Ch1_ + "] imgb=[" + Ch2_ + "] pearson cytofluo");
		Plot.getValues(x1, y2);
		Table.setColumn(Ch1_st, x1);
		Table.setColumn(Ch2_st, y2);
		run("JACoP ", "imga=[" + Ch1_ + "] imgb=[" + Ch3_ + "] pearson cytofluo");
		Plot.getValues(x1, y3);
		Table.setColumn(Ch3_st, y3);
		
		run("JACoP ", "imga=[" + Ch1_ + "] imgb=[" + Mask_2 + "] pearson cytofluo");
		Plot.getValues(x1, y5);
		Table.setColumn(Mask_, y5);
		
		if (Ch4) {
			selectWindow("C4-"+list[i] + "_dup_max");
			rename(Ch4_st + list[i]);
			Ch4_=getTitle();
			run("JACoP ", "imga=[" + Ch1_ + "] imgb=[" + Ch4_ + "] pearson cytofluo");
			Plot.getValues(x1, y4);
			Table.setColumn(Ch4_st, y4);
		}
		Table.save(dir2 + "/" + Results_name + "_" + i + ".csv");
		run("Close");
		selectWindow(Mask_2);
		saveAs("Tif", dir2 + "/" + Results_name + "_" + i + "_mask");
		run("Close All");
		
}}






