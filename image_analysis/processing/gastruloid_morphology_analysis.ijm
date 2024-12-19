waitForUser("Activate the Shape descriptors. Analyze > Set Measurements");
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");

list = getFileList(dir1);

for (i=0; i<list.length; i++) {
	if (endsWith(list[i], ".czi")){
		path = dir1 + list[i];
		run("Bio-Formats Importer", "open=[" + path  + "]"); 
		run("Split Channels");
		close("C2-"+list[i]);
		close("C3-"+list[i]);
		selectWindow("C1-"+list[i]);
		original=getTitle();
		run("Grays");
		run("8-bit");
		run("Duplicate...", "title=duplicate"+i);
		selectWindow("duplicate"+ i);
		run("Gaussian Blur...", "sigma=3");
		selectWindow("duplicate"+ i);
		run("Normalize Local Contrast", "block_radius_x=200 block_radius_y=200 standard_deviations=3 center");
		run("Threshold...", "Otsu dark");
		waitForUser("Adjust the proper threshold");
		run("Convert to Mask");
		run("Invert");
		run("Fill Holes");
		run("Analyze Particles...", "size=10000-Infinity show=Outlines display exclude include add");
		selectWindow(original);
		run("Duplicate...", "title=duplicate2");
		run("Scale Bar...", "width=300 height=8 font=28 color=White background=None location=[Lower Right] bold overlay");
		saveAs("TIFF", dir2+list[i]);
		selectWindow(original);
		roiManager("Select", i);
		run("Flatten");
		run("Scale Bar...", "width=300 height=8 font=28 color=White background=None location=[Lower Right] bold overlay");
		saveAs("Jpeg", dir2+list[i]+"_segmented");
	}
}
selectWindow("Results");
saveAs("Results",dir2+"Results.csv");
run("Close");
roiManager("deselect");
roiManager("Save",dir2+"RoiSet.zip");
roiManager("reset");
run("Close All");

