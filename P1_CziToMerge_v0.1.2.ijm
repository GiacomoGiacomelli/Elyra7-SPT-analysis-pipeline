//need to have the filled maximum projection folder
//Select the folder containing the multiple strains (it will do everything at once)


dir = getDirectory("Choose a directory");  // Open file dialog to select a directory
fileList = getFileList(dir);  // Get list of files and folders in the selected directory

for (i = 0; i < fileList.length; i++) {
	print(fileList[i]);
	path_in = replace(dir, "\\","/") + fileList[i]+"phase_czi/";
	path_in_max = replace(dir, "\\","/") + fileList[i]+"max/";  
	path_out = replace(dir, "\\","/") + fileList[i]+"phase/"; // Define the path where the folder should be created
	path_out_bin = replace(dir, "\\","/") + fileList[i]+"binary/"; // Define the path where the folder should be created
	path_out_merge = replace(dir, "\\","/") + fileList[i]+"merge/";
	path_in_bin = replace(dir, "\\","/") + fileList[i]+"binary/";
	path_in_phase = replace(dir, "\\","/") + fileList[i]+"phase/";
	
	//////////// Check if the folder already exists ///////////
	if (File.exists(replace(dir, "\\","/") + fileList[i]+"phase")) {
    	print("Folder already exists: " + replace(dir, "\\","/") + fileList[i]+"phase");
	} else {
    	// Create the folder
    	File.makeDirectory(replace(dir, "\\","/") + fileList[i]+"phase");  	
	}
	
	if (File.exists(replace(dir, "\\","/") + fileList[i]+"binary")) {
    	print("Folder already exists: " + replace(dir, "\\","/") + fileList[i]+"binary");
	} else {
    	// Create the folder
    	File.makeDirectory(replace(dir, "\\","/") + fileList[i]+"binary");
	}
	
	if (File.exists(replace(dir, "\\","/") + fileList[i]+"merge")) {
    	print("Folder already exists: " + replace(dir, "\\","/") + fileList[i]+"merge");
	} else {
    	// Create the folder
    	File.makeDirectory(replace(dir, "\\","/") + fileList[i]+"merge");  	
	}
	//////////////////////////////////////////////////////////
	
	list = getFileList(path_in);
	list1 = getFileList(path_in_max);
	for (d=0; d<list.length; d++) {
		print(list[d]);
		print(list1[d]);
	run("Bio-Formats (Windowless)", "open=["+path_in_max+list1[d]+"] use_virtual_stack");
 	run("Bio-Formats (Windowless)", "open=["+path_in+list[d]+ "] use_virtual_stack");
 	title = "series00";
 	f=d+1;
 	if (f>9)
    	  title = "series0";
 	saveAs("tiff",path_out+title+f);
 	
 	imgnum = nImages(); // Get the number of open images
	if (imgnum == 0) {
	    print("No images are currently open.");
	} else {
	// Preprocessing	
	selectImage(1);
	rename("max");
	run("Smooth");
	run("Enhance Contrast...", "saturated=0.3");
	run("Subtract Background...", "rolling=50 dark");
	// Scaling
	run("Scale...", "x=10 y=10 width=12800 height=12800 interpolation=Bicubic average create title=Results");
 	run("8-bit");
 	rename("max_8");
 	close("max");
 	 	
 	// Preprocessing	
	selectImage(1);
	rename("bright");
	run("Smooth");
	run("Enhance Contrast...", "saturated=0.3");
	run("Subtract Background...", "rolling=10 light");
	// Scaling
	run("Scale...", "x=10 y=10 width=12800 height=12800 interpolation=Bicubic average create title=Results");
 	run("8-bit");
 	rename("bright_8");
 	close("bright");
 	
 	run("Image Calculator...", "image1=[bright_8] operation=Subtract image2=[max_8] create");
 	close("bright_8");
 	close("max_8");
	setAutoThreshold("Default");
	run("Convert to Mask");
	run("Despeckle");
	run("Fill Holes");
	for (l = 0; l < 5; l++) {
	    run("Dilate");
	}
	saveAs("Tiff", path_out_bin + title + f);
	close();
	}
	}
 	
list_tiff = getFileList(path_in_bin);
	list = newArray();
	for (w = 0; w < list_tiff.length; w++) {
    	if (endsWith(list_tiff[w], ".tif") || endsWith(list_tiff[w], ".tiff")) {
        	list = Array.concat(list, list_tiff[w]); // Add the TIFF file to the array
    	}
	}
	
	list1 = getFileList(path_in_max);
	list2 = getFileList(path_in_phase);
	for (d=0; d<list.length; d++) {
		print(list[d]);
		print(list1[d]);
		print(list2[d]);
	run("Bio-Formats (Windowless)", "open=["+path_in_max+list1[d]+"] use_virtual_stack");
	rename("ori");
	run("Scale...", "x=10 y=10 width=12800 height=12800 interpolation=Bicubic average create title=Results");
 	run("8-bit");
	rename("max");
 	close("ori");
 	run("Bio-Formats (Windowless)", "open=["+path_in_bin+list[d]+ "] use_virtual_stack");
 	rename("bin");
 	run("Bio-Formats (Windowless)", "open=["+path_in_phase+list2[d]+"] use_virtual_stack");
	rename("p");
	run("Scale...", "x=10 y=10 width=12800 height=12800 interpolation=Bicubic average create title=Results");
 	run("8-bit");
	rename("phase");
 	close("p");
 	
 	title = "series00";
 	f=d+1;
 	if (f>9)
    	  title = "series0";
    	  
 	// Merge multiple channels into a single composite image
	// Ensure the images for each channel are open

	// Assign each channel to variables by name
	gray = "phase"; // Replace with the name of your first channel image
	red = "bin"; // Replace with the name of your second channel image
	green = "max"; // Replace with the name of your second channel image
	// Check if the images are open
	if (isOpen(gray) && isOpen(red)) {
    	// Merge the channels into a single composite image
    	run("Merge Channels...", "c1=[" + red + "] c2=[" + green + "] c4=[" + gray + "] create keep");
    	saveAs("Tiff", path_out_merge + title + f);
    	close("max");
    	close("bin");
    	close("phase");
    	close();
	} else {
    	print("One or more channel images are not open. Please check the names.");
	}
	}
}

print("end")
