dir = getDirectory("Choose a directory");  // Open file dialog to select a directory
fileList = getFileList(dir);  // Get list of files and folders in the selected directory
membrane = 1;
psize=0.0096782;

for (i = 0; i < fileList.length; i++) {
	print(fileList[i]);
	path_in_bin = replace(dir, "\\","/") + fileList[i]+"merge/";
	path_in_bin2 = replace(dir, "\\","/") + fileList[i]+"binary/";
	path_in_roi = replace(dir, "\\","/") + fileList[i]+"rois/";
		
	list1 = getFileList(path_in_bin);
	list2 = getFileList(path_in_roi);
	
	for (d=0; d<list1.length; d++) {
		//print(list1[d]);
		//print(list2[d]);
		print(list1[d]);
		name = replace(list1[d], ".tif", "");
		roiFile = name + ".zip";

		//if (File.exists(path_in_roi + roiFile)) {
   		// roiManager("Open", path_in_roi + roiFile);
		//} else {
    	//	print("No matching ROI file for " + list1[d]);
		//}
		print(roiFile);
		
		
		path_out_cells = replace(dir, "\\","/") + fileList[i]+"binary/cells"+(d+1)+"/";
		path_out_cells1 = replace(dir, "\\","/") + fileList[i]+"binary/smallcells"+(d+1)+"/";
	
	// Check if the folder already exists
	if (File.exists(replace(dir, "\\","/") + fileList[i]+"binary/cells"+(d+1)+"/")) {
    	print("Folder already exists: " + replace(dir, "\\","/") + fileList[i]+"binary/cells"+(d+1)+"/");
	} else {
    	// Create the folder
    	File.makeDirectory(replace(dir, "\\","/") + fileList[i]+"binary/cells"+(d+1)+"/");  	
	}
	
	// Check if the folder already exists
	if (File.exists(replace(dir, "\\","/") + fileList[i]+"binary/smallcells"+(d+1)+"/")) {
    	print("Folder already exists: " + replace(dir, "\\","/") + fileList[i]+"binary/smallcells"+(d+1)+"/");
	} else {
    	// Create the folder
    	File.makeDirectory(replace(dir, "\\","/") + fileList[i]+"binary/smallcells"+(d+1)+"/");  	
	}
	
		run("Bio-Formats (Windowless)", "open=["+path_in_bin2+list1[d]+"] use_virtual_stack");
		run("Select All");
		run("Fill", "slice");
		run("ROI Manager...");
		roiManager("Open", path_in_roi + roiFile);
		n = roiManager("count");
 			for (m=0; m<2*n; m++) {
 				if (m<n) {
      		roiManager("select", m);
      		run("Clear", "slice");
      		roiManager("rename", m);
      		run("Properties... ", "name=0157-0027 position=none stroke=none width=0 fill=none list");
      		wait(25);
      		if (m<10) {
	  		saveAs("Results", path_out_cells + "Cell00" + m + ".txt");
	  		if (isOpen("Cell00"+m+".txt")) { 
       			selectWindow("Cell00"+m+".txt"); 
       		run("Close"); 
	  		}
	  		if (membrane==1){
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
				run("Set Measurements...", "centroid decimal=3");
				run("Measure");
				x = getResult("X", nResults-1);
				y = getResult("Y", nResults-1);
				selectWindow("Results"); run("Close");
				// round coordinates for doWand
				xInt = round(x/psize);
				yInt = round(y/psize);
				roiManager("Deselect");
				doWand(xInt, yInt);
				roiManager("Add");
      		}
      		}
	  		if (m > 9 && m < 100) {
	  		saveAs("Results", path_out_cells + "Cell0" + m + ".txt");
	  		if (isOpen("Cell0"+m+".txt")) { 
       			selectWindow("Cell0"+m+".txt"); 
       		run("Close"); 
	  		}
	  		if (membrane==1){
	  		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
				run("Set Measurements...", "centroid decimal=3");
				run("Measure");
				x = getResult("X", nResults-1);
				y = getResult("Y", nResults-1);
				selectWindow("Results"); run("Close");
				// round coordinates for doWand
				xInt = round(x/psize);
				yInt = round(y/psize);
				roiManager("Deselect");
				doWand(xInt, yInt);
				roiManager("Add");
      		}
 			}
	  		if (m > 99) {
	  		saveAs("Results", path_out_cells + "Cell" + m + ".txt");
	  		if (isOpen("Cell"+m+".txt")) { 
       			selectWindow("Cell"+m+".txt"); 
       		run("Close"); 
	  		}
	  		if (membrane==1){
	  		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
      		run("Erode");
				run("Set Measurements...", "centroid decimal=3");
				run("Measure");
				x = getResult("X", nResults-1);
				y = getResult("Y", nResults-1);
				selectWindow("Results"); run("Close");
				// round coordinates for doWand
				xInt = round(x/psize);
				yInt = round(y/psize);
				roiManager("Deselect");
				doWand(xInt, yInt);
				roiManager("Add");
	  		}
 			}
			}
			
			if (m>(n-1)) {
				roiManager("select", m);
				mn=m-n;
				roiManager("rename", mn);
      			run("Properties... ", "name=0157-0027 position=none stroke=none width=0 fill=none list");
      			wait(25);
      			
	      			if (mn<10) {
				  		saveAs("Results", path_out_cells1 + "Cell00" + mn + ".txt");
				  		if (isOpen("Cell00"+mn+".txt")) { 
			       			selectWindow("Cell00"+mn+".txt"); 
			       		run("Close"); 
				  		}
			      			}
			      	if (mn>9 && mn < 100) {
				  		saveAs("Results", path_out_cells1 + "Cell0" + mn + ".txt");
				  		if (isOpen("Cell0"+mn+".txt")) { 
			       			selectWindow("Cell0"+mn+".txt"); 
			       		run("Close"); 
				  		}
			      			}
			      	if (mn>99) {
				  		saveAs("Results", path_out_cells1 + "Cell" + mn + ".txt");
				  		if (isOpen("Cell"+mn+".txt")) { 
			       			selectWindow("Cell"+mn+".txt"); 
			       		run("Close"); 
				  		}
			}
			}
			}
 			roiManager("reset"); // Clear all ROIs
    		run("Close"); // Close the image 	
		
	}
}
print("end");
