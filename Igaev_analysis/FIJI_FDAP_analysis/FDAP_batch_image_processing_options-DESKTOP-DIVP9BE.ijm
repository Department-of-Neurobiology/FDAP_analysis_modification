/*
 * Macro template to process multiple images in a folder
 */


// THERE SHOULDN'T BE ANY FILES WITH THE CHOSEN EXTENSION IN SUBFOLDERS

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".nd2") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
// roiManager("Delete");
run("Set Measurements...", "integrated redirect=None decimal=3");
roiManager("reset");
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	print("Processing: " + input + File.separator + file);
	open(input + File.separator + file);
	selectWindow(file +" - C=0"); 
	Stack.setFrame(2);
	resetMinAndMax();
	name = File.nameWithoutExtension();
	id=getImageID();
	waitForUser("Draw ROI and check, close file if ROI is inappropriate");
	if (isOpen(id) == false) {
	 	f = File.open(output + File.separator + name + "_inappr_ROI" + ".txt");
	 	File.close(f);
        // Run("Close);
	}
	else {
	//all the analysis steps
	roiManager("Add");
	roiManager("Multi Measure");
	saveAs("Results", output + File.separator + name + "_RawIntDen1.csv");
	roiManager("Save", output + File.separator + name + ".zip");
	roiManager("Delete");
	close();
	// run("Close");
	}
	// doCommand("Start Animation [\\]");
	// waitForUser("Check area");
	
}
// saveAs("Results", output + File.separator + "test_batch.csv"); //is possible to save batch 