// Focal adhesion quantification ImageJ script from Cook, Ueno, Lodoen, 2017

// requires CLAHE plugin
// should be used to analyze individual cells at a time

run("Subtract Background...", "rolling=10 sliding");

run("CLAHE ", "blocksize=19 histogram=256 maximum=6");

run("Exp");

run("Enhance Contrast", "saturated=0.35");
run("Apply LUT");


setAutoThreshold("Default dark");
setOption("BlackBackground", false);
run("Convert to Mask");

run("Analyze Particles...", "size=0.15-Infinity circularity=0.00-0.99 show=Masks display exclude");

