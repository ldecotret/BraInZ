/* This macro will try to determine the distance cells have migrated from the central cell mass
 *  It will do the following things:
 *  	1) Find the central mass of cells
 *  	2) Fit an ellipse of the central mass
 *  	3) Find the properties of the ellipse
 *  	4) Determine the distance of particle away from the closest point on the ellipse
 */

input = getTitle();

run("Duplicate...", "title=DUP");		
run("Gaussian Blur...", "sigma=2"); 	//process the image so it does not have speckles when thresholded
run("Options...", "iterations=1 count=1 black do=Nothing");
run("Auto Threshold", "method=Huang white"); 	//thresholding the image, "Huang" is a very generous thresholding algorithm that will probably capture the central mass of cells

run("Set Measurements...", "area center fit limit redirect=None decimal=4");
run("Analyze Particles...", "  show=Nothing display add");

//This part allows to find the particle with the largest area - which most likely is the central mass
selectWindow("Results");
area = Table.getColumn("Area");
indices_max = Array.findMaxima(area, 1);
imax = indices_max[0];

roiManager("Show None");
roiManager("select", imax);
run("Fit Ellipse");
roiManager("reset");

selectWindow("DUP");
waitForUser("Adjust the ellipse...");
roiManager("add");
run("Clear Results");

roiManager("select", 0);
roiManager("Measure");

//This part stores the value of the fitted ellipse properties
centre_x = getResult("XM", 0);
centre_y = getResult("YM", 0);
e_angle = getResult("Angle", 0);
major_axis = getResult("Major", 0);
major_r = major_axis/2;
minor_axis = getResult("Minor", 0);
minor_r = minor_axis/2;
error = 1e-5;

//Fit the ellipse

run("Create Mask");

//Create an image with the central ellipse subtracted from the image and wathershed the remaining "cells"
imageCalculator("Subtract create", "DUP","Mask");
selectWindow("Result of DUP");
run("Watershed");


selectWindow("Results");
p_count = nResults;
Table.deleteRows(0, nResults);	//Clear results

roiManager("reset");	//Reset the ROI manager

//Go back to the image of the remaining particles and find their position
selectWindow("Result of DUP");
run("Set Measurements...", "area center redirect=None decimal=4");
run("Analyze Particles...", "size=10-Infinity show=Masks display add");

//For each particle identified, calculate the distance from the particle
selectWindow("Results");
for (row = 0; row<nResults; row++) {
	particle_x = getResult("XM", row);
	particle_y = getResult("YM", row);
	distance_centre = sqrt((Math.sqr(particle_x - centre_x)) + (Math.sqr(particle_y - centre_y)));
	distance_ellipse = estimate_distance(particle_x, particle_y, major_r, minor_r, centre_x, centre_x, e_angle, error);
	setResult("Distance from ellipse centre", row, distance_centre);
	setResult("Distance from closest point on ellipse", row, distance_ellipse);
}
selectWindow("Results");
run("Summarize");


/* The following scripts attempts to label the analysed particles on the duplicated image - this is NOT COMPLETE - DO NOT INCLUDE THIS SECTION INTO THE SCRIPT -----
 *  
selectWindow(input);
		run("Select All");
		run("Copy");
		
		selectWindow("DUP");
		run("RGB Color");
		run("Add Slice");
		run("Paste");

				selectWindow(input);
		setSlice(1);
		run("Duplicate...", "title=Positions");
		run("8-bit");
		run("Invert");
		run("Duplicate...", "title=Numbers");
		selectWindow("Positions");
		setThreshold(250, 255);
		run("Convert to Mask");
		run("Red");
		
		selectWindow("Numbers");
		setThreshold(150, 250);
		run("Convert to Mask");
		run("Green");
				
		selectWindow(input);
		setSlice(2);
		run("Add Slice");
		setSlice(2);		
		run("Duplicate...", "title=Raw");

		imageCalculator("Add", "Raw","Positions");
		imageCalculator("Add", "Raw","Numbers");


		selectWindow("Raw");
		run("Select All");
		run("Copy");
	
		selectWindow(input);
		setSlice(3);
		run("Paste");
		run("Reverse");
----------------------------------------------------------------------------------- END OF THE SECTION		
*/











//The following functions allows the calculation for the distance away from the closest point on the ellipse
//These functions are derived from this paper: https://www.ma.ic.ac.uk/~rn/distance2ellipse.pdf
  
 
function ellipe_tan_dot(rx, ry, px, py, theta) {
	return ( ((Math.sqr(rx) - Math.sqr(ry)) * cos(theta) * sin(theta)) - (px * rx * sin(theta)) + (py * ry * cos(theta)));
}
//  Dot product of the equation of the line formed by the point with another point on the ellipse's boundary and the tangent of the ellipse
//  at that point on the boundary.


function ellipe_tan_dot_derivative(rx, ry, px, py, theta) {    //The derivative of ellipe_tan_dot.
    return ( ((Math.sqr(rx) - Math.sqr(ry)) * (Math.sqr(cos(theta)) - Math.sqr(sin(theta)))) - (px * rx * cos(theta)) - (py * ry * sin(theta)));
}

function estimate_distance(x, y, rx, ry, x0, y0, angle, error) {
	/*
	 * Given a point (x, y), and an ellipse with major - minor semi-axis (rx, ry),
    its center at (x0, y0), and with a counter clockwise rotation of
    `angle` degrees, will return the distance between the ellipse and the
    closest point on the ellipses boundary.
	 */
	x -= x0;
    y -= y0;
    
    if (angle != 0) {
        // rotate the points onto an ellipse whose rx, and ry lay on the x, y axis
        
        angle = (PI / 180.0 * angle);
        new_x = (x * cos(angle)) - (y * sin(angle));
        new_y = (x * sin(angle)) + (y * cos(angle));
        x = new_x;
        y = new_y;

    }

    theta = atan2(rx * y, ry * x);
    while (Math.abs(ellipe_tan_dot(rx, ry, x, y, theta)) > error) {
    	if (ellipe_tan_dot_derivative(rx, ry, x, y, theta) == 0) {
    		theta -= (ellipe_tan_dot(rx, ry, x, y, theta)) / 1e-5;
    	}else{
    	
        theta -= (ellipe_tan_dot(rx, ry, x, y, theta)) / (ellipe_tan_dot_derivative(rx, ry, x, y, theta));
        
    }    
        px = rx * cos(theta);
    	py = ry * sin(theta);

    	

    	
    return sqrt((Math.sqr(x - px) + Math.sqr(y - py)));
}
