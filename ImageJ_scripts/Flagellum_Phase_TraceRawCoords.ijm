fileSuffix=".tif";
baseInertia=0.1;

path=getDirectory("");
files=getFileList(path);

var prefilterRadius=1;
var unsharpRadius1=20;
var unsharpRadius2=15;
var unsharpWeight=0.6;
var backgroundSubRadius=25;
var autoThreshold="Triangle"; //"Yen", "Triangle";
var refineIterations=1;
var refineRadius=2;

print("==Thresholding==");
if (File.exists(path+"threshold.txt")==true) {
	print("Using custom threshold parameters");
	data=File.openAsString(path+"threshold.txt");
	data=split(data, "\r\n");
	prefilterRadius=parseInt(data[0]);
	unsharpRadius1=parseInt(data[1]);
	unsharpRadius2=parseInt(data[2]);
	unsharpWeight=parseFloat(data[3]);
	backgroundSubRadius=parseInt(data[4]);
	autoThreshold=data[5];
	refineIterations=parseInt(data[6]);
	refineRadius=parseInt(data[7]);
} else {
	print("Using default threshold parameters");
}
print("Pre-filter radius", prefilterRadius, "px");
print("Unsharp filter radius 1", unsharpRadius1, "px");
print("Unsharp filter radius 2", unsharpRadius2, "px");
print("Unsharp filter weight", unsharpWeight, "(0.1 to 0.9)");
print("Background subtraction radius", backgroundSubRadius, "px");
print("Auto threshold method", autoThreshold, "");
print("Refinement filter iterations", refineIterations, "");
print("Refinement filter radius", refineRadius, "px");

for (i=0; i<lengthOf(files); i++) {
	if (endsWith(files[i], ".tif")==true) {
		print("==Analysing TIFF file==");
		print(files[i]);
		t=0;
		//Analyse every tiff, and loop through any extra metadata files
		while (File.exists(path+files[i]+"."+t+".txt")==true || t==0) {
			//Determine output directory and metadata path
			//Output data will be placed in an output directory
			if (t==0) {
				outdir=path+files[i]+"_analysis"+File.separator();
				metadata=files[i]+".txt";
				if (File.exists(path+metadata)==true) {
					print("Analysing primary metadata");
					print(metadata);
				} else {
					print("Analysing TIFF");
					print(files[i]);
				}
			} else {
				outdir=path+files[i]+"."+t+"_analysis"+File.separator();
				metadata=files[i]+"."+t+".txt";
				print("Analysing metadata index "+t);
				print(metadata);
			}
			//Only analyse the .tif if the output file does not exists or is empty
			//Delete the output file or directory to force re-analysis
			direxists=File.exists(outdir);
			outexists=false;
			if (direxists==true) {
				outexists=File.exists(outdir+File.separator()+"rawcoordinates.txt");
				if (outexists==true) {
					outstr=File.openAsString(outdir+File.separator()+"rawcoordinates.txt");
					if (lengthOf(outstr)==0) {
						outexists=false;
					}
				}
			}
			if (direxists==false || outexists==false) {
				//Make the output directory
				if (direxists==false) {
					File.makeDirectory(outdir);
				}
				//Write an empty output file
				outfile=File.open(outdir+File.separator()+"rawcoordinates.txt");
				File.close(outfile);
				//Open the image
				open(path+files[i]);
				ori=getImageID();
				if (is("hyperstack")==true) {
					//Trace from first channel for hyperstacks
					run("Duplicate...", "duplicate channels=1");
					src=getImageID();
					selectImage(ori);
					close();
				} else {
					src=getImageID();
				}
				//Determine analysis parameters
				//Check for saved point selection
				if (selectionType()==10) {
					getSelectionCoordinates(sx, sy);
					print("Base reference point", sx[0], sy[0], "px");
					flabasex=sx[0];
					flabasey=sy[0];
				} else {
					flabasex=getWidth()/2;
					flabasey=getWidth()/2;
				}
				vidcrop=newArray(4);
				vidcrop[0]=0;
				vidcrop[1]=0;
				vidcrop[2]=getWidth();
				vidcrop[3]=getHeight();
				vidend=nSlices();
				//Look for metadata text file
				if (File.exists(path+metadata)==true) {
					data=File.openAsString(path+metadata);
					data=split(data, "\n");
					if (lengthOf(data)>=1) {
						coords=split(data[0], " ");
						flabasex=parseInt(coords[0]);
						flabasey=parseInt(coords[1]);
					}
					if (lengthOf(data)>=2) {
						vidcrop=split(data[1], " ");
						for (j=0; j<lengthOf(vidcrop); j++) {
							vidcrop[j]=parseInt(vidcrop[j]);
						}
					}
					if (lengthOf(data)>=3) {	
						vidend=parseInt(data[2]);
					}
				}
				makeRectangle(vidcrop[0], vidcrop[1], vidcrop[2], vidcrop[3]);
				run("Crop");
	
				if (nSlices()>1) {
					doFlagellumTracing(src, outdir, flabasex, flabasey, vidcrop, vidend);
				}
				selectImage(src);
				close();
			} else {
				print("File skipped, output directory already exists");
			}
			t++;
		}
	}
}

function doFlagellumTracing(src, path, flabasex, flabasey, vidcrop, vidend) {
	coordinateLog="";
	flabasex=flabasex-vidcrop[0];
	flabasey=flabasey-vidcrop[1];
	numTracedPoints=newArray(minOf(nSlices(), vidend));

	selectImage(src);
	for (i=0; i<minOf(nSlices(), vidend); i++) {
		setSlice(i+1);
		setBatchMode(true);
		run("Select None");
		run("Duplicate...", " ");
		thr=getImageID();
		imageThresholdingPhase(thr, prefilterRadius, unsharpRadius1, unsharpRadius2, unsharpWeight, backgroundSubRadius, autoThreshold, refineIterations, refineRadius);
		imageThresholdingRemoveBright(src, thr);
		binaryLargestParticle(thr);
		selectImage(thr);
		run("Duplicate...", " ");
		ske=getImageID();
		success=binaryGenerateTraceSkeleton(ske, src, flabasex, flabasey);
		if (success==true) {
			selectImage(src);
			getSelectionCoordinates(sx, sy);
			selectImage(thr);
			run("Distance Map");
			widths=newArray(lengthOf(sx));
			for (j=0; j<lengthOf(sx); j++) {
				widths[j]=getPixel(sx[j], sy[j]);
			}
			flabasex=(sx[0]-flabasex)*baseInertia+flabasex;
			flabasey=(sy[0]-flabasey)*baseInertia+flabasey;
		} else {
			sx=newArray(1);
			sx[0]=flabasex;
			sy=newArray(1);
			sy[0]=flabasey;
			widths=newArray(1);
			widths[0]=0;
		}
		selectImage(thr);
		close();
		for (j=0; j<lengthOf(sx); j++) {
			coordinateLog+=""+(sx[j]+vidcrop[0])+","+(sy[j]+vidcrop[1])+","+widths[j]+" ";
		}
		numTracedPoints[i]=lengthOf(sx);
		coordinateLog+="\n";
		setBatchMode(false);
	}

	Plot.create("Trace profile", "Frame", "Points traced", numTracedPoints);
	Plot.show();
	saveAs("Jpeg",path+File.separator()+"traceProfile.jpg");
	close();

	outfile=File.open(path+File.separator()+"rawCoordinates.txt");
	print(outfile, coordinateLog);
	File.close(outfile);
}

function imageThresholdingPhase(imageID, prefilterRadius, unsharpRadius1, unsharpRadius2, unsharpWeight, backgroundSubRadius, autoThreshold, refineIterations, refineRadius) {
	selectImage(imageID);
	run("32-bit");
	run("Gaussian Blur...", "sigma="+prefilterRadius);
	run("Unsharp Mask...", "radius="+unsharpRadius1+" mask="+unsharpWeight);
	run("Unsharp Mask...", "radius="+unsharpRadius2+" mask="+unsharpWeight);
	run("Subtract Background...", "rolling="+backgroundSubRadius+" light");
	run("Unsharp Mask...", "radius="+unsharpRadius2+" mask="+unsharpWeight);
	setAutoThreshold(autoThreshold);
	getThreshold(autoMin, autoMax);
	run("Convert to Mask");
	for (i=0; i<refineIterations; i++) {
		run("Gaussian Blur...", "sigma="+refineRadius);
		changeValues(0, 128, 0);
		changeValues(1, 255, 255);
	}
	run("8-bit");
	run("Fill Holes");
}

function imageThresholdingRemoveBright(src, thr) {
	selectImage(src);
	getRawStatistics(area, mean);
	limit=mean*0.9;
	selectImage(thr);
	run("Find Maxima...", "noise=1 output=[Point Selection]");
	getSelectionCoordinates(sx, sy);
	count=0;
	remove=newArray(lengthOf(sx));
	for (i=0; i<lengthOf(sx); i++) {
		selectImage(thr);
		doWand(sx[i], sy[i]);
		getSelectionCoordinates(cx, cy);
		selectImage(src);
		makeSelection("polygon", cx, cy);
		getRawStatistics(area, mean);
		if (mean>limit) {
			remove[i]=true;
		} else {
			count++;
			remove[i]=false;
		}
	}
	selectImage(thr);
	if (count>0) {
		for (i=0; i<lengthOf(sx); i++) {
			if (remove[i]==true) {
				doWand(sx[i], sy[i]);
				setColor(0);
				fill();
			}
		}
	}
	selectImage(src);
	run("Select None");
	selectImage(thr);
	run("Select None");
}

function binaryLargestParticle(imageID) {
	selectImage(imageID);
	run("Find Maxima...", "noise=10 output=[Point Selection]");
	getSelectionCoordinates(x, y);
	maxArea=0;
	maxAreaIndex=0;
	for (i=0; i<lengthOf(x); i++) {
		doWand(x[i], y[i]);
		getRawStatistics(area);
		if (area>maxArea) {
			maxArea=area;
			maxAreaIndex=i;
		}
	}
	for (i=0; i<lengthOf(x); i++) {
		if (i!=maxAreaIndex) {
			doWand(x[i], y[i]);
			setColor(0);
			fill();
		}
	}
	run("Select None");
}

function binaryGenerateTraceSkeleton(thrImageID, targetImageID, baseX, baseY) {
	ter=makeSkeletonTerminus(thrImageID);
	run("Duplicate...", " ");
	ter2=getImageID();
	changeValues(2, 255, 0);
	run("Find Maxima...", "noise=0.5 output=[Point Selection]");
	if (selectionType()!=-1) {
		getSelectionCoordinates(terx, tery);
	} else {
		selectImage(thrImageID);
		close();
		selectImage(ter2);
		close();
		selectImage(ter);
		close();
		return false;
	}
	selectImage(ter2);
	close();
	selectImage(ter);

	mind=getWidth()*getHeight();
	minind=-1;
	for (k=0; k<lengthOf(terx); k++) {
		curd=pow(pow(terx[k]-baseX, 2)+pow(tery[k]-baseY, 2), 0.5);
		if (curd<mind) {
			mind=curd;
			minind=k;
		}
	}

	traceSkeleton(terx[minind], tery[minind]);
	if (selectionType()!=-1) {
		getSelectionCoordinates(matx, maty);
	} else {
		selectImage(thrImageID);
		close();
		selectImage(ter);
		close();
		return false;
	}

	selectImage(ter);
	close();

	selectImage(thrImageID);
	close();

	selectImage(targetImageID);
	if (lengthOf(matx)>0) {
		makeSelection("polyline", matx, maty);
	} else {
		return false;
	}
	return true;
}

function makeMedialAxisTransform(thr) {
	id=makeSkeletonFilter(thr, "medial");
	return id;
}

function makeSkeletonTerminus(thr) {
	id=makeSkeletonFilter(thr, "terminus");
	return id;
}

function makeSkeletonFilter(thr, type) {
	selectImage(thr);
	run("Select None");

	run("Duplicate...", "slice");
	mat=getImageID();
	changeValues(0, 254, 0);
	run("Make Binary");

	run("Duplicate...", "slice");
	tmp=getImageID();
	run("Skeletonize");
trimSkeletonBranches();
	changeValues(1, 255, 1);
	run("Select All");
	setPasteMode("Copy");
	run("Copy");
	selectImage(tmp);
	close();

	selectImage(mat);
	if (type=="medial") {
		run("Distance Map");
		setPasteMode("Multiply");
		run("Paste");
		run("Select None");
	} else if (type=="terminus") {
		run("Paste");
		run("Select None");
		changeValues(1, 255, 1);
		//run("Convolve...", "text1=[1 1 1\n1 0 1\n1 1 1\n] slice");
		for (x=0; x<getWidth(); x++) {
			for (y=0; y<getHeight(); y++) {
				if (getPixel(x, y)!=0) {
					//Connected neighbours
					dx=newArray(0, 1, 1, 1, 0, -1, -1, -1);
					dy=newArray(-1, -1, 0, 1, 1, 1, 0, -1);
					count=0;
					for (i=0; i<lengthOf(dx); i++) {
						if (getPixel(x+dx[i], y+dy[i])!=0) {
							count++;
						}
					}
					setPixel(x, y, count);
				}
			}
		}
	}

	return mat;
}

function traceSkeleton(x, y) {
	//Traces a skeleton from a start point (should be a terminus) until it hits a branch
	ox=newArray(1);
	oy=newArray(1);
	ox[0]=x;
	oy[0]=y;
	cont=true;
	while (cont==true) {
		setPixel(x, y, 0);
		dx=0;
		dy=0;
		for (a=-1; a<=1; a++) {
			for (b=-1; b<=1; b++) {
				v=getPixel(x+a, y+b);
				if (v>0 && v<=2) {
					dx=a;
					dy=b;
				}
			}
		}
		if (dx!=0 || dy!=0) {
			x+=dx;
			y+=dy;
			setPixel(x, y, 0);
			ox=Array.concat(ox, x);
			oy=Array.concat(oy, y);
		} else {
			cont=false;
		}
	}
	if (lengthOf(ox)>1) {
		makeSelection("polyline", ox, oy);
	}
}

function simplifySelection(detail) {
	if (selectionType()!=-1) {
		getSelectionCoordinates(sx, sy);
		osx=newArray(0);
		osy=newArray(0);
		osx=Array.concat(osx, sx[0]);
		osy=Array.concat(osy, sy[0]);
		dist=0;
		i=0;
		while (i<lengthOf(sx)-1) {
			ddist=pow(pow(sx[i]-sx[i+1], 2)+pow(sy[i]-sy[i+1], 2), 0.5);
			dist+=ddist;
			i++;
			if (dist>detail*lengthOf(osx)) {
				osx=Array.concat(osx, sx[i]);
				osy=Array.concat(osy, sy[i]);
			}
		}
		osx=Array.concat(osx, sx[lengthOf(sx)-1]);
		osy=Array.concat(osy, sy[lengthOf(sx)-1]);
		makeSelection("Polyline", osx, osy);
	}
}
		

function simplifySelectionOld(detail) {
	if (selectionType()!=-1) {
		getSelectionCoordinates(sx, sy);
		osx=newArray(0);
		osy=newArray(0);
		for (i=0; i<lengthOf(sx)-detail/2; i+=detail) {
			osx=Array.concat(osx, sx[i]);
			osy=Array.concat(osy, sy[i]);
		}
		osx=Array.concat(osx, sx[lengthOf(sx)-1]);
		osy=Array.concat(osy, sy[lengthOf(sx)-1]);
		makeSelection("Polyline", osx, osy);
	}
}

function trimSkeletonBranches() {
	branchLength=20;
	src=getImageID();
	run("Select None");
	run("Duplicate...", " ");
	tmp=getImageID();
	selectImage(src);
	run("Divide...", "value=255");
	run("Convolve...", "text1=[1 1 1\n1 0 1\n1 1 1\n]");
	selectImage(tmp);
	run("Copy");
	setPasteMode("Transparent-white");
	close();
	selectImage(src);
	run("Paste");
	run("Select None");
	setPasteMode("Copy");
	run("Duplicate...", " ");
	tmp=getImageID();
	changeValues(2, 255, 0);
	run("Find Maxima...", "noise=0.5 output=[Point Selection]");
	if (selectionType!=-1) {
		getSelectionCoordinates(sx, sy);
		selectImage(src);
		for (i=0; i<lengthOf(sx); i++) {
			l=0;
			cont=true;
			vx=newArray(branchLength+1);
			vy=newArray(branchLength+1);
			vx[0]=sx[i];
			vy[0]=sy[i];
			while (cont==true) {
				na=0;
				nb=0;
				for (a=-1; a<2; a++) {
					for (b=-1; b<2; b++) {
						if ((a!=0 || b!=0) && getPixel(vx[l]+a, vy[l]+b)>0 && getPixel(vx[l]+a, vy[l]+b)<=2) {
							if (l==0) {
								na=a;
								nb=b;
							} else if (vx[l-1]!=vx[l]+a || vy[l-1]!=vy[l]+b) {
								na=a;
								nb=b;
							}
						}
					}
				}
				if (na==0 && nb==0) {
					cont=false;
				} else if (l<branchLength) {
					l++;
					vx[l]=vx[l-1]+na;
					vy[l]=vy[l-1]+nb;
				} else {
					cont=false;
				}
			}
			if (l<branchLength) {
				for (j=0; j<=l; j++) {
					setPixel(vx[j], vy[j], 0);
				}
			}
		}
	}
	selectImage(tmp);
	close();
	selectImage(src);
	changeValues(1, 255, 255);
	run("Skeletonize");
}

