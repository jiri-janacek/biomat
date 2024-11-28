package biomat;

import ij.IJ;
import ij.plugin.*;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;

/**
 * 3D Euclidean distance map
 * BINARY images.
 *
 * @author Jiri Janacek
 */

public class Distance_Map_3D implements PlugIn {

	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
        ImageProcessor ip = imp.getProcessor();
//        if (! (ip.isBinary() && (imp.getNSlices() > 1) && (imp.getNSlices() == imp.getStackSize()))) {
//            IJ.error("single z-stack of 8-bit binary images (0 and 255) required.");
        if (! ip.isBinary()) {
            IJ.error("8-bit binary images (0 and 255) required.");
			return;
		}

		imp.startTiming();  
		run(imp);		
		imp.updateAndDraw();
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private void run(ImagePlus imp) {
		int dims[] = imp.getDimensions();
		Calibration cal = imp.getCalibration();
		double dx = cal.pixelWidth;
		double dy = cal.pixelHeight;
		double dz = cal.pixelDepth;
		IJ.showStatus("Calculating Distance Map 3D.");
		ImageStack res;
		if (imp.isHyperStack())	
			res = filterHyperstack(imp, dx, dy, dz);
		else res = filter(imp.getStack(), dims[3], dims[4], dx, dy, dz);
//find maximum
		double maxv = 0.;
		int nums = res.getSize();
		for (int i = 1; i <= nums; i ++) {
			ImageProcessor iprc = res.getProcessor(i);
			ImageStatistics stat = ImageStatistics.getStatistics(iprc, ImageStatistics.MIN_MAX, null);
			if (maxv < stat.max) maxv = stat.max;
		}	
		ImagePlus outp = new ImagePlus("Distance Map 3D", res); 
		outp.setDimensions(dims[2], dims[3], dims[4]); // xyCzt
		outp.setCalibration(cal);
		outp.setC(1); //channel position
		outp.setDisplayRange(0., maxv);
		outp.setOpenAsHyperStack(imp.isHyperStack());
		outp.show();
		outp.updateAndDraw();
	}
	
	public static ImageStack filter(ImageStack ims, int nz, int nt, double dx, double dy, double dz) {
		
		
		class Aac1 {
			class Vert1 {
				public double x, v, stp;
			}
			public float[] data;
			public int nx, ny, nz;
			private Vert1[] stack;
			Aac1(int nx0, int ny0, int nz0, float[] data0) {
				nx = nx0; ny = ny0; nz = nz0; data = data0;
				int mxyz = (nx > ny) ? nx : ny;
				if (nz > mxyz) mxyz= nz;
				stack = new Vert1[mxyz];
				for (int i = 0; i < mxyz; i ++)
					stack[i]= new Vert1();				
			}
			public void setif(float val) {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++)
					if (data[i] > 0.)
						data[i]= val;
			}
			public void sqrt() {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++)
					if (data[i] > 0.)
						data[i]= (float) Math.sqrt(data[i]);
			}
			public void erode_by_parabola(int dim, double dx, double maxd2) {
				int i, j, k, l, m, top;
				double x, x1, v1, stp1;
				int nData, nInter, nLine, nStep;
				if (dim == 1) {
					nData= ny * nz; nInter = 1; nLine= nx;
				}
				else if (dim == 2) {
					nData= nx * nz; nInter = nx;	nLine= ny;
				}
				else {
					nData= nx * ny; nInter = nx * ny;	nLine= nz;
				}
				nStep = nInter * (nLine - 1);

				for(i= 0, l= 0; i < nData; i+= nInter) {
					for(j= 0; j < nInter; j++) {
						x= 0.0; top= -1; stp1= 0.0;
			 			m= l++;
						for(k= 0; k < nLine ; k++, m+= nInter, x+= dx) 
							if (data[m] < maxd2) {
								x1= x; v1= data[m];
								if (top >= 0) 
									stp1= (0.5 * (stack[top].x + x1 + (v1 - stack[top].v) / (x1 - stack[top].x)));
								while ((top >= 0) && (stp1 <= stack[top].stp)) {
									top--;
									if (top >= 0)  
										stp1= (0.5 * (stack[top].x + x1 + (v1 - stack[top].v) / (x1 - stack[top].x)));
								}
								top++;
								stack[top].stp= stp1;
								stack[top].x= x1;
								stack[top].v= v1;
							}
						m-= nInter; x-= dx;
						if (top >= 0) {
							while (k > 0) {
								while ((top > 0)&&(stack[top].stp > x)) top--;
								data[m]= (float) (stack[top].v + Math.pow(stack[top].x - x, 2.));
								k--; m-= nInter; x-= dx;
							}
						}
					}
					l += nStep;
				}			
			}
		}

		
		int width = ims.getWidth();
		int height = ims.getHeight();
		int size = ims.getSize();
		
//output
		ImageStack imout = ImageStack.create(width, height, size, 32); //float
		
		IJ.showStatus("Calculating Distance Map 3D.");
//utils			
		int nx = width;
		int ny = height;
		
		float img[] = new float [ny * nx * nz];//buffer
		Aac1 pimg = new Aac1(nx, ny, nz, img);
		
		for (int i = 0; i < nt; i++) {
			ims.getVoxels(0, 0, i * nz, nx, ny, nz, img);
			pimg.setif(Float.MAX_VALUE);
			
//opening by paraboloid			
	//squared ED map = erosion by paraboloid
			pimg.erode_by_parabola(1, dx, Float.MAX_VALUE);			
			pimg.erode_by_parabola(2, dy, Float.MAX_VALUE);
			pimg.erode_by_parabola(3, dz, Float.MAX_VALUE);
			
			pimg.sqrt();
			imout.setVoxels(0, 0, i * nz, nx, ny, nz, img);
		}

		return imout;
	}

	private static ImageStack filterHyperstack(ImagePlus imp, double dx, double dy, double dz) {
		int dims[] = imp.getDimensions();
		if (dims[2] == 1) {
			ImageStack stack = filter(imp.getStack(), dims[3], dims[4] , dx, dy, dz);
			return stack;
	    }
        ImagePlus[] channels = ChannelSplitter.split(imp);
        int nch = channels.length;
        for (int i = 0; i < nch; i++) {
			ImageStack stack = filter(channels[i].getStack(), imp.getNSlices(), imp.getNFrames(), dx, dy, dz);
			channels[i].setStack(stack);
		}
		ImagePlus outp = RGBStackMerge.mergeChannels(channels, false);
		ImageStack stack =  outp.getStack();
		return stack;
	}

	public void showAbout() {
		IJ.showMessage("Distance map",
			"of 3D objects"
		);
	}
	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads
	 * an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Distance_Map_3D.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the sample
		ImagePlus image = IJ.openImage//("C:/src/ImageJ/Biomat/hyperball.tif");
				("https://raw.githubusercontent.com/jiri-janacek/biomat/93f0fbe74646db4588feb99e86cb3ba35288dcc2/media/testball.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
