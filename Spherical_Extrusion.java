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
 * Spherical extrusion
 * BINARY images.
 *
 * @author Jiri Janacek
 */

public class Spherical_Extrusion implements PlugIn {

	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
        ImageProcessor ip = imp.getProcessor();
        if (!ip.isBinary()) {
            IJ.error("8-bit binary image (0 and 255) required.");
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
		IJ.showStatus("Calculating Spherical extrusion");
		ImageStack res = filter(imp.getStack(), dx, dy);
//find maximum
		double maxv = 0.;
		int nums = res.getSize();
		for (int i = 1; i <= nums; i ++) {
			ImageProcessor iprc = res.getProcessor(i);
			ImageStatistics stat = ImageStatistics.getStatistics(iprc, ImageStatistics.MIN_MAX, null);
			if (maxv < stat.max) maxv = stat.max;
		}	
		ImagePlus outp = new ImagePlus("Spherical extrusion", res); 
		outp.setDimensions(dims[2], dims[3], dims[4]);
		outp.setCalibration(cal);
		outp.setDisplayRange(0., maxv);
		outp.show();
		outp.updateAndDraw();
	}

	public static ImageStack filter(ImageStack ims, double dx, double dy) {
		
		
		class Aac1 {
			class Vert1 {
				public double x, v, stp;
			}
			public float[] data;
			public int nx, ny;
			private Vert1[] stack;
			Aac1(int nx0, int ny0, float[] data0) {
				nx = nx0; ny = ny0; data = data0;
				int mxy = (nx > ny) ? nx : ny;
				stack = new Vert1[mxy];
				for (int i = 0; i < mxy; i ++)
					stack[i]= new Vert1();				
			}
			public void setif(float val) {
				int size = nx * ny;
				for (int i = 0; i < size; i ++)
					if (data[i] > 0.)
						data[i]= val;
			}
			public void sqrt() {
				int size = nx * ny;
				for (int i = 0; i < size; i ++)
					if (data[i] > 0.)
						data[i]= (float) Math.sqrt(data[i]);
			}
			public void inv() {
				int size = nx * ny;
				for (int i = 0; i < size; i ++)
					data[i] = -	data[i];
			}
			public void erode_by_parabola(int dim, double dx, double maxd2) {
				int i, j, k, l, m, top;
				double x, x1, v1, stp1;
				int nData, nInter, nLine, nStep;
				if (dim == 1) {
					nData= ny; nInter = 1; nLine= nx;
				}
				else {
					nData= nx; nInter = nx;	nLine= ny;
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
		
		IJ.showStatus("Calculating Spherical extrusion.");
//utils			
		int nx = width;
		int ny = height;
		
		float img[] = new float [ny * nx];//buffer
		Aac1 pimg = new Aac1(nx, ny, img);
		
		for (int slice = 0; slice < size; slice ++) {
			ims.getVoxels(0, 0, slice, nx, ny, 1, img);
			pimg.setif(Float.MAX_VALUE);
			
//opening by paraboloid			
		//squared ED map = erosion by paraboloid
			pimg.erode_by_parabola(1, dx, Float.MAX_VALUE);			
			pimg.erode_by_parabola(2, dy, Float.MAX_VALUE);
			
		//dilation by paraboloid
			pimg.inv();
			pimg.erode_by_parabola(1, dx, Float.MAX_VALUE);			
			pimg.erode_by_parabola(2, dy, Float.MAX_VALUE);
			pimg.setif((float)0.);
			pimg.inv();
			
			pimg.sqrt();
			imout.setVoxels(0, 0, slice, nx, ny, 1, img);
			IJ.showProgress(slice, size);
		}
		return imout;
	}

	public void showAbout() {
		IJ.showMessage("Spherical extrusion",
			"for calculation of volume of objects from 2D projections"
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
		Class<?> clazz = Spherical_Extrusion.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the sample
		ImagePlus image = IJ.openImage("https://raw.githubusercontent.com/jiri-janacek/biomat/c3f75436ccf4b863dbdf6267a352b129b28a89a7/media/simobjinv.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
