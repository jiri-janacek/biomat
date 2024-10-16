package biomat;

import ij.IJ;
import ij.plugin.*;
//import ij.process.ImageProcessor;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;

/**
 * tensor of second order moments of heat kernel
 * GRAY images.
 *
 * @author Jiri Janacek
 */

public class Heat_Kernel_Tensor implements PlugIn {

	// plugin parameters
	public static double sigma;
	public static double time;
	public static double step;
	
	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		int type = imp.getType();
		
		if (! ((type == ImagePlus.GRAY8) || (type == ImagePlus.GRAY16)
				|| (type == ImagePlus.GRAY32))) {
			IJ.error("Heat Kernel Tensor", "unsupported image type");
			return;
		}

		if (!showDialog())
			return;

		imp.startTiming();  
		run(imp);		
		imp.updateAndDraw();
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Heat Kernel Tensor filter");

		gd.addNumericField("sigma", 10.0, 2);
//TODO parametry vypoctu
		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		sigma = gd.getNextNumber();
		
		return true;
	}

	private void run(ImagePlus imp) {
		int dims[] = imp.getDimensions();
		Calibration cal = imp.getCalibration();
		double dx = cal.pixelWidth;
		double dy = cal.pixelHeight;
		time = sigma * sigma / 2.0;
		step = 0.5 * Math.pow(dx * dy, 0.25);
		IJ.showStatus("Calculating Heat Kernel Tensor.");
		ImageStack res = filter(imp.getStack(), dx, dy, sigma);
		ImagePlus outp = new ImagePlus("Heat Kernel Tensor", res); 
		outp.setDimensions(3 * dims[2], dims[3], dims[4]);
		outp.setCalibration(cal);
		outp.setC(1); //channel position
		outp.show();
		outp.updateAndDraw();
	}
	
	public static ImageStack filter(ImageStack ims, double dx, double dy, double sigma) {
		
		int width = ims.getWidth();
		int height = ims.getHeight();
		int size = ims.getSize();
		
		ImageStack imout = ImageStack.create(width, height, 3 * size, 32); //float 3 channels

		for (int slice = 0; slice < size; slice ++) {
			
			int[] lab = new int[width * height];
			int l = 0;
			int q = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, q ++)
					if (ims.getVoxel(j, i, slice) > 0.) lab[q] = ++ l;
					else lab[q] = 0;
			
			int maxl = l;
			int nums[] = new int[maxl];
			int totnum = 0;
			l = 0;
			int p = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0) {
						int num = 0;
						if ((i > 0) && (lab[p - width] > 0)) num ++; 
						if ((i < (height - 1)) && (lab[p + width] > 0)) num ++; 
						if ((j > 0) && (lab[p - 1] > 0)) num ++; 
						if ((j < (width - 1)) && (lab[p + 1] > 0)) num ++;
						nums[l ++] = num;
						totnum += num;
					}
			double wght[] = new double[2];
			wght[0] = 1. / (dx * dx);
			wght[1] = 1. / (dy * dy);

			int sparse[] = new int[totnum];
			short wi[] = new short[totnum];
			int r = 0;
			p = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0) {
						if ((i > 0) && (lab[p - width] > 0)) {
							sparse[r] = lab[p - width] - 1;
							wi[r] = 1;
							r ++;
						}
						if ((i < (height - 1)) && (lab[p + width] > 0)) {
							sparse[r]= lab[p + width] - 1; 
							wi[r] = 1;
							r ++;
						}
						if ((j > 0) && (lab[p - 1] > 0)) {
							sparse[r] = lab[p - 1] - 1; 
							wi[r] =  0;
							r ++;
						}
						if ((j < (width - 1)) && (lab[p + 1] > 0)) {
							sparse[r] = lab[p + 1] - 1; 
							wi[r] = 0;
							r ++;
						}
					}
	
			double x0 = width / 2.;
			double y0 = height / 2.;
	
			double ux[] = new double[maxl];
			p = 0; l = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0)
						ux[l ++] = dx * (j - x0);
	
			double uy[] = new double[maxl];
			p = 0; l = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0)
						uy[l ++] = dy * (i - y0);
	
			double uxx[] = new double[maxl];
			p = 0; l = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0)
						uxx[l ++] = dx * (j - x0) * dx * (j - x0);
	
			double uxy[] = new double[maxl];
			p = 0; l = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0)
						uxy[l ++] = dx * (j - x0) * dy * (i - y0);
	
			double uyy[] = new double[maxl];
			p = 0; l = 0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++, p ++)
					if (lab[p] > 0)
						uyy[l ++] = dy * (i - y0) * dy * (i - y0);		
	
			IJ.showProgress(slice * 5, size * 5);
			solveu(maxl, nums, sparse, wi, wght, ux);
			IJ.showProgress(slice * 5 + 1, size * 5);
			solveu(maxl, nums, sparse, wi, wght, uy);
			IJ.showProgress(slice * 5 + 2, size * 5);
			solveu(maxl, nums, sparse, wi, wght, uxx);
			IJ.showProgress(slice * 5 + 3, size * 5);
			solveu(maxl, nums, sparse, wi, wght, uxy);
			IJ.showProgress(slice * 5 + 4, size * 5);
			solveu(maxl, nums, sparse, wi, wght, uyy);
			IJ.showProgress(slice * 5 + 5, size * 5);
	
			p = 0;
			for (int i = 0; i < height; i ++) 
				for (int j = 0; j < width; j ++, p ++)
				if (lab[p] > 0) {
					l = lab[p] - 1;
					imout.setVoxel(j, i, slice * 3, uxx[l] - ux[l] * ux[l]);
					imout.setVoxel(j, i, slice * 3 + 1, uxy[l] - ux[l] * uy[l]);
					imout.setVoxel(j, i, slice * 3 + 2, uyy[l] - uy[l] * uy[l]);
				}
			
		}
		return imout;
	}
	
	public static void solveu(int maxl, int[] nums, int[] sparse, short wi[], double wght[], double[] u) {
		//Crank-Nicolson solution of heat equation
		double u1[] = new double[maxl];
		double lapu[] = new double[maxl];
		double sumw[] = new double[maxl];
		for (double t = 0.; t < time; t += step) {
			double norm = 0.;
			int a = 0; //sparse;
			for (int i = 0; i < maxl; i ++) {
				norm += Math.abs(u[i]);
				lapu[i] = 0.;
				sumw[i] = 0.;
				for (int j = 0; j < nums[i]; j ++, a ++) {
					double weight = wght[wi[a]];
					lapu[i] += weight * u[sparse[a]];
					sumw[i] += weight;
				}
				lapu[i] -= sumw[i] * u[i];
			}
			
			for (int i = 0; i < maxl; i ++) 
				u1[i] = u[i];
			double dist;
			//Gauss-Seidel iterations
			do {
				dist = 0.;
				a = 0; //sparse;
				for (int i = 0; i < maxl; i ++) {
					double sum = 0.;
					for (int j = 0; j < nums[i]; j ++, a ++) 
						sum += wght[wi[a]] * u[sparse[a]];
					double pom = u[i];
					u[i] = (sum + lapu[i] + 2. * u1[i] / step) / (2. / step + sumw[i]);
					dist += Math.abs(u[i] - pom);
				}
			}
			while (dist > (1e-6 * norm));

		}
	}

	public void showAbout() {
		IJ.showMessage("Heat Kernel Tensor",
			"for structure tensor calculation on image"
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
		Class<?> clazz = Heat_Kernel_Tensor.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://raw.githubusercontent.com/jiri-janacek/biomat/master/media/MAX_2_4cortexa1.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
