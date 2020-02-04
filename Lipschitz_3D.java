/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package biomat;

import ij.IJ;
import ij.plugin.*;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;

/**
 * for background subtraction in Z stack
 * GRAY8, GRAY16, GRAY32, COLOR_RGB hyperstacks.
 *
 * @author Jiri Janacek
 */
public class Lipschitz_3D implements PlugIn {

	// plugin parameters
	public double slope;
	public boolean tophat;
	
	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		if (!showDialog())
			return;

		imp.startTiming();  
		run(imp);		
		imp.updateAndDraw();
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Lipschitz 3D filter");

		// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("slope", 1.00, 2);
		gd.addCheckbox("top hat", true);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		slope = gd.getNextNumber();
		tophat = gd.getNextBoolean();

		return true;
	}

	private void run(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		double dx = cal.pixelWidth;
		double dy = cal.pixelHeight;
		double dz = cal.pixelDepth;
		IJ.showStatus("Calculating Lipschitz 3D.");
		if (imp.isHyperStack()) {
			filterHyperstack(imp, dx, dy, dz, slope, tophat);
			return;
		}
		ImageStack res = filter(imp.getStack(), imp.getStackSize(), dx, dy, dz, slope, tophat);
		imp.setStack(res);
	}
	
	public static ImageStack filter(ImageStack ims, int depth, double dx, double dy, double dz, double slope, boolean tophat) {
		if (ims.getBitDepth() == 24)
			return filterRGB(ims, depth, dx, dy, dz, slope, tophat);
		
		// get width and height
		int width = ims.getWidth();
		int height = ims.getHeight();
		int nums = ims.getSize();
		int numz = nums / depth;

		double ddx = dx * slope;
		double ddy = dy * slope;
		double ddz = dz * slope;
		double ddxy = Math.sqrt(dx * dx + dy * dy) * slope;
		double ddxz = Math.sqrt(dx * dx + dz * dz) * slope;
		double ddyz = Math.sqrt(dy * dy + dz * dz) * slope;

		ImageStack imf;
		if (ims.getBitDepth() == 32) imf = ims.duplicate();
		else imf = ims.convertToFloat();

		double p, q;
		for (int l = 0; l < numz; l++ ) {
			for (int k = l * depth; k < (l + 1) * depth; k++) {
				for (int j = 0; j < height; j++)
					for (int i = 0; i < width; i++) {
						p = ims.getVoxel(i, j, k);					
						if (i > 0) {
							q = imf.getVoxel(i - 1, j, k) + ddx;
							if (q < p) 
								p = q;
							if (j > 0) {
								q = imf.getVoxel(i - 1, j - 1, k) + ddxy;
								if (q < p) p = q;
							}
							if (k > 0) {
								q = imf.getVoxel(i - 1, j, k - 1) + ddxz;
								if (q < p) p = q;
							}
						}
						if (j > 0) {
							q = imf.getVoxel(i, j - 1, k) + ddy;
							if (q < p) 
								p = q;
							if (k > 0) {
								q = imf.getVoxel(i, j - 1, k - 1) + ddyz;
								if (q < p) p = q;
							}
						}
						if (k > 0) {
							q = imf.getVoxel(i, j, k - 1) + ddz;
							if (q < p) 
								p = q;
						}	
						imf.setVoxel(i, j, k, p);
					}
				IJ.showProgress(2 * l * depth + k, 2 * nums - 1);
			}
			for (int k = (l + 1) * depth - 1; k >= l * depth; k --) {
				for (int j = height - 1; j >= 0; j --)
					for (int i = width - 1; i >= 0; i --) {
						p = imf.getVoxel(i, j, k);					
						if (i < width - 1) {
							q = imf.getVoxel(i + 1, j, k) + ddx;
							if (q < p) 
								p = q;
							if (j < height - 1) {
								q = imf.getVoxel(i + 1, j + 1, k) + ddxy;
								if (q < p) p = q;
							}
							if (k < depth - 1) {
								q = imf.getVoxel(i + 1, j, k + 1) + ddxz;
								if (q < p) p = q;
							}
						}
						if (j < height - 1) {
							q = imf.getVoxel(i, j + 1, k) + ddy;
							if (q < p) 
								p = q;
							if (k < depth - 1) {
								q = imf.getVoxel(i, j + 1, k + 1) + ddyz;
								if (q < p) p = q;
							}
						}
						if (k < depth - 1) {
							q = imf.getVoxel(i, j, k + 1) + ddz;
							if (q < p) 
								p = q;
						}					
						imf.setVoxel(i, j, k, p);
					}
				IJ.showProgress(2 * l * depth + 2 * depth - k - 1, 2 * nums - 1);
			}
		}
		
		if (tophat) {
			for (int k = 0; k < nums; k ++)
				for (int j = 0; j < height; j ++)
					for (int i = 0; i < width; i ++)
						ims.setVoxel(i, j, k, ims.getVoxel(i, j, k) - imf.getVoxel(i, j, k));
		}
		else {
			for (int k = 0; k < nums; k++)
				for (int j = 0; j < height; j++)
					for (int i = 0; i < width; i++)
						ims.setVoxel(i, j, k, imf.getVoxel(i, j, k));
		}
		return ims;
}
	
	private static void filterHyperstack(ImagePlus imp, double dx, double dy, double dz, double slope, boolean tophat) {
		if (imp.getNChannels() == 1) {
			ImageStack stack = filter(imp.getStack(), imp.getNSlices(), dx, dy, dz, slope, tophat);
			imp.setStack(stack);
			return;
	    }
        ImagePlus[] channels = ChannelSplitter.split(imp);
        int nch = channels.length;
        for (int i = 0; i < nch; i++) {
			ImageStack stack = filter(channels[i].getStack(), imp.getNSlices(), dx, dy, dz, slope, tophat);
			channels[i].setStack(stack);
		}
		ImagePlus imp2 = RGBStackMerge.mergeChannels(channels, false);
		int dims[] = imp.getDimensions();
		imp.setImage(imp2);
		imp.setDimensions(dims[2], dims[3], dims[4]);
		imp.setC(1); //channel position
	}

	private static ImageStack filterRGB(ImageStack rgb_in, int nz, double dx, double dy, double dz, double slope, boolean tophat) {
		ImageStack[] channels = ChannelSplitter.splitRGB(rgb_in, false);
		ImageStack red = filter(channels[0], nz, dx, dy, dz, slope, tophat);
		ImageStack green = filter(channels[1], nz, dx, dy, dz, slope, tophat);
		ImageStack blue = filter(channels[2], nz, dx, dy, dz, slope, tophat);
        return RGBStackMerge.mergeStacks(red, green, blue, false);
	}

	public void showAbout() {
		IJ.showMessage("Lipschitz 3D",
			"for background subtraction in stack"
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
		Class<?> clazz = Lipschitz_3D.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://imagej.net/_images/2/2e/Capillaries_brain.zip");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
