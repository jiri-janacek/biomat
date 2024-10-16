
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
import ij.gui.GenericDialog;

/**
 * linear stack contrast enhancement of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Jiri Janacek
 */
public class StackLin_Contrast implements PlugIn {
	// plugin parameters
	public double first;
	public double last;

	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus image = IJ.getImage();
		int type = image.getType();
		
		if ((type != ImagePlus.GRAY8)&&(type != ImagePlus.GRAY16)&&
				(type != ImagePlus.GRAY32)&&(type != ImagePlus.COLOR_RGB)) {
			IJ.error("StackLin Contrast", "image type not supported");
			return;
		}
		if (! showDialog())
			return;

		image.startTiming();  
		// slice numbers start with 1 for historical reasons
		for (int i = 1; i <= image.getStackSize(); i ++) {
			int pos[] = image.convertIndexToPosition(i);
			double coeff;
			if (image.getNSlices() == 1) 
				coeff = first;
			else
				coeff = ((double)(image.getNSlices() - pos[1]) * first 
						+ (double)(pos[1] - 1) * last) / (double)(image.getNSlices() - 1);
				image.getStack().getProcessor(i).multiply(coeff);
		}
		
		image.updateAndDraw();
		IJ.showTime(image, image.getStartTime(), "", image.getStackSize());
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("StackLin Contrast Adjustment");

		// default value is 1.00, 2 digits right of the decimal point
		gd.addNumericField("first", 1.00, 2);
		gd.addNumericField("last", 1.00, 2);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		first = gd.getNextNumber();
		last = gd.getNextNumber();

		return true;
	}

	/**
	 * Process an image.
	 * <p>
	 * Please provide this method even if {@link ij.plugin.filter.PlugInFilter} does require it;
	 * the method {@link ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)} can only
	 * handle 2-dimensional data.
	 * </p>
	 * <p>
	 * If your plugin does not change the pixels in-place, make this method return the results and
	 * change the {@link #setup(java.lang.String, ij.ImagePlus)} method to return also the
	 * <i>DOES_NOTHING</i> flag.
	 * </p>
	 *
	 * @param image the image (possible multi-dimensional)
	 */
	public void process(ImagePlus image) {
		int type = image.getType();
		
		if ((type == ImagePlus.GRAY8) || (type == ImagePlus.GRAY16)
				|| (type == ImagePlus.GRAY32) || (type == ImagePlus.COLOR_RGB))
		// slice numbers start with 1 for historical reasons
			for (int i = 1; i <= image.getStackSize(); i ++) {
				int pos[]= image.convertIndexToPosition(i);
				double coeff;
				if (image.getNSlices() == 1) 
					coeff = first;
				else
					coeff = ((double)(image.getNSlices()-pos[1]) * first + (double)(pos[1] - 1) * last) / (double)(image.getNSlices() - 1);
				if ((type == ImagePlus.GRAY8)||(type == ImagePlus.GRAY16)
						|| (type == ImagePlus.GRAY32) || (type == ImagePlus.COLOR_RGB))
					image.getStack().getProcessor(i).multiply(coeff);
			}
		else {
			throw new RuntimeException("image type not supported");
		}
	}

	public void showAbout() {
		IJ.showMessage("StackLin Contrast Adjustment",
			"for improvement of declining contrast in stack"
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
		Class<?> clazz = StackLin_Contrast.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://imagej.github.io/media/plugins/capillaries-brain.zip");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
