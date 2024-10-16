package biomat;

import ij.IJ;
import ij.plugin.*;
import ij.process.*;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
//import ij.gui.GenericDialog;
//import ij.measure.Calibration;
import ij.measure.Calibration;

/**
 * color coding of 2D tensor
 * 3 channel GRAY images.
 *
 * @author Jiri Janacek
 */

public class Tensor_Color_Coding implements PlugIn {

	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		int type = imp.getType();
		int dims[] = imp.getDimensions();
		
		if (! ((dims[2] % 3) == 0) //xyCzt 
				|| ! ((type == ImagePlus.GRAY8) || (type == ImagePlus.GRAY16)
				|| (type == ImagePlus.GRAY32))) {
			IJ.error("Tensor_Color_Coding", "unsupported image type");
			return;
		}

		imp.startTiming();  
		run(imp);		
		imp.updateAndDraw();
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private void run(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		int dims[] = imp.getDimensions();
		IJ.showStatus("Calculating Tensor Color Coding.");
		ImageStack res = filter(imp.getStack());
		ImagePlus outp = new ImagePlus("Color Coded Tensor", res); 
		outp.setDimensions(dims[2] / 3, dims[3], dims[4]);
		outp.setCalibration(cal);
		outp.setC(1); //channel position
		outp.show();
		outp.updateAndDraw();
	}
	
	public static ImageStack filter(ImageStack ims) {
		
		int width = ims.getWidth();
		int height = ims.getHeight();
		int size = ims.getSize() / 3;
		
		ImageStack imout = ImageStack.create(width, height, size, 24); //RGB

		for (int slice = 0; slice < size; slice ++) {
			
			ImageProcessor iop= imout.getProcessor(slice + 1);
			
			double maxval= 0.0;
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++) {
					double a11 = ims.getVoxel(j, i, 3 * slice);
					double a22 = ims.getVoxel(j, i, 3 * slice + 2);
					double val= a11 + a22;
					if (val > maxval) maxval= val;
				}
			
			if (maxval > 0.) {
				double fact = 255. / maxval;		
				for (int i = 0; i < height; i ++)
					for (int j = 0; j < width; j ++) {
						double hue, sat, val;
						float rgb[] = new float[3];
						double a11 = ims.getVoxel(j, i, 3 * slice);
						double a12 = ims.getVoxel(j, i, 3 * slice + 1);
						double a22 = ims.getVoxel(j, i, 3 * slice + 2);
//						if ((a11 != 0) || (a12 != 0) || (a22 != 0)) {
							hue = Math.atan2(2.0 * a12, a22 - a11);
							if (hue < 0.0) hue += 2 * Math.PI;
							val = a11 + a22;
							double pom = (a22 - a11) * (a22 - a11) + 4 * a12 * a12;
							if ((pom > 0.0) && (val > 0.0)) sat = Math.sqrt(pom) / val;
							else sat = 0.0;
							if (sat > 1.) sat = 1.;
							hsv_to_rgb((float)(180.0 * hue / Math.PI), (float)sat, (float)val, rgb);
							int irgb[] = new int[3];
							for (int c = 0; c < 3; c++) irgb[c] = (int)Math.floor(fact * rgb[c]); 
							iop.putPixel(j, i, irgb);
//						}
					}			
			}
		}
		return imout;
	}
	
	static void hsv_to_rgb(float h, float s, float v, float[] rgb)
	{
		if(s ==  0)	{
			rgb[0] = v; rgb[1] = v;	rgb[2] = v;
		}
		else
		{
			if(h == 360.0f)
				h = 0.0f;
			
			h /= 60.0f;

			int i = (int)Math.floor(h);
			
			float f = h - (float) i;

			float p = v * (1.0f -  s);
			float q = v * (1.0f - (s * f));
			float t = v * (1.0f - (s * (1.0f - f)));

			switch(i) {
				case 0:
					rgb[0] = v; rgb[1] = t; rgb[2] = p;
					break;
				case 1:
					rgb[0] = q; rgb[1] = v; rgb[2] = p;
					break;
				case 2:
					rgb[0] = p; rgb[1] = v; rgb[2] = t;
					break;
				case 3:
					rgb[0] = p; rgb[1] = q; rgb[2] = v;
					break;
				case 4:
					rgb[0] = t; rgb[1] = p; rgb[2] = v;
					break;
				case 5:
					rgb[0] = v; rgb[1] = p; rgb[2] = q;
					break;
			}
		}
	}

	public void showAbout() {
		IJ.showMessage("Tensor Color Coding",
			"for presenting structure tensor image"
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
		Class<?> clazz = Tensor_Color_Coding.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://github.com/jiri-janacek/biomat/raw/0381683d9a39c3bc723f8e9013d543fa50c712f0/media/MAX_2_4cortexa1_tens.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
