package biomat;

import java.awt.Point;

import ij.IJ;
import ij.plugin.*;
import ij.plugin.filter.Analyzer;
import ij.process.*;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.ResultsTable;

/**
 * average of 2D tensor
 * 3 channel GRAY images.
 *
 * @author Jiri Janacek
 */

public class Tensor_Average implements PlugIn {

	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		int type = imp.getType();
		int dims[] = imp.getDimensions();
		
		if (! ((dims[2] % 3) == 0) //xyCzt 
				|| ! ((type == ImagePlus.GRAY8) || (type == ImagePlus.GRAY16)
				|| (type == ImagePlus.GRAY32))) {
			IJ.error("Tensor_Average", "unsupported image type");
			return;
		}

		imp.startTiming();  
		run(imp);		
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private void run(ImagePlus imp) {
		IJ.showStatus("Calculating Tensor Average.");
		int width = imp.getWidth();
		int height = imp.getHeight();
		Roi contour= imp.getRoi();
		int pos= imp.getCurrentSlice();
		ImageStack ims = imp.getStack();
		ImageProcessor ip1= ims.getProcessor(3 * ((pos - 1) / 3) + 1);
		ImageProcessor ip2= ims.getProcessor(3 * ((pos - 1) / 3) + 2);
		ImageProcessor ip3= ims.getProcessor(3 * ((pos - 1) / 3) + 3);
        int num= 0; 
        double a[][]= new double[2][2];
        a[0][0] = 0.;
  	  	a[0][1] = 0.;
  	  	a[1][1] = 0.;

		if ((contour != null) && contour.isArea()) {
          for (Point p : contour) {
        	  num++;
        	  a[0][0] += ip1.getPixelValue(p.x, p.y);
        	  a[0][1] += ip2.getPixelValue(p.x, p.y);
        	  a[1][1] += ip3.getPixelValue(p.x, p.y);
          }
        }
		else {
			for (int i = 0; i < height; i ++)
				for (int j = 0; j < width; j ++) {
					 num++;
					 a[0][0] += ip1.getPixelValue(j, i);
					 a[0][1] += ip2.getPixelValue(j, i);
					 a[1][1] += ip3.getPixelValue(j, i);				
				}
		}
		if (num > 0) {
			a[0][0] /= num;
			a[0][1] /= num;
			a[1][1] /= num;
		}

		a[1][0] = a[0][1];
		double angle = Math.atan2(2.0 * a[0][1], a[1][1] - a[0][0]);
		angle += Math.PI;
		angle *= 90.0 / Math.PI;
		double val = a[0][0] + a[1][1];
		double pom = (a[1][1] - a[0][0]) * (a[1][1] - a[0][0]) + 4 * a[0][1] * a[0][1];
		double shape = 0.;
		if (val > 0.0) {
			double coh = Math.sqrt(pom) / val;
			shape = (1. + coh) / (1. - coh); 
		}
		
	    ResultsTable rt=Analyzer.getResultsTable();
	    rt.incrementCounter();
	    rt.addValue("Slice", (pos + 2) / 3);
	    rt.addValue("Pixels", num);
//	    rt.addValue("A11", a[0][0]);
//	    rt.addValue("A12", a[0][1]);
//	    rt.addValue("A22", a[1][1]);
	    rt.addValue("Value", val);
	    rt.addValue("Elongation", shape);
	    rt.addValue("Angle", angle);
	    rt.show("Results");
		
	}

	public void showAbout() {
		IJ.showMessage("Tensor Average",
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
		Class<?> clazz = Tensor_Average.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("D:/data/tif/MAX_2_4cortexa1_tens.tif");
				//("https://imagej.net/_images/2/2e/Capillaries_brain.zip");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
