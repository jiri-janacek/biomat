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
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;

/**
 * for detection of fibers in stack
 * GRAY8, GRAY16, GRAY32 images.
 *
 * @author Jiri Janacek
 */
public class Tensor_Line_3D implements PlugIn {

	// plugin parameters
	public double sigma;
	
	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		int type = imp.getType();
		
		if (!((type == ImagePlus.GRAY8)||(type == ImagePlus.GRAY16)||
				(type == ImagePlus.GRAY32))) {
			IJ.error("Tensor Line 3D", "unsupported image type");
			return;
		}
		if (!((imp.getNSlices() > 1)&&(imp.getNSlices()==imp.getStackSize()))) {
			IJ.error("Tensor Line 3D", "only single z-stack is supported");
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
		GenericDialog gd = new GenericDialog("Tensor Line 3D filter");

		// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("sigma (units)", 3.00, 2);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		sigma = gd.getNextNumber();

		return true;
	}

	private void run(ImagePlus imp) {
		Calibration cal= imp.getCalibration();
		double dx= cal.pixelWidth;
		double dy= cal.pixelHeight;
		double dz= cal.pixelDepth;
		ImageStack res = filter(imp.getStack(), dx, dy, dz, sigma);
//		imp.setStack(res);  //in place
//find maximum
		double maxv = 0.;
		int nums = res.getSize();
		for (int i = 1; i <= nums; i++) {
			ImageProcessor iprc= res.getProcessor(i);
			ImageStatistics stat= ImageStatistics.getStatistics(iprc, ImageStatistics.MIN_MAX, null);
			if (maxv < stat.max) maxv = stat.max;
		}
//new image		
		ImagePlus outp = new ImagePlus("Tensor Line 3D Filtered", res); 
		outp.setDisplayRange(0., maxv);
		outp.show();
		outp.updateAndDraw();
//new image
	}
	
	public static ImageStack filter(ImageStack ims, double dx, double dy, double dz, double sigma) {
		
		class Aac1 {
			public float[] data;
			public int nx, ny, nz;
			Aac1(int nx0, int ny0, int nz0, float[] data0) {
				nx= nx0; ny= ny0; nz= nz0; data= data0;
			}
			public float get(int i, int j, int k) {
				if ((i>=0)&&(i<nx)&&(j>=0)&&(j<ny)&&(k>=0)&&(k<nz))
					return data[i+nx*(j+ny*k)];
				else return (float) 0.;
			}
		}

		class Aac6 {
			public float[] data;
			public int nx, ny, nz;
			Aac6(int nx0, int ny0, int nz0, float[] data0) {
				nx= nx0; ny= ny0; nz= nz0; data= data0;
			}
			public void get(int i, int j, int k, float[] val) {
				if ((i>=0)&&(i<nx)&&(j>=0)&&(j<ny)&&(k>=0)&&(k<nz))
					for (int c= 0; c < 6; c++)
						val[c]= data[c+6*(i+nx*(j+ny*k))];
			}
			public void add(int i, int j, int k, float[] val, double coeff) {
				if ((i>=0)&&(i<nx)&&(j>=0)&&(j<ny)&&(k>=0)&&(k<nz))
					for (int c= 0; c < 6; c++)
						data[c+6*(i+nx*(j+ny*k))]+= val[c]*coeff;
			}
			public void zero() {
				int size= 6 * nx * ny * nz;
				for (int i= 0; i < size; i++) data[i]= (float) 0.;
			}
		}

		// get width and height
		int width = ims.getWidth();
		int height = ims.getHeight();
		int depth= ims.getSize();

		double sx= sigma / dx;
		double sy= sigma / dy;
		double sz= sigma / dz;

//algorithmus TODO kosticky
		ImageStack imf;
		if (ims.getBitDepth() == 32) imf = ims.duplicate();
		else imf = ims.convertToFloat();
		
		ImagePlus pimf= new ImagePlus("", imf);
		GaussianBlur3D.blur(pimf, sx, sy, sz);

		IJ.showStatus("Calculating Tensor Line 3D.");
		double sigma2= sigma * sigma;
		int nx= width;
		int ny= height;
		int nz= depth;
		float img[]= new float [nz * ny * nx];
		imf.getVoxels(0, 0, 0, nx, ny, nz, img);
		Aac1 pimg= new Aac1(nx, ny, nz, img);
		float tens[]= new float [nz * ny * nx * 6];
		Aac6 pten= new Aac6(nx, ny, nz, tens);
		pten.zero();
		for (int k = 0; k < nz; k++) {
			int k0= k - 1; if (k0 < 0) k0= 0;
			int k2= k + 1; if (k2 == nz) k2= nz - 1;
			for (int j = 0; j < ny; j++) {
				int j0= j - 1; if (j0 < 0) j0= 0;
				int j2= j + 1; if (j2 == ny) j2= ny - 1;
				for (int i = 0; i < nx; i++) {
					int i0= i - 1; if (i0 < 0) i0= 0;
					int i2= i + 1; if (i2 == nx) i2= nx - 1;
					
					double x000= pimg.get(i0, j0, k0);
					double x001= pimg.get(i, j0, k0);
					double x002= pimg.get(i2, j0, k0);
					double x010= pimg.get(i0, j, k0);
					double x011= pimg.get(i, j, k0);
					double x012= pimg.get(i2, j, k0);
					double x020= pimg.get(i0, j2, k0);
					double x021= pimg.get(i, j2, k0);
					double x022= pimg.get(i2, j2, k0);
					double x100= pimg.get(i0, j0, k);
					double x101= pimg.get(i, j0, k);
					double x102= pimg.get(i2, j0, k);
					double x110= pimg.get(i0, j, k);
					double x111= pimg.get(i, j, k);
					double x112= pimg.get(i2, j, k);
					double x120= pimg.get(i0, j2, k);
					double x121= pimg.get(i, j2, k);
					double x122= pimg.get(i2, j2, k);
					double x200= pimg.get(i0, j0, k2);
					double x201= pimg.get(i, j0, k2);
					double x202= pimg.get(i2, j0, k2);
					double x210= pimg.get(i0, j, k2);
					double x211= pimg.get(i, j, k2);
					double x212= pimg.get(i2, j, k2);
					double x220= pimg.get(i0, j2, k2);
					double x221= pimg.get(i, j2, k2);
					double x222= pimg.get(i2, j2, k2);
					
					double xs00= x000+x100+x200;
//					double xs01= x001+x101+x201;
					double xs02= x002+x102+x202;
//					double xs10= x010+x110+x210;
//					double xs11= x011+x111+x211;
//					double xs12= x012+x112+x212;
					double xs20= x020+x120+x220;
//					double xs21= x021+x121+x221;
					double xs22= x022+x122+x222;

					double x0s0= x000+x010+x020;
					double x0s1= x001+x011+x021;
					double x0s2= x002+x012+x022;
					double x1s0= x100+x110+x120;
					double x1s1= x101+x111+x121;
					double x1s2= x102+x112+x122;
					double x2s0= x200+x210+x220;
					double x2s1= x201+x211+x221;
					double x2s2= x202+x212+x222;
					
					double x00s= x000+x001+x002;
					double x01s= x010+x011+x012;
					double x02s= x020+x021+x022;
					double x10s= x100+x101+x102;
					double x11s= x110+x111+x112;
					double x12s= x120+x121+x122;
					double x20s= x200+x201+x202;
					double x21s= x210+x211+x212;
					double x22s= x220+x221+x222;
					
					double ddx= (x0s2+x1s2+x2s2-x0s0-x1s0-x2s0)/(18.*dx);
					double ddy= (x02s+x12s+x22s-x00s-x10s-x20s)/(18.*dy);
					double ddz= (x20s+x21s+x22s-x00s-x01s-x02s)/(18.*dz);
					double dxx= (x0s2+x1s2+x2s2-2.*(x0s1+x1s1+x2s1)+x0s0+x1s0+x2s0)/(9.*dx*dx);
					double dxy= (xs22-xs20-xs02+xs00)/(12.*dx*dy);
					double dxz= (x2s2-x2s0-x0s2+x0s0)/(12.*dx*dz);
					double dyy= (x02s+x12s+x22s-2.*(x01s+x11s+x21s)+x00s+x10s+x20s)/(9.*dy*dy);
					double dyz= (x22s-x20s-x02s+x00s)/(12.*dy*dz);
					double dzz= (x20s+x21s+x22s-2.*(x10s+x11s+x12s)+x00s+x01s+x02s)/(9.*dz*dz);
					
					float qten[]= new float [6];
					if (x111 > 0.01) { //TODO adapt
						qten[0]= (float) (dxx + x111 / sigma2 - ddx * ddx / x111);  
						qten[1]= (float) (dxy - ddx * ddy / x111);
						qten[2]= (float) (dxz - ddx * ddz / x111);
						qten[3]= (float) (dyy + x111 / sigma2 - ddy * ddy / x111);  
						qten[4]= (float) (dyz - ddy * ddz / x111);
						qten[5]= (float) (dzz + x111 / sigma2 - ddz * ddz / x111);  

						double si= sigma2 * ddx / x111;
						double sj= sigma2 * ddy / x111;
						double sk= sigma2 * ddz / x111;
						int di, dj, dk;
						di= (int) Math.floor(si);
						dj= (int) Math.floor(sj);
						dk= (int) Math.floor(sk);
						double ai= si - di;
						double aj= sj - dj;
						double ak= sk - dk;
						pten.add(i+di, j+dj, k+dk, qten, (1.-ai)*(1.-aj)*(1.-ak));
						pten.add(i+di+1, j+dj, k+dk, qten, ai*(1.-aj)*(1.-ak));
						pten.add(i+di, j+dj+1, k+dk, qten, (1.-ai)*aj*(1.-ak));
						pten.add(i+di+1, j+dj+1, k+dk, qten, ai*aj*(1.-ak));
						pten.add(i+di, j+dj, k+dk+1, qten, (1.-ai)*(1.-aj)*ak);
						pten.add(i+di+1, j+dj, k+dk+1, qten, ai*(1.-aj)*ak);
						pten.add(i+di, j+dj+1, k+dk+1, qten, (1.-ai)*aj*ak);
						pten.add(i+di+1, j+dj+1, k+dk+1, qten, ai*aj*ak);

//						pten.add(i, j, k, qten, 1.);//bez shiftu
					}
				}
			}
			IJ.showProgress(k, 2*nz-1);
		}
		
		double sqrt3= Math.sqrt(3.);
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					float qten[]= new float [6];
					pten.get(i, j, k, qten);

					double b= - qten[0] - qten[3] - qten[5];
					double c= qten[0] * qten[3] + qten[0] * qten[5] + qten[3] * qten[5]
						- qten[1] * qten[1] - qten[2] * qten[2] - qten[4] * qten[4];
					double d= - 2.0 * qten[1] * qten[2] * qten[4] - qten[0] * qten[3] * qten[5]
						+ qten[0] * qten[4] * qten[4] + qten[3] * qten[2] * qten[2] + qten[5] * qten[1] * qten[1];
					double p= c - b * b / 3.0;
					double q= d - b * c / 3.0 + 2.0 * b * b * b / 27.0;
					double disc= - 27. * q * q - 4. * p * p * p;
					double lambda1;
					double lambda2;
					double lambda3;
					if (disc > 0.0) {					
						double r= Math.sqrt(- p / 3.0);
						double s= - 0.5 * q / (r * r * r);
						double fi3;
						if (s <= -1.0) fi3= Math.PI / 3.0;
						else if (s >= 1.0) fi3= 0.0;
						else fi3= Math.acos(s) / 3.0;
						double cf3= Math.cos(fi3);
						double s3sf3= sqrt3 * Math.sin(fi3);
						lambda1= 2. * r * cf3 - b / 3.0;
						lambda2= - r * (cf3 - s3sf3) - b / 3.0;
						lambda3= - r * (cf3 + s3sf3) - b / 3.0;
					}
					else if (b < 0.) {
						lambda1= lambda2= lambda3= Math.sqrt(-b);
					}
					else {
						lambda1= lambda2= lambda3= 0.;
					}

					//double val= lambda1;
					double lambda= (Math.sqrt(lambda1)+Math.sqrt(lambda2)+Math.sqrt(lambda3))/3.;
					double val= Math.sqrt(lambda1*((lambda1+lambda2+lambda3)/3.-lambda*lambda));
					if (!(val >= 0.)) val= 0.;
					
					imf.setVoxel(i, j, k, val);
				}
		IJ.showProgress(nz+k, 2*nz-1);
	}

		return imf;
}
	
	public void showAbout() {
		IJ.showMessage("Tensor Line 3D",
			"for detection of fibers in stack"
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
		Class<?> clazz = Tensor_Line_3D.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://imagej.net/_images/e/e8/Capillaries_adipose.zip");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
