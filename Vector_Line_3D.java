
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
import ij.process.*;
import ij.gui.GenericDialog;
import ij.measure.Calibration;

class Vec3d {
//3d vectors utilities
	public static void substvect3d(double x, double y, double z, double v[]) {
		v[0] = x; v[1] = y; v[2] = z;
	}
	public static void conv2vect3d(double x1[], double x2[], double a, double x[]) {
		int j;
		
		for (j = 0; j < 3; j ++) 
		x[j] = a * x1[j] + (1.0 - a) * x2[j];
	}

	public static void copyvect3d(double u[], double v[]) {
		int j;
		
		for (j = 0; j < 3; j ++) 
			v[j] = u[j];
	}

	public static void normalize3d(double x[]) {
		int j;
		double s;

		for (j = 0, s = 0.0; j < 3; j ++) s += x[j] * x[j];
		s = Math.sqrt(s);
		for (j = 0; j < 3; j ++) x[j] /= s;
	}
	public static void axis_to_matrix(double a[], double phi, double m[][])
	{
	    double sp = Math.sin(phi);
	    double cp = Math.cos(phi);
		int i, j, k;
		for (i = 0; i < 3; i ++) {
			for (j = 0; j < 3; j ++) {
				m[i][j] = (1 - cp) * a[i] * a[j];
				for (k = 0; k < 3; k++) 
					m[i][j] += 0.5 * (j - k) * (k - i) * (i - j) * a[k] * sp;
			}
			m[i][i] += cp;
		}
	}

	public static void mult_mat_vec(double mat[][],double v[]) {
		double[] temp = new double[3];
		int i, j;
		for (i = 0; i < 3; i ++) 
			for (j = 0, temp[i] = 0; j < 3; j ++)
				temp[i] += v[j] * mat[j][i];
		for (i = 0; i < 3; i ++) v[i] = temp[i];
	}
}

class IcosaDirs {
	public double[][] dirs;
	IcosaDirs(int divs) {
		int num_nodes, vrch;
		double dp, dq;
		int i, j, k, p1, r1, p2, r2, p3, r3, 
			t1, t2, t3;
		double[][] x;
		int[] connec0 = 
			{ 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5,
			  0, 5, 1, 1, 6, 2, 2, 6, 7, 2, 7, 3,
			  3, 7, 8, 3, 8, 4, 4, 8, 9, 4, 9, 5,
			  5, 9,10, 5,10, 1, 1,10, 6, 6,11, 7,
			  7,11, 8, 8,11, 9, 9,11,10,10,11, 6 };
		int[] edges =
			{ 0, 1, 0, 2, 0, 3, 0, 4, 0, 5,
			  1, 2, 2, 3, 3, 4, 4, 5, 5, 1,
			  1,10, 1, 6, 2, 6, 2, 7, 3, 7,
			  3, 8, 4, 8, 4, 9, 5, 9, 5,10,
			  6, 7, 7, 8, 8, 9, 9,10,10, 6,
			  6,11, 7,11, 8,11, 9,11,10,11 };

		num_nodes = 10 * divs * divs + 2;

		x = new double[num_nodes][3];
		dirs = new double[num_nodes / 2][3];

	//vertices
		Vec3d.substvect3d(0.0, 0.0, 1.0, x[0]);
		if (divs <= 0) 
			Vec3d.substvect3d(0.0, 0.0, -1.0, x[1]);
		else {
			Vec3d.substvect3d(0.0, 0.0, -1.0, x[11]);
			dp = 1.0 / Math.sqrt(5.0); dq = Math.sqrt(1.0 - dp * dp);
			for (i = 0; i < 5; i ++) {
				Vec3d.substvect3d(dq * Math.cos(2.0 * Math.PI * i / 5.0), 
					dq * Math.sin(2.0 * Math.PI * i / 5.0),
					dp, x[1 + i]);
				Vec3d.substvect3d(dq * Math.cos(2.0 * Math.PI * (i + 0.5) / 5.0), 
					dq * Math.sin(2.0 * Math.PI * (i + 0.5) / 5.0),
					- dp, x[6 + i]);
			}

	//edges
			if (divs > 1) for (i = 0; i < 30; i ++) {
				vrch = 11 + i * (divs - 1);
				for (j = 1; j < divs; j ++) 
					Vec3d.conv2vect3d(x[edges[2 * i + 1]], x[edges[2 * i]], (float)j / (float)divs, x[vrch + j]);
			}

			vrch = 12 + 30 * (divs - 1);

	//planes
			if (divs > 2) for (i = 0; i < 20; i ++) {
				t1 = connec0[3 * i];
				t2 = connec0[3 * i + 1];
				t3 = connec0[3 * i + 2];
				j = 0;
				while (!((edges[2 * j] == t1) && (edges[2 * j + 1] == t2))
				&& ! ((edges[2 * j] == t2) && (edges[2 * j + 1] == t1))) j++;
				p1 = j;
				if (edges[2 * p1] == t1) {
					p3 = 1; p2 = 12 +  p1 * (divs - 1);
				}
				else {
					p3 = -1; p2 = 12 + (p1 + 1) * (divs - 1) - 1;
				}
				j = 0;
				while (!((edges[2 * j] == t1) && (edges[2 * j + 1] == t3))
				&& ! ((edges[2 * j] == t3) && (edges[2 * j + 1] == t1))) j ++;
				r1 = j;
				if (edges[2 * r1] == t1) {
					r3 = 1; r2 = 12 +  r1 * (divs - 1);
				}
				else {
					r3 = -1; r2 = 12 + (r1 + 1) * (divs - 1) - 1;
				}
				j = 0;
				while (!((edges[2 * j] == t2) && (edges[2 * j + 1] == t3))
				&& !((edges[2 * j] == t3) && (edges[2 * j + 1] == t2))) j ++;

				for (j = 2; j < divs; j ++) {
					p2 += p3;
					r2 += r3;
					Vec3d.conv2vect3d(x[r2], x[p2], 1.0/(float)j, x[vrch]);
					vrch ++;
					for (k = 2; k < j; k ++) {
						Vec3d.conv2vect3d(x[r2], x[p2], (float)k/(float)j, x[vrch]);
						vrch ++;
					}
				}
			}
		}

	//normalize
		for (i = 0; i < num_nodes; i++) Vec3d.normalize3d(x[i]);

	//symmetrize
		double eps = 1e-6;
	    for (i= 0, j = 0; i < num_nodes; i ++)
			if ((x[i][2]>eps) || ((Math.abs(x[i][2]) <= eps) && (x[i][1] > eps))
				|| ((Math.abs(x[i][2]) <= eps) && (Math.abs(x[i][1]) <= eps) && (x[i][0] > 0.0))) {

				Vec3d.copyvect3d(x[i],dirs[j]);
				j ++;
			}
	
	}
}

/**
 * for detection of fibers in stack
 * GRAY8, GRAY16, GRAY32 images.
 *
 * @author Jiri Janacek
 */
public class Vector_Line_3D implements PlugIn {

	// plugin parameters
	public double sigma;
	public static int num_scales;
	
	public void run(String arg) {
		if (IJ.versionLessThan("1.46j"))
			return;
		ImagePlus imp = IJ.getImage();
		int type = imp.getType();
		
		if (! ((type == ImagePlus.GRAY8)||(type == ImagePlus.GRAY16)||
				(type == ImagePlus.GRAY32))) {
			IJ.error("Vector Line 3D", "unsupported image type");
			return;
		}
		if (! ((imp.getNSlices() > 1) && (imp.getNSlices() == imp.getStackSize()))) {
			IJ.error("Vector Line 3D", "only single z-stack is supported");
			return;
		}
		if (! showDialog())
			return;

		imp.startTiming();  
		run(imp);		
		imp.updateAndDraw();
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Vector Line 3D filter");

// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("sigma (units)", 4.00, 2);
		gd.addNumericField("scale number", 2, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

// get entered values
		sigma = gd.getNextNumber();
		num_scales = (int) gd.getNextNumber();

		return true;
	}

	private void run(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		double dx = cal.pixelWidth;
		double dy = cal.pixelHeight;
		double dz = cal.pixelDepth;
		ImageStack res = filter(imp.getStack(), dx, dy, dz, sigma);
// find maximum
		double maxv = 0.;
		int nums = res.getSize();
		for (int i = 1; i <= nums; i++) {
			ImageProcessor iprc = res.getProcessor(i);
			ImageStatistics stat = ImageStatistics.getStatistics(iprc, ImageStatistics.MIN_MAX, null);
			if (maxv < stat.max) maxv = stat.max;
		}
// new image		
		ImagePlus outp = new ImagePlus("Vector Line 3D Filtered", res); 
		outp.setDisplayRange(0., maxv);
		outp.show();
		outp.updateAndDraw();
	}
	
	public static ImageStack filter(ImageStack ims, double dx, double dy, double dz, double sigma) {
		
		class Aac1 {
			public float[] data;
			public int nx, ny, nz;
			Aac1(int nx0, int ny0, int nz0, float[] data0) {
				nx = nx0; ny = ny0; nz = nz0; data = data0;
			}
			public float get(int i, int j, int k) {
				if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz))
					return data[i + nx * (j + ny * k)];
				else return (float) 0.;
			}
			public void set(int i, int j, int k, float val) {
				if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz))
					data[i + nx * (j + ny * k)] = val;
			}
			public void zero() {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++) data[i] = (float) 0.;
			}
			public void maxval() {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++) data[i] = Float.MAX_VALUE;
			}
			public void max(float[] val) {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++) 
					if (val[i] > data[i]) data[i] = val[i];
			}
			public void min(float[] val) {
				int size = nx * ny * nz;
				for (int i = 0; i < size; i ++) 
					if (val[i] < data[i]) data[i] = val[i];
			}
		}

		class Aac {
			public float[] data;
			public int nc, nx, ny, nz, nv;
			Aac(int nc0, int nx0, int ny0, int nz0, int nv0, float[] data0) {
				nc = nc0; nx = nx0; ny = ny0; nz = nz0; nv = nv0; data = data0;
			}
			public float get(int c, int i, int j, int k, int v) {
				if ((c >= 0) && (c < nc) && (i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz))
					return data[c + nc * (i + nx * (j + ny * (k + nz * v)))];
				else return (float) 0.;
			}
			public void set(int c, int i, int j, int k, int v, float val) {
				if ((c >= 0) && (c < nc) && (i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz) && (v >= 0) && (v < nv))
					data[c + nc * (i + nx * (j + ny * (k + nz * v)))] = val;
			}
		}

// get width and height
		int width = ims.getWidth();
		int height = ims.getHeight();
		int depth = ims.getSize();
		int[] dims0 = {width, height, depth};

// result
		ImageStack imr = ImageStack.create(dims0[0], dims0[1], dims0[2], 32);
// algorithmus
		// boxes size - max sigma
		double sx = sigma / dx;
		double sy = sigma / dy;
		double sz = sigma / dz;

		int[] sigmai = {(int) (2. * sx + 0.5), (int) (2. * sy + 0.5), (int) (2. * sz + 0.5)};
		int[] dims = new int[3];
		int[] nums = new int[3];
		double[] ldims = new double[3];
		int size = 256; //max size boxes
		for (int i = 0; i < 3; i ++) 
			if (dims0[i] > size) {
				nums[i] = (int) ((dims0[i] - 2. * sigmai[i]) / (size - 2. * sigmai[i]) + 0.5);
				ldims[i] = (dims0[i] + 2. * (nums[i] - 1) * sigmai[i]) / nums[i];
				dims[i] = (int) (ldims[i] + 0.5);	
			}
			else {
				nums[i] = 1;
				ldims[i] = dims[i] = dims0[i];
			}
		double mins = sx; if (sy < mins) mins = sy; if (sz < mins) mins = sz;
		int n_sc = 1;
		while ((mins > 1.) && (n_sc < num_scales)) {
			n_sc ++;
			mins = mins / 2.;
		}
//directions
		IcosaDirs Dirs = new IcosaDirs(2);
		int num_nodes = Dirs.dirs.length;
		int ndir = 6;
		double sqrt3 = Math.sqrt(3.);
		
		ImageStack imf = ImageStack.create(dims[0], dims[1], dims[2], 32);
		ImagePlus pimf = new ImagePlus("", imf);		
		
		float[] img = new float[dims[0]*dims[1]*dims[2]];
		Aac1 pimg = new Aac1(dims[0], dims[1], dims[2], img);

		float[] imh = new float[3*dims[0]*dims[1]*dims[2]*n_sc];
		Aac pimh = new Aac(3,dims[0], dims[1], dims[2], n_sc, imh);

		float[] imv = new float[dims[0]*dims[1]*dims[2]];
		Aac1 pimv = new Aac1(dims[0], dims[1], dims[2], imv);

		float[] imu = new float[dims[0]*dims[1]*dims[2]];
		Aac1 pimu = new Aac1(dims[0], dims[1], dims[2], imu);

		float[] imt = new float[dims[0]*dims[1]*dims[2]];
		Aac1 pimt = new Aac1(dims[0], dims[1], dims[2], imt);

		//boxec
		int numbox = nums[0] * nums[1] * nums[2];

		for (int ix = 0; ix < nums[0]; ix ++) {
			int x1, x2, x3, x4;
			if (nums[0] == 1) {
				x2 = dims0[0];
				x1 = x2 - dims[0];
				x4 = x2;
				x3 = x1;	
			}
			else if (ix == (nums[0] - 1)) {
				x2 = dims0[0];
				x1 = x2 - dims[0];
				x4 = dims0[0];
				x3 = (int) (x4 - ldims[0] + 0.5 + sigmai[0]);
			}
			else {
				x1 = (int) (ix * (ldims[0] - 2 * sigmai[0]));
				x2 = x1 + dims[0];
				x3 = x1 + ((ix > 0)?sigmai[0]:0);
				x4 = (int) (x2 - sigmai[0] + 0.5);
			}
			for (int iy = 0; iy < nums[1]; iy ++)	{
				int y1, y2, y3, y4;
				if (nums[1] == 1) {
					y2 = dims0[1];
					y1 = y2 - dims[1];
					y4 = y2;
					y3 = y1;				
				}
				else if (iy == (nums[1] - 1)) {
					y2 = dims0[1];
					y1 = y2 - dims[1];
					y4 = dims0[1];
					y3 = (int) (y4 - ldims[1] + 0.5 + sigmai[1]);
				}
				else {
					y1 = (int) (iy * (ldims[1] - 2 * sigmai[1]));
					y2 = y1 + dims[1];
					y3 = y1 + ((iy > 0)?sigmai[1]:0);
					y4 = (int) (y2 - sigmai[1] + 0.5);
				}
				for (int iz = 0; iz < nums[2]; iz ++) {
					int z1, z2, z3, z4;
					if (nums[2] == 1) {
						z2 = dims0[2];
						z1 = z2 - dims[2];
						z4 = z2;
						z3 = z1;					
					}
					else if (iz == (nums[2] - 1)) {
						z2 = dims0[2];
						z1 = z2 - dims[2];
						z4 = dims0[2];
						z3 = (int) (z4 - ldims[2] + 0.5 + sigmai[2]);
					}
					else {
						z1 = (int) (iz * (ldims[2] - 2 * sigmai[2]));
						z2 = z1 + dims[2];
						z3 = z1 + ((iz > 0)?sigmai[2]:0);
						z4 = (int) (z2 - sigmai[2] + 0.5);
					}
										
					int nx = dims[0];
					int ny = dims[1];
					int nz = dims[2];
						
// scales				
					double sigmav = sigma;
					for (int v = 0; v < n_sc; v ++) {
						//	get source voxels
						ims.getVoxels(x1, y1, z1, dims[0], dims[1], dims[2], img);
						sx = sigmav / dx;
						sy = sigmav / dy;
						sz = sigmav / dz;
						//	gaussian
						imf.setVoxels(0, 0, 0, dims[0], dims[1], dims[2], img);
						GaussianBlur3D.blur(pimf, sx, sy, sz);
						imf.getVoxels(0, 0, 0, dims[0], dims[1], dims[2], img);
						//	algoritmus
						for (int k = 0; k < nz; k ++) {
							int k0 = k - 1; if (k0 < 0) k0 = 0;
							int k2 = k + 1; if (k2 == nz) k2 = nz - 1;
							for (int j = 0; j < ny; j ++) {
								int j0 = j - 1; if (j0 < 0) j0 = 0;
								int j2 = j + 1; if (j2 == ny) j2 = ny - 1;
								for (int i = 0; i < nx; i ++) {
									int i0 = i - 1; if (i0 < 0) i0 = 0;
									int i2 = i + 1; if (i2 == nx) i2 = nx - 1;
									double x000 = pimg.get(i0, j0, k0);
									double x001 = pimg.get(i, j0, k0);
									double x002 = pimg.get(i2, j0, k0);
									double x010 = pimg.get(i0, j, k0);
									double x011 = pimg.get(i, j, k0);
									double x012 = pimg.get(i2, j, k0);
									double x020 = pimg.get(i0, j2, k0);
									double x021 = pimg.get(i, j2, k0);
									double x022 = pimg.get(i2, j2, k0);
									double x100 = pimg.get(i0, j0, k);
									double x101 = pimg.get(i, j0, k);
									double x102 = pimg.get(i2, j0, k);
									double x110 = pimg.get(i0, j, k);
	//								double x111 = pimg.get(i, j, k);
									double x112 = pimg.get(i2, j, k);
									double x120 = pimg.get(i0, j2, k);
									double x121 = pimg.get(i, j2, k);
									double x122 = pimg.get(i2, j2, k);
									double x200 = pimg.get(i0, j0, k2);
									double x201 = pimg.get(i, j0, k2);
									double x202 = pimg.get(i2, j0, k2);
									double x210 = pimg.get(i0, j, k2);
									double x211 = pimg.get(i, j, k2);
									double x212 = pimg.get(i2, j, k2);
									double x220 = pimg.get(i0, j2, k2);
									double x221 = pimg.get(i, j2, k2);
									double x222 = pimg.get(i2, j2, k2);
									
									double ddx = (x002 + x012 + x022 + x102 + x112 + x122 + x202 + x212 + x222
											 - x000 - x010 - x020 - x100 - x110 - x120 - x200 - x210 - x220) / (18. * dx);
									double ddy = (x020 + x021 + x022 + x120 + x121 + x122 + x220 + x221 + x222
											 - x000 - x001 - x002 - x100 - x101 - x102 - x200 - x201 - x202) / (18. * dy);
									double ddz = (x200 + x201 + x202 + x210 + x211 + x212 + x220 + x221 + x222
											 - x000 - x001 - x002 - x010 - x011 - x012 - x020 - x021 - x022) / (18. * dz);
										
									pimh.set(0, i, j, k, v, (float) ddx);
									pimh.set(1, i, j, k, v, (float) ddy);
									pimh.set(2, i, j, k, v, (float) ddz);
									
								}
							}
						}
						sigmav = sigmav/2.;
					}
// end scales				
					pimv.zero();
					for(int n = 0; n < num_nodes; n ++) {
/*						double[] d1 = new double[3];
						double[] d2 = new double[3];
						double pangle;

						d1[0] = Dirs.dirs[n][0];
						pangle = 0.;
						d1[1] = Math.cos(pangle) * Dirs.dirs[n][1] - Math.sin(pangle) * Dirs.dirs[n][2];
						d1[2] = Math.sin(pangle) * Dirs.dirs[n][1] + Math.cos(pangle) * Dirs.dirs[n][2];
						pangle = 0.;
						d2[0] = Math.cos(pangle) * d1[0] - Math.sin(pangle) * d1[1];
						d2[1] = Math.sin(pangle) * d1[0] + Math.cos(pangle) * d1[1];
						d2[2] = d1[2];
*/
						double[] d2 = {Dirs.dirs[n][0], Dirs.dirs[n][1], Dirs.dirs[n][2]};
						
						double dangle =  Math.PI / ndir;
						double rangle = Math.acos(d2[2]);
						double[] axe = {-d2[1], d2[0], 0};
						double[][] mat = new double[3][3];
						Vec3d.axis_to_matrix(axe, rangle, mat); 

						double angle = 0;
						pimu.maxval();
						for (int m = 0; m < ndir; m ++, angle += dangle) {
							double[] e1 = {Math.cos(angle),Math.sin(angle), 0.0};
							Vec3d.mult_mat_vec(mat,e1);
							for (int k = 0; k < nz; k ++) {
								for (int j = 0; j < ny; j ++) {
									for (int i = 0; i < nx; i ++) {
										double max = 0.;
// scales										
										sigmav = sigma;
										for (int v = 0; v < n_sc; v ++) {
											double min = Float.MAX_VALUE; 
											int dx1 = Math.round((float)(e1[0] * sigmav * sqrt3 / dx)); 
											int dy1 = Math.round((float)(e1[1] * sigmav * sqrt3 / dy)); 
											int dz1 = Math.round((float)(e1[2] * sigmav * sqrt3 / dz)); 
											int ii, jj, kk;
											ii = i + dx1;
											jj = j + dy1;
											kk = k + dz1;
											if ((ii >= 0) && (ii < nx) && (jj >= 0) && (jj < ny) && (kk >= 0) && (kk < nz)) {
												double val= pimh.get(0, ii, jj, kk, v) * e1[0]
														+ pimh.get(1, ii, jj, kk, v) * e1[1]
														+ pimh.get(0, ii, jj, kk, v) * e1[2];
												val = - val * sigmav;
												if (val < min) min = val;
											}
											
											ii = i - dx1;
											jj = j - dy1;
											kk = k - dz1;
											if ((ii >= 0) && (ii < nx) && (jj >= 0) && (jj < ny) && (kk >= 0) && (kk < nz)) {
												double val= pimh.get(0, ii, jj, kk, v) * e1[0]
														+ pimh.get(1, ii, jj, kk, v) * e1[1]
														+ pimh.get(0, ii, jj, kk, v) * e1[2];
												val = val * sigmav;
												if (val < min) min = val;											
											}
											if (min == Float.MAX_VALUE) min = 0.; 
											if (min > max) max = min;
											sigmav = sigmav / 2.;
										}
// end scales										
										pimt.set(i, j, k, (float) (max));
									}
								}
							}

							pimu.min(imt);
							angle += dangle;

						}

						pimv.max(imu);
						IJ.showProgress(n, num_nodes - 1);					
					}
					imf.setVoxels(0, 0, 0, dims[0], dims[1], dims[2], imv);
					//	set destination voxels
					float[] img1 = new float[(x4 - x3) * (y4 - y3) * (z4 - z3)]; //length must be exact

					imf.getVoxels(x3 - x1, y3 - y1, z3 - z1, x4 - x3, y4 - y3, z4 - z3, img1); //by getvoxels
					imr.setVoxels(x3, y3, z3, x4 - x3, y4 - y3, z4 - z3, img1);
					// counter
					if (numbox > 1) IJ.showStatus("Done "+String.format("%.2g",
							(iz + nums[2] * (iy + ix * nums[1])) / (double)(numbox - 1)));
				}
			}
		}

		return imr;
}
	
	public void showAbout() {
		IJ.showMessage("Vector Line 3D",
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
		Class<?> clazz = Vector_Line_3D.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the capillaries sample
		ImagePlus image = IJ.openImage("https://imagej.net/_images/e/e9/Capillaries_heart.zip");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
