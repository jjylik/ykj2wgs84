public class YkjToWgs84Converter {

	public double[] convert(double E, double N) {

		double f = 1.0D / 297.0;
		double el = f / (2 - f);
		double a = 6378388.0;
		double A1 = (a / (1 + el)) * (1 + (Math.pow(el, 2) / 4) + (Math.pow(el, 4) / 64));
		double E0 = 3500000.0;
		double lambda_0 = 27 * (Math.PI / 180);
		double k_0 = 1.0D;

		double eps = N / (A1 * k_0);
		double n = (E - E0) / (A1 * k_0);

		double h1 = ((((1.0D / 2.0) * el) - ((2.0 / 3.0) * Math.pow(el, 2))) + ((37.0 / 96.0) * Math.pow(el, 3))) - ((1.0D / 360.0) * Math.pow(el, 4));
		double h2 = (((1.0D / 48.0) * Math.pow(el, 2)) + ((1.0D / 15.0) * Math.pow(el, 3))) - ((437.0 / 1140.0) * Math.pow(el, 4));
		double h3 = ((17.0 / 480.0) * Math.pow(el, 3)) - ((37.0 / 840.0) * Math.pow(el, 4));
		double h4 = (4397.0 / 161280.0) * Math.pow(el, 4);

		double r1 = h1 * Math.sin(2 * eps) * Math.cosh(2 * n);
		double r2 = h2 * Math.sin(4 * eps) * Math.cosh(4 * n);
		double r3 = h3 * Math.sin(6 * eps) * Math.cosh(6 * n);
		double r4 = h4 * Math.sin(8 * eps) * Math.cosh(8 * n);

		double s1 = h1 * Math.sinh(2 * n) * Math.cos(2 * eps);
		double s2 = h2 * Math.sinh(4 * n) * Math.cos(4 * eps);
		double s3 = h3 * Math.sinh(6 * n) * Math.cos(6 * eps);
		double s4 = h4 * Math.sinh(8 * n) * Math.cos(8 * eps);

		double eps_f = eps - (r1 + r2 + r3 + r4);
		double n_f = n - (s1 + s2 + s3 + s4);
		double e = Math.sqrt((2 * f) - Math.pow(f, 2));

		double beta = Math.asin(sech(n_f) * Math.sin(eps_f));
		double l = Math.asin(Math.tanh(n_f) / (Math.cos(beta)));

		double Q1 = asinh(Math.tan(beta));
		double Q = Q1 + (e * atanh(Math.tanh(Q1) * e));
		Q = Q1 + (e * atanh(Math.tanh(Q) * e));
		Q = Q1 + (e * atanh(Math.tanh(Q) * e));
		Q = Q1 + (e * atanh(Math.tanh(Q) * e));

		double phi = Math.atan(Math.sinh(Q));
		double lambda_r = lambda_0 + l;
		double lambda = lambda_r;

		double NN = a * Math.pow((1 - (Math.pow(e, 2) * Math.pow(Math.sin(phi), 2))), -1.0D / 2.0);
		double X0 = (NN + 50) * Math.cos(phi) * Math.cos(lambda);
		double Y0 = (NN + 50) * Math.cos(phi) * Math.sin(lambda);
		double Z0 = ((NN * (1 - Math.pow(e, 2))) + 50) * Math.sin(phi);

		double DX = -96.062;
		double DY = -82.428;
		double DZ = -121.754;
		double Rx = -4.801 * (1.0D / 3600.0) * (Math.PI / 180.0);
		double Ry = -0.345 * (1.0D / 3600.0) * (Math.PI / 180.0);
		double Rz = 1.376 * (1.0D / 3600.0) * (Math.PI / 180.0);
		double m = 1.496;

		double X1 = DX + ((1 + (m / 1000000)) * (((1 * X0) + (Rz * Y0)) - (Ry * Z0)));
		double Y1 = DY + ((1 + (m / 1000000)) * ((-Rz * X0) + (1 * Y0) + (Rx * Z0)));
		double Z1 = DZ + ((1 + (m / 1000000)) * (((Ry * X0) - (Rx * Y0)) + (1 * Z0)));

		lambda = Math.atan(Y1 / X1);
		double a_grs = 6378137.0;
		double f_grs = 1 / 298.2572;
		double e_grs = Math.pow(((2 * f_grs) - Math.pow(f_grs, 2)), 0.5);
		double phi_0 = Math.atan(Z1 / ((1 - Math.pow(e_grs, 2)) * Math.pow(
		(Math.pow(X1, 2) + Math.pow(Y1, 2)), (0.5))));
		double phi_i = phi_0;
		double N_i = 0.0;
		double h_i = 0.0;
		for (int i = 0; i < 3; i++) {
			N_i = a_grs * Math.pow((1 - (Math.pow(e_grs, 2) * Math.pow(
			Math.sin(phi_i), 2))), (-0.5));
			if (Math.abs(phi_0) < (45.0 * (Math.PI / 180.0))) {
				h_i = ((Math.sqrt((Math.pow(X1, 2) + Math.pow(Y1, 2)))) / (Math.cos(phi_i))) - N_i;
			} else {
				h_i = (Z1 / (Math.sin(phi_i))) - ((1.0 - Math.pow(e_grs, 2)) * N_i);
			}
			phi_i = Math.atan((Z1 / Math.sqrt(Math.pow(X1, 2) + Math.pow(Y1, 2))) * (1.0D / (1.0D - ((Math.pow(e_grs, 2) * N_i) / (N_i + h_i)))));
		}

		return new double[] {
			phi_i * (180.0 / Math.PI),
			lambda * (180.0 / Math.PI)
		};
	}

}
