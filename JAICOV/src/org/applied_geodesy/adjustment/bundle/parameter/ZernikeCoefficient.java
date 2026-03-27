/***********************************************************************
* Copyright by Michael Loesler, https://software.applied-geodesy.org   *
*                                                                      *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 3 of the License, or    *
* at your option any later version.                                    *
*                                                                      *
* This program is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with this program; if not, see <http://www.gnu.org/licenses/>  *
* or write to the                                                      *
* Free Software Foundation, Inc.,                                      *
* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.            *
*                                                                      *
***********************************************************************/


package org.applied_geodesy.adjustment.bundle.parameter;

import org.applied_geodesy.adjustment.MathExtension;
import org.applied_geodesy.adjustment.bundle.camera.distortion.ZernikeDistortionModel;

public class ZernikeCoefficient extends PolynomialCoefficient<ZernikeDistortionModel> {
	
	public class ZernikePolynomial {
		private final int n, m;
		private final long p[];
		private final long c[];
		private final double length;

		private ZernikePolynomial(int order) {
			this(order, Boolean.FALSE);
		}

		private ZernikePolynomial(int order, boolean normalise) {
			// Schwiegerling (2014, p. 81)
			this.n = (int)Math.ceil((-3 + Math.sqrt(9 + 8*order)) / 2); // Eq. 2:100
			this.m = 2 * order - this.n * (this.n + 2);                 // Eq. 2:101

			int halfnm = (this.n - Math.abs(this.m)) / 2;

			this.p = new long[halfnm + 1];
			this.c = new long[halfnm + 1];

			for (int k = 0; k <= halfnm; k++) {
				this.p[k] = this.n - 2 * k;
				this.c[k] = ((k % 2 == 0) ? 1 : -1) * MathExtension.binomial(this.n - k, k) * MathExtension.binomial(this.n - 2 * k, halfnm - k);
			}

			this.length = normalise ? Math.sqrt((1 + ((this.m != 0) ? 1 : 0)) * (this.n + 1)/Math.PI) : 1.0;
		}

		/**
		 * Returns radial order <code>n</code>
		 * @return n
		 */
		public final int getRadialOrder() {
			return this.n;
		}

		/**
		 * Returns azimuthal frequency <code>m</code>
		 * @return m
		 */
		public final int getAzimuthalFrequency() {
			return this.m;
		}

		/**
		 * Returns exponent <code>p</code> of idx-th radial polynomial
		 * 
		 *    m               p
		 *  R(r) = len * c * r
		 *    n
		 * 
		 * @return p
		 */
		public long getRadialExponent(int idx) {
			return this.p[idx];
		}

		/**
		 * Returns coefficient <code>cn = len * c</code> of idx-th radial polynomial
		 * 
		 *    m               p
		 *  R(r) = len * c * r
		 *    n
		 * 
		 * @return cn
		 */
		public double getRadialCoefficient(int idx) {
			return this.length * this.c[idx];
		}

		/**
		 * Returns the number of radial polynomials
		 * @return num
		 */
		public int getNumberOfRadialTerms() {
			return this.c.length;
		}

		/**
		 * Returns the value of the radial polynomial defined by
		 * 
		 *    m               p
		 *  R(r) = len * c * r
		 *    n
		 *    
		 *  where r is the distance  
		 * 
		 * @param r
		 * @return R
		 */
		public double getRadial(double r) {
			double R = 0;
			for (int i = 0; i < this.c.length; i++)
				R += this.c[i] * Math.pow(r, this.p[i]);
			return this.length * R;
		}

		/**
		 * Returns the value of the azimuthal polynomial defined by
		 *
		 *     m     cos( |m| * phi )   if m >= 0
		 *  G(phi) = sin( |m| * phi )    else
		 * 
		 * where phi is the azimuth
		 * 
		 * @param phi
		 * @return G
		 */
		public double getAzimuthal(double phi) {
			if (this.m < 0)
				return Math.sin(-this.m * phi);

			return Math.cos(this.m * phi);
		}

		/**
		 * Returns the value of the Zernike polynomial Z(r, phi) defined by
		 * 
		 *   +/-m        m       m
		 * Z(r, phi) = R(r) * G(phi)
		 *     n         n 
		 * 
		 * @param r
		 * @param phi
		 * @return Z
		 */
		public double getValue(double r, double phi) {
			double R = this.getRadial(r);
			double G = this.getAzimuthal(phi);
			return R * G;
		}

		@Override
		public String toString() {
			return "ZernikePolynomial [n=" + n + ", m=" + m + "]";
		}
	}

	private final ZernikePolynomial zernikePolynomial; 
	
	public ZernikeCoefficient(ParameterType parameterType, ZernikeDistortionModel reference, int order) {
		super(parameterType, reference, order);
		
		this.zernikePolynomial = new ZernikePolynomial(order);
	}
	
	public final ZernikePolynomial getZernikePolynomial() {
		return this.zernikePolynomial;
	}

}
