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

package org.applied_geodesy.adjustment.bundle.derivation;

import java.util.Set;

import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.ZernikeDistortionModel;
import org.applied_geodesy.adjustment.bundle.derivation.PartialDerivativeFactory.CollinearityEquationFactory;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ZernikeCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.ZernikeCoefficient.ZernikePolynomial;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

final class ZernikeDistortionModelFactory {
	private ZernikeDistortionModelFactory() {}
	
	static void apply(ZernikeDistortionModel distortionModel, CollinearityEquationFactory collinearityEquation, Set<Integer> columns, Matrix A, Vector w) {
		double r0 = distortionModel.getR0();
		double r02 = r0 * r0;
		
		double xxs = collinearityEquation.xs * collinearityEquation.xs;
		double yys = collinearityEquation.ys * collinearityEquation.ys;
		double xys = collinearityEquation.xs * collinearityEquation.ys;

		double phi = Math.atan2(collinearityEquation.ys, collinearityEquation.xs);
		double r2 = xxs + yys;
		double rn2 = r2/r02;
		double const2rnr0 = 2.0/rn2/r02;
		
		for (UnknownParameter<? extends DistortionModel> unknownParameter : distortionModel) {
			if (unknownParameter.getParameterType() != ParameterType.ZERNIKE_POLYNOMIAL_Z)
				continue;
			
			ZernikeCoefficient Zi = (ZernikeCoefficient)unknownParameter;
			ZernikePolynomial zernikePoly = Zi.getZernikePolynomial();
			
			double zi = Zi.getValue();
			double m = zernikePoly.getAzimuthalFrequency();
			
			double sinmphi = Math.sin(m*phi);
			double cosmphi = Math.cos(m*phi);
			
			double par_xs_Zi = 0;
			double par_ys_Zi = 0;
			
			int numOfRadialTerms = zernikePoly.getNumberOfRadialTerms();
			
			for (int j = 0; j < numOfRadialTerms; j++) {
				long pj = zernikePoly.getRadialExponent(j);
				long constExp = (pj/2 - 1);
				
				double cj = zernikePoly.getRadialCoefficient(j);
				double constC = cj/r02 * Math.pow(rn2, constExp);
				
				double constXsin = (-pj * collinearityEquation.xs * sinmphi + m * collinearityEquation.ys * cosmphi);
				double constYsin = (-pj * collinearityEquation.ys * sinmphi - m * collinearityEquation.xs * cosmphi);

				double constXcos = ( pj * collinearityEquation.xs * cosmphi + m * collinearityEquation.ys * sinmphi);
				double constYcos = ( pj * collinearityEquation.ys * cosmphi - m * collinearityEquation.xs * sinmphi);

				if (m < 0) {
					// distortion coefficients in x and y
					double deltaX_sin = zi * constC * constXsin;
					double deltaY_sin = zi * constC * constYsin;

					// chain rule coefficients
					double par_deltaXsin_xs = zi * constC * (constExp * collinearityEquation.xs*const2rnr0 * constXsin - pj*sinmphi + m/r2*(pj*xys*cosmphi + m*yys*sinmphi));
					double par_deltaXsin_ys = zi * constC * (constExp * collinearityEquation.ys*const2rnr0 * constXsin +  m*cosmphi - m/r2*(pj*xxs*cosmphi + m*xys*sinmphi));

					double par_deltaYsin_xs = zi * constC * (constExp * collinearityEquation.xs*const2rnr0 * constYsin -  m*cosmphi + m/r2*(pj*yys*cosmphi - m*xys*sinmphi));
					double par_deltaYsin_ys = zi * constC * (constExp * collinearityEquation.ys*const2rnr0 * constYsin - pj*sinmphi - m/r2*(pj*xys*cosmphi - m*xxs*sinmphi));

					// apply distortion model to object point, exterior parameters and principle distance
					DistortionModelFactory.apply(collinearityEquation, A, w, deltaX_sin, deltaY_sin, par_deltaXsin_xs, par_deltaXsin_ys, par_deltaYsin_xs, par_deltaYsin_ys);

					par_xs_Zi += constC * constXsin;
					par_ys_Zi += constC * constYsin;
				}
				else {
					// distortion coefficients in x and y
					double deltaX_cos = zi * constC * constXcos;
					double deltaY_cos = zi * constC * constYcos;

					// chain rule coefficients
					double par_deltaXcos_xs = zi * constC * (constExp * collinearityEquation.xs*const2rnr0 * constXcos + pj*cosmphi + m/r2*(pj*xys*sinmphi - m*yys*cosmphi));
					double par_deltaXcos_ys = zi * constC * (constExp * collinearityEquation.ys*const2rnr0 * constXcos +  m*sinmphi - m/r2*(pj*xxs*sinmphi - m*xys*cosmphi));

					double par_deltaYcos_xs = zi * constC * (constExp * collinearityEquation.xs*const2rnr0 * constYcos -  m*sinmphi + m/r2*(pj*yys*sinmphi + m*xys*cosmphi));
					double par_deltaYcos_ys = zi * constC * (constExp * collinearityEquation.ys*const2rnr0 * constYcos + pj*cosmphi - m/r2*(pj*xys*sinmphi + m*xxs*cosmphi));

					// apply distortion model to object point, exterior parameters and principle distance
					DistortionModelFactory.apply(collinearityEquation, A, w, deltaX_cos, deltaY_cos, par_deltaXcos_xs, par_deltaXcos_ys, par_deltaYcos_xs, par_deltaYcos_ys);
					
					par_xs_Zi += constC * constXcos;
					par_ys_Zi += constC * constYcos;
				}
				
				//  dZernX: par_deltaX[sin|cos]_xs * par_xs_<PARAM>  +  par_deltaX[sin|cos]_ys * par_ys_<PARAM>
				//  dZernY: par_deltaX[sin|cos]_xs * par_xs_<PARAM>  +  par_deltaX[sin|cos]_ys * par_ys_<PARAM>
			}
			
			/** disortion model **/
			int column = Zi.getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				columns.add(column);
				A.set(0, column, par_xs_Zi);
				A.set(1, column, par_ys_Zi);
			}
		}
	}
}
