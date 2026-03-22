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

package org.applied_geodesy.adjustment.bundle;

import java.util.Set;

import org.applied_geodesy.adjustment.bundle.PartialDerivativeFactory.CollinearityEquationFactory;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadialDistanceDistortionModel;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

final class RadialDistanceDistortionModelFactory {
	private RadialDistanceDistortionModelFactory() {}

	static void apply(RadialDistanceDistortionModel distortionModel, CollinearityEquationFactory collinearityEquation, Set<Integer> columns, Matrix A, Vector w) {
		double r0 = distortionModel.getR0();
		
		double r2 = collinearityEquation.xs*collinearityEquation.xs + collinearityEquation.ys*collinearityEquation.ys;
		double r02 = r0 * r0;
		
		double xxs2 = 2.0 * collinearityEquation.xs * collinearityEquation.xs;
		double yys2 = 2.0 * collinearityEquation.ys * collinearityEquation.ys;
		double xys2 = 2.0 * collinearityEquation.xs * collinearityEquation.ys;
		
		for (UnknownParameter<? extends DistortionModel> unknownParameter : distortionModel) {
			if (unknownParameter.getParameterType() != ParameterType.DISTANCE_POLYNOMIAL_D)
				continue;
			
			PolynomialCoefficient<?> Di = (PolynomialCoefficient<?>)unknownParameter;
			
			double di = Di.getValue();
			int expi = Di.getOrder();

			double dRi = Math.pow(r2, expi) - Math.pow(r02, expi);
			double dDisti = (di * dRi) / collinearityEquation.N;

			// distortion coefficients in x and y
			double deltaX = collinearityEquation.xs * dDisti;
			double deltaY = collinearityEquation.ys * dDisti;


			// chain rule coefficients
			double par_N_X = collinearityEquation.r13;
			double par_N_Y = collinearityEquation.r23;
			double par_N_Z = collinearityEquation.r33;

			double par_N_X0 = -collinearityEquation.r13;
			double par_N_Y0 = -collinearityEquation.r23;
			double par_N_Z0 = -collinearityEquation.r33;
					
			double par_N_omega = -collinearityEquation.r33 * collinearityEquation.dY + collinearityEquation.r23 * collinearityEquation.dZ;              // -cosOmega*cosPhi * dY - sinOmega*cosPhi * dZ;
			double par_N_phi   =  collinearityEquation.kx * collinearityEquation.cosKappa - collinearityEquation.ky * collinearityEquation.sinKappa;    //  cosPhi*dX + sinOmega*sinPhi*dY - cosOmega*sinPhi*dZ
			double par_N_kappa =  0.0;
			
			double constRadi = (di * expi * Math.pow(r2, expi - 1)) / collinearityEquation.N;
			double par_deltaX_xs = xxs2 * constRadi + dDisti;
			double par_deltaX_ys = xys2 * constRadi;
			double par_dDistX_N  = -deltaX / collinearityEquation.N;

			double par_deltaY_xs = xys2 * constRadi;
			double par_deltaY_ys = yys2 * constRadi + dDisti;
			double par_dDistY_N  = -deltaY / collinearityEquation.N;

			double par_xs_Di = (collinearityEquation.xs * dRi) / collinearityEquation.N;
			double par_ys_Di = (collinearityEquation.ys * dRi) / collinearityEquation.N;
			
			// apply distortion model to object point, exterior parameters and principle distance
			DistortionModelFactory.apply(collinearityEquation, A, w, deltaX, deltaY, par_deltaX_xs, par_deltaX_ys, par_deltaY_xs, par_deltaY_ys);

			//  dDistX: par_dDistX_xs * par_xs_<PARAM>  +  par_dDistX_ys * par_ys_<PARAM>  +  par_dDistX_N * par_N_<PARAM>;
			//  dDistY: par_dDistY_xs * par_xs_<PARAM>  +  par_dDistY_ys * par_ys_<PARAM>  +  par_dDistY_N * par_N_<PARAM>;

			/** disortion model **/
			int column = Di.getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				columns.add(column);
				A.set(0, column, par_xs_Di);
				A.set(1, column, par_ys_Di);
			}
			
			// Object point coordinates
			column = collinearityEquation.objectCoordinate.getX().getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_X * par_dDistX_N);
				A.add(1, column, par_N_X * par_dDistY_N);
			}

			column = collinearityEquation.objectCoordinate.getY().getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_Y * par_dDistX_N);
				A.add(1, column, par_N_Y * par_dDistY_N);
			}

			column = collinearityEquation.objectCoordinate.getZ().getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_Z * par_dDistX_N);
				A.add(1, column, par_N_Z * par_dDistY_N);
			}

			// exterior orientation of the image
			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_X0 * par_dDistX_N);
				A.add(1, column, par_N_X0 * par_dDistY_N);
			}

			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_Y0 * par_dDistX_N);
				A.add(1, column, par_N_Y0 * par_dDistY_N);
			}

			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_Z0 * par_dDistX_N);
				A.add(1, column, par_N_Z0 * par_dDistY_N);
			}

			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_omega * par_dDistX_N);
				A.add(1, column, par_N_omega * par_dDistY_N);
			}

			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_phi * par_dDistX_N);
				A.add(1, column, par_N_phi * par_dDistY_N);
			}

			column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				A.add(0, column, par_N_kappa * par_dDistX_N);
				A.add(1, column, par_N_kappa * par_dDistY_N);
			}
		}
	}
}
