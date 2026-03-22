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

import org.applied_geodesy.adjustment.bundle.PartialDerivativeFactory.CollinearityEquationFactory;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

class DistortionModelFactory {
	private DistortionModelFactory() {}

	static void apply(CollinearityEquationFactory collinearityEquation, Matrix A, Vector w, double deltaX, double deltaY, double par_deltaX_xs, double par_deltaX_ys, double par_deltaY_xs, double par_deltaY_ys) {
		w.add(0, -deltaX); // xobs - (x0 + xs + deltaX)
		w.add(1, -deltaY); // yobs - (y0 + ys + deltaY)
		
		/** Object point coordinates **/
		int column = collinearityEquation.objectCoordinate.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_X + par_deltaX_ys * collinearityEquation.par_ys_X);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_X + par_deltaY_ys * collinearityEquation.par_ys_X);
		}

		column = collinearityEquation.objectCoordinate.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_Y + par_deltaX_ys * collinearityEquation.par_ys_Y);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_Y + par_deltaY_ys * collinearityEquation.par_ys_Y);
		}

		column = collinearityEquation.objectCoordinate.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_Z + par_deltaX_ys * collinearityEquation.par_ys_Z);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_Z + par_deltaY_ys * collinearityEquation.par_ys_Z);
		}


		/**  interior orientation of the camera **/
		column = collinearityEquation.interiorOrientation.getPrincipleDistance().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_c + par_deltaX_ys * collinearityEquation.par_ys_c);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_c + par_deltaY_ys * collinearityEquation.par_ys_c);
		}
		
		
		/**  exterior orientation of the image **/
		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_X0 + par_deltaX_ys * collinearityEquation.par_ys_X0);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_X0 + par_deltaY_ys * collinearityEquation.par_ys_X0);
		}

		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_Y0 + par_deltaX_ys * collinearityEquation.par_ys_Y0);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_Y0 + par_deltaY_ys * collinearityEquation.par_ys_Y0);
		}

		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_Z0 + par_deltaX_ys * collinearityEquation.par_ys_Z0);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_Z0 + par_deltaY_ys * collinearityEquation.par_ys_Z0);
		}

		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_omega + par_deltaX_ys * collinearityEquation.par_ys_omega);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_omega + par_deltaY_ys * collinearityEquation.par_ys_omega);
		}

		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_phi + par_deltaX_ys * collinearityEquation.par_ys_phi);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_phi + par_deltaY_ys * collinearityEquation.par_ys_phi);
		}

		column = collinearityEquation.exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			A.add(0, column, par_deltaX_xs * collinearityEquation.par_xs_kappa + par_deltaX_ys * collinearityEquation.par_ys_kappa);
			A.add(1, column, par_deltaY_xs * collinearityEquation.par_xs_kappa + par_deltaY_ys * collinearityEquation.par_ys_kappa);
		}
	}
}
