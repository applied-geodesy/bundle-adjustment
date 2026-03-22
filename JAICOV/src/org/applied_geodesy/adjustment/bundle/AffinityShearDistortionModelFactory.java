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
import org.applied_geodesy.adjustment.bundle.camera.distortion.AffinityShearDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

final class AffinityShearDistortionModelFactory {
	private AffinityShearDistortionModelFactory() {}

	static void apply(AffinityShearDistortionModel distortionModel, CollinearityEquationFactory collinearityEquation, Set<Integer> columns, Matrix A, Vector w) {
		UnknownParameter<? extends DistortionModel> Cx = distortionModel.getCx();
		UnknownParameter<? extends DistortionModel> Cy = distortionModel.getCy();
		
		double cx = Cx.getValue();
		double cy = Cy.getValue();
		
		// distortion coefficients in x and y
		double deltaX = cx*collinearityEquation.xs + cy*collinearityEquation.ys;
		double deltaY = 0.0;
		
		// chain rule coefficients
		double par_deltaX_xs = cx;
		double par_deltaX_ys = cy;

		double par_deltaY_xs = 0.0;
		double par_deltaY_ys = 0.0;
		
		double par_xs_C1 = collinearityEquation.xs;
		double par_xs_C2 = collinearityEquation.ys;
		
		double par_ys_C1 = 0.0;
		double par_ys_C2 = 0.0;
		
		// apply distortion model to object point, exterior parameters and principle distance
		DistortionModelFactory.apply(collinearityEquation, A, w, deltaX, deltaY, par_deltaX_xs, par_deltaX_ys, par_deltaY_xs, par_deltaY_ys);
		
		//  dAffX: par_dAffX_xs * par_xs_<PARAM>  +  par_dAffX_ys * par_ys_<PARAM>
		//  dAffY: par_dAffY_xs * par_xs_<PARAM>  +  par_dAffY_ys * par_ys_<PARAM>
		
		/** disortion model **/
		int column = Cx.getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_C1);
			A.set(1, column, par_ys_C1);
		}

		column = Cy.getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_C2);
			A.set(1, column, par_ys_C2);
		}
	}
}
