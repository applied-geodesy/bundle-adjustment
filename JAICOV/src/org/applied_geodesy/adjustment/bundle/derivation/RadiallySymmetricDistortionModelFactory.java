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
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadiallySymmetricDistortionModel;
import org.applied_geodesy.adjustment.bundle.derivation.PartialDerivativeFactory.CollinearityEquationFactory;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

final class RadiallySymmetricDistortionModelFactory {
	private RadiallySymmetricDistortionModelFactory() {}

	static void apply(RadiallySymmetricDistortionModel distortionModel, CollinearityEquationFactory collinearityEquation, Set<Integer> columns, Matrix A, Vector w) {
		double r0 = distortionModel.getR0();
		
		double r2 = collinearityEquation.xs*collinearityEquation.xs + collinearityEquation.ys*collinearityEquation.ys;
		double r02 = r0 * r0;
		
		double xxs2 = 2.0 * collinearityEquation.xs * collinearityEquation.xs;
		double yys2 = 2.0 * collinearityEquation.ys * collinearityEquation.ys;
		double xys2 = 2.0 * collinearityEquation.xs * collinearityEquation.ys;
		
		for (UnknownParameter<? extends DistortionModel> unknownParameter : distortionModel) {
			if (unknownParameter.getParameterType() != ParameterType.RADIAL_POLYNOMIAL_A)
				continue;
			
			PolynomialCoefficient<?> Ai = (PolynomialCoefficient<?>)unknownParameter;
			
			double ai = Ai.getValue();
			int expi = Ai.getOrder();
			
			double dRi = Math.pow(r2, expi) - Math.pow(r02, expi);
			double dRadi = ai * dRi;
			
			// distortion coefficients in x and y
			double deltaX = collinearityEquation.xs * dRadi;
			double deltaY = collinearityEquation.ys * dRadi;
			
			// chain rule coefficients
			double constRadi = ai * expi * Math.pow(r2, expi - 1);
			double par_deltaX_xs = xxs2 * constRadi + dRadi;
			double par_deltaX_ys = xys2 * constRadi;
			
			double par_deltaY_xs = xys2 * constRadi;
			double par_deltaY_ys = yys2 * constRadi + dRadi;
			
			double par_xs_Ai = collinearityEquation.xs * dRi;
			double par_ys_Ai = collinearityEquation.ys * dRi;
			
			// apply distortion model to object point, exterior parameters and principle distance
			DistortionModelFactory.apply(collinearityEquation, A, w, deltaX, deltaY, par_deltaX_xs, par_deltaX_ys, par_deltaY_xs, par_deltaY_ys);
			
			//  dRadX: par_dRadX_xs * par_xs_<PARAM>  +  par_dRadX_ys * par_ys_<PARAM>
			//  dRadY: par_dRadY_xs * par_xs_<PARAM>  +  par_dRadY_ys * par_ys_<PARAM>
			
			/** disortion model **/
			int column = Ai.getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				columns.add(column);
				A.set(0, column, par_xs_Ai);
				A.set(1, column, par_ys_Ai);
			}
		}
	}
}
