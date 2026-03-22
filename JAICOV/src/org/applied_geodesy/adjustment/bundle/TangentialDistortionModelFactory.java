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
import org.applied_geodesy.adjustment.bundle.camera.distortion.TangentialDistortionModel;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

final class TangentialDistortionModelFactory {
	private TangentialDistortionModelFactory() {}

	static void apply(TangentialDistortionModel distortionModel, CollinearityEquationFactory collinearityEquation, Set<Integer> columns, Matrix A, Vector w) {
		UnknownParameter<? extends DistortionModel> Bx = distortionModel.getBx();
		UnknownParameter<? extends DistortionModel> By = distortionModel.getBy();
		
		double bx = Bx.getValue();
		double by = By.getValue();
		
		double r2 = collinearityEquation.xs*collinearityEquation.xs + collinearityEquation.ys*collinearityEquation.ys;
		
		double xxs2 = 2.0 * collinearityEquation.xs * collinearityEquation.xs;
		double yys2 = 2.0 * collinearityEquation.ys * collinearityEquation.ys;
		double xys2 = 2.0 * collinearityEquation.xs * collinearityEquation.ys;
		
		// distortion coefficients in x and y
		double sum = 1.0;
		
		double deltaX = bx * (r2 + xxs2) + by * xys2;
		double deltaY = by * (r2 + yys2) + bx * xys2;
		
		double par_deltaX_xs = 2.0 * (3.0*bx*collinearityEquation.xs + by*collinearityEquation.ys); 
		double par_deltaX_ys = 2.0 * (    by*collinearityEquation.xs + bx*collinearityEquation.ys);
		
		double par_deltaY_xs = 2.0 * (by*collinearityEquation.xs +     bx*collinearityEquation.ys);
		double par_deltaY_ys = 2.0 * (bx*collinearityEquation.xs + 3.0*by*collinearityEquation.ys);
		
		// apply distortion model to object point, exterior parameters and principle distance
		DistortionModelFactory.apply(collinearityEquation, A, w, deltaX, deltaY, par_deltaX_xs, par_deltaX_ys, par_deltaY_xs, par_deltaY_ys);
		
		for (UnknownParameter<? extends DistortionModel> unknownParameter : distortionModel) {
			if (unknownParameter.getParameterType() != ParameterType.TANGENTIAL_POLYNOMIAL_B)
				continue;
			
			PolynomialCoefficient<?> Bi = (PolynomialCoefficient<?>)unknownParameter;
			
			double bi = Bi.getValue();
			int expi = Bi.getOrder();
			
			double ri = Math.pow(r2, expi);
			double dTani = bi * ri;
			sum += dTani;
			
			// distortion coefficients in x and y
			double deltaXi = deltaX * dTani;
			double deltaYi = deltaY * dTani;
			
			// chain rule coefficients
			double par_xs_Bi = deltaX * ri;
			double par_ys_Bi = deltaY * ri;
			
			double constTani = 2.0 * bi * expi * Math.pow(r2, expi - 1);
			double constTanXi = deltaX * constTani;
			double constTanYi = deltaY * constTani;	
			
			double par_deltaX_xsi = dTani * par_deltaX_xs + collinearityEquation.xs * constTanXi; 
			double par_deltaX_ysi = dTani * par_deltaX_ys + collinearityEquation.ys * constTanXi;
					
			double par_deltaY_xsi = dTani * par_deltaY_xs + collinearityEquation.xs * constTanYi; 
			double par_deltaY_ysi = dTani * par_deltaY_ys + collinearityEquation.ys * constTanYi;
			
			// apply distortion model to object point, exterior parameters and principle distance
			DistortionModelFactory.apply(collinearityEquation, A, w, deltaXi, deltaYi, par_deltaX_xsi, par_deltaX_ysi, par_deltaY_xsi, par_deltaY_ysi);
			
			//  dTanX: par_dTanX_xs * par_xs_<PARAM>  +  par_dTanX_ys * par_ys_<PARAM>
			//  dTanY: par_dTanY_xs * par_xs_<PARAM>  +  par_dTanY_ys * par_ys_<PARAM>
			
			/** disortion model **/
			int column = Bi.getColumn();
			if (column >= 0 && column != Integer.MAX_VALUE) {
				columns.add(column);
				A.set(0, column, par_xs_Bi);
				A.set(1, column, par_ys_Bi);
			}
		}
				
		// chain rule coefficients
		double par_xs_B1 = sum * (r2 + xxs2);
		double par_ys_B1 = sum * xys2;
		
		double par_xs_B2 = sum * xys2;
		double par_ys_B2 = sum * (r2 + yys2);
		
		/** disortion model **/
		int column = Bx.getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B1);
			A.set(1, column, par_ys_B1);
		}

		column = By.getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B2);
			A.set(1, column, par_ys_B2);
		}
	}
}