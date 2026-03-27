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

package org.applied_geodesy.adjustment.bundle.camera.distortion;

import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class TangentialDistortionModel extends PolynomialDistortionModel {
	private final UnknownParameter<TangentialDistortionModel> Bx = new UnknownParameter<TangentialDistortionModel>(ParameterType.TANGENTIAL_DISTORTION_Bx, this);
	private final UnknownParameter<TangentialDistortionModel> By = new UnknownParameter<TangentialDistortionModel>(ParameterType.TANGENTIAL_DISTORTION_By, this);
	
	
	public TangentialDistortionModel(Camera camera) {
		super(camera, 0);
		
		this.Bx.setColumn(Integer.MAX_VALUE);
		this.By.setColumn(Integer.MAX_VALUE);
				
		super.add(-1, this.Bx);
		super.add(-2, this.By);
	}
	
	/**
	 * Returns the x coefficient Bx
	 * 
	 *                 2       2                              2          4               (2*i)
	 * dTanX = [Bx * (r  +  2*x)  +  By * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * 
	 *                 2       2                              2          4               (2*i)
	 * dTanY = [By * (r  +  2*y)  +  Bx * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * 
	 * @return Bx
	 */
	public UnknownParameter<TangentialDistortionModel> getBx() {
		return this.Bx;
	}

	/**
	 * Returns the y coefficient By
	 * 
	 *                 2       2                              2          4               (2*i)
	 * dTanX = [Bx * (r  +  2*x)  +  By * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * 
	 *                 2       2                              2          4               (2*i)
	 * dTanY = [By * (r  +  2*y)  +  Bx * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * 
	 * @return By
	 */
	public UnknownParameter<TangentialDistortionModel> getBy() {
		return this.By;
	}
	
	/**
	 * Add coefficient to distortion model 
	 *         2       2                              2          4               (2*i)
	 * [Bx * (r  +  2*x)  +  By * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * 
	 *         2       2                              2          4               (2*i)
	 * [By * (r  +  2*y)  +  Bx * 2*x*y] * (1 + B1 * r  +  B2 * r   + ... + Bi * r)
	 * @param order
	 * @return coefficient
	 */
	@Override
	public PolynomialCoefficient<TangentialDistortionModel> add(int order) {
		if (order <= 0)
			throw new IllegalArgumentException("Error, polynomial coefficient order must be a real positive integer. " + order);

		PolynomialCoefficient<TangentialDistortionModel> coefficient = new PolynomialCoefficient<TangentialDistortionModel>(ParameterType.TANGENTIAL_POLYNOMIAL_B, this, order);
		super.add(order, coefficient);

		return coefficient;
	}
	
	@Override
	public Type getType() {
		return Type.TANGENTIAL_DISTORTION;
	}
}
