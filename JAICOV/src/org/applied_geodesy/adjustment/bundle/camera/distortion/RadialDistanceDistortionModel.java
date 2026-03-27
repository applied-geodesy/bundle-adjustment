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

public class RadialDistanceDistortionModel extends RadiallySymmetricDistortionModel {

	public RadialDistanceDistortionModel(Camera camera, double r0) {
		super(camera, r0);
	}
	
	/**
	 * Add coefficient to distortion model 
	 *        2          4          6                (2*i)           2           4           6                (2*i) 
	 * (D1 * r  +  D2 * r  +  D3 * r  + ... +  Di * r      - (D1 * r0  +  D2 * r0  +  D3 * r0  + ... +  Di * r0) ) / N
	 * @param order
	 * @return coefficient
	 */
	@Override
	public PolynomialCoefficient<? extends DistortionModel> add(int order) {
		if (order <= 0)
			throw new IllegalArgumentException("Error, polynomial coefficient order must be a real positive integer. " + order);

		PolynomialCoefficient<RadialDistanceDistortionModel> coefficient = new PolynomialCoefficient<RadialDistanceDistortionModel>(ParameterType.DISTANCE_POLYNOMIAL_D, this, order);
		super.add(order, coefficient);

		return coefficient;
	}
	
	@Override
	public Type getType() {
		return Type.DISTANCE_DISTORTION;
	}
}
