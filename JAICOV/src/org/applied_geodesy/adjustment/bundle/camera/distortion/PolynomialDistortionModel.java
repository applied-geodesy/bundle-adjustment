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
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;

abstract class PolynomialDistortionModel extends DistortionModel {
	private final double r0;
	
	PolynomialDistortionModel(Camera camera, double r0) {
		super(camera);
		this.r0 = r0;
	}
	
	public final double getR0() {
		return this.r0;
	}
	
	/**
	 * Add polynomial coefficient to distortion model 
	 * @param order
	 * @return coefficient
	 */
	public abstract PolynomialCoefficient<?> add(int order);
	
	/**
	 * Returns the polynomial coefficient of the distortion model for the specified order  
	 * @param order
	 * @return coefficient
	 */
	public abstract PolynomialCoefficient<?> get(int order);
}
