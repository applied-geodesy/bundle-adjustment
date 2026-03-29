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
import org.applied_geodesy.adjustment.bundle.parameter.ZernikeCoefficient;

public abstract class ZernikeDistortionModel extends PolynomialDistortionModel {
	
	public static class X extends ZernikeDistortionModel {
		public X(Camera camera, double r0) {
			super(camera, r0, Type.ZERNIKE_X);
		}
		
		@Override
		public ZernikeCoefficient add(int order) {
			return this.add(order, ParameterType.ZERNIKE_POLYNOMIAL_X);
		}
	}
	
	public static class Y extends ZernikeDistortionModel {
		public Y(Camera camera, double r0) {
			super(camera, r0, Type.ZERNIKE_Y);
		}
		
		@Override
		public ZernikeCoefficient add(int order) {
			return this.add(order, ParameterType.ZERNIKE_POLYNOMIAL_Y);
		}
	}
	
	public static class Gradient extends ZernikeDistortionModel {
		public Gradient(Camera camera, double r0) {
			super(camera, r0, Type.ZERNIKE_GRADIENT);
		}
		
		@Override
		public ZernikeCoefficient add(int order) {
			return this.add(order, ParameterType.ZERNIKE_POLYNOMIAL_Z);
		}
	}

	private final Type type;
	ZernikeDistortionModel(Camera camera, double r0, Type type) {
		super(camera, r0);
		this.type = type;
	}
	
	ZernikeCoefficient add(int order, ParameterType parameterType) {
		if (order <= 0)
			throw new IllegalArgumentException("Error, polynomial coefficient order must be a real positive integer. " + order);
		
		ZernikeCoefficient coefficient = new ZernikeCoefficient(parameterType, this, order);
		super.add(order, coefficient);

		return coefficient;
	}
	
	@Override
	public final Type getType() {
		return this.type; 
	}
}
