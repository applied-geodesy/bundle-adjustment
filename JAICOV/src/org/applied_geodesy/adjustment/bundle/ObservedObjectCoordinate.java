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

import java.util.Iterator;

import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

public class ObservedObjectCoordinate implements Referenceable<ObjectCoordinate>, ObservationParameterGroup<ObservedObjectCoordinate> {
	private final ObjectCoordinate objectCoordinate;
	private final double corrCoefXY;
	private final double corrCoefXZ;
	private final double corrCoefYZ;
	
	private ObservationParameter<ObservedObjectCoordinate> x = new ObservationParameter<ObservedObjectCoordinate>(ParameterType.OBJECT_COORDINATE_X, this);
	private ObservationParameter<ObservedObjectCoordinate> y = new ObservationParameter<ObservedObjectCoordinate>(ParameterType.OBJECT_COORDINATE_Y, this);
	private ObservationParameter<ObservedObjectCoordinate> z = new ObservationParameter<ObservedObjectCoordinate>(ParameterType.OBJECT_COORDINATE_Z, this);
			
	
	ObservedObjectCoordinate(ObjectCoordinate objectCoordinate, double x, double y, double z, double sigmax, double sigmay, double sigmaz, double corrCoefXY, double corrCoefXZ, double corrCoefYZ) {
		if (Math.abs(corrCoefXY) >= 1)
			throw new IllegalArgumentException("Error, correlation coefficient rho(x,y) must be in the open interval (-1 1): " + corrCoefXY);
		
		if (Math.abs(corrCoefXZ) >= 1)
			throw new IllegalArgumentException("Error, correlation coefficient rho(x,z) must be in the open interval (-1 1): " + corrCoefXZ);
		
		if (Math.abs(corrCoefYZ) >= 1)
			throw new IllegalArgumentException("Error, correlation coefficient rho(y,z) must be in the open interval (-1 1): " + corrCoefYZ);
		
		this.objectCoordinate = objectCoordinate;
		
		this.x.setValue(x);
		this.y.setValue(y);
		this.z.setValue(z);
		
		this.x.setVariance(sigmax * sigmax);
		this.y.setVariance(sigmay * sigmay);
		this.z.setVariance(sigmaz * sigmaz);
		
		this.corrCoefXY = corrCoefXY;
		this.corrCoefXZ = corrCoefXZ;
		this.corrCoefYZ = corrCoefYZ;
	}
	
	public ObjectCoordinate getObjectCoordinate() {
		return this.objectCoordinate;
	}
	
	public ObservationParameter<ObservedObjectCoordinate> getX() {
		return this.x;
	}
	
	public ObservationParameter<ObservedObjectCoordinate> getY() {
		return this.y;
	}
	
	public ObservationParameter<ObservedObjectCoordinate> getZ() {
		return this.z;
	}
	
	public double getCorrelationCoefficientXY() {
		return this.corrCoefXY;
	}
	
	public double getCorrelationCoefficientXZ() {
		return this.corrCoefXZ;
	}
	
	public double getCorrelationCoefficientYZ() {
		return this.corrCoefYZ;
	}

	@Override
	public Iterator<ObservationParameter<ObservedObjectCoordinate>> iterator() {
		return new Iterator<ObservationParameter<ObservedObjectCoordinate>>() {
			private byte component = 0;
			@Override
			public boolean hasNext() {
				return this.component < 3;
			}

			@Override
			public ObservationParameter<ObservedObjectCoordinate> next() {			
				switch(this.component++) {
				case 0:
					return x;
				case 1:
					return y;
				case 2:
					return z;
				default:
					return null;
				}
			}
		};
	}
	
	@Override
	public int getNumberOfParameters() {
		return 3;
	}

	@Override
	public ObjectCoordinate getReference() {
		return this.objectCoordinate;
	}
	
	@Override
	public String toString() {
		return "ObservedObjectCoordinate [objectCoordinate=" + objectCoordinate + ", corrCoefXY=" + corrCoefXY
				+ ", corrCoefXZ=" + corrCoefXZ + ", corrCoefYZ=" + corrCoefYZ + ", x=" + x + ", y=" + y + ", z=" + z
				+ "]";
	}
}
