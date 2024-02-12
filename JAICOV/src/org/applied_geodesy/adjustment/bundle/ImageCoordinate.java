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

import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

public class ImageCoordinate implements Referenceable<Image> {
	private final ObjectCoordinate objectCoordinate;
	private final Image image;
	
	private ObservationParameter<ImageCoordinate> x = new ObservationParameter<ImageCoordinate>(ParameterType.IMAGE_COORDINATE_X, this);
	private ObservationParameter<ImageCoordinate> y = new ObservationParameter<ImageCoordinate>(ParameterType.IMAGE_COORDINATE_Y, this);
	
	ImageCoordinate(ObjectCoordinate objectCoordinate, Image image, double xp, double yp, double sigmax, double sigmay) {
		this.objectCoordinate = objectCoordinate;
		this.image = image;
		
		this.x.setValue(xp);
		this.y.setValue(yp);
		
		this.x.setVariance(sigmax * sigmax);
		this.y.setVariance(sigmay * sigmay);
	}
	
	public ObjectCoordinate getObjectCoordinate() {
		return this.objectCoordinate;
	}
	
	public ObservationParameter<ImageCoordinate> getX() {
		return this.x;
	}
	
	public ObservationParameter<ImageCoordinate> getY() {
		return this.y;
	}

	@Override
	public Image getReference() {
		return this.image;
	}

	@Override
	public String toString() {
		return "ImageCoordinate [name=" + this.objectCoordinate.getName() + ", x=" + this.x.getValue() + ", y=" + this.y.getValue() + ", image=" + this.image.getId() + ", ]";
	}
}
