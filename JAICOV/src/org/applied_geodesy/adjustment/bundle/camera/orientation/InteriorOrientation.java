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

package org.applied_geodesy.adjustment.bundle.camera.orientation;

import java.util.Iterator;

import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.camera.Modelable;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class InteriorOrientation implements Modelable, Iterable<UnknownParameter<InteriorOrientation>> {
	private final Camera camera;

	private final UnknownParameter<InteriorOrientation> principlePointX   = new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_POINT_X, this);
	private final UnknownParameter<InteriorOrientation> principlePointY   = new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_POINT_Y, this);
	private final UnknownParameter<InteriorOrientation> principleDistance = new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_DISTANCE, this);
	
	public InteriorOrientation(Camera camera) {
		this.camera = camera;
	}
	
	public UnknownParameter<InteriorOrientation> getPrinciplePointX() {
		return this.principlePointX;
	}
	
	public UnknownParameter<InteriorOrientation> getPrinciplePointY() {
		return this.principlePointY;
	}
	
	public UnknownParameter<InteriorOrientation> getPrincipleDistance() {
		return this.principleDistance;
	}
	
	@Override
	public final int getNumberOfParameters() {
		return 3;
	}

	@Override
	public Iterator<UnknownParameter<InteriorOrientation>> iterator() {
		return new Iterator<UnknownParameter<InteriorOrientation>>() {
			private byte component = 0;
			@Override
			public boolean hasNext() {
				return this.component < getNumberOfParameters();
			}

			@Override
			public UnknownParameter<InteriorOrientation> next() {
				switch (this.component++) {
				case 0:
					return principlePointX;
				case 1:
					return principlePointY;
				case 2:
					return principleDistance;
				default:
					return null;	
				}
			}
		};
	}
	
	@Override
	public Camera getReference() {
		return this.camera;
	}

	@Override
	public String toString() {
		return "InteriorOrientation [camera=" + camera + ", principlePointX=" + principlePointX + ", principlePointY="
				+ principlePointY + ", principleDistance=" + principleDistance + "]";
	}

}
