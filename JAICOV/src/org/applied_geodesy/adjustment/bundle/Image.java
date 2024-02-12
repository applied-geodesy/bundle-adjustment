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
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;

public class Image implements Referenceable<Camera>, Iterable<ImageCoordinate> {
	private final long id;
	private final Camera camera;
	
	private ExteriorOrientation exteriorOrientation = new ExteriorOrientation();
	private Map<ObjectCoordinate, ImageCoordinate> imageCoordinates = new LinkedHashMap<ObjectCoordinate, ImageCoordinate>();
	
	Image(long id, Camera camera) {
		this.id = id;
		this.camera = camera;
	}
	
	public ExteriorOrientation getExteriorOrientation() {
		return this.exteriorOrientation;
	}
	
	public final long getId() {
		return this.id;
	}
	
	public ImageCoordinate add(ObjectCoordinate objectCoordinate, double xp, double yp, double sigmax, double sigmay) {
		if (this.imageCoordinates.containsKey(objectCoordinate))
			return this.imageCoordinates.get(objectCoordinate);
		
		ImageCoordinate imgCoord = new ImageCoordinate(objectCoordinate, this, xp, yp, sigmax, sigmay);
		this.imageCoordinates.put(objectCoordinate, imgCoord);
		
		objectCoordinate.add(this);
		return imgCoord;
	}
	
	public ImageCoordinate get(ObjectCoordinate objectCoordinate) {
		return this.imageCoordinates.get(objectCoordinate);
	}
	
	public int getNumberOfImageCoordinates() {
		return this.imageCoordinates.size();
	}

	@Override
	public Camera getReference() {
		return this.camera;
	}

	@Override
	public Iterator<ImageCoordinate> iterator() {
		return this.imageCoordinates.values().iterator();
	}

	@Override
	public String toString() {
		return "Image [id=" + id + ", camera=" + camera.getId() + "]";
	}
}
