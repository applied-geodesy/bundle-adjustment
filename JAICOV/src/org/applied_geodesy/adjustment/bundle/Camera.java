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

import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;

public class Camera implements Iterable<Image> {
	private final long id;
	private final double r0;
	private InteriorOrientation interiorOrientation = new InteriorOrientation();
	private Map<Long, Image> images = new LinkedHashMap<Long, Image>();
	
	public Camera(long id, double r0) {
		this.id = id;
		this.r0 = r0;
	}
	
	public InteriorOrientation getInteriorOrientation() {
		return this.interiorOrientation;
	}
	
	public final long getId() {
		return this.id;
	}
	
	public final double getR0() {
		return this.r0;
	}
	
	public int getNumberOfImages() {
		return this.images.size();
	}
	
	public void removeAll() {
		this.images.clear();
	}
	
	public boolean remove(Image image) {
		if (!this.images.containsKey(image.getId()))
			return false;
		this.images.remove(image.getId());
		return true;
	}
	
	public Image add(long imageId) {
		if (this.images.containsKey(imageId))
			return this.images.get(imageId);
		Image image = new Image(imageId, this);
		this.images.put(imageId, image);
		return image;
	}
	
//	public Image get(int index) {
//		return this.images.get(index);
//	}

	@Override
	public Iterator<Image> iterator() {
		return this.images.values().iterator();
	}

	@Override
	public String toString() {
		return "Camera [id=" + id + ", r0=" + r0 + ", images=" + images.size() + "]";
	}
}
