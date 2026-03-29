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

package org.applied_geodesy.adjustment.bundle.camera;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.camera.distortion.AffinityShearDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadialDistanceDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadiallySymmetricDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.TangentialDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.ZernikeDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.orientation.InteriorOrientation;

public class Camera implements Iterable<Image> {
	private final long id;
	private InteriorOrientation interiorOrientation = new InteriorOrientation(this);
	private Map<DistortionModel.Type, DistortionModel> distortionModels;
	private Map<Long, Image> images = new LinkedHashMap<Long, Image>();
	
	public Camera(long id, double r0, DistortionModel.Type... distortionModelTypes) {
		this.id = id;
		
		if (distortionModelTypes.length > 0) {
			Map<DistortionModel.Type, DistortionModel> models = new LinkedHashMap<DistortionModel.Type, DistortionModel>(distortionModelTypes.length);
			Arrays.sort(distortionModelTypes); // ensuring always identical order
			for (DistortionModel.Type type : distortionModelTypes) {
				if (models.containsKey(type))
					throw new IllegalArgumentException("Error, duplicate type of distortion model detected. " + type);
				
				switch (type) {
				case AFFINITY_AND_SHEAR:
					models.put(type, new AffinityShearDistortionModel(this));
					break;
				case TANGENTIAL_DISTORTION:
					models.put(type, new TangentialDistortionModel(this));
					break;
				case RADIAL_DISTORTION:
					models.put(type, new RadiallySymmetricDistortionModel(this, r0));
					break;
				case DISTANCE_DISTORTION:
					models.put(type, new RadialDistanceDistortionModel(this, r0));
					break;
				case ZERNIKE_GRADIENT:
					models.put(type, new ZernikeDistortionModel.Gradient(this, r0));
					break;
				case ZERNIKE_X:
					models.put(type, new ZernikeDistortionModel.X(this, r0));
					break;
				case ZERNIKE_Y:
					models.put(type, new ZernikeDistortionModel.Y(this, r0));
					break;
				}
			}
			this.distortionModels = Collections.unmodifiableMap(models);
		}
		else
			this.distortionModels = Collections.emptyMap();
	}
	
	public InteriorOrientation getInteriorOrientation() {
		return this.interiorOrientation;
	}
	
	public final long getId() {
		return this.id;
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
	
	public DistortionModel getDistortionModel(DistortionModel.Type type) {
		return this.distortionModels.get(type);
	}
	
	public Collection<DistortionModel> getDistortionModels() {
		return this.distortionModels.values();
	}
	
	@Override
	public Iterator<Image> iterator() {
		return this.images.values().iterator();
	}

	@Override
	public String toString() {
		return "Camera [id=" + id + ", images=" + images.size() + "]";
	}
}
