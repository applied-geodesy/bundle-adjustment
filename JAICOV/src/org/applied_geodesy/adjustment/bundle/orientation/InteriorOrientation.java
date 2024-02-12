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

package org.applied_geodesy.adjustment.bundle.orientation;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class InteriorOrientation implements Iterable<UnknownParameter<InteriorOrientation>> {
	private Map<ParameterType, UnknownParameter<InteriorOrientation>> params = new LinkedHashMap<ParameterType, UnknownParameter<InteriorOrientation>>(10);
	
	public InteriorOrientation() {
		this.params.put(ParameterType.PRINCIPAL_POINT_X,  new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_POINT_X, this));
		this.params.put(ParameterType.PRINCIPAL_POINT_Y,  new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_POINT_Y, this));
		
		this.params.put(ParameterType.PRINCIPAL_DISTANCE, new UnknownParameter<InteriorOrientation>(ParameterType.PRINCIPAL_DISTANCE, this));
		
		this.params.put(ParameterType.RADIAL_DISTORTION_A1, new UnknownParameter<InteriorOrientation>(ParameterType.RADIAL_DISTORTION_A1, this));
		this.params.put(ParameterType.RADIAL_DISTORTION_A2, new UnknownParameter<InteriorOrientation>(ParameterType.RADIAL_DISTORTION_A2, this));
		this.params.put(ParameterType.RADIAL_DISTORTION_A3, new UnknownParameter<InteriorOrientation>(ParameterType.RADIAL_DISTORTION_A3, this));
		
		this.params.put(ParameterType.TANGENTIAL_DISTORTION_B1, new UnknownParameter<InteriorOrientation>(ParameterType.TANGENTIAL_DISTORTION_B1, this));
		this.params.put(ParameterType.TANGENTIAL_DISTORTION_B2, new UnknownParameter<InteriorOrientation>(ParameterType.TANGENTIAL_DISTORTION_B2, this));
		
		this.params.put(ParameterType.AFFINITY_AND_SHEAR_C1, new UnknownParameter<InteriorOrientation>(ParameterType.AFFINITY_AND_SHEAR_C1, this));
		this.params.put(ParameterType.AFFINITY_AND_SHEAR_C2, new UnknownParameter<InteriorOrientation>(ParameterType.AFFINITY_AND_SHEAR_C2, this));
	
		this.params.put(ParameterType.DISTANCE_DISTORTION_D1, new UnknownParameter<InteriorOrientation>(ParameterType.DISTANCE_DISTORTION_D1, this));
		this.params.put(ParameterType.DISTANCE_DISTORTION_D2, new UnknownParameter<InteriorOrientation>(ParameterType.DISTANCE_DISTORTION_D2, this));
		this.params.put(ParameterType.DISTANCE_DISTORTION_D3, new UnknownParameter<InteriorOrientation>(ParameterType.DISTANCE_DISTORTION_D3, this));	
	}
	
	public UnknownParameter<InteriorOrientation> get(ParameterType parameterType) {
		return this.params.get(parameterType);
	}

	@Override
	public Iterator<UnknownParameter<InteriorOrientation>> iterator() {
		return params.values().iterator();
	}

	@Override
	public String toString() {
		return "InteriorOrientation [params=" + params + "]";
	}
}
