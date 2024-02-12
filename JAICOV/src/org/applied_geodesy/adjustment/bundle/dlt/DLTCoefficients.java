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

package org.applied_geodesy.adjustment.bundle.dlt;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.Referenceable;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class DLTCoefficients implements Referenceable<Image>, Iterable<UnknownParameter<DLTCoefficients>> {
	private Map<ParameterType, UnknownParameter<DLTCoefficients>> params = new LinkedHashMap<ParameterType, UnknownParameter<DLTCoefficients>>(20);
	private final Image image;
	public DLTCoefficients(Image image) {
		this.image = image;
		
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14, this));
		
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24, this));
		
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32, this));
		this.params.put(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33, new UnknownParameter<DLTCoefficients>(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33, this));
		
		this.params.put(ParameterType.PRINCIPAL_POINT_X, new UnknownParameter<DLTCoefficients>(ParameterType.PRINCIPAL_POINT_X, this));
		this.params.put(ParameterType.PRINCIPAL_POINT_Y, new UnknownParameter<DLTCoefficients>(ParameterType.PRINCIPAL_POINT_Y, this));
		this.params.put(ParameterType.PRINCIPAL_DISTANCE, new UnknownParameter<DLTCoefficients>(ParameterType.PRINCIPAL_DISTANCE, this));
		
		this.params.put(ParameterType.CAMERA_COORDINATE_X, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_COORDINATE_X, this));
		this.params.put(ParameterType.CAMERA_COORDINATE_Y, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_COORDINATE_Y, this));
		this.params.put(ParameterType.CAMERA_COORDINATE_Z, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_COORDINATE_Z, this));
		
		this.params.put(ParameterType.CAMERA_OMEGA, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_OMEGA, this));
		this.params.put(ParameterType.CAMERA_PHI, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_PHI, this));
		this.params.put(ParameterType.CAMERA_KAPPA, new UnknownParameter<DLTCoefficients>(ParameterType.CAMERA_KAPPA, this));
	}
	
	public UnknownParameter<DLTCoefficients> get(ParameterType parameterType) {
		return this.params.get(parameterType);
	}

	@Override
	public Iterator<UnknownParameter<DLTCoefficients>> iterator() {
		return this.params.values().iterator();
	}

	@Override
	public Image getReference() {
		return this.image;
	}

	@Override
	public String toString() {
		return "DLTCoefficients [params=" + this.params + "]";
	}
}
