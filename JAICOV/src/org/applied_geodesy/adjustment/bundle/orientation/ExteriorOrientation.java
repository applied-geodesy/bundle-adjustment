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

public class ExteriorOrientation implements Iterable<UnknownParameter<ExteriorOrientation>> {
	private Map<ParameterType, UnknownParameter<ExteriorOrientation>> params = new LinkedHashMap<ParameterType, UnknownParameter<ExteriorOrientation>>(6);

	public ExteriorOrientation() {
		this.params.put(ParameterType.CAMERA_COORDINATE_X, new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_COORDINATE_X, this));
		this.params.put(ParameterType.CAMERA_COORDINATE_Y, new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_COORDINATE_Y, this));
		this.params.put(ParameterType.CAMERA_COORDINATE_Z, new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_COORDINATE_Z, this));

		this.params.put(ParameterType.CAMERA_OMEGA, new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_OMEGA, this));
		this.params.put(ParameterType.CAMERA_PHI,   new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_PHI, this));
		this.params.put(ParameterType.CAMERA_KAPPA, new UnknownParameter<ExteriorOrientation>(ParameterType.CAMERA_KAPPA, this));
	}

	public UnknownParameter<ExteriorOrientation> get(ParameterType parameterType) {
		return this.params.get(parameterType);
	}
	
	public double[][] getRotationMatrix() {
		
		double omega = this.params.get(ParameterType.CAMERA_OMEGA).getValue();
		double phi   = this.params.get(ParameterType.CAMERA_PHI).getValue();
		double kappa = this.params.get(ParameterType.CAMERA_KAPPA).getValue();
		
		double cosOmega = Math.cos(omega);
		double sinOmega = Math.sin(omega);
		
		double cosPhi = Math.cos(phi);
		double sinPhi = Math.sin(phi);
		
		double cosKappa = Math.cos(kappa);
		double sinKappa = Math.sin(kappa);
		
		// Rotationsmatrix (Gl 2.30, S. 61)
        double r11 =  cosPhi * cosKappa;
        double r12 = -cosPhi * sinKappa;
        double r13 =  sinPhi;
        
        double r21 =  cosOmega * sinKappa + sinOmega * sinPhi * cosKappa;
        double r22 =  cosOmega * cosKappa - sinOmega * sinPhi * sinKappa;
        double r23 = -sinOmega * cosPhi;
        
        double r31 = sinOmega * sinKappa - cosOmega * sinPhi * cosKappa;
        double r32 = sinOmega * cosKappa + cosOmega * sinPhi * sinKappa;
        double r33 = cosOmega * cosPhi;
        
        return new double[][] {
        	{r11, r12, r13},
        	{r21, r22, r23},
        	{r31, r32, r33}
        };
	}
	
	@Override
	public Iterator<UnknownParameter<ExteriorOrientation>> iterator() {
		return params.values().iterator();
	}
}
