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

import java.util.Iterator;

import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class AffinityShearDistortionModel extends DistortionModel {

	private final UnknownParameter<AffinityShearDistortionModel> Cx = new UnknownParameter<AffinityShearDistortionModel>(ParameterType.AFFINITY_AND_SHEAR_Cx, this);
	private final UnknownParameter<AffinityShearDistortionModel> Cy = new UnknownParameter<AffinityShearDistortionModel>(ParameterType.AFFINITY_AND_SHEAR_Cy, this);
	
	public AffinityShearDistortionModel(Camera camera) {
		super(camera);
		
		this.Cx.setColumn(Integer.MAX_VALUE);
		this.Cy.setColumn(Integer.MAX_VALUE);
	}

	/**
	 * Returns the x coefficient Cx
	 * dAffX = Cx * xs + Cy * ys;
	 * dAffY = 0; 
	 * @return Cx
	 */
	public UnknownParameter<AffinityShearDistortionModel> getCx() {
		return this.Cx;
	}

	/**
	 * Returns the y coefficient Cy
	 * dAffX = Cx * xs + Cy * ys;
	 * dAffY = 0; 
	 * @return Cy
	 */
	public UnknownParameter<AffinityShearDistortionModel> getCy() {
		return this.Cy;
	}

	@Override
	public Iterator<UnknownParameter<? extends DistortionModel>> iterator() {
		return new Iterator<UnknownParameter<? extends DistortionModel>>() {
			private byte component = 0;
			
			@Override
			public boolean hasNext() {
				return this.component < 2;
			}

			@Override
			public UnknownParameter<? extends DistortionModel> next() {
				if (this.component > 1)
					return null;
				
				return this.component++ == 0 ? Cx : Cy;
			}
		};
	}
	
	@Override
	public final int getNumberOfParameters() {
		return 2;
	}
	
	@Override
	public Type getType() {
		return Type.AFFINITY_AND_SHEAR;
	}

	@Override
	public String toString() {
		return "AffinityShearDistortionModel [Cx=" + Cx + ", Cy=" + Cy + "]";
	}
}
