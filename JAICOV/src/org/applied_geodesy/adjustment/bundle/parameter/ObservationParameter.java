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

package org.applied_geodesy.adjustment.bundle.parameter;

import org.applied_geodesy.adjustment.bundle.Referenceable;

public class ObservationParameter<T> extends Parameter implements Referenceable<T> {
	private int row = -1;
	private double variance;
	private final T reference;
	
	public ObservationParameter(ParameterType parameterType, T reference) {
		super(parameterType);
		this.reference = reference;
	}
	
	public ObservationParameter(T parameter) { // Java 25
		if (!(parameter instanceof UnknownParameter<?>))
				throw new IllegalArgumentException("Error, parameter must be an UnknownParameter");
		this(((UnknownParameter<?>)parameter).getParameterType(), parameter);
		this.setValue(((UnknownParameter<?>)parameter).getValue());
	}
		
	public void setRow(int row) {
		this.row = row;
	}
	
	public int getRow() {
		return this.row;
	}
	
	public double getVariance() {
		return this.variance;
	}
	
	public void setVariance(double variance) {
		if (variance <= 0)
			throw new IllegalArgumentException("Error, variance must be positive: " + variance);
		this.variance = variance;
	}

	@Override
	public T getReference() {
		return this.reference;
	}
}