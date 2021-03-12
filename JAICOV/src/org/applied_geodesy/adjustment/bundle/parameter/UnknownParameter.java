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

public class UnknownParameter<T> extends Parameter implements Referenceable<T> {
	/** Value: -1 == not set, Integer.MAX_VALUE == fixed, else column in normal equation system **/
	private int column = -1;
	private final T reference;
	
	
	public UnknownParameter(ParameterType parameterType, T reference) {
		super(parameterType);
		this.reference = reference;
	}
	
	public void setColumn(int column) {
		this.column = column;
	}
	
	public int getColumn() {
		return this.column;
	}
	
	public T getReference() {
		return this.reference;
	}
}
