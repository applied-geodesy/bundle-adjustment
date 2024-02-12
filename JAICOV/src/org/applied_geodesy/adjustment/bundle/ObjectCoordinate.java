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
import java.util.LinkedHashSet;
import java.util.Set;

import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

public class ObjectCoordinate implements Iterable<Image> {
	private final String name;
	private boolean datum = Boolean.TRUE;
	
	private UnknownParameter<ObjectCoordinate> x = new UnknownParameter<ObjectCoordinate>(ParameterType.OBJECT_COORDINATE_X, this);
	private UnknownParameter<ObjectCoordinate> y = new UnknownParameter<ObjectCoordinate>(ParameterType.OBJECT_COORDINATE_Y, this);
	private UnknownParameter<ObjectCoordinate> z = new UnknownParameter<ObjectCoordinate>(ParameterType.OBJECT_COORDINATE_Z, this);
	
	private Set<Image> images = new LinkedHashSet<Image>();
	
	public ObjectCoordinate(String name, double x, double y, double z) {
		this.name = name;
		this.x.setValue(x);
		this.y.setValue(y);
		this.z.setValue(z);
	}
	
	public String getName() {
		return this.name;
	}
	
	public UnknownParameter<ObjectCoordinate> getX() {
		return this.x;
	}
	
	public UnknownParameter<ObjectCoordinate> getY() {
		return this.y;
	}
	
	public UnknownParameter<ObjectCoordinate> getZ() {
		return this.z;
	}
	
	public void removeAll() {
		this.images.clear();
	}
	
	public boolean remove(Image image) {
		if (!this.images.contains(image))
			return false;
		this.images.remove(image);
		return true;
	}
	
	public boolean add(Image image) {
		if (this.images.contains(image) || image.get(this) == null) // ImageCoordinate#add
			return false;
		
		this.images.add(image);
		return true;
	}
	
	public boolean isDatum() {
		return this.datum;
	}
	
	public void setDatum(boolean datum) {
		this.datum = datum;
	}

	@Override
	public Iterator<Image> iterator() {
		return this.images.iterator();
	}

	@Override
	public String toString() {
		return "ObjectCoordinate [name=" + this.name + ", datum=" + this.datum + ", x=" + this.x.getValue() + ", y=" + this.y.getValue() + ", z=" + this.z.getValue() + ", images=" + images.size() + "]";
	}
	
}

