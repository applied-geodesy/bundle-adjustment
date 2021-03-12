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

package org.applied_geodesy.adjustment.defect;

public class RankDefect {
	private DefectType  rx   = DefectType.NOT_SET, 
						ry   = DefectType.NOT_SET,
						rz   = DefectType.NOT_SET,
						tx   = DefectType.NOT_SET,
						ty   = DefectType.NOT_SET,
						tz   = DefectType.NOT_SET,
						mxyz = DefectType.NOT_SET;
	
	public RankDefect() {}
	
	public void setScale(DefectType defectType) {
		this.mxyz = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
		
	public void setRotationX(DefectType defectType) {
		this.rx = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
	
	public void setRotationY(DefectType defectType) {
		this.ry = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
	
	public void setRotationZ(DefectType defectType) {
		this.rz = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
	
	public void setTranslationX(DefectType defectType) {
		this.tx = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
	
	public void setTranslationY(DefectType defectType) {
		this.ty = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
	
	public void setTranslationZ(DefectType defectType) {
		this.tz = defectType == DefectType.FIXED?DefectType.FIXED:DefectType.FREE;
	}
		
	public boolean estimateScale() {
		return this.mxyz == DefectType.FREE;
	}
	
	public boolean estimateRotationX() {
		return this.rx == DefectType.FREE;
	}
	
	public boolean estimateRotationY() {
		return this.ry == DefectType.FREE;
	}
	
	public boolean estimateRotationZ() {
		return this.rz == DefectType.FREE;
	}
	
	public boolean estimateTranslationX() {
		return this.tx == DefectType.FREE;
	}
	
	public boolean estimateTranslationY() {
		return this.ty == DefectType.FREE;
	}
	
	public boolean estimateTranslationZ() {
		return this.tz == DefectType.FREE;
	}
	
	public DefectType getScale() {
		return this.mxyz;
	}
		
	public DefectType getRotationX() {
		return this.rx;
	}
	
	public DefectType getRotationY() {
		return this.ry;
	}
	
	public DefectType getRotationZ() {
		return this.rz;
	}
	
	public DefectType getTranslationX() {
		return this.tx;
	}
	
	public DefectType getTranslationY() {
		return this.ty;
	}
	
	public DefectType getTranslationZ() {
		return this.tz;
	}
	
	public int getDefect() {
		int d = 0;
		d += this.estimateScale()?1:0;
		d += this.estimateRotationX()?1:0;
		d += this.estimateRotationY()?1:0;
		d += this.estimateRotationZ()?1:0;
		d += this.estimateTranslationX()?1:0;
		d += this.estimateTranslationY()?1:0;
		d += this.estimateTranslationZ()?1:0;
		
		return d;
	}
	
	public void reset() {
		this.tx = this.ty = this.tz = DefectType.NOT_SET;
		this.rx = this.ry = this.rz = DefectType.NOT_SET;
		this.mxyz = DefectType.NOT_SET;
	}

	@Override
	public String toString() {
		return "RankDefect [rx=" + rx + ", ry=" + ry + ", rz=" + rz + ", tx=" + tx + ", ty=" + ty + ", tz=" + tz
				+ ", mxyz=" + mxyz + ", defect=" + getDefect() + "]";
	}
}
