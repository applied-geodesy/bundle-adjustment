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

public enum ParameterType {
	
	/** Bildhauptpunkt x **/
	PRINCIPAL_POINT_X(111),
	/** Bildhauptpunkt y **/
	PRINCIPAL_POINT_Y(112),
	/** Brennweite c **/
	PRINCIPAL_DISTANCE(113),
	
	/** Polynomkoeffizient der radial-symmetrischen Verzeichnung A1 **/
	RADIAL_DISTORTION_A1(121),
	/** Polynomkoeffizient der radial-symmetrischen Verzeichnung A2 **/
	RADIAL_DISTORTION_A2(122),
	/** Polynomkoeffizient der radial-symmetrischen Verzeichnung A3 **/
	RADIAL_DISTORTION_A3(123),
	
	/** Radial-asymmetrische und tangentiale Verzeichnung B1 **/
	TANGENTIAL_DISTORTION_B1(131),
	/** Radial-asymmetrische und tangentiale Verzeichnung B2 **/
	TANGENTIAL_DISTORTION_B2(132),
	
	/** Affinitaet und Scherung C1 **/
	AFFINITY_AND_SHEAR_C1(141),
	/** Affinitaet und Scherung C2 **/
	AFFINITY_AND_SHEAR_C2(142),
	
	/** Polynomkoeffizient der streckenabhaenigen Verzeichnung D1 **/
	DISTANCE_DISTORTION_D1(151),
	/** Polynomkoeffizient der streckenabhaenigen Verzeichnung D2 **/
	DISTANCE_DISTORTION_D2(152),
	/** Polynomkoeffizient der streckenabhaenigen Verzeichnung D3 **/
	DISTANCE_DISTORTION_D3(153),
	
	/** Koordinate der Kamera X0 **/
	CAMERA_COORDINATE_X(251),
	/** Koordinate der Kamera Y0 **/
	CAMERA_COORDINATE_Y(252),
	/** Koordinate der Kamera X0 **/
	CAMERA_COORDINATE_Z(253),
	
	/** Drehwinkel der Kamera Omega **/
	CAMERA_OMEGA(261),
	/** Drehwinkel der Kamera Phi **/
	CAMERA_PHI(262),
	/** Drehwinkel der Kamera Kappa **/
	CAMERA_KAPPA(263),
	
	/** Koordinate des Objektpunktes X **/
	OBJECT_COORDINATE_X(311),
	/** Koordinate des Objektpunktes Y **/
	OBJECT_COORDINATE_Y(312),
	/** Koordinate des Objektpunktes Z **/
	OBJECT_COORDINATE_Z(313),
	
	/** Koordinate des Bildpunktes xP **/
	IMAGE_COORDINATE_X(411),
	/** Koordinate des Bildpunktes yP **/
	IMAGE_COORDINATE_Y(412),
	
	/** Massstab s **/
	SCALE_BAR_LENGTH(511),
	
	/** DLT Parameter **/
	DIRECT_LINEAR_TRANSFORMATION_B11(611),
	DIRECT_LINEAR_TRANSFORMATION_B12(612),
	DIRECT_LINEAR_TRANSFORMATION_B13(613),
	DIRECT_LINEAR_TRANSFORMATION_B14(614),
	
	DIRECT_LINEAR_TRANSFORMATION_B21(621),
	DIRECT_LINEAR_TRANSFORMATION_B22(622),
	DIRECT_LINEAR_TRANSFORMATION_B23(623),
	DIRECT_LINEAR_TRANSFORMATION_B24(624),
	
	DIRECT_LINEAR_TRANSFORMATION_B31(631),
	DIRECT_LINEAR_TRANSFORMATION_B32(632),
	DIRECT_LINEAR_TRANSFORMATION_B33(633)
	;
	
	private int id;
	private ParameterType(int id) {
		this.id = id;
	}

	public final int getId() {
		return id;
	}

	public static ParameterType getEnumByValue(int value) {
		for(ParameterType element : ParameterType.values()) {
			if(element.id == value)
				return element;
		}
		return null;
	}  
}	
