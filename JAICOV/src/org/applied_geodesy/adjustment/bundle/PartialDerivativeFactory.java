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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.UpperSymmBandMatrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;


class PartialDerivativeFactory {

	private PartialDerivativeFactory() {}
	
	static double getMisclosure(ObservationParameter<?> observation) {
		ParameterType observationParamType = observation.getParameterType();
		if (observationParamType == ParameterType.IMAGE_COORDINATE_X || observationParamType == ParameterType.IMAGE_COORDINATE_Y)
			return observation.getValue() - getCollinearityEquationValue( (ImageCoordinate)observation.getReference(), observationParamType);
		else if (observationParamType == ParameterType.SCALE_BAR_LENGTH) 
			return observation.getValue() - getDistanceValue((ScaleBar)observation.getReference());
		return 0;
	}
	
	public static double getDistanceValue(ScaleBar scaleBar) {
		ObjectCoordinate objectCoordinateA = scaleBar.getObjectCoordinateA();
		ObjectCoordinate objectCoordinateB = scaleBar.getObjectCoordinateB();
		
		double XA = objectCoordinateA.getX().getValue();
		double YA = objectCoordinateA.getY().getValue();
		double ZA = objectCoordinateA.getZ().getValue();
		
		double XB = objectCoordinateB.getX().getValue();
		double YB = objectCoordinateB.getY().getValue();
		double ZB = objectCoordinateB.getZ().getValue();
		
		double dX = XB-XA;
        double dY = YB-YA;
        double dZ = ZB-ZA;
        
        return Math.sqrt( dX*dX + dY*dY + dZ*dZ );
	}
	
	public static double getCollinearityEquationValue(ImageCoordinate imageCoordinate, ParameterType observationParamType) {
		Image image = imageCoordinate.getReference();
		Camera camera = image.getReference();
		ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
		
		double X = objectCoordinate.getX().getValue();
		double Y = objectCoordinate.getY().getValue();
		double Z = objectCoordinate.getZ().getValue();
		
		double r0 = camera.getR0();
		
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
		
		double c  = interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE).getValue();
		double x0 = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X).getValue();
		double y0 = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y).getValue();
	
		double A1 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A1).getValue();
		double A2 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A2).getValue();
		double A3 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).getValue();
		
		double B1 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B1).getValue();
		double B2 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B2).getValue();
		
		double C1 = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).getValue();
		double C2 = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).getValue();
		
		double X0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getValue();
		double Y0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getValue();
		double Z0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getValue();
		
		double omega = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getValue();
		double phi   = exteriorOrientation.get(ParameterType.CAMERA_PHI).getValue();
		double kappa = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getValue();

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
		
		double dX = X - X0;
		double dY = Y - Y0;
		double dZ = Z - Z0;
		
		double kx = r11*dX + r21*dY + r31*dZ;
		double ky = r12*dX + r22*dY + r32*dZ;
		double N  = r13*dX + r23*dY + r33*dZ;
		
		double kxN = kx/N;
		double kyN = ky/N;
		
		double xs = -c*kxN;
		double ys = -c*kyN;
        
		double r = Math.hypot(xs, ys);
		double dRad  = -(A1*r0*r0 + A2*r0*r0*r0*r0 + A3*r0*r0*r0*r0*r0*r0)*r + A1*r*r*r + A2*r*r*r*r*r + A3*r*r*r*r*r*r*r;
		
		double dRadX = xs * dRad/r;
		double dRadY = ys * dRad/r;
		
		double dTanX = B1 * (r*r + 2.0*xs*xs) + 2.0 * B2 * xs * ys;
		double dTanY = B2 * (r*r + 2.0*ys*ys) + 2.0 * B1 * xs * ys;
		
		double dAffX = C1*xs + C2*ys;
		double dAffY = 0;

		double deltaX = dRadX + dTanX + dAffX;
		double deltaY = dRadY + dTanY + dAffY;

		if (observationParamType == ParameterType.IMAGE_COORDINATE_X)
			return x0 + xs + deltaX;
		else if (observationParamType == ParameterType.IMAGE_COORDINATE_Y)
			return y0 + ys + deltaY;
		return 0;
	}
	
	static void getPartialDerivative(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ObservationParameter<?> observation) {
		ParameterType observationParamType = observation.getParameterType();

		if (observationParamType == ParameterType.IMAGE_COORDINATE_X) // ParameterType.IMAGE_COORDINATE_Y
			getPartialDerivativeImageCoordinate(sigma2apriori, NEQ, neq, (ImageCoordinate)observation.getReference());
		else if (observationParamType == ParameterType.SCALE_BAR_LENGTH)
			getPartialDerivativeScaleBar(sigma2apriori, NEQ, neq, (ScaleBar)observation.getReference());
	}
	
	private static void getPartialDerivativeScaleBar(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ScaleBar scaleBar) {
		ObjectCoordinate objectCoordinateA = scaleBar.getObjectCoordinateA();
		ObjectCoordinate objectCoordinateB = scaleBar.getObjectCoordinateB();
		
		double XA = objectCoordinateA.getX().getValue();
		double YA = objectCoordinateA.getY().getValue();
		double ZA = objectCoordinateA.getZ().getValue();
		
		double XB = objectCoordinateB.getX().getValue();
		double YB = objectCoordinateB.getY().getValue();
		double ZB = objectCoordinateB.getZ().getValue();
		
		double dX = XB - XA;
        double dY = YB - YA;
        double dZ = ZB - ZA;
        
        double lengthAB = Math.sqrt( dX*dX + dY*dY + dZ*dZ );
        
        double ax = dX/lengthAB;
        double ay = dY/lengthAB;
        double az = dZ/lengthAB;
        
        DenseMatrix A = new DenseMatrix(1, NEQ.numColumns());
        DenseVector w = new DenseVector(1);
        UpperSymmBandMatrix P = new UpperSymmBandMatrix(1, 0);
		
        P.set(0, 0, sigma2apriori/scaleBar.getLength().getVariance());
        w.set(0, scaleBar.getLength().getValue() - lengthAB);
        
        List<Integer> columns = new ArrayList<Integer>();
        int column = -1;
		
		// point coordinate A
		column = objectCoordinateA.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -ax);
		}
			
		column = objectCoordinateA.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -ay);
		}
		
		column = objectCoordinateA.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -az);
		}
		
		// point coordinate B
		column = objectCoordinateB.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) { 
			columns.add(column);
			A.set(0, column, +ax);
		}

		column = objectCoordinateB.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, +ay);
		}

		column = objectCoordinateB.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, +az);
		}
		
		Collections.sort(columns);
		for (int row = 0; row < A.numRows(); row++) {
			for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
				int colAT = columns.get(columnATIdx);
				double aT = A.get(row, colAT);
				neq.add(colAT, aT * P.get(row, row) * w.get(row));
				
				for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
					int colA = columns.get(columnAIdx);	
					double a  = A.get(row, colA);
					NEQ.add(colAT, colA, aT * P.get(row, row) * a);
				}
			}	
		}
	}
	
	private static void getPartialDerivativeImageCoordinate(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ImageCoordinate imageCoordinate) {
		
		ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
		Image image = imageCoordinate.getReference();
		Camera camera = image.getReference();
		
		double X = objectCoordinate.getX().getValue();
		double Y = objectCoordinate.getY().getValue();
		double Z = objectCoordinate.getZ().getValue();
		
		double r0 = camera.getR0();
		
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
		
		double c  = interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE).getValue();
		double x0 = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X).getValue();
		double y0 = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y).getValue();
	
		double A1 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A1).getValue();
		double A2 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A2).getValue();
		double A3 = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).getValue();
		
		double B1 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B1).getValue();
		double B2 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B2).getValue();
		
		double C1 = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).getValue();
		double C2 = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).getValue();
		
		double X0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getValue();
		double Y0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getValue();
		double Z0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getValue();
		
		double omega = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getValue();
		double phi   = exteriorOrientation.get(ParameterType.CAMERA_PHI).getValue();
		double kappa = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getValue();
		
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
		
		double dX = X - X0;
		double dY = Y - Y0;
		double dZ = Z - Z0;
		
		double kx = r11*dX + r21*dY + r31*dZ;
		double ky = r12*dX + r22*dY + r32*dZ;
		double N  = r13*dX + r23*dY + r33*dZ;
		
		double kxN = kx/N;
		double kyN = ky/N;
		
		double xs = -c*kxN;
		double ys = -c*kyN;
        
		double r = Math.hypot(xs, ys);
		double dRad = -(A1*r0*r0 + A2*r0*r0*r0*r0 + A3*r0*r0*r0*r0*r0*r0)*r + A1*r*r*r + A2*r*r*r*r*r + A3*r*r*r*r*r*r*r;
		
		double dRadX = xs * dRad/r;
		double dRadY = ys * dRad/r;
		
		double dTanX = B1 * (r*r + 2.0*xs*xs) + 2.0 * B2 * xs * ys;
		double dTanY = B2 * (r*r + 2.0*ys*ys) + 2.0 * B1 * xs * ys;
		
		double dAffX = C1*xs + C2*ys;
		double dAffY = 0;

		double deltaX = dRadX + dTanX + dAffX;
		double deltaY = dRadY + dTanY + dAffY;
		
		double varianceX = imageCoordinate.getX().getVariance();
		double varianceY = imageCoordinate.getY().getVariance();
				
		// 3 ... object point, 10 ... interior orientation, 6 ... exterior orientation, datum defect
		DenseVector w = new DenseVector(2);
		DenseMatrix A = new DenseMatrix(2, NEQ.numColumns());
		
		UpperSymmBandMatrix P = new UpperSymmBandMatrix(2, 0);
		P.set(0, 0, sigma2apriori/varianceX);
		P.set(1, 1, sigma2apriori/varianceY);
		
		w.set(0, imageCoordinate.getX().getValue() - (x0 + xs + deltaX));
		w.set(1, imageCoordinate.getY().getValue() - (y0 + ys + deltaY));
		
		List<Integer> columns = new ArrayList<Integer>();
		int column = -1;
				
		// Object point coordinates
		column = objectCoordinate.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -((c*r11)/N+(r13*xs)/N+(C1*c*r11)/N+(C2*c*r12)/N+(C1*r13*xs)/N+(C2*r13*ys)/N+(B1*(c*c)*(kxN*r11*-3.0+(kxN*kxN)*r13*3.0+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*2.0)/N+(c*r11*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)+(r13*xs*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)-(B2*(c*c)*kyN*r11*2.0)/N+(B2*(c*c)*cosPhi*kxN*sinKappa*2.0)/N+(B2*(c*c)*kxN*kyN*r13*4.0)/N+((c*c)*1.0/(r*r)*xs*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N+((c*c*c)*kxN*1.0/(r*r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/N));
			A.set(1, column, -((c*r12)/N+(r13*ys)/N+(B2*(c*c)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13*3.0+cosPhi*kyN*sinKappa*3.0)*2.0)/N+(c*r12*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)+(r13*ys*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)-(B1*(c*c)*kyN*r11*2.0)/N+(B1*(c*c)*cosPhi*kxN*sinKappa*2.0)/N+(B1*(c*c)*kxN*kyN*r13*4.0)/N+((c*c)*1.0/(r*r)*ys*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N+((c*c*c)*kyN*1.0/(r*r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/N));
		}
		
		column = objectCoordinate.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -((c*(cosOmega*sinKappa+C2*cosOmega*cosKappa+C1*cosOmega*sinKappa+cosPhi*kxN*sinOmega+cosKappa*r13*sinOmega+A1*cosOmega*(r*r)*sinKappa-A1*cosOmega*(r0*r0)*sinKappa+A2*cosOmega*(r*r*r*r)*sinKappa-A2*cosOmega*(r0*r0*r0*r0)*sinKappa+A3*cosOmega*(r*r*r*r*r*r)*sinKappa-A3*cosOmega*(r0*r0*r0*r0*r0*r0)*sinKappa+C1*cosPhi*kxN*sinOmega+C2*cosPhi*kyN*sinOmega+C1*cosKappa*r13*sinOmega-C2*r13*sinOmega*sinKappa+A1*(c*c)*cosPhi*(kxN*kxN*kxN)*sinOmega*2.0+A1*(c*c)*cosOmega*(kxN*kxN)*sinKappa*2.0-B2*c*cosOmega*cosKappa*kxN*2.0-B1*c*cosOmega*cosKappa*kyN*2.0-B1*c*cosOmega*kxN*sinKappa*6.0-B2*c*cosOmega*kyN*sinKappa*2.0-B1*c*cosPhi*(kxN*kxN)*sinOmega*6.0-B1*c*cosPhi*(kyN*kyN)*sinOmega*2.0+A1*cosPhi*kxN*(r*r)*sinOmega-A1*cosPhi*kxN*(r0*r0)*sinOmega+A2*cosPhi*kxN*(r*r*r*r)*sinOmega-A2*cosPhi*kxN*(r0*r0*r0*r0)*sinOmega+A3*cosPhi*kxN*(r*r*r*r*r*r)*sinOmega-A3*cosPhi*kxN*(r0*r0*r0*r0*r0*r0)*sinOmega+A1*cosKappa*(r*r)*r13*sinOmega-A1*cosKappa*(r0*r0)*r13*sinOmega+A2*cosKappa*(r*r*r*r)*r13*sinOmega-A2*cosKappa*(r0*r0*r0*r0)*r13*sinOmega+A3*cosKappa*(r*r*r*r*r*r)*r13*sinOmega-A3*cosKappa*(r0*r0*r0*r0*r0*r0)*r13*sinOmega+A1*(c*c)*cosOmega*cosKappa*kxN*kyN*2.0+A1*(c*c)*cosPhi*kxN*(kyN*kyN)*sinOmega*2.0+A1*(c*c)*cosKappa*(kxN*kxN)*r13*sinOmega*2.0-B2*c*cosPhi*kxN*kyN*sinOmega*4.0-B1*c*cosKappa*kxN*r13*sinOmega*6.0-B2*c*cosKappa*kyN*r13*sinOmega*2.0+B2*c*kxN*r13*sinOmega*sinKappa*2.0+B1*c*kyN*r13*sinOmega*sinKappa*2.0+A2*(c*c)*cosPhi*(kxN*kxN*kxN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kxN*kxN*kxN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosOmega*(kxN*kxN)*(r*r)*sinKappa*4.0+A3*(c*c)*cosOmega*(kxN*kxN)*(r*r*r*r)*sinKappa*6.0+A2*(c*c)*cosPhi*kxN*(kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*kxN*(kyN*kyN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosKappa*(kxN*kxN)*(r*r)*r13*sinOmega*4.0+A3*(c*c)*cosKappa*(kxN*kxN)*(r*r*r*r)*r13*sinOmega*6.0-A1*(c*c)*kxN*kyN*r13*sinOmega*sinKappa*2.0+A2*(c*c)*cosOmega*cosKappa*kxN*kyN*(r*r)*4.0+A3*(c*c)*cosOmega*cosKappa*kxN*kyN*(r*r*r*r)*6.0-A2*(c*c)*kxN*kyN*(r*r)*r13*sinOmega*sinKappa*4.0-A3*(c*c)*kxN*kyN*(r*r*r*r)*r13*sinOmega*sinKappa*6.0))/N));
			A.set(1, column, -((c*(cosOmega*cosKappa+cosPhi*kyN*sinOmega-r13*sinOmega*sinKappa+A1*cosOmega*cosKappa*(r*r)-A1*cosOmega*cosKappa*(r0*r0)+A2*cosOmega*cosKappa*(r*r*r*r)-A2*cosOmega*cosKappa*(r0*r0*r0*r0)+A3*cosOmega*cosKappa*(r*r*r*r*r*r)-A3*cosOmega*cosKappa*(r0*r0*r0*r0*r0*r0)+A1*(c*c)*cosPhi*(kyN*kyN*kyN)*sinOmega*2.0-B1*c*cosOmega*cosKappa*kxN*2.0-B2*c*cosOmega*cosKappa*kyN*6.0-B2*c*cosOmega*kxN*sinKappa*2.0-B1*c*cosOmega*kyN*sinKappa*2.0-B2*c*cosPhi*(kxN*kxN)*sinOmega*2.0-B2*c*cosPhi*(kyN*kyN)*sinOmega*6.0+A1*cosPhi*kyN*(r*r)*sinOmega-A1*cosPhi*kyN*(r0*r0)*sinOmega+A2*cosPhi*kyN*(r*r*r*r)*sinOmega-A2*cosPhi*kyN*(r0*r0*r0*r0)*sinOmega+A3*cosPhi*kyN*(r*r*r*r*r*r)*sinOmega-A3*cosPhi*kyN*(r0*r0*r0*r0*r0*r0)*sinOmega-A1*(r*r)*r13*sinOmega*sinKappa+A1*(r0*r0)*r13*sinOmega*sinKappa-A2*(r*r*r*r)*r13*sinOmega*sinKappa+A2*(r0*r0*r0*r0)*r13*sinOmega*sinKappa-A3*(r*r*r*r*r*r)*r13*sinOmega*sinKappa+A3*(r0*r0*r0*r0*r0*r0)*r13*sinOmega*sinKappa+A1*(c*c)*cosOmega*cosKappa*(kyN*kyN)*2.0+A1*(c*c)*cosOmega*kxN*kyN*sinKappa*2.0+A1*(c*c)*cosPhi*(kxN*kxN)*kyN*sinOmega*2.0-A1*(c*c)*(kyN*kyN)*r13*sinOmega*sinKappa*2.0-B1*c*cosPhi*kxN*kyN*sinOmega*4.0-B2*c*cosKappa*kxN*r13*sinOmega*2.0-B1*c*cosKappa*kyN*r13*sinOmega*2.0+A2*(c*c)*cosOmega*cosKappa*(kyN*kyN)*(r*r)*4.0+A3*(c*c)*cosOmega*cosKappa*(kyN*kyN)*(r*r*r*r)*6.0+B1*c*kxN*r13*sinOmega*sinKappa*2.0+B2*c*kyN*r13*sinOmega*sinKappa*6.0+A2*(c*c)*cosPhi*(kyN*kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kyN*kyN*kyN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosPhi*(kxN*kxN)*kyN*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kxN*kxN)*kyN*(r*r*r*r)*sinOmega*6.0-A2*(c*c)*(kyN*kyN)*(r*r)*r13*sinOmega*sinKappa*4.0-A3*(c*c)*(kyN*kyN)*(r*r*r*r)*r13*sinOmega*sinKappa*6.0+A1*(c*c)*cosKappa*kxN*kyN*r13*sinOmega*2.0+A2*(c*c)*cosOmega*kxN*kyN*(r*r)*sinKappa*4.0+A3*(c*c)*cosOmega*kxN*kyN*(r*r*r*r)*sinKappa*6.0+A2*(c*c)*cosKappa*kxN*kyN*(r*r)*r13*sinOmega*4.0+A3*(c*c)*cosKappa*kxN*kyN*(r*r*r*r)*r13*sinOmega*6.0))/N));
		}
		
		column = objectCoordinate.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -((r33*xs+C1*r33*xs+C2*r33*ys+c*sinOmega*sinKappa+A1*(r*r)*r33*xs-A1*(r0*r0)*r33*xs+A2*(r*r*r*r)*r33*xs-A2*(r0*r0*r0*r0)*r33*xs+A3*(r*r*r*r*r*r)*r33*xs-A3*(r0*r0*r0*r0*r0*r0)*r33*xs-A1*(c*c*c)*(kxN*kxN*kxN)*r33*2.0+B1*(c*c)*(kxN*kxN)*r33*6.0+B1*(c*c)*(kyN*kyN)*r33*2.0+C2*c*cosKappa*sinOmega+C1*c*sinOmega*sinKappa-c*cosOmega*cosKappa*r13-A1*(c*c*c)*kxN*(kyN*kyN)*r33*2.0+A1*(c*c*c)*(kxN*kxN)*sinOmega*sinKappa*2.0-C1*c*cosOmega*cosKappa*r13+C2*c*cosOmega*r13*sinKappa-A2*(c*c*c)*(kxN*kxN*kxN)*(r*r)*r33*4.0-A3*(c*c*c)*(kxN*kxN*kxN)*(r*r*r*r)*r33*6.0-B2*(c*c)*cosKappa*kxN*sinOmega*2.0-B1*(c*c)*cosKappa*kyN*sinOmega*2.0+B2*(c*c)*kxN*kyN*r33*4.0-B1*(c*c)*kxN*sinOmega*sinKappa*6.0-B2*(c*c)*kyN*sinOmega*sinKappa*2.0+A1*c*(r*r)*sinOmega*sinKappa-A1*c*(r0*r0)*sinOmega*sinKappa+A2*c*(r*r*r*r)*sinOmega*sinKappa-A2*c*(r0*r0*r0*r0)*sinOmega*sinKappa+A3*c*(r*r*r*r*r*r)*sinOmega*sinKappa-A3*c*(r0*r0*r0*r0*r0*r0)*sinOmega*sinKappa+B1*(c*c)*cosOmega*cosKappa*kxN*r13*6.0+B2*(c*c)*cosOmega*cosKappa*kyN*r13*2.0-A1*c*cosOmega*cosKappa*(r*r)*r13+A1*c*cosOmega*cosKappa*(r0*r0)*r13-A2*c*cosOmega*cosKappa*(r*r*r*r)*r13+A2*c*cosOmega*cosKappa*(r0*r0*r0*r0)*r13-A3*c*cosOmega*cosKappa*(r*r*r*r*r*r)*r13+A3*c*cosOmega*cosKappa*(r0*r0*r0*r0*r0*r0)*r13+A1*(c*c*c)*cosKappa*kxN*kyN*sinOmega*2.0-B2*(c*c)*cosOmega*kxN*r13*sinKappa*2.0-B1*(c*c)*cosOmega*kyN*r13*sinKappa*2.0-A1*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*r13*2.0-A2*(c*c*c)*kxN*(kyN*kyN)*(r*r)*r33*4.0-A3*(c*c*c)*kxN*(kyN*kyN)*(r*r*r*r)*r33*6.0+A2*(c*c*c)*(kxN*kxN)*(r*r)*sinOmega*sinKappa*4.0+A3*(c*c*c)*(kxN*kxN)*(r*r*r*r)*sinOmega*sinKappa*6.0-A2*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*(r*r)*r13*4.0-A3*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*(r*r*r*r)*r13*6.0+A1*(c*c*c)*cosOmega*kxN*kyN*r13*sinKappa*2.0+A2*(c*c*c)*cosKappa*kxN*kyN*(r*r)*sinOmega*4.0+A3*(c*c*c)*cosKappa*kxN*kyN*(r*r*r*r)*sinOmega*6.0+A2*(c*c*c)*cosOmega*kxN*kyN*(r*r)*r13*sinKappa*4.0+A3*(c*c*c)*cosOmega*kxN*kyN*(r*r*r*r)*r13*sinKappa*6.0)/N));
			A.set(1, column, -((r33*ys+c*cosKappa*sinOmega+c*cosOmega*r13*sinKappa+A1*(r*r)*r33*ys-A1*(r0*r0)*r33*ys+A2*(r*r*r*r)*r33*ys-A2*(r0*r0*r0*r0)*r33*ys+A3*(r*r*r*r*r*r)*r33*ys-A3*(r0*r0*r0*r0*r0*r0)*r33*ys-A1*(c*c*c)*(kyN*kyN*kyN)*r33*2.0+B2*(c*c)*(kxN*kxN)*r33*2.0+B2*(c*c)*(kyN*kyN)*r33*6.0+A1*(c*c*c)*cosKappa*(kyN*kyN)*sinOmega*2.0-A1*(c*c*c)*(kxN*kxN)*kyN*r33*2.0-A2*(c*c*c)*(kyN*kyN*kyN)*(r*r)*r33*4.0-A3*(c*c*c)*(kyN*kyN*kyN)*(r*r*r*r)*r33*6.0-B1*(c*c)*cosKappa*kxN*sinOmega*2.0-B2*(c*c)*cosKappa*kyN*sinOmega*6.0+A1*c*cosKappa*(r*r)*sinOmega-A1*c*cosKappa*(r0*r0)*sinOmega+A2*c*cosKappa*(r*r*r*r)*sinOmega-A2*c*cosKappa*(r0*r0*r0*r0)*sinOmega+A3*c*cosKappa*(r*r*r*r*r*r)*sinOmega-A3*c*cosKappa*(r0*r0*r0*r0*r0*r0)*sinOmega+B1*(c*c)*kxN*kyN*r33*4.0-B2*(c*c)*kxN*sinOmega*sinKappa*2.0-B1*(c*c)*kyN*sinOmega*sinKappa*2.0+B2*(c*c)*cosOmega*cosKappa*kxN*r13*2.0+B1*(c*c)*cosOmega*cosKappa*kyN*r13*2.0-B1*(c*c)*cosOmega*kxN*r13*sinKappa*2.0-B2*(c*c)*cosOmega*kyN*r13*sinKappa*6.0+A1*c*cosOmega*(r*r)*r13*sinKappa-A1*c*cosOmega*(r0*r0)*r13*sinKappa+A2*c*cosOmega*(r*r*r*r)*r13*sinKappa-A2*c*cosOmega*(r0*r0*r0*r0)*r13*sinKappa+A3*c*cosOmega*(r*r*r*r*r*r)*r13*sinKappa-A3*c*cosOmega*(r0*r0*r0*r0*r0*r0)*r13*sinKappa+A1*(c*c*c)*kxN*kyN*sinOmega*sinKappa*2.0+A1*(c*c*c)*cosOmega*(kyN*kyN)*r13*sinKappa*2.0+A2*(c*c*c)*cosKappa*(kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c*c)*cosKappa*(kyN*kyN)*(r*r*r*r)*sinOmega*6.0-A2*(c*c*c)*(kxN*kxN)*kyN*(r*r)*r33*4.0-A3*(c*c*c)*(kxN*kxN)*kyN*(r*r*r*r)*r33*6.0+A2*(c*c*c)*cosOmega*(kyN*kyN)*(r*r)*r13*sinKappa*4.0+A3*(c*c*c)*cosOmega*(kyN*kyN)*(r*r*r*r)*r13*sinKappa*6.0-A1*(c*c*c)*cosOmega*cosKappa*kxN*kyN*r13*2.0+A2*(c*c*c)*kxN*kyN*(r*r)*sinOmega*sinKappa*4.0+A3*(c*c*c)*kxN*kyN*(r*r*r*r)*sinOmega*sinKappa*6.0-A2*(c*c*c)*cosOmega*cosKappa*kxN*kyN*(r*r)*r13*4.0-A3*(c*c*c)*cosOmega*cosKappa*kxN*kyN*(r*r*r*r)*r13*6.0)/N));
		}
		
	
		// interior orientation of the camera
		column = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, 1.0);
			A.set(1, column, 0.0);
		}
		
		column = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, 0.0);
			A.set(1, column, 1.0);
		}
		
		column = interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -kxN-C1*kxN-C2*kyN+B1*c*((kxN*kxN)*3.0+kyN*kyN)*2.0-kxN*(r*r-r0*r0)*(A1+A2*(r*r)+A2*(r0*r0)+A3*(r*r*r*r)+A3*(r0*r0*r0*r0)+A3*(r*r)*(r0*r0))+B2*c*kxN*kyN*4.0+c*1.0/(r*r)*xs*(kxN*kxN+kyN*kyN)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))+(c*c)*kxN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)));
			A.set(1, column, -kyN+B2*c*(kxN*kxN+(kyN*kyN)*3.0)*2.0-kyN*(r*r-r0*r0)*(A1+A2*(r*r)+A2*(r0*r0)+A3*(r*r*r*r)+A3*(r0*r0*r0*r0)+A3*(r*r)*(r0*r0))+B1*c*kxN*kyN*4.0+c*1.0/(r*r)*ys*(kxN*kxN+kyN*kyN)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))+(c*c)*kyN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)));
		}
		
		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs*(r*r - r0*r0));
			A.set(1, column, ys*(r*r - r0*r0));
		}
		
		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs*(r*r*r*r - r0*r0*r0*r0));
			A.set(1, column, ys*(r*r*r*r - r0*r0*r0*r0));
		}
		
		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs*(r*r*r*r*r*r - r0*r0*r0*r0*r0*r0));	
			A.set(1, column, ys*(r*r*r*r*r*r - r0*r0*r0*r0*r0*r0));
		}
		
		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, r*r + 2.0*xs*xs); // c*c*(3.0*kxN*kxN + kyN*kyN);
			A.set(1, column, 2.0*xs*ys); 
		}
		
		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, 2.0*xs*ys);
			A.set(1, column, r*r + 2.0*ys*ys); //c*c*(3.0*kyN*kyN + kxN*kxN);
		}
		
		column = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs);
			A.set(1, column, 0.0);
		}

		column = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, ys);
			A.set(1, column, 0.0);
		}
				
		// exterior orientation of the image
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (c*r11)/N+(r13*xs)/N+(C1*c*r11)/N+(C2*c*r12)/N+(C1*r13*xs)/N+(C2*r13*ys)/N+(B1*(c*c)*(kxN*r11*-3.0+(kxN*kxN)*r13*3.0+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*2.0)/N+(c*r11*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)+(r13*xs*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)-(B2*(c*c)*kyN*r11*2.0)/N+(B2*(c*c)*cosPhi*kxN*sinKappa*2.0)/N+(B2*(c*c)*kxN*kyN*r13*4.0)/N+((c*c)*1.0/(r*r)*xs*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N+((c*c*c)*kxN*1.0/(r*r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/N);
			A.set(1, column, (c*r12)/N+(r13*ys)/N+(B2*(c*c)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13*3.0+cosPhi*kyN*sinKappa*3.0)*2.0)/N+(c*r12*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)+(r13*ys*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)-(B1*(c*c)*kyN*r11*2.0)/N+(B1*(c*c)*cosPhi*kxN*sinKappa*2.0)/N+(B1*(c*c)*kxN*kyN*r13*4.0)/N+((c*c)*1.0/(r*r)*ys*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N+((c*c*c)*kyN*1.0/(r*r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/N);
		}
		
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (c*(cosOmega*sinKappa+C2*cosOmega*cosKappa+C1*cosOmega*sinKappa+cosPhi*kxN*sinOmega+cosKappa*r13*sinOmega+A1*cosOmega*(r*r)*sinKappa-A1*cosOmega*(r0*r0)*sinKappa+A2*cosOmega*(r*r*r*r)*sinKappa-A2*cosOmega*(r0*r0*r0*r0)*sinKappa+A3*cosOmega*(r*r*r*r*r*r)*sinKappa-A3*cosOmega*(r0*r0*r0*r0*r0*r0)*sinKappa+C1*cosPhi*kxN*sinOmega+C2*cosPhi*kyN*sinOmega+C1*cosKappa*r13*sinOmega-C2*r13*sinOmega*sinKappa+A1*(c*c)*cosPhi*(kxN*kxN*kxN)*sinOmega*2.0+A1*(c*c)*cosOmega*(kxN*kxN)*sinKappa*2.0-B2*c*cosOmega*cosKappa*kxN*2.0-B1*c*cosOmega*cosKappa*kyN*2.0-B1*c*cosOmega*kxN*sinKappa*6.0-B2*c*cosOmega*kyN*sinKappa*2.0-B1*c*cosPhi*(kxN*kxN)*sinOmega*6.0-B1*c*cosPhi*(kyN*kyN)*sinOmega*2.0+A1*cosPhi*kxN*(r*r)*sinOmega-A1*cosPhi*kxN*(r0*r0)*sinOmega+A2*cosPhi*kxN*(r*r*r*r)*sinOmega-A2*cosPhi*kxN*(r0*r0*r0*r0)*sinOmega+A3*cosPhi*kxN*(r*r*r*r*r*r)*sinOmega-A3*cosPhi*kxN*(r0*r0*r0*r0*r0*r0)*sinOmega+A1*cosKappa*(r*r)*r13*sinOmega-A1*cosKappa*(r0*r0)*r13*sinOmega+A2*cosKappa*(r*r*r*r)*r13*sinOmega-A2*cosKappa*(r0*r0*r0*r0)*r13*sinOmega+A3*cosKappa*(r*r*r*r*r*r)*r13*sinOmega-A3*cosKappa*(r0*r0*r0*r0*r0*r0)*r13*sinOmega+A1*(c*c)*cosOmega*cosKappa*kxN*kyN*2.0+A1*(c*c)*cosPhi*kxN*(kyN*kyN)*sinOmega*2.0+A1*(c*c)*cosKappa*(kxN*kxN)*r13*sinOmega*2.0-B2*c*cosPhi*kxN*kyN*sinOmega*4.0-B1*c*cosKappa*kxN*r13*sinOmega*6.0-B2*c*cosKappa*kyN*r13*sinOmega*2.0+B2*c*kxN*r13*sinOmega*sinKappa*2.0+B1*c*kyN*r13*sinOmega*sinKappa*2.0+A2*(c*c)*cosPhi*(kxN*kxN*kxN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kxN*kxN*kxN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosOmega*(kxN*kxN)*(r*r)*sinKappa*4.0+A3*(c*c)*cosOmega*(kxN*kxN)*(r*r*r*r)*sinKappa*6.0+A2*(c*c)*cosPhi*kxN*(kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*kxN*(kyN*kyN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosKappa*(kxN*kxN)*(r*r)*r13*sinOmega*4.0+A3*(c*c)*cosKappa*(kxN*kxN)*(r*r*r*r)*r13*sinOmega*6.0-A1*(c*c)*kxN*kyN*r13*sinOmega*sinKappa*2.0+A2*(c*c)*cosOmega*cosKappa*kxN*kyN*(r*r)*4.0+A3*(c*c)*cosOmega*cosKappa*kxN*kyN*(r*r*r*r)*6.0-A2*(c*c)*kxN*kyN*(r*r)*r13*sinOmega*sinKappa*4.0-A3*(c*c)*kxN*kyN*(r*r*r*r)*r13*sinOmega*sinKappa*6.0))/N);
			A.set(1, column, (c*(cosOmega*cosKappa+cosPhi*kyN*sinOmega-r13*sinOmega*sinKappa+A1*cosOmega*cosKappa*(r*r)-A1*cosOmega*cosKappa*(r0*r0)+A2*cosOmega*cosKappa*(r*r*r*r)-A2*cosOmega*cosKappa*(r0*r0*r0*r0)+A3*cosOmega*cosKappa*(r*r*r*r*r*r)-A3*cosOmega*cosKappa*(r0*r0*r0*r0*r0*r0)+A1*(c*c)*cosPhi*(kyN*kyN*kyN)*sinOmega*2.0-B1*c*cosOmega*cosKappa*kxN*2.0-B2*c*cosOmega*cosKappa*kyN*6.0-B2*c*cosOmega*kxN*sinKappa*2.0-B1*c*cosOmega*kyN*sinKappa*2.0-B2*c*cosPhi*(kxN*kxN)*sinOmega*2.0-B2*c*cosPhi*(kyN*kyN)*sinOmega*6.0+A1*cosPhi*kyN*(r*r)*sinOmega-A1*cosPhi*kyN*(r0*r0)*sinOmega+A2*cosPhi*kyN*(r*r*r*r)*sinOmega-A2*cosPhi*kyN*(r0*r0*r0*r0)*sinOmega+A3*cosPhi*kyN*(r*r*r*r*r*r)*sinOmega-A3*cosPhi*kyN*(r0*r0*r0*r0*r0*r0)*sinOmega-A1*(r*r)*r13*sinOmega*sinKappa+A1*(r0*r0)*r13*sinOmega*sinKappa-A2*(r*r*r*r)*r13*sinOmega*sinKappa+A2*(r0*r0*r0*r0)*r13*sinOmega*sinKappa-A3*(r*r*r*r*r*r)*r13*sinOmega*sinKappa+A3*(r0*r0*r0*r0*r0*r0)*r13*sinOmega*sinKappa+A1*(c*c)*cosOmega*cosKappa*(kyN*kyN)*2.0+A1*(c*c)*cosOmega*kxN*kyN*sinKappa*2.0+A1*(c*c)*cosPhi*(kxN*kxN)*kyN*sinOmega*2.0-A1*(c*c)*(kyN*kyN)*r13*sinOmega*sinKappa*2.0-B1*c*cosPhi*kxN*kyN*sinOmega*4.0-B2*c*cosKappa*kxN*r13*sinOmega*2.0-B1*c*cosKappa*kyN*r13*sinOmega*2.0+A2*(c*c)*cosOmega*cosKappa*(kyN*kyN)*(r*r)*4.0+A3*(c*c)*cosOmega*cosKappa*(kyN*kyN)*(r*r*r*r)*6.0+B1*c*kxN*r13*sinOmega*sinKappa*2.0+B2*c*kyN*r13*sinOmega*sinKappa*6.0+A2*(c*c)*cosPhi*(kyN*kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kyN*kyN*kyN)*(r*r*r*r)*sinOmega*6.0+A2*(c*c)*cosPhi*(kxN*kxN)*kyN*(r*r)*sinOmega*4.0+A3*(c*c)*cosPhi*(kxN*kxN)*kyN*(r*r*r*r)*sinOmega*6.0-A2*(c*c)*(kyN*kyN)*(r*r)*r13*sinOmega*sinKappa*4.0-A3*(c*c)*(kyN*kyN)*(r*r*r*r)*r13*sinOmega*sinKappa*6.0+A1*(c*c)*cosKappa*kxN*kyN*r13*sinOmega*2.0+A2*(c*c)*cosOmega*kxN*kyN*(r*r)*sinKappa*4.0+A3*(c*c)*cosOmega*kxN*kyN*(r*r*r*r)*sinKappa*6.0+A2*(c*c)*cosKappa*kxN*kyN*(r*r)*r13*sinOmega*4.0+A3*(c*c)*cosKappa*kxN*kyN*(r*r*r*r)*r13*sinOmega*6.0))/N);
		}
		
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (r33*xs+C1*r33*xs+C2*r33*ys+c*sinOmega*sinKappa+A1*(r*r)*r33*xs-A1*(r0*r0)*r33*xs+A2*(r*r*r*r)*r33*xs-A2*(r0*r0*r0*r0)*r33*xs+A3*(r*r*r*r*r*r)*r33*xs-A3*(r0*r0*r0*r0*r0*r0)*r33*xs-A1*(c*c*c)*(kxN*kxN*kxN)*r33*2.0+B1*(c*c)*(kxN*kxN)*r33*6.0+B1*(c*c)*(kyN*kyN)*r33*2.0+C2*c*cosKappa*sinOmega+C1*c*sinOmega*sinKappa-c*cosOmega*cosKappa*r13-A1*(c*c*c)*kxN*(kyN*kyN)*r33*2.0+A1*(c*c*c)*(kxN*kxN)*sinOmega*sinKappa*2.0-C1*c*cosOmega*cosKappa*r13+C2*c*cosOmega*r13*sinKappa-A2*(c*c*c)*(kxN*kxN*kxN)*(r*r)*r33*4.0-A3*(c*c*c)*(kxN*kxN*kxN)*(r*r*r*r)*r33*6.0-B2*(c*c)*cosKappa*kxN*sinOmega*2.0-B1*(c*c)*cosKappa*kyN*sinOmega*2.0+B2*(c*c)*kxN*kyN*r33*4.0-B1*(c*c)*kxN*sinOmega*sinKappa*6.0-B2*(c*c)*kyN*sinOmega*sinKappa*2.0+A1*c*(r*r)*sinOmega*sinKappa-A1*c*(r0*r0)*sinOmega*sinKappa+A2*c*(r*r*r*r)*sinOmega*sinKappa-A2*c*(r0*r0*r0*r0)*sinOmega*sinKappa+A3*c*(r*r*r*r*r*r)*sinOmega*sinKappa-A3*c*(r0*r0*r0*r0*r0*r0)*sinOmega*sinKappa+B1*(c*c)*cosOmega*cosKappa*kxN*r13*6.0+B2*(c*c)*cosOmega*cosKappa*kyN*r13*2.0-A1*c*cosOmega*cosKappa*(r*r)*r13+A1*c*cosOmega*cosKappa*(r0*r0)*r13-A2*c*cosOmega*cosKappa*(r*r*r*r)*r13+A2*c*cosOmega*cosKappa*(r0*r0*r0*r0)*r13-A3*c*cosOmega*cosKappa*(r*r*r*r*r*r)*r13+A3*c*cosOmega*cosKappa*(r0*r0*r0*r0*r0*r0)*r13+A1*(c*c*c)*cosKappa*kxN*kyN*sinOmega*2.0-B2*(c*c)*cosOmega*kxN*r13*sinKappa*2.0-B1*(c*c)*cosOmega*kyN*r13*sinKappa*2.0-A1*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*r13*2.0-A2*(c*c*c)*kxN*(kyN*kyN)*(r*r)*r33*4.0-A3*(c*c*c)*kxN*(kyN*kyN)*(r*r*r*r)*r33*6.0+A2*(c*c*c)*(kxN*kxN)*(r*r)*sinOmega*sinKappa*4.0+A3*(c*c*c)*(kxN*kxN)*(r*r*r*r)*sinOmega*sinKappa*6.0-A2*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*(r*r)*r13*4.0-A3*(c*c*c)*cosOmega*cosKappa*(kxN*kxN)*(r*r*r*r)*r13*6.0+A1*(c*c*c)*cosOmega*kxN*kyN*r13*sinKappa*2.0+A2*(c*c*c)*cosKappa*kxN*kyN*(r*r)*sinOmega*4.0+A3*(c*c*c)*cosKappa*kxN*kyN*(r*r*r*r)*sinOmega*6.0+A2*(c*c*c)*cosOmega*kxN*kyN*(r*r)*r13*sinKappa*4.0+A3*(c*c*c)*cosOmega*kxN*kyN*(r*r*r*r)*r13*sinKappa*6.0)/N);
			A.set(1, column, (r33*ys+c*cosKappa*sinOmega+c*cosOmega*r13*sinKappa+A1*(r*r)*r33*ys-A1*(r0*r0)*r33*ys+A2*(r*r*r*r)*r33*ys-A2*(r0*r0*r0*r0)*r33*ys+A3*(r*r*r*r*r*r)*r33*ys-A3*(r0*r0*r0*r0*r0*r0)*r33*ys-A1*(c*c*c)*(kyN*kyN*kyN)*r33*2.0+B2*(c*c)*(kxN*kxN)*r33*2.0+B2*(c*c)*(kyN*kyN)*r33*6.0+A1*(c*c*c)*cosKappa*(kyN*kyN)*sinOmega*2.0-A1*(c*c*c)*(kxN*kxN)*kyN*r33*2.0-A2*(c*c*c)*(kyN*kyN*kyN)*(r*r)*r33*4.0-A3*(c*c*c)*(kyN*kyN*kyN)*(r*r*r*r)*r33*6.0-B1*(c*c)*cosKappa*kxN*sinOmega*2.0-B2*(c*c)*cosKappa*kyN*sinOmega*6.0+A1*c*cosKappa*(r*r)*sinOmega-A1*c*cosKappa*(r0*r0)*sinOmega+A2*c*cosKappa*(r*r*r*r)*sinOmega-A2*c*cosKappa*(r0*r0*r0*r0)*sinOmega+A3*c*cosKappa*(r*r*r*r*r*r)*sinOmega-A3*c*cosKappa*(r0*r0*r0*r0*r0*r0)*sinOmega+B1*(c*c)*kxN*kyN*r33*4.0-B2*(c*c)*kxN*sinOmega*sinKappa*2.0-B1*(c*c)*kyN*sinOmega*sinKappa*2.0+B2*(c*c)*cosOmega*cosKappa*kxN*r13*2.0+B1*(c*c)*cosOmega*cosKappa*kyN*r13*2.0-B1*(c*c)*cosOmega*kxN*r13*sinKappa*2.0-B2*(c*c)*cosOmega*kyN*r13*sinKappa*6.0+A1*c*cosOmega*(r*r)*r13*sinKappa-A1*c*cosOmega*(r0*r0)*r13*sinKappa+A2*c*cosOmega*(r*r*r*r)*r13*sinKappa-A2*c*cosOmega*(r0*r0*r0*r0)*r13*sinKappa+A3*c*cosOmega*(r*r*r*r*r*r)*r13*sinKappa-A3*c*cosOmega*(r0*r0*r0*r0*r0*r0)*r13*sinKappa+A1*(c*c*c)*kxN*kyN*sinOmega*sinKappa*2.0+A1*(c*c*c)*cosOmega*(kyN*kyN)*r13*sinKappa*2.0+A2*(c*c*c)*cosKappa*(kyN*kyN)*(r*r)*sinOmega*4.0+A3*(c*c*c)*cosKappa*(kyN*kyN)*(r*r*r*r)*sinOmega*6.0-A2*(c*c*c)*(kxN*kxN)*kyN*(r*r)*r33*4.0-A3*(c*c*c)*(kxN*kxN)*kyN*(r*r*r*r)*r33*6.0+A2*(c*c*c)*cosOmega*(kyN*kyN)*(r*r)*r13*sinKappa*4.0+A3*(c*c*c)*cosOmega*(kyN*kyN)*(r*r*r*r)*r13*sinKappa*6.0-A1*(c*c*c)*cosOmega*cosKappa*kxN*kyN*r13*2.0+A2*(c*c*c)*kxN*kyN*(r*r)*sinOmega*sinKappa*4.0+A3*(c*c*c)*kxN*kyN*(r*r*r*r)*sinOmega*sinKappa*6.0-A2*(c*c*c)*cosOmega*cosKappa*kxN*kyN*(r*r)*r13*4.0-A3*(c*c*c)*cosOmega*cosKappa*kxN*kyN*(r*r*r*r)*r13*6.0)/N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, B1*(((c*c)*(kxN*kxN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*6.0)/N+((c*c)*(kyN*kyN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N-((c*c)*kxN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*6.0)/N-((c*c)*kyN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*2.0)/N)+(xs*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/N+(c*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0)))/N+(C1*c*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0)))/N+(C2*c*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0)))/N+(C1*xs*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/N+(C2*ys*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/N+(c*dRad*kxN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N-((c*c)*kxN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*2.0)/N-((c*c)*kyN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*2.0)/N))/2.0+(c*dRad*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0)))/(N*r)+(dRad*xs*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/(N*r)-(B2*(c*c)*kxN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*2.0)/N-(B2*(c*c)*kyN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*2.0)/N+((c*c)*1.0/(r*r)*xs*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(Y*(kxN*kxN)*r33-Y0*(kxN*kxN)*r33+Y*(kyN*kyN)*r33-Y0*(kyN*kyN)*r33-Y*kxN*sinOmega*sinKappa+Y0*kxN*sinOmega*sinKappa+Z*cosPhi*(kxN*kxN)*sinOmega-Z0*cosPhi*(kxN*kxN)*sinOmega+Z*cosPhi*(kyN*kyN)*sinOmega-Z0*cosPhi*(kyN*kyN)*sinOmega+Z*cosOmega*cosKappa*kyN-Z0*cosOmega*cosKappa*kyN-Y*cosKappa*kyN*sinOmega+Y0*cosKappa*kyN*sinOmega+Z*cosOmega*kxN*sinKappa-Z0*cosOmega*kxN*sinKappa+Y*cosOmega*cosKappa*kxN*r13-Y0*cosOmega*cosKappa*kxN*r13-Y*cosOmega*kyN*r13*sinKappa+Y0*cosOmega*kyN*r13*sinKappa+Z*cosKappa*kxN*r13*sinOmega-Z0*cosKappa*kxN*r13*sinOmega-Z*kyN*r13*sinOmega*sinKappa+Z0*kyN*r13*sinOmega*sinKappa))/N+(B2*(c*c)*kxN*kyN*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*4.0)/N);
			A.set(1, column, B2*(((c*c)*(kxN*kxN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*6.0)/N-((c*c)*kxN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*2.0)/N-((c*c)*kyN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*6.0)/N)+(ys*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/N+(c*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0)))/N+(c*dRad*kyN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*2.0)/N-((c*c)*kxN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*2.0)/N-((c*c)*kyN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*2.0)/N))/2.0+(c*dRad*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0)))/(N*r)+(dRad*ys*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega))/(N*r)-(B1*(c*c)*kxN*((cosKappa*sinOmega+cosOmega*r13*sinKappa)*(Y-Y0)-(cosOmega*cosKappa-r13*sinOmega*sinKappa)*(Z-Z0))*2.0)/N-(B1*(c*c)*kyN*((sinOmega*sinKappa-cosOmega*cosKappa*r13)*(Y-Y0)-(cosOmega*sinKappa+cosKappa*r13*sinOmega)*(Z-Z0))*2.0)/N+((c*c)*1.0/(r*r)*ys*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(Y*(kxN*kxN)*r33-Y0*(kxN*kxN)*r33+Y*(kyN*kyN)*r33-Y0*(kyN*kyN)*r33-Y*kxN*sinOmega*sinKappa+Y0*kxN*sinOmega*sinKappa+Z*cosPhi*(kxN*kxN)*sinOmega-Z0*cosPhi*(kxN*kxN)*sinOmega+Z*cosPhi*(kyN*kyN)*sinOmega-Z0*cosPhi*(kyN*kyN)*sinOmega+Z*cosOmega*cosKappa*kyN-Z0*cosOmega*cosKappa*kyN-Y*cosKappa*kyN*sinOmega+Y0*cosKappa*kyN*sinOmega+Z*cosOmega*kxN*sinKappa-Z0*cosOmega*kxN*sinKappa+Y*cosOmega*cosKappa*kxN*r13-Y0*cosOmega*cosKappa*kxN*r13-Y*cosOmega*kyN*r13*sinKappa+Y0*cosOmega*kyN*r13*sinKappa+Z*cosKappa*kxN*r13*sinOmega-Z0*cosKappa*kxN*r13*sinOmega-Z*kyN*r13*sinOmega*sinKappa+Z0*kyN*r13*sinOmega*sinKappa))/N+(B1*(c*c)*kxN*kyN*(Y*r33-Y0*r33+Z*cosPhi*sinOmega-Z0*cosPhi*sinOmega)*4.0)/N);		}
		
		column = exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column,  c*cosKappa-B1*(((c*c)*kxN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*6.0)/N-((c*c)*kyN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*2.0)/N+((c*c)*(kxN*kxN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*6.0)/N+((c*c)*(kyN*kyN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N)+C1*c*cosKappa-C2*c*sinKappa+(c*kxN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/N+(C1*c*kxN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/N+(C2*c*kyN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/N-(c*kxN*1.0/(r*r*r)*(((c*c)*kxN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*2.0)/N-((c*c)*kyN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*2.0)/N+((c*c)*(kxN*kxN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/2.0-(B2*(c*c)*kyN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*2.0)/N+(c*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0))*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega))/(N*r)+(B2*(c*c)*kxN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*2.0)/N+((c*c*c)*kxN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(X*cosPhi*(kxN*kxN)-X0*cosPhi*(kxN*kxN)+X*cosPhi*(kyN*kyN)-X0*cosPhi*(kyN*kyN)-X*kyN*r13*sinKappa+X0*kyN*r13*sinKappa-Y*kxN*r11*sinOmega+Y0*kxN*r11*sinOmega-Y*kyN*r12*sinOmega+Y0*kyN*r12*sinOmega-Z*kyN*r33*sinKappa+Z0*kyN*r33*sinKappa-Z*cosOmega*(kxN*kxN)*r13+Z0*cosOmega*(kxN*kxN)*r13-Z*cosOmega*(kyN*kyN)*r13+Z0*cosOmega*(kyN*kyN)*r13+Y*(kxN*kxN)*r13*sinOmega-Y0*(kxN*kxN)*r13*sinOmega+Y*(kyN*kyN)*r13*sinOmega-Y0*(kyN*kyN)*r13*sinOmega+X*cosKappa*kxN*r13-X0*cosKappa*kxN*r13+Z*cosOmega*kxN*r11-Z0*cosOmega*kxN*r11))/N-(B2*(c*c)*kxN*kyN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*4.0)/N+(c*kxN*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0))*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/(N*r));
			A.set(1, column, -c*sinKappa-B2*(((c*c)*kxN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*2.0)/N-((c*c)*kyN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*6.0)/N+((c*c)*(kxN*kxN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*6.0)/N)+(c*kyN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/N-(c*kyN*1.0/(r*r*r)*(((c*c)*kxN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*2.0)/N-((c*c)*kyN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*2.0)/N+((c*c)*(kxN*kxN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*2.0)/N)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/2.0-(B1*(c*c)*kyN*(X*cosKappa*r13-X0*cosKappa*r13+Z*cosOmega*r11-Z0*cosOmega*r11-Y*r11*sinOmega+Y0*r11*sinOmega)*2.0)/N+(B1*(c*c)*kxN*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*2.0)/N-(c*(X*r13*sinKappa-X0*r13*sinKappa+Y*r12*sinOmega-Y0*r12*sinOmega+Z*r33*sinKappa-Z0*r33*sinKappa)*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0)))/(N*r)+((c*c*c)*kyN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(X*cosPhi*(kxN*kxN)-X0*cosPhi*(kxN*kxN)+X*cosPhi*(kyN*kyN)-X0*cosPhi*(kyN*kyN)-X*kyN*r13*sinKappa+X0*kyN*r13*sinKappa-Y*kxN*r11*sinOmega+Y0*kxN*r11*sinOmega-Y*kyN*r12*sinOmega+Y0*kyN*r12*sinOmega-Z*kyN*r33*sinKappa+Z0*kyN*r33*sinKappa-Z*cosOmega*(kxN*kxN)*r13+Z0*cosOmega*(kxN*kxN)*r13-Z*cosOmega*(kyN*kyN)*r13+Z0*cosOmega*(kyN*kyN)*r13+Y*(kxN*kxN)*r13*sinOmega-Y0*(kxN*kxN)*r13*sinOmega+Y*(kyN*kyN)*r13*sinOmega-Y0*(kyN*kyN)*r13*sinOmega+X*cosKappa*kxN*r13-X0*cosKappa*kxN*r13+Z*cosOmega*kxN*r11-Z0*cosOmega*kxN*r11))/N-(B1*(c*c)*kxN*kyN*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega)*4.0)/N+(c*kyN*(A1*(r*r*r)+A2*(r*r*r*r*r)+A3*(r*r*r*r*r*r*r)-A1*r*(r0*r0)-A2*r*(r0*r0*r0*r0)-A3*r*(r0*r0*r0*r0*r0*r0))*(X*cosPhi-X0*cosPhi-Z*cosOmega*r13+Z0*cosOmega*r13+Y*r13*sinOmega-Y0*r13*sinOmega))/(N*r));
		}
			
		column = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, ys+C1*ys-B2*(c*c)*(kxN*kxN)*2.0+B2*(c*c)*(kyN*kyN)*2.0+C2*c*kxN+A1*(r*r)*ys-A1*(r0*r0)*ys+A2*(r*r*r*r)*ys-A2*(r0*r0*r0*r0)*ys+A3*(r*r*r*r*r*r)*ys-A3*(r0*r0*r0*r0*r0*r0)*ys+B1*(c*c)*kxN*kyN*4.0);
			A.set(1, column, c*(kxN-B1*c*(kxN*kxN)*2.0+B1*c*(kyN*kyN)*2.0+A1*kxN*(r*r)-A1*kxN*(r0*r0)+A2*kxN*(r*r*r*r)-A2*kxN*(r0*r0*r0*r0)+A3*kxN*(r*r*r*r*r*r)-A3*kxN*(r0*r0*r0*r0*r0*r0)-B2*c*kxN*kyN*4.0));
		}
		
		Collections.sort(columns);
		for (int row = 0; row < A.numRows(); row++) {
			for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
				int colAT = columns.get(columnATIdx);
				double aT = A.get(row, colAT);
				neq.add(colAT, aT * P.get(row, row) * w.get(row));
				
				for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
					int colA = columns.get(columnAIdx);	
					double a  = A.get(row, colA);
					NEQ.add(colAT, colA, aT * P.get(row, row) * a);
				}
			}	
		}
	}

}
