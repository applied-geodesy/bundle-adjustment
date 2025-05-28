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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmBandMatrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;

class PartialDerivativeFactory {
	private PartialDerivativeFactory() {}
	
	static double getWeightedSumOfSquaredResiduals(double sigma2apriori, ObservationParameterGroup<?> observations) {
		double omega = 0.0;
		if (observations instanceof ImageCoordinate) 
			omega = getWeightedSumOfSquaredResiduals(sigma2apriori, (ImageCoordinate)observations);
		else if (observations instanceof ScaleBar)
			omega = getWeightedSumOfSquaredResiduals(sigma2apriori, (ScaleBar)observations);
		return omega;
	}
	
	private static double getWeightedSumOfSquaredResiduals(double sigma2apriori, ScaleBar scaleBar) {
		Map<ParameterType, Double> residuals = PartialDerivativeFactory.getMisclosures(scaleBar);
		double residuum = residuals.get(ParameterType.SCALE_BAR_LENGTH);
		return residuum * residuum * sigma2apriori / scaleBar.getLength().getVariance();
	}
	
	private static double getWeightedSumOfSquaredResiduals(double sigma2apriori, ImageCoordinate imageCoordinate) {
		Map<ParameterType, Double> residuals = PartialDerivativeFactory.getMisclosures(imageCoordinate);
		double residuumX = residuals.get(ParameterType.IMAGE_COORDINATE_X);
		double residuumY = residuals.get(ParameterType.IMAGE_COORDINATE_Y);

		double varianceX  = imageCoordinate.getX().getVariance();
		double varianceY  = imageCoordinate.getY().getVariance();
		double corrCoefXY = imageCoordinate.getCorrelationCoefficient();

		if (corrCoefXY == 0)
			return residuumX * residuumX * sigma2apriori/varianceX + residuumY * residuumY * sigma2apriori/varianceY;
		else {
			double invDet = sigma2apriori / ((corrCoefXY + 1.0) * varianceX * varianceY);
			double qxx =  invDet * varianceY;
			double qyy =  invDet * varianceX;
			double qxy = -invDet * corrCoefXY * Math.sqrt(varianceX * varianceY);
			
			return residuumX * (qxx * residuumX + qxy * residuumY) + residuumY * (qxy * residuumX + qyy * residuumY);
		}
	}

	static Map<ParameterType, Double> getMisclosures(ObservationParameterGroup<?> observations) {
		Map<ParameterType, Double> values = Collections.emptyMap();
		
		if (observations instanceof ImageCoordinate) 
			values = getCollinearityEquationValues((ImageCoordinate)observations);			
		else if (observations instanceof ScaleBar)
			values = getDistanceValue((ScaleBar)observations);
		
		Map<ParameterType, Double> misclosures = new HashMap<ParameterType, Double>(observations.getNumberOfParameters());
		for (ObservationParameter<?> observation : observations) {
			ParameterType parameterType = observation.getParameterType();
			double value = values.get(parameterType);
			// estimate residuals (== observed - calculated)
			misclosures.put(parameterType, observation.getValue() - value);
		}
		
		return misclosures;
	}

	private static Map<ParameterType, Double> getDistanceValue(ScaleBar scaleBar) {
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

		return Map.ofEntries( 
				Map.entry(ParameterType.SCALE_BAR_LENGTH, Math.sqrt( dX*dX + dY*dY + dZ*dZ ))
		);
	}

	private static Map<ParameterType, Double> getCollinearityEquationValues(ImageCoordinate imageCoordinate) {
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
		
		double D1 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).getValue();
		double D2 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).getValue();
		double D3 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).getValue();

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

		// Rotation (Luhmann (2023) Eq. 2.31, p. 62)
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
		double dRad  = A1*r*r*r + A2*r*r*r*r*r + A3*r*r*r*r*r*r*r - (A1*r0*r0 + A2*r0*r0*r0*r0 + A3*r0*r0*r0*r0*r0*r0)*r;

		double dRadX = xs * dRad/r;
		double dRadY = ys * dRad/r;

		double dTanX = B1 * (r*r + 2.0*xs*xs) + 2.0 * B2 * xs * ys;
		double dTanY = B2 * (r*r + 2.0*ys*ys) + 2.0 * B1 * xs * ys;

		double dAffX = C1*xs + C2*ys;
		double dAffY = 0;

		double dDist = 1.0/N * (D1*r*r*r + D2*r*r*r*r*r + D3*r*r*r*r*r*r*r - (D1*r0*r0 + D2*r0*r0*r0*r0 + D3*r0*r0*r0*r0*r0*r0)*r);
		double dDistX = xs*dDist/r;
		double dDistY = ys*dDist/r;
		
		double deltaX = dRadX + dTanX + dAffX + dDistX;
		double deltaY = dRadY + dTanY + dAffY + dDistY;
		
		return Map.ofEntries( 
				Map.entry(ParameterType.IMAGE_COORDINATE_X, x0 + xs + deltaX),
				Map.entry(ParameterType.IMAGE_COORDINATE_Y, y0 + ys + deltaY) 
		);
	}
	
	static GaussMarkovEquations getPartialDerivative(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ObservationParameterGroup<?> observations) {
		if (observations instanceof ImageCoordinate) 
			return getPartialDerivativeImageCoordinate(sigma2apriori, NEQ, neq, (ImageCoordinate)observations);
		else if (observations instanceof ScaleBar)
			return getPartialDerivativeScaleBar(sigma2apriori, NEQ, neq, (ScaleBar)observations);
		
		return null;
	}
	
	private static GaussMarkovEquations getPartialDerivativeScaleBar(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ScaleBar scaleBar) {
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
		
		int numberOfRows = scaleBar.getNumberOfParameters();

		DenseMatrix A = new DenseMatrix(numberOfRows, neq.size());
		DenseVector w = new DenseVector(numberOfRows);
		UpperSymmBandMatrix P = new UpperSymmBandMatrix(numberOfRows, 0);

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
		for (int row = 0; row < numberOfRows; row++) {
			for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
				int colAT = columns.get(columnATIdx);
				double aT = A.get(row, colAT);
				neq.add(colAT, aT * P.get(row, row) * w.get(row));
				if (NEQ != null) {
					for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
						int colA = columns.get(columnAIdx);	
						double a = A.get(row, colA);
						NEQ.add(colAT, colA, aT * P.get(row, row) * a);
					}
				}
			}	
		}
		return new GaussMarkovEquations(A, P, w);
	}

	private static GaussMarkovEquations getPartialDerivativeImageCoordinate(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ImageCoordinate imageCoordinate) {
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
		
		double D1 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).getValue();
		double D2 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).getValue();
		double D3 = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).getValue();
		
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

		// Rotation (Luhmann (2023) Eq. 2.31, p. 62)
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
		double dRad = A1*r*r*r + A2*r*r*r*r*r + A3*r*r*r*r*r*r*r - (A1*r0*r0 + A2*r0*r0*r0*r0 + A3*r0*r0*r0*r0*r0*r0)*r;

		double dRadX = xs * dRad/r;
		double dRadY = ys * dRad/r;

		double dTanX = B1 * (r*r + 2.0*xs*xs) + 2.0 * B2 * xs * ys;
		double dTanY = B2 * (r*r + 2.0*ys*ys) + 2.0 * B1 * xs * ys;

		double dAffX = C1*xs + C2*ys;
		double dAffY = 0;

		double dDist = 1.0/N * (D1*r*r*r + D2*r*r*r*r*r + D3*r*r*r*r*r*r*r - (D1*r0*r0 + D2*r0*r0*r0*r0 + D3*r0*r0*r0*r0*r0*r0)*r);
		double dDistX = xs*dDist/r;
		double dDistY = ys*dDist/r;
		
		double deltaX = dRadX + dTanX + dAffX + dDistX;
		double deltaY = dRadY + dTanY + dAffY + dDistY;

		double varianceX  = imageCoordinate.getX().getVariance();
		double varianceY  = imageCoordinate.getY().getVariance();
		double corrCoefXY = imageCoordinate.getCorrelationCoefficient();
		
		int numberOfRows = imageCoordinate.getNumberOfParameters();

		// 3 ... object point, 13 ... interior orientation, 6 ... exterior orientation, datum defect
		DenseVector w = new DenseVector(numberOfRows);
		DenseMatrix A = new DenseMatrix(numberOfRows, neq.size());

		Matrix P = null;
		if (corrCoefXY == 0) {
			P = new UpperSymmBandMatrix(numberOfRows, 0);
			P.set(0, 0, sigma2apriori/varianceX);
			P.set(1, 1, sigma2apriori/varianceY);
		}
		else {
			double invDet = sigma2apriori / ((corrCoefXY + 1.0) * varianceX * varianceY);
			P = new UpperSymmPackMatrix(numberOfRows);
			P.set(0, 0,  invDet * varianceY);
			P.set(1, 1,  invDet * varianceX);
			P.set(0, 1, -invDet * corrCoefXY * Math.sqrt(varianceX * varianceY));
		}

		w.set(0, imageCoordinate.getX().getValue() - (x0 + xs + deltaX));
		w.set(1, imageCoordinate.getY().getValue() - (y0 + ys + deltaY));

		List<Integer> columns = new ArrayList<Integer>();
		int column = -1;

		// Object point coordinates
		column = objectCoordinate.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -(c*1.0/(r*r*r)*((r*r*r)*r11+dDist*(r*r)*r11+dRad*(r*r)*r11-kxN*(r*r*r)*r13+C1*(r*r*r)*r11-C2*cosPhi*(r*r*r)*sinKappa-C1*kxN*(r*r*r)*r13-C2*kyN*(r*r*r)*r13-dDist*kxN*(r*r)*r13*2.0-dRad*kxN*(r*r)*r13-(c*c)*dDist*(kxN*kxN)*r11+(c*c)*dDist*(kxN*kxN*kxN)*r13-(c*c)*dRad*(kxN*kxN)*r11+(c*c)*dRad*(kxN*kxN*kxN)*r13+B1*c*(kxN*kxN)*(r*r*r)*r13*6.0+B1*c*(kyN*kyN)*(r*r*r)*r13*2.0+(c*c)*dDist*kxN*(kyN*kyN)*r13+(c*c)*dRad*kxN*(kyN*kyN)*r13+A1*(c*c)*(kxN*kxN)*(r*r*r)*r11*3.0-A1*(c*c)*(kxN*kxN*kxN)*(r*r*r)*r13*3.0+A2*(c*c)*(kxN*kxN)*(r*r*r*r*r)*r11*5.0-A2*(c*c)*(kxN*kxN*kxN)*(r*r*r*r*r)*r13*5.0+A3*(c*c)*(kxN*kxN)*(r*r*r*r*r*r*r)*r11*7.0-A3*(c*c)*(kxN*kxN*kxN)*(r*r*r*r*r*r*r)*r13*7.0-B1*c*kxN*(r*r*r)*r11*6.0-B2*c*kyN*(r*r*r)*r11*2.0+B2*c*cosPhi*kxN*(r*r*r)*sinKappa*2.0+B1*c*cosPhi*kyN*(r*r*r)*sinKappa*2.0+B2*c*kxN*kyN*(r*r*r)*r13*4.0+(c*c)*cosPhi*dDist*kxN*kyN*sinKappa+(c*c)*cosPhi*dRad*kxN*kyN*sinKappa-A1*(c*c)*kxN*(kyN*kyN)*(r*r*r)*r13*3.0-A2*(c*c)*kxN*(kyN*kyN)*(r*r*r*r*r)*r13*5.0-A3*(c*c)*kxN*(kyN*kyN)*(r*r*r*r*r*r*r)*r13*7.0-A1*(c*c)*(kxN*kxN)*r*(r0*r0)*r11+A1*(c*c)*(kxN*kxN*kxN)*r*(r0*r0)*r13-A2*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0)*r11+A2*(c*c)*(kxN*kxN*kxN)*r*(r0*r0*r0*r0)*r13-A3*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0*r0*r0)*r11+A3*(c*c)*(kxN*kxN*kxN)*r*(r0*r0*r0*r0*r0*r0)*r13+A1*(c*c)*kxN*(kyN*kyN)*r*(r0*r0)*r13+A2*(c*c)*kxN*(kyN*kyN)*r*(r0*r0*r0*r0)*r13+A3*(c*c)*kxN*(kyN*kyN)*r*(r0*r0*r0*r0*r0*r0)*r13-A1*(c*c)*cosPhi*kxN*kyN*(r*r*r)*sinKappa*3.0-A2*(c*c)*cosPhi*kxN*kyN*(r*r*r*r*r)*sinKappa*5.0-A3*(c*c)*cosPhi*kxN*kyN*(r*r*r*r*r*r*r)*sinKappa*7.0+A1*(c*c)*cosPhi*kxN*kyN*r*(r0*r0)*sinKappa+A2*(c*c)*cosPhi*kxN*kyN*r*(r0*r0*r0*r0)*sinKappa+A3*(c*c)*cosPhi*kxN*kyN*r*(r0*r0*r0*r0*r0*r0)*sinKappa))/N+1.0/(N*N)*(c*c*c)*kxN*1.0/(r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)));
			A.set(1, column, -(c*1.0/(r*r*r)*(-cosPhi*(r*r*r)*sinKappa-kyN*(r*r*r)*r13-cosPhi*dDist*(r*r)*sinKappa-cosPhi*dRad*(r*r)*sinKappa-dDist*kyN*(r*r)*r13*2.0-dRad*kyN*(r*r)*r13+(c*c)*dDist*(kyN*kyN*kyN)*r13+(c*c)*dRad*(kyN*kyN*kyN)*r13+B2*c*(kxN*kxN)*(r*r*r)*r13*2.0+B2*c*(kyN*kyN)*(r*r*r)*r13*6.0+(c*c)*cosPhi*dDist*(kyN*kyN)*sinKappa+(c*c)*cosPhi*dRad*(kyN*kyN)*sinKappa+(c*c)*dDist*(kxN*kxN)*kyN*r13+(c*c)*dRad*(kxN*kxN)*kyN*r13-A1*(c*c)*(kyN*kyN*kyN)*(r*r*r)*r13*3.0-A2*(c*c)*(kyN*kyN*kyN)*(r*r*r*r*r)*r13*5.0-A3*(c*c)*(kyN*kyN*kyN)*(r*r*r*r*r*r*r)*r13*7.0-B2*c*kxN*(r*r*r)*r11*2.0-B1*c*kyN*(r*r*r)*r11*2.0-(c*c)*dDist*kxN*kyN*r11-(c*c)*dRad*kxN*kyN*r11+B1*c*cosPhi*kxN*(r*r*r)*sinKappa*2.0+B2*c*cosPhi*kyN*(r*r*r)*sinKappa*6.0+B1*c*kxN*kyN*(r*r*r)*r13*4.0+A1*(c*c)*kxN*kyN*(r*r*r)*r11*3.0+A2*(c*c)*kxN*kyN*(r*r*r*r*r)*r11*5.0+A3*(c*c)*kxN*kyN*(r*r*r*r*r*r*r)*r11*7.0-A1*(c*c)*cosPhi*(kyN*kyN)*(r*r*r)*sinKappa*3.0-A2*(c*c)*cosPhi*(kyN*kyN)*(r*r*r*r*r)*sinKappa*5.0-A3*(c*c)*cosPhi*(kyN*kyN)*(r*r*r*r*r*r*r)*sinKappa*7.0-A1*(c*c)*(kxN*kxN)*kyN*(r*r*r)*r13*3.0-A2*(c*c)*(kxN*kxN)*kyN*(r*r*r*r*r)*r13*5.0-A3*(c*c)*(kxN*kxN)*kyN*(r*r*r*r*r*r*r)*r13*7.0+A1*(c*c)*(kyN*kyN*kyN)*r*(r0*r0)*r13+A2*(c*c)*(kyN*kyN*kyN)*r*(r0*r0*r0*r0)*r13+A3*(c*c)*(kyN*kyN*kyN)*r*(r0*r0*r0*r0*r0*r0)*r13+A1*(c*c)*cosPhi*(kyN*kyN)*r*(r0*r0)*sinKappa+A2*(c*c)*cosPhi*(kyN*kyN)*r*(r0*r0*r0*r0)*sinKappa+A3*(c*c)*cosPhi*(kyN*kyN)*r*(r0*r0*r0*r0*r0*r0)*sinKappa+A1*(c*c)*(kxN*kxN)*kyN*r*(r0*r0)*r13+A2*(c*c)*(kxN*kxN)*kyN*r*(r0*r0*r0*r0)*r13+A3*(c*c)*(kxN*kxN)*kyN*r*(r0*r0*r0*r0*r0*r0)*r13-A1*(c*c)*kxN*kyN*r*(r0*r0)*r11-A2*(c*c)*kxN*kyN*r*(r0*r0*r0*r0)*r11-A3*(c*c)*kxN*kyN*r*(r0*r0*r0*r0*r0*r0)*r11))/N+1.0/(N*N)*(c*c*c)*kyN*1.0/(r*r)*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)));
		}

		column = objectCoordinate.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -(c*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/N+(cosPhi*sinOmega*xs)/N+(B1*(c*c)*(cosPhi*(kxN*kxN)*sinOmega*3.0+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa*3.0+cosKappa*kxN*r13*sinOmega*3.0-kyN*r13*sinOmega*sinKappa)*2.0)/N-(C1*c*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/N-(C2*c*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/N+(C1*cosPhi*sinOmega*xs)/N+(C2*cosPhi*sinOmega*ys)/N+(B2*(c*c)*kxN*(cosOmega*cosKappa-r13*sinOmega*sinKappa)*2.0)/N+(B2*(c*c)*kyN*(cosOmega*sinKappa+cosKappa*r13*sinOmega)*2.0)/N-(c*dDist*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/(N*r)-(c*dRad*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/(N*r)+1.0/(N*N)*(c*c)*1.0/(r*r)*xs*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa)+((c*c*c)*dDist*kxN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*dRad*kxN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+(cosPhi*dRad*sinOmega*xs)/(N*r)+((c*c)*1.0/(r*r)*xs*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+(B2*(c*c)*cosPhi*kxN*kyN*sinOmega*4.0)/N-(c*cosPhi*dDist*kxN*sinOmega*2.0)/(N*r));
			A.set(1, column, -(c*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/N+(cosPhi*sinOmega*ys)/N+(B2*(c*c)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega*3.0+cosOmega*cosKappa*kyN*3.0+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa*3.0)*2.0)/N+(B1*(c*c)*kxN*(cosOmega*cosKappa-r13*sinOmega*sinKappa)*2.0)/N+(B1*(c*c)*kyN*(cosOmega*sinKappa+cosKappa*r13*sinOmega)*2.0)/N-(c*dDist*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/(N*r)-(c*dRad*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/(N*r)+1.0/(N*N)*(c*c)*1.0/(r*r)*ys*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa)+((c*c*c)*dDist*kyN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*dRad*kyN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+(cosPhi*dRad*sinOmega*ys)/(N*r)+((c*c)*1.0/(r*r)*ys*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+(B1*(c*c)*cosPhi*kxN*kyN*sinOmega*4.0)/N-(c*cosPhi*dDist*kyN*sinOmega*2.0)/(N*r));
		}

		column = objectCoordinate.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -(c*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/N-(B1*(c*c)*((kxN*kxN)*r33*3.0+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa*3.0+cosOmega*cosKappa*kxN*r13*3.0-cosOmega*kyN*r13*sinKappa)*2.0)/N+(c*kxN*r33)/N-(C1*c*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/N-(C2*c*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/N+(C1*c*kxN*r33)/N+(C2*c*kyN*r33)/N+(B2*(c*c)*kxN*(cosKappa*sinOmega+cosOmega*r13*sinKappa)*2.0)/N+(B2*(c*c)*kyN*(sinOmega*sinKappa-cosOmega*cosKappa*r13)*2.0)/N-(c*dDist*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/(N*r)-(c*dRad*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/(N*r)-((c*c)*1.0/(r*r)*xs*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N-1.0/(N*N)*(c*c)*1.0/(r*r)*xs*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))-(B2*(c*c)*kxN*kyN*r33*4.0)/N-((c*c*c)*dDist*kxN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N-((c*c*c)*dRad*kxN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N+(c*dDist*kxN*r33*2.0)/(N*r)+(c*dRad*kxN*r33)/(N*r));
			A.set(1, column, -(c*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/N-(B2*(c*c)*((kxN*kxN)*r33+(kyN*kyN)*r33*3.0-cosKappa*kyN*sinOmega*3.0-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa*3.0)*2.0)/N+(c*kyN*r33)/N+(B1*(c*c)*kxN*(cosKappa*sinOmega+cosOmega*r13*sinKappa)*2.0)/N+(B1*(c*c)*kyN*(sinOmega*sinKappa-cosOmega*cosKappa*r13)*2.0)/N-(c*dDist*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/(N*r)-(c*dRad*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/(N*r)-((c*c)*1.0/(r*r)*ys*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N-1.0/(N*N)*(c*c)*1.0/(r*r)*ys*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))-(B1*(c*c)*kxN*kyN*r33*4.0)/N-((c*c*c)*dDist*kyN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N-((c*c*c)*dRad*kyN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N+(c*dDist*kyN*r33*2.0)/(N*r)+(c*dRad*kyN*r33)/(N*r));
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
			A.set(0, column, -kxN-C1*kxN-C2*kyN-(dDist*kxN)/r-(dRad*kxN)/r+B1*c*((kxN*kxN)*3.0+kyN*kyN)*2.0+B2*c*kxN*kyN*4.0+c*1.0/(r*r)*xs*(kxN*kxN+kyN*kyN)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))+(c*c)*dDist*kxN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)+(c*c)*dRad*kxN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)+(c*1.0/(r*r)*xs*(kxN*kxN+kyN*kyN)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)))/N);
			A.set(1, column, -kyN-(dDist*kyN)/r-(dRad*kyN)/r+B2*c*(kxN*kxN+(kyN*kyN)*3.0)*2.0+B1*c*kxN*kyN*4.0+c*1.0/(r*r)*ys*(kxN*kxN+kyN*kyN)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))+(c*c)*dDist*kyN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)+(c*c)*dRad*kyN*1.0/(r*r*r)*(kxN*kxN+kyN*kyN)+(c*1.0/(r*r)*ys*(kxN*kxN+kyN*kyN)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)))/N);
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
			A.set(0, column, r*r + 2.0*xs*xs); //c*c*(3.0*kxN*kxN + kyN*kyN);
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

		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs/N*(r*r - r0*r0)); 
			A.set(1, column, ys/N*(r*r - r0*r0));
		}
		
		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs/N*(r*r*r*r - r0*r0*r0*r0));
			A.set(1, column, ys/N*(r*r*r*r - r0*r0*r0*r0));
		}
		
		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, xs/N*(r*r*r*r*r*r - r0*r0*r0*r0*r0*r0));
			A.set(1, column, ys/N*(r*r*r*r*r*r - r0*r0*r0*r0*r0*r0));
		}

		
		// exterior orientation of the image
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (1.0/(r*r*r)*(c*(r*r*r)*r11+(r*r*r)*r13*xs+C1*c*(r*r*r)*r11+C2*c*(r*r*r)*r12+C1*(r*r*r)*r13*xs+C2*(r*r*r)*r13*ys+c*dDist*(r*r)*r11+c*dRad*(r*r)*r11+dRad*(r*r)*r13*xs-(c*c*c)*dDist*(kxN*kxN)*r11+(c*c*c)*dDist*(kxN*kxN*kxN)*r13-(c*c*c)*dRad*(kxN*kxN)*r11+(c*c*c)*dRad*(kxN*kxN*kxN)*r13-B1*(c*c)*kxN*(r*r*r)*r11*6.0-B2*(c*c)*kyN*(r*r*r)*r11*2.0+(c*c*c)*dDist*kxN*(kyN*kyN)*r13+(c*c*c)*dRad*kxN*(kyN*kyN)*r13+B1*(c*c)*(kxN*kxN)*(r*r*r)*r13*6.0+B1*(c*c)*(kyN*kyN)*(r*r*r)*r13*2.0-c*dDist*kxN*(r*r)*r13*2.0+(c*c*c)*cosPhi*dDist*kxN*kyN*sinKappa+(c*c*c)*cosPhi*dRad*kxN*kyN*sinKappa+B2*(c*c)*cosPhi*kxN*(r*r*r)*sinKappa*2.0+B1*(c*c)*cosPhi*kyN*(r*r*r)*sinKappa*2.0+B2*(c*c)*kxN*kyN*(r*r*r)*r13*4.0-A1*(c*c)*kxN*(r*r*r)*r11*xs*3.0-A2*(c*c)*kxN*(r*r*r*r*r)*r11*xs*5.0-A3*(c*c)*kxN*(r*r*r*r*r*r*r)*r11*xs*7.0+A1*(c*c)*(kxN*kxN)*(r*r*r)*r13*xs*3.0+A2*(c*c)*(kxN*kxN)*(r*r*r*r*r)*r13*xs*5.0+A3*(c*c)*(kxN*kxN)*(r*r*r*r*r*r*r)*r13*xs*7.0+A1*(c*c)*(kyN*kyN)*(r*r*r)*r13*xs*3.0+A2*(c*c)*(kyN*kyN)*(r*r*r*r*r)*r13*xs*5.0+A3*(c*c)*(kyN*kyN)*(r*r*r*r*r*r*r)*r13*xs*7.0-A1*(c*c)*(kxN*kxN)*r*(r0*r0)*r13*xs-A2*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0)*r13*xs-A3*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0*r0*r0)*r13*xs-A1*(c*c)*(kyN*kyN)*r*(r0*r0)*r13*xs-A2*(c*c)*(kyN*kyN)*r*(r0*r0*r0*r0)*r13*xs-A3*(c*c)*(kyN*kyN)*r*(r0*r0*r0*r0*r0*r0)*r13*xs+A1*(c*c)*cosPhi*kyN*(r*r*r)*sinKappa*xs*3.0+A2*(c*c)*cosPhi*kyN*(r*r*r*r*r)*sinKappa*xs*5.0+A3*(c*c)*cosPhi*kyN*(r*r*r*r*r*r*r)*sinKappa*xs*7.0+A1*(c*c)*kxN*r*(r0*r0)*r11*xs+A2*(c*c)*kxN*r*(r0*r0*r0*r0)*r11*xs+A3*(c*c)*kxN*r*(r0*r0*r0*r0*r0*r0)*r11*xs-A1*(c*c)*cosPhi*kyN*r*(r0*r0)*sinKappa*xs-A2*(c*c)*cosPhi*kyN*r*(r0*r0*r0*r0)*sinKappa*xs-A3*(c*c)*cosPhi*kyN*r*(r0*r0*r0*r0*r0*r0)*sinKappa*xs))/N+1.0/(N*N)*(c*c)*1.0/(r*r)*xs*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)));
			A.set(1, column, (1.0/(r*r*r)*(c*(r*r*r)*r12+(r*r*r)*r13*ys+c*dDist*(r*r)*r12+c*dRad*(r*r)*r12+dRad*(r*r)*r13*ys+(c*c*c)*dDist*(kyN*kyN*kyN)*r13+(c*c*c)*dRad*(kyN*kyN*kyN)*r13-B2*(c*c)*kxN*(r*r*r)*r11*2.0-B1*(c*c)*kyN*(r*r*r)*r11*2.0+(c*c*c)*cosPhi*dDist*(kyN*kyN)*sinKappa+(c*c*c)*cosPhi*dRad*(kyN*kyN)*sinKappa+(c*c*c)*dDist*(kxN*kxN)*kyN*r13+(c*c*c)*dRad*(kxN*kxN)*kyN*r13+B2*(c*c)*(kxN*kxN)*(r*r*r)*r13*2.0+B2*(c*c)*(kyN*kyN)*(r*r*r)*r13*6.0-(c*c*c)*dDist*kxN*kyN*r11-(c*c*c)*dRad*kxN*kyN*r11-c*dDist*kyN*(r*r)*r13*2.0+B1*(c*c)*cosPhi*kxN*(r*r*r)*sinKappa*2.0+B2*(c*c)*cosPhi*kyN*(r*r*r)*sinKappa*6.0+B1*(c*c)*kxN*kyN*(r*r*r)*r13*4.0-A1*(c*c)*kxN*(r*r*r)*r11*ys*3.0-A2*(c*c)*kxN*(r*r*r*r*r)*r11*ys*5.0-A3*(c*c)*kxN*(r*r*r*r*r*r*r)*r11*ys*7.0+A1*(c*c)*(kxN*kxN)*(r*r*r)*r13*ys*3.0+A2*(c*c)*(kxN*kxN)*(r*r*r*r*r)*r13*ys*5.0+A3*(c*c)*(kxN*kxN)*(r*r*r*r*r*r*r)*r13*ys*7.0+A1*(c*c)*(kyN*kyN)*(r*r*r)*r13*ys*3.0+A2*(c*c)*(kyN*kyN)*(r*r*r*r*r)*r13*ys*5.0+A3*(c*c)*(kyN*kyN)*(r*r*r*r*r*r*r)*r13*ys*7.0-A1*(c*c)*(kxN*kxN)*r*(r0*r0)*r13*ys-A2*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0)*r13*ys-A3*(c*c)*(kxN*kxN)*r*(r0*r0*r0*r0*r0*r0)*r13*ys-A1*(c*c)*(kyN*kyN)*r*(r0*r0)*r13*ys-A2*(c*c)*(kyN*kyN)*r*(r0*r0*r0*r0)*r13*ys-A3*(c*c)*(kyN*kyN)*r*(r0*r0*r0*r0*r0*r0)*r13*ys+A1*(c*c)*cosPhi*kyN*(r*r*r)*sinKappa*ys*3.0+A2*(c*c)*cosPhi*kyN*(r*r*r*r*r)*sinKappa*ys*5.0+A3*(c*c)*cosPhi*kyN*(r*r*r*r*r*r*r)*sinKappa*ys*7.0+A1*(c*c)*kxN*r*(r0*r0)*r11*ys+A2*(c*c)*kxN*r*(r0*r0*r0*r0)*r11*ys+A3*(c*c)*kxN*r*(r0*r0*r0*r0*r0*r0)*r11*ys-A1*(c*c)*cosPhi*kyN*r*(r0*r0)*sinKappa*ys-A2*(c*c)*cosPhi*kyN*r*(r0*r0*r0*r0)*sinKappa*ys-A3*(c*c)*cosPhi*kyN*r*(r0*r0*r0*r0*r0*r0)*sinKappa*ys))/N+1.0/(N*N)*(c*c)*1.0/(r*r)*ys*(-kxN*r11+(kxN*kxN)*r13+(kyN*kyN)*r13+cosPhi*kyN*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0)));
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (c*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/N-(B1*(c*c)*(cosPhi*(kxN*kxN)*sinOmega*3.0+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa*3.0+cosKappa*kxN*r13*sinOmega*3.0-kyN*r13*sinOmega*sinKappa)*2.0)/N+(C1*c*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/N+(C2*c*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/N+(c*cosPhi*kxN*sinOmega)/N-(B2*(c*c)*kxN*(cosOmega*cosKappa-r13*sinOmega*sinKappa)*2.0)/N-(B2*(c*c)*kyN*(cosOmega*sinKappa+cosKappa*r13*sinOmega)*2.0)/N+(c*dDist*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/(N*r)+(c*dRad*(cosOmega*sinKappa+cosKappa*r13*sinOmega))/(N*r)-((c*c*c)*dDist*kxN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N-((c*c*c)*dRad*kxN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*kxN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+(C1*c*cosPhi*kxN*sinOmega)/N+(C2*c*cosPhi*kyN*sinOmega)/N+1.0/(N*N)*(c*c*c)*kxN*1.0/(r*r)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa)-(B2*(c*c)*cosPhi*kxN*kyN*sinOmega*4.0)/N+(c*cosPhi*dDist*kxN*sinOmega*2.0)/(N*r)+(c*cosPhi*dRad*kxN*sinOmega)/(N*r));
			A.set(1, column, (c*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/N-(B2*(c*c)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega*3.0+cosOmega*cosKappa*kyN*3.0+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa*3.0)*2.0)/N+(c*cosPhi*kyN*sinOmega)/N-(B1*(c*c)*kxN*(cosOmega*cosKappa-r13*sinOmega*sinKappa)*2.0)/N-(B1*(c*c)*kyN*(cosOmega*sinKappa+cosKappa*r13*sinOmega)*2.0)/N+(c*dDist*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/(N*r)+(c*dRad*(cosOmega*cosKappa-r13*sinOmega*sinKappa))/(N*r)-((c*c*c)*dDist*kyN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N-((c*c*c)*dRad*kyN*1.0/(r*r*r)*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*kyN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa))/N+1.0/(N*N)*(c*c*c)*kyN*1.0/(r*r)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*(kxN*kxN)*sinOmega+cosPhi*(kyN*kyN)*sinOmega+cosOmega*cosKappa*kyN+cosOmega*kxN*sinKappa+cosKappa*kxN*r13*sinOmega-kyN*r13*sinOmega*sinKappa)-(B1*(c*c)*cosPhi*kxN*kyN*sinOmega*4.0)/N+(c*cosPhi*dDist*kyN*sinOmega*2.0)/(N*r)+(c*cosPhi*dRad*kyN*sinOmega)/(N*r));
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, (c*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/N+(r33*xs)/N+(B1*(c*c)*((kxN*kxN)*r33*3.0+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa*3.0+cosOmega*cosKappa*kxN*r13*3.0-cosOmega*kyN*r13*sinKappa)*2.0)/N+(C1*r33*xs)/N+(C2*r33*ys)/N+(C1*c*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/N+(C2*c*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/N+(dRad*r33*xs)/(N*r)-(B2*(c*c)*kxN*(cosKappa*sinOmega+cosOmega*r13*sinKappa)*2.0)/N-(B2*(c*c)*kyN*(sinOmega*sinKappa-cosOmega*cosKappa*r13)*2.0)/N+(c*dDist*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/(N*r)+(c*dRad*(sinOmega*sinKappa-cosOmega*cosKappa*r13))/(N*r)-((c*c*c)*kxN*1.0/(r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N-1.0/(N*N)*(c*c*c)*kxN*1.0/(r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))+(B2*(c*c)*kxN*kyN*r33*4.0)/N+((c*c*c)*dDist*kxN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N+((c*c*c)*dRad*kxN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N-(c*dDist*kxN*r33*2.0)/(N*r));
			A.set(1, column, (c*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/N+(r33*ys)/N+(B2*(c*c)*((kxN*kxN)*r33+(kyN*kyN)*r33*3.0-cosKappa*kyN*sinOmega*3.0-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa*3.0)*2.0)/N+(dRad*r33*ys)/(N*r)-(B1*(c*c)*kxN*(cosKappa*sinOmega+cosOmega*r13*sinKappa)*2.0)/N-(B1*(c*c)*kyN*(sinOmega*sinKappa-cosOmega*cosKappa*r13)*2.0)/N+(c*dDist*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/(N*r)+(c*dRad*(cosKappa*sinOmega+cosOmega*r13*sinKappa))/(N*r)-((c*c*c)*kyN*1.0/(r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0)))/N-1.0/(N*N)*(c*c*c)*kyN*1.0/(r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))+(B1*(c*c)*kxN*kyN*r33*4.0)/N+((c*c*c)*dDist*kyN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N+((c*c*c)*dRad*kyN*1.0/(r*r*r)*((kxN*kxN)*r33+(kyN*kyN)*r33-cosKappa*kyN*sinOmega-kxN*sinOmega*sinKappa+cosOmega*cosKappa*kxN*r13-cosOmega*kyN*r13*sinKappa))/N-(c*dDist*kyN*r33*2.0)/(N*r));
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column); 
			A.set(0, column, -(c*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13))/N+(xs*(dY*r33+cosPhi*dZ*sinOmega))/N+(B1*(c*c)*(dY*(kxN*kxN)*r33*3.0+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa*3.0-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa*3.0+cosPhi*dZ*(kxN*kxN)*sinOmega*3.0+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13*3.0-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega*3.0-dZ*kyN*r13*sinOmega*sinKappa)*2.0)/N-(C1*c*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13))/N+(C2*c*(-cosOmega*cosKappa*dZ+cosKappa*dY*sinOmega+cosOmega*dY*r13*sinKappa+dZ*r13*sinOmega*sinKappa))/N+(C1*xs*(dY*r33+cosPhi*dZ*sinOmega))/N+(C2*ys*(dY*r33+cosPhi*dZ*sinOmega))/N-(B2*(c*c)*kxN*(-cosOmega*cosKappa*dZ+cosKappa*dY*sinOmega+cosOmega*dY*r13*sinKappa+dZ*r13*sinOmega*sinKappa)*2.0)/N+(B2*(c*c)*kyN*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13)*2.0)/N-(c*dDist*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13))/(N*r)-(c*dRad*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13))/(N*r)+(dRad*xs*(dY*r33+cosPhi*dZ*sinOmega))/(N*r)+(B2*(c*c)*kxN*kyN*(dY*r33+cosPhi*dZ*sinOmega)*4.0)/N-(c*dDist*kxN*(dY*r33+cosPhi*dZ*sinOmega)*2.0)/(N*r)+((c*c)*1.0/(r*r)*xs*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N+1.0/(N*N)*(c*c)*1.0/(r*r)*xs*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa)+((c*c*c)*dDist*kxN*1.0/(r*r*r)*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*dRad*kxN*1.0/(r*r*r)*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N);
			A.set(1, column, -(c*(cosOmega*cosKappa*dZ-cosKappa*dY*sinOmega-cosOmega*dY*r13*sinKappa-dZ*r13*sinOmega*sinKappa))/N+(ys*(dY*r33+cosPhi*dZ*sinOmega))/N+(B2*(c*c)*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33*3.0+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega*3.0-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega*3.0+cosOmega*cosKappa*dZ*kyN*3.0+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa*3.0+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa*3.0)*2.0)/N-(B1*(c*c)*kxN*(-cosOmega*cosKappa*dZ+cosKappa*dY*sinOmega+cosOmega*dY*r13*sinKappa+dZ*r13*sinOmega*sinKappa)*2.0)/N+(B1*(c*c)*kyN*(cosOmega*dZ*sinKappa-dY*sinOmega*sinKappa+cosKappa*dZ*r13*sinOmega+cosOmega*cosKappa*dY*r13)*2.0)/N+(c*dDist*(-cosOmega*cosKappa*dZ+cosKappa*dY*sinOmega+cosOmega*dY*r13*sinKappa+dZ*r13*sinOmega*sinKappa))/(N*r)+(c*dRad*(-cosOmega*cosKappa*dZ+cosKappa*dY*sinOmega+cosOmega*dY*r13*sinKappa+dZ*r13*sinOmega*sinKappa))/(N*r)+(dRad*ys*(dY*r33+cosPhi*dZ*sinOmega))/(N*r)+(B1*(c*c)*kxN*kyN*(dY*r33+cosPhi*dZ*sinOmega)*4.0)/N-(c*dDist*kyN*(dY*r33+cosPhi*dZ*sinOmega)*2.0)/(N*r)+((c*c)*1.0/(r*r)*ys*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N+1.0/(N*N)*(c*c)*1.0/(r*r)*ys*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa)+((c*c*c)*dDist*kyN*1.0/(r*r*r)*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N+((c*c*c)*dRad*kyN*1.0/(r*r*r)*(dY*(kxN*kxN)*r33+dY*(kyN*kyN)*r33+cosOmega*dZ*kxN*sinKappa-cosKappa*dY*kyN*sinOmega-dY*kxN*sinOmega*sinKappa+cosPhi*dZ*(kxN*kxN)*sinOmega+cosPhi*dZ*(kyN*kyN)*sinOmega+cosOmega*cosKappa*dZ*kyN+cosOmega*cosKappa*dY*kxN*r13-cosOmega*dY*kyN*r13*sinKappa+cosKappa*dZ*kxN*r13*sinOmega-dZ*kyN*r13*sinOmega*sinKappa))/N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, -B1*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*6.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*6.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N)+(c*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega))/N-(C2*c*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa))/N+(c*kxN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/N+(C1*c*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega))/N-(c*dDist*kxN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N))/2.0-(c*dRad*kxN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N))/2.0-(B2*(c*c)*kyN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N+(B2*(c*c)*kxN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N+(c*dDist*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega))/(N*r)+(c*dRad*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega))/(N*r)+(C1*c*kxN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/N+(C2*c*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/N+((c*c*c)*kxN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*dX*(kxN*kxN)+cosPhi*dX*(kyN*kyN)+cosKappa*dX*kxN*r13+cosOmega*dZ*kxN*r11-dY*kxN*r11*sinOmega-dX*kyN*r13*sinKappa-dY*kyN*r12*sinOmega-dZ*kyN*r33*sinKappa-cosOmega*dZ*(kxN*kxN)*r13-cosOmega*dZ*(kyN*kyN)*r13+dY*(kxN*kxN)*r13*sinOmega+dY*(kyN*kyN)*r13*sinOmega))/N+1.0/(N*N)*(c*c*c)*kxN*1.0/(r*r)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*dX*(kxN*kxN)+cosPhi*dX*(kyN*kyN)+cosKappa*dX*kxN*r13+cosOmega*dZ*kxN*r11-dY*kxN*r11*sinOmega-dX*kyN*r13*sinKappa-dY*kyN*r12*sinOmega-dZ*kyN*r33*sinKappa-cosOmega*dZ*(kxN*kxN)*r13-cosOmega*dZ*(kyN*kyN)*r13+dY*(kxN*kxN)*r13*sinOmega+dY*(kyN*kyN)*r13*sinOmega)-(B2*(c*c)*kxN*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*4.0)/N+(c*dDist*kxN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/(N*r)+(c*dRad*kxN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/(N*r));
			A.set(1, column, -B2*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*6.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*6.0)/N)-(c*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa))/N+(c*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/N-(c*dDist*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa))/(N*r)-(c*dRad*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa))/(N*r)-(c*dDist*kyN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N))/2.0-(c*dRad*kyN*1.0/(r*r*r)*(((c*c)*(kxN*kxN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*(kyN*kyN)*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/N+((c*c)*kxN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N-((c*c)*kyN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N))/2.0-(B1*(c*c)*kyN*(cosKappa*dX*r13+cosOmega*dZ*r11-dY*r11*sinOmega)*2.0)/N+(B1*(c*c)*kxN*(dX*r13*sinKappa+dY*r12*sinOmega+dZ*r33*sinKappa)*2.0)/N+((c*c*c)*kyN*1.0/(r*r)*(A1*(r*r)*3.0-A1*(r0*r0)+A2*(r*r*r*r)*5.0-A2*(r0*r0*r0*r0)+A3*(r*r*r*r*r*r)*7.0-A3*(r0*r0*r0*r0*r0*r0))*(cosPhi*dX*(kxN*kxN)+cosPhi*dX*(kyN*kyN)+cosKappa*dX*kxN*r13+cosOmega*dZ*kxN*r11-dY*kxN*r11*sinOmega-dX*kyN*r13*sinKappa-dY*kyN*r12*sinOmega-dZ*kyN*r33*sinKappa-cosOmega*dZ*(kxN*kxN)*r13-cosOmega*dZ*(kyN*kyN)*r13+dY*(kxN*kxN)*r13*sinOmega+dY*(kyN*kyN)*r13*sinOmega))/N+1.0/(N*N)*(c*c*c)*kyN*1.0/(r*r)*(D1*(r*r)*3.0-D1*(r0*r0)+D2*(r*r*r*r)*5.0-D2*(r0*r0*r0*r0)+D3*(r*r*r*r*r*r)*7.0-D3*(r0*r0*r0*r0*r0*r0))*(cosPhi*dX*(kxN*kxN)+cosPhi*dX*(kyN*kyN)+cosKappa*dX*kxN*r13+cosOmega*dZ*kxN*r11-dY*kxN*r11*sinOmega-dX*kyN*r13*sinKappa-dY*kyN*r12*sinOmega-dZ*kyN*r33*sinKappa-cosOmega*dZ*(kxN*kxN)*r13-cosOmega*dZ*(kyN*kyN)*r13+dY*(kxN*kxN)*r13*sinOmega+dY*(kyN*kyN)*r13*sinOmega)-(B1*(c*c)*kxN*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*4.0)/N+(c*dDist*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega)*2.0)/(N*r)+(c*dRad*kyN*(cosPhi*dX-cosOmega*dZ*r13+dY*r13*sinOmega))/(N*r));
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, ys+C1*ys-B2*(c*c)*(kxN*kxN)*2.0+B2*(c*c)*(kyN*kyN)*2.0+C2*c*kxN+(ys*(dDist+dRad))/r+B1*(c*c)*kxN*kyN*4.0);
			A.set(1, column, c*(kxN-B1*c*(kxN*kxN)*2.0+B1*c*(kyN*kyN)*2.0-B2*c*kxN*kyN*4.0)+(c*kxN*(dDist+dRad))/r);
		}
		
		Collections.sort(columns);
		for (int row = 0; row < numberOfRows; row++) {
			for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
				int colAT = columns.get(columnATIdx);
				double aT = A.get(row, colAT);

				if (corrCoefXY == 0)
					neq.add(colAT, aT * P.get(row, row) * w.get(row));
				else
					for (int colP = 0; colP < numberOfRows; colP++)
						neq.add(colAT, aT * P.get(row, colP) * w.get(colP));

				if (NEQ != null) {
					for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
						int colA = columns.get(columnAIdx);
						if (corrCoefXY == 0)
							NEQ.add(colAT, colA, aT * P.get(row, row) * A.get(row, colA));
						else
							for (int colP = 0; colP < numberOfRows; colP++)
								NEQ.add(colAT, colA, aT * P.get(row, colP) * A.get(colP, colA));
					}
				}
			}	
		}
		return new GaussMarkovEquations(A, P, w);
	}
}
