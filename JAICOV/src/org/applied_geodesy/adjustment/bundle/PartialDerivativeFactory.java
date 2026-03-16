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
import org.applied_geodesy.adjustment.bundle.parameter.DirectlyObservedParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmBandMatrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;

class PartialDerivativeFactory {
	private PartialDerivativeFactory() {}
	
	static GaussMarkovEquations getPartialDerivative(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ObservationParameterGroup<?> observations) throws UnsupportedOperationException {
		if (observations instanceof ImageCoordinate) 
			return getPartialDerivativeImageCoordinate(sigma2apriori, NEQ, neq, (ImageCoordinate)observations);
		else if (observations instanceof ScaleBar)
			return getPartialDerivativeScaleBar(sigma2apriori, NEQ, neq, (ScaleBar)observations);
		else if (observations instanceof DirectlyObservedParameterGroup)
			return getPartialDerivativeDirectlyObservedParameters(sigma2apriori, NEQ, neq, (DirectlyObservedParameterGroup)observations);
		
		throw new UnsupportedOperationException("Error, unsupported or unknown observation type.");
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

		return stackNormalEquationSystem(NEQ, neq, A, P, w, columns, true);
	}

	private static GaussMarkovEquations getPartialDerivativeImageCoordinate(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ImageCoordinate imageCoordinate) {
		ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
		Image image = imageCoordinate.getReference();
		Camera camera = image.getReference();

		double X = objectCoordinate.getX().getValue();
		double Y = objectCoordinate.getY().getValue();
		double Z = objectCoordinate.getZ().getValue();

		double r0  = camera.getR0();
		double r02 = r0*r0;
		double r04 = r02*r02;
		double r06 = r04*r02;

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
		double B3 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B3).getValue();
		double B4 = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B4).getValue();

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

		// Rotation (Luhmann (2023, Eq. 2.31, p. 62)
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
		
		double xxs2 = 2.0 * xs * xs;
		double yys2 = 2.0 * ys * ys;
		double xys2 = 2.0 * xs * ys;

		//double r  = Math.hypot(xs, ys);
		double r2 = xs*xs + ys*ys;
		double r4 = r2*r2;
		double r6 = r4*r2;
		
		
		double dRad  = A1*r2 + A2*r4 + A3*r6 - (A1*r02 + A2*r04 + A3*r06);
		double dRadX = xs * dRad;
		double dRadY = ys * dRad;
		
		double constTanSym = (1.0 + B3*r2 + B4*r4);
		
		double dTanX = (B1 * (r2 + xxs2) + B2 * xys2) * constTanSym;
		double dTanY = (B2 * (r2 + yys2) + B1 * xys2) * constTanSym;
		
		double dAffX = C1*xs + C2*ys;
		double dAffY = 0.0;

		double dDist  = (D1*r2 + D2*r4 + D3*r6 - (D1*r02 + D2*r04 + D3*r06)) / N;
		double dDistX = xs * dDist;
		double dDistY = ys * dDist;
		
		// total corrections
		double deltaX = dRadX + dTanX + dAffX + dDistX;
		double deltaY = dRadY + dTanY + dAffY + dDistY;
		
		// partial derivatives
		// collinearity x-equation
		double par_xs_X = -(r13*xs + c*r11) / N;
		double par_xs_Y = -(r23*xs + c*r21) / N;
		double par_xs_Z = -(r33*xs + c*r31) / N;

		double par_xs_x0 = 1.0;
		double par_xs_y0 = 0.0;
		double par_xs_c = -kxN;

		double par_xs_X0 = -par_xs_X; // -c/N2 * (r13*kx - r11*N);
		double par_xs_Y0 = -par_xs_Y; // -c/N2 * (r23*kx - r21*N);
		double par_xs_Z0 = -par_xs_Z; // -c/N2 * (r33*kx - r31*N);

		double par_xs_omega = (xs * (r33*dY - r23*dZ) + c * (r31*dY - r21*dZ)) / N;                          // -c/N * ( kx/N * (r33*(Y - Y0) - r23*(Z - Z0)) - r31*(Y - Y0) + r21*(Z - Z0) )
		double par_xs_phi   = (xs * (ky * sinKappa - kx * cosKappa) + c * N * cosKappa) / N;                 // -c/N * ( kx/N * (ky * sinKappa - kx * cosKappa) - N * cosKappa )
		double par_xs_kappa = ys;                                                                            // -c/N * ky;

		double par_xs_A1 = xs*(r2 - r02);
		double par_xs_A2 = xs*(r4 - r04);
		double par_xs_A3 = xs*(r6 - r06);
		
		double par_xs_B1 = constTanSym * (xxs2 + r2);
		double par_xs_B2 = constTanSym * xys2;
		double par_xs_B3 = r2*(B1 * (xxs2 + r2) + B2 * xys2);
		double par_xs_B4 = r4*(B1 * (xxs2 + r2) + B2 * xys2);
		
		double par_xs_C1 = xs;
		double par_xs_C2 = ys;
		
		double par_xs_D1 = par_xs_A1 / N;
		double par_xs_D2 = par_xs_A2 / N;
		double par_xs_D3 = par_xs_A3 / N;
		
		// partial derivatives
		// collinearity y-equation
		double par_ys_X = -(r13*ys + c*r12) / N;
		double par_ys_Y = -(r23*ys + c*r22) / N;
		double par_ys_Z = -(r33*ys + c*r32) / N;

		double par_ys_x0 = 0.0;
		double par_ys_y0 = 1.0;
		double par_ys_c = -kyN;

		double par_ys_X0 = -par_ys_X; // -c/N2 * (r13*ky - r12*N)
		double par_ys_Y0 = -par_ys_Y; // -c/N2 * (r23*ky - r22*N)
		double par_ys_Z0 = -par_ys_Z; // -c/N2 * (r33*ky - r32*N)

		double par_ys_omega = (ys * (r33*dY - r23*dZ) + c * (r32*dY - r22*dZ)) / N;                         // -c/N * ( ky/N * (r33*(Y - Y0) - r23*(Z - Z0)) - r32*(Y - Y0) + r22*(Z - Z0) )
		double par_ys_phi   = (ys * (ky * sinKappa - kx * cosKappa) - c * N * sinKappa) / N;                // -c/N * ( ky/N * (ky * sinKappa - kx * cosKappa) + N * sinKappa )
		double par_ys_kappa = -xs;                                                                          //  c/N * kx
		
		double par_ys_A1 = ys*(r2 - r02);
		double par_ys_A2 = ys*(r4 - r04);
		double par_ys_A3 = ys*(r6 - r06);
		
		double par_ys_B1 = constTanSym * xys2;
		double par_ys_B2 = constTanSym * (yys2 + r2);
		double par_ys_B3 = r2*(B2 * (yys2 + r2) + B1 * xys2);     
		double par_ys_B4 = r4*(B2 * (yys2 + r2) + B1 * xys2);
		
		double par_ys_C1 = 0.0;
		double par_ys_C2 = 0.0;
		
		double par_ys_D1 = par_ys_A1 / N;
		double par_ys_D2 = par_ys_A2 / N;
		double par_ys_D3 = par_ys_A3 / N;
		

		// chain rule coefficients: dRad
		double constRad = A1 + 2.0*A2*r2 + 3.0*A3*r4;
		double par_dRadX_xs = xxs2 * constRad + dRad;
		double par_dRadX_ys = xys2 * constRad;

		double par_dRadY_xs = xys2 * constRad;
		double par_dRadY_ys = yys2 * constRad + dRad;
		
		// chain rule coefficients: dTan
		double constTan = 2.0*(B3 + 2.0*B4*r2);
		double par_dTanX_xs = xs*constTan*(B1*(xxs2 + r2) + B2*xys2) + 2.0*(3*B1*xs + B2*ys)*constTanSym;
		double par_dTanX_ys = ys*constTan*(B1*(xxs2 + r2) + B2*xys2) + 2.0*(  B2*xs + B1*ys)*constTanSym;
		
		double par_dTanY_xs = xs*constTan*(B2*(yys2 + r2) + B1*xys2) + 2.0*(  B1*ys + B2*xs)*constTanSym;
		double par_dTanY_ys = ys*constTan*(B2*(yys2 + r2) + B1*xys2) + 2.0*(3*B2*ys + B1*xs)*constTanSym;
		
		// chain rule coefficients: dAff
		double par_dAffX_xs = C1;
		double par_dAffX_ys = C2;
		
		double par_dAffY_xs = 0.0;
		double par_dAffY_ys = 0.0;
		
		// chain rule coefficients: dDist
		double par_N_X = r13;
		double par_N_Y = r23;
		double par_N_Z = r33;

		double par_N_X0 = -r13;
		double par_N_Y0 = -r23;
		double par_N_Z0 = -r33;
				
		double par_N_omega = -r33 * dY + r23 * dZ;              // -cosOmega*cosPhi * dY - sinOmega*cosPhi * dZ;
		double par_N_phi   =  kx * cosKappa - ky * sinKappa;    //  cosPhi*dX + sinOmega*sinPhi*dY - cosOmega*sinPhi*dZ
		double par_N_kappa =  0.0;
		
		double constDist = (D1 + 2*D2*r2 + 3*D3*r4) / N;
		double par_dDistX_xs = xxs2 * constDist + dDist;
		double par_dDistX_ys = xys2 * constDist;
		double par_dDistX_N  = -dDistX / N; 
		
		double par_dDistY_xs = xys2 * constDist;
		double par_dDistY_ys = yys2 * constDist + dDist;
		double par_dDistY_N  = -dDistY / N; 
		
		// partial derivatives dRad:
		//  dRadX: par_dRadX_xs * par_xs_<PARAM>  +  par_dRadX_ys * par_ys_<PARAM>
		//  dRadY: par_dRadY_xs * par_xs_<PARAM>  +  par_dRadY_ys * par_ys_<PARAM>
		
		// partial derivatives dTan:
		//  dTanX: par_dTanX_xs * par_xs_<PARAM>  +  par_dTanX_ys * par_ys_<PARAM>
		//  dTanY: par_dTanY_xs * par_xs_<PARAM>  +  par_dTanY_ys * par_ys_<PARAM>
		
		// partial derivatives dAff:
		//  dAffX: par_dAffX_xs * par_xs_<PARAM>  +  par_dAffX_ys * par_ys_<PARAM>
		//  dAffY: par_dAffY_xs * par_xs_<PARAM>  +  par_dAffY_ys * par_ys_<PARAM>
		
		// partial derivatives dDist:
		//  dDistX: par_dDistX_xs * par_xs_<PARAM>  +  par_dDistX_ys * par_ys_<PARAM>  +  par_dDistX_N * par_N_<PARAM>;
		//  dDistY: par_dDistY_xs * par_xs_<PARAM>  +  par_dDistY_ys * par_ys_<PARAM>  +  par_dDistY_N * par_N_<PARAM>;
		
		// sum of coefficients in x-equation 
		// depending on par_xs_<PARAM> 
		double par_corrX_xs = par_dRadX_xs + par_dTanX_xs + par_dAffX_xs + par_dDistX_xs;
		// depending on par_ys_<PARAM>
		double par_corrX_ys = par_dRadX_ys + par_dTanX_ys + par_dAffX_ys + par_dDistX_ys;
		
		// sum of coefficients in y-equation 
		// depending on par_xs_<PARAM> 
		double par_corrY_xs = par_dRadY_xs + par_dTanY_xs + par_dAffY_xs + par_dDistY_xs;
		// depending on par_ys_<PARAM> 
		double par_corrY_ys = par_dRadY_ys + par_dTanY_ys + par_dAffY_ys + par_dDistY_ys;
		
		// stochastic model
		double varianceX  = imageCoordinate.getX().getVariance();
		double varianceY  = imageCoordinate.getY().getVariance();
		double corrCoefXY = imageCoordinate.getCorrelationCoefficientXY();
		
		boolean diagonalWeighting = corrCoefXY == 0;
		
		int numberOfRows = imageCoordinate.getNumberOfParameters();
		
		DenseVector w = new DenseVector(numberOfRows);
		DenseMatrix A = new DenseMatrix(numberOfRows, neq.size());

		Matrix P = null;
		if (diagonalWeighting) {
			P = new UpperSymmBandMatrix(numberOfRows, 0);
			P.set(0, 0, sigma2apriori/varianceX);
			P.set(1, 1, sigma2apriori/varianceY);
		}
		else {
			double invDet = sigma2apriori / ((1.0 - corrCoefXY*corrCoefXY) * varianceX * varianceY);
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
			A.set(0, column, par_xs_X * (1.0 + par_corrX_xs) + par_ys_X * par_corrX_ys + par_N_X * par_dDistX_N);
			A.set(1, column, par_ys_X * (1.0 + par_corrY_ys) + par_xs_X * par_corrY_xs + par_N_X * par_dDistY_N);
		}

		column = objectCoordinate.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_Y * (1.0 + par_corrX_xs) + par_ys_Y * par_corrX_ys + par_N_Y * par_dDistX_N);
			A.set(1, column, par_ys_Y * (1.0 + par_corrY_ys) + par_xs_Y * par_corrY_xs + par_N_Y * par_dDistY_N);
		}

		column = objectCoordinate.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_Z * (1.0 + par_corrX_xs) + par_ys_Z * par_corrX_ys + par_N_Z * par_dDistX_N);
			A.set(1, column, par_ys_Z * (1.0 + par_corrY_ys) + par_xs_Z * par_corrY_xs + par_N_Z * par_dDistY_N);
		}


		// interior orientation of the camera
		column = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_x0);
			A.set(1, column, par_ys_x0);
		}

		column = interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_y0);
			A.set(1, column, par_ys_y0);
		}

		column = interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_c * (1.0 + par_corrX_xs) + par_ys_c * par_corrX_ys); 
			A.set(1, column, par_ys_c * (1.0 + par_corrY_ys) + par_xs_c * par_corrY_xs); 
		}

		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_A1);
			A.set(1, column, par_ys_A1);
		}

		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_A2);
			A.set(1, column, par_ys_A2);
		}

		column = interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_A3);
			A.set(1, column, par_ys_A3);
		}

		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B1); 
			A.set(1, column, par_ys_B1);
		}

		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B2);
			A.set(1, column, par_ys_B2);
		}
		
		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B3).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B3);
			A.set(1, column, par_ys_B3);
		}
		
		column = interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B4).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_B4);
			A.set(1, column, par_ys_B4);
		}

		column = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_C1);
			A.set(1, column, par_ys_C1);
		}

		column = interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_C2);
			A.set(1, column, par_ys_C2);
		}

		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_D1); 
			A.set(1, column, par_ys_D1);
		}
		
		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_D2);
			A.set(1, column, par_ys_D2);
		}
		
		column = interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_D3);
			A.set(1, column, par_ys_D3);
		}

		
		// exterior orientation of the image
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_X0 * (1.0 + par_corrX_xs) + par_ys_X0 * par_corrX_ys + par_N_X0 * par_dDistX_N);
			A.set(1, column, par_ys_X0 * (1.0 + par_corrY_ys) + par_xs_X0 * par_corrY_xs + par_N_X0 * par_dDistY_N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_Y0 * (1.0 + par_corrX_xs) + par_ys_Y0 * par_corrX_ys + par_N_Y0 * par_dDistX_N);
			A.set(1, column, par_ys_Y0 * (1.0 + par_corrY_ys) + par_xs_Y0 * par_corrY_xs + par_N_Y0 * par_dDistY_N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_Z0 * (1.0 + par_corrX_xs) + par_ys_Z0 * par_corrX_ys + par_N_Z0 * par_dDistX_N);
			A.set(1, column, par_ys_Z0 * (1.0 + par_corrY_ys) + par_xs_Z0 * par_corrY_xs + par_N_Z0 * par_dDistY_N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column); 
			A.set(0, column, par_xs_omega * (1.0 + par_corrX_xs) + par_ys_omega * par_corrX_ys + par_N_omega * par_dDistX_N);
			A.set(1, column, par_ys_omega * (1.0 + par_corrY_ys) + par_xs_omega * par_corrY_xs + par_N_omega * par_dDistY_N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_phi * (1.0 + par_corrX_xs) + par_ys_phi * par_corrX_ys + par_N_phi * par_dDistX_N);
			A.set(1, column, par_ys_phi * (1.0 + par_corrY_ys) + par_xs_phi * par_corrY_xs + par_N_phi * par_dDistY_N);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, par_xs_kappa * (1.0 + par_corrX_xs) + par_ys_kappa * par_corrX_ys + par_N_kappa * par_dDistX_N);
			A.set(1, column, par_ys_kappa * (1.0 + par_corrY_ys) + par_xs_kappa * par_corrY_xs + par_N_kappa * par_dDistY_N);
		}
		
		return stackNormalEquationSystem(NEQ, neq, A, P, w, columns, diagonalWeighting);
	}
	
	private static GaussMarkovEquations getPartialDerivativeDirectlyObservedParameters(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, DirectlyObservedParameterGroup observedParameterGroup) {
		boolean diagonalWeighting = !observedParameterGroup.hasFullyPopulatedWeightMatrix();

		int numberOfRows = observedParameterGroup.getNumberOfParameters();
		
		Matrix P = observedParameterGroup.getWeightMatrix(sigma2apriori);
		DenseVector w = new DenseVector(numberOfRows);
		DenseMatrix A = new DenseMatrix(numberOfRows, neq.size());
		
		List<Integer> columns = new ArrayList<Integer>();
		int row = 0;
		
		for (ObservationParameter<? extends UnknownParameter<?>> observedParameter : observedParameterGroup) {
			UnknownParameter<?> unknownParameter = observedParameter.getReference();
			int column = unknownParameter.getColumn();
			
			if (column >= 0 && column != Integer.MAX_VALUE) {
				columns.add(column);
				A.set(row, column, 1.0);
			}
			
			double parameterValue = unknownParameter.getValue();
			double observation    = observedParameter.getValue();
			w.set(row++, observation - parameterValue);
		}
		return stackNormalEquationSystem(NEQ, neq, A, P, w, columns, diagonalWeighting);
	}
	
	private static GaussMarkovEquations stackNormalEquationSystem(UpperSymmPackMatrix NEQ, DenseVector neq, Matrix A, Matrix P, DenseVector w, List<Integer> columns, boolean diagonalWeighting) {
		Collections.sort(columns);
		
		int numberOfRows = w.size();
		for (int row = 0; row < numberOfRows; row++) {
			for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
				int colAT = columns.get(columnATIdx);
				double aT = A.get(row, colAT);
				
				if (diagonalWeighting)
					neq.add(colAT, aT * P.get(row, row) * w.get(row));
				else
					for (int colP = 0; colP < numberOfRows; colP++)
						neq.add(colAT, aT * P.get(row, colP) * w.get(colP));

				if (NEQ != null) {
					for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
						int colA = columns.get(columnAIdx);
						
						if (diagonalWeighting)
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
