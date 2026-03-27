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

package org.applied_geodesy.adjustment.bundle.derivation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.ScaleBar;
import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.camera.Image;
import org.applied_geodesy.adjustment.bundle.camera.ImageCoordinate;
import org.applied_geodesy.adjustment.bundle.camera.distortion.AffinityShearDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadialDistanceDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadiallySymmetricDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.TangentialDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.ZernikeDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.camera.orientation.InteriorOrientation;
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

final public class PartialDerivativeFactory {
	static int cnt = 0;
	final static class CollinearityEquationFactory {
		// substitute parameters
	    // rotations angles
	    final double cosOmega, sinOmega;
	    final double cosPhi, sinPhi;
	    final double cosKappa, sinKappa;

	    // rotation elements
	    final double r11, r12, r13;
	    final double r21, r22, r23;
	    final double r31, r32, r33;

	    // coordinates
	    final double xs, ys, x, y;
	    final double dX, dY, dZ;
	    
	    // numerators and denominator
	    final double kx, ky, N, kxN, kyN;
	    
	    // partial derivatives
		final double par_xs_X, par_xs_Y, par_xs_Z;
		final double par_xs_c, par_xs_x0, par_xs_y0;
		
		final double par_xs_X0, par_xs_Y0, par_xs_Z0;
		final double par_xs_omega, par_xs_phi, par_xs_kappa;
		
		final double par_ys_X, par_ys_Y, par_ys_Z;
		final double par_ys_c, par_ys_x0, par_ys_y0;
		
		final double par_ys_X0, par_ys_Y0, par_ys_Z0;
		final double par_ys_omega, par_ys_phi, par_ys_kappa;
		
		final InteriorOrientation interiorOrientation;
		final ExteriorOrientation exteriorOrientation;
		final ObjectCoordinate objectCoordinate;
	    
		private CollinearityEquationFactory(InteriorOrientation interiorOrientation, ExteriorOrientation exteriorOrientation, ObjectCoordinate objectCoordinate) {
			this.interiorOrientation = interiorOrientation;
			this.exteriorOrientation = exteriorOrientation;
			this.objectCoordinate    = objectCoordinate;
			
	    	double c  = interiorOrientation.getPrincipleDistance().getValue();
	    	double x0 = interiorOrientation.getPrinciplePointX().getValue();
	    	double y0 = interiorOrientation.getPrinciplePointY().getValue();
			
	    	double X0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getValue();
	    	double Y0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getValue();
	    	double Z0 = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getValue();

	    	double omega = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getValue();
	    	double phi   = exteriorOrientation.get(ParameterType.CAMERA_PHI).getValue();
	    	double kappa = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getValue();
	    	
	    	double X = objectCoordinate.getX().getValue();
			double Y = objectCoordinate.getY().getValue();
			double Z = objectCoordinate.getZ().getValue();

			this.cosOmega = Math.cos(omega);
			this.sinOmega = Math.sin(omega);

			this.cosPhi = Math.cos(phi);
			this.sinPhi = Math.sin(phi);

			this.cosKappa = Math.cos(kappa);
			this.sinKappa = Math.sin(kappa);

			// Rotation (Luhmann 2023, Eq. 2.31, p. 62)
			this.r11 =  this.cosPhi * this.cosKappa;
			this.r12 = -this.cosPhi * this.sinKappa;
			this.r13 =  this.sinPhi;

			this.r21 =  this.cosOmega * this.sinKappa + this.sinOmega * this.sinPhi * this.cosKappa;
			this.r22 =  this.cosOmega * this.cosKappa - this.sinOmega * this.sinPhi * this.sinKappa;
			this.r23 = -this.sinOmega * this.cosPhi;

			this.r31 = this.sinOmega * this.sinKappa - this.cosOmega * this.sinPhi * this.cosKappa;
			this.r32 = this.sinOmega * this.cosKappa + this.cosOmega * this.sinPhi * this.sinKappa;
			this.r33 = this.cosOmega * this.cosPhi;

			this.dX = X - X0;
			this.dY = Y - Y0;
			this.dZ = Z - Z0;
			
			this.kx = this.r11*this.dX + this.r21*this.dY + this.r31*this.dZ;
			this.ky = this.r12*this.dX + this.r22*this.dY + this.r32*this.dZ;
			this.N  = this.r13*this.dX + this.r23*this.dY + this.r33*this.dZ;

			this.kxN = this.kx/this.N;
			this.kyN = this.ky/this.N;

			this.xs = -c*this.kxN;
			this.ys = -c*this.kyN;

			this.x = x0 + this.xs;
			this.y = y0 + this.ys;


			// partial derivatives
			// collinearity x-equation
			this.par_xs_X = -(this.r13*this.xs + c*this.r11) / this.N;
			this.par_xs_Y = -(this.r23*this.xs + c*this.r21) / this.N;
			this.par_xs_Z = -(this.r33*this.xs + c*this.r31) / this.N;

			this.par_xs_x0 = 1.0;
			this.par_xs_y0 = 0.0;
			this.par_xs_c = -this.kxN;

			this.par_xs_X0 = -this.par_xs_X; // -c/N2 * (r13*kx - r11*N);
			this.par_xs_Y0 = -this.par_xs_Y; // -c/N2 * (r23*kx - r21*N);
			this.par_xs_Z0 = -this.par_xs_Z; // -c/N2 * (r33*kx - r31*N);

			this.par_xs_omega = (this.xs * (this.r33*this.dY - this.r23*this.dZ) + c * (this.r31*this.dY - this.r21*this.dZ)) / this.N;  // -c/N * ( kx/N * (r33*(Y - Y0) - r23*(Z - Z0)) - r31*(Y - Y0) + r21*(Z - Z0) )
			this.par_xs_phi   = (this.xs * (this.ky * this.sinKappa - this.kx * this.cosKappa) + c * this.N * this.cosKappa) / this.N;   // -c/N * ( kx/N * (ky * sinKappa - kx * cosKappa) - N * cosKappa )
			this.par_xs_kappa = this.ys;                                                                                                 // -c/N * ky;

			// partial derivatives
			// collinearity y-equation
			this.par_ys_X = -(this.r13*this.ys + c*this.r12) / this.N;
			this.par_ys_Y = -(this.r23*this.ys + c*this.r22) / this.N;
			this.par_ys_Z = -(this.r33*this.ys + c*this.r32) / this.N;

			this.par_ys_x0 = 0.0;
			this.par_ys_y0 = 1.0;
			this.par_ys_c = -this.kyN;

			this.par_ys_X0 = -this.par_ys_X; // -c/N2 * (r13*ky - r12*N)
			this.par_ys_Y0 = -this.par_ys_Y; // -c/N2 * (r23*ky - r22*N)
			this.par_ys_Z0 = -this.par_ys_Z; // -c/N2 * (r33*ky - r32*N)

			this.par_ys_omega = (this.ys * (this.r33*this.dY - this.r23*this.dZ) + c * (this.r32*this.dY - this.r22*this.dZ)) / this.N; // -c/N * ( ky/N * (r33*(Y - Y0) - r23*(Z - Z0)) - r32*(Y - Y0) + r22*(Z - Z0) )
			this.par_ys_phi   = (this.ys * (this.ky * this.sinKappa - this.kx * this.cosKappa) - c * this.N * this.sinKappa) / this.N;  // -c/N * ( ky/N * (ky * sinKappa - kx * cosKappa) + N * sinKappa )
			this.par_ys_kappa = -this.xs;                                                                                               //  c/N * kx
	    }
		
		private static CollinearityEquationFactory getInstance(InteriorOrientation interiorOrientation, ExteriorOrientation exteriorOrientation, ObjectCoordinate objectCoordinate) {
			return new CollinearityEquationFactory(interiorOrientation, exteriorOrientation, objectCoordinate);
		}
	}
	
	private PartialDerivativeFactory() {}
	
	public static GaussMarkovEquations getPartialDerivative(double sigma2apriori, UpperSymmPackMatrix NEQ, DenseVector neq, ObservationParameterGroup<?> observations) throws UnsupportedOperationException {
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
		
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
		
		CollinearityEquationFactory collinearityEquation = CollinearityEquationFactory.getInstance(interiorOrientation, exteriorOrientation, objectCoordinate);
		
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

		w.set(0, imageCoordinate.getX().getValue() - collinearityEquation.x);
		w.set(1, imageCoordinate.getY().getValue() - collinearityEquation.y);
		
		Set<Integer> columns = new HashSet<Integer>();
		int column = -1;

		// Object point coordinates
		column = objectCoordinate.getX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_X);
			A.set(1, column, collinearityEquation.par_ys_X);
		}

		column = objectCoordinate.getY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_Y);
			A.set(1, column, collinearityEquation.par_ys_Y);
		}

		column = objectCoordinate.getZ().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_Z);
			A.set(1, column, collinearityEquation.par_ys_Z);
		}


		// interior orientation of the camera
		column = interiorOrientation.getPrinciplePointX().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_x0);
			A.set(1, column, collinearityEquation.par_ys_x0);
		}

		column = interiorOrientation.getPrinciplePointY().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_y0);
			A.set(1, column, collinearityEquation.par_ys_y0);
		}

		column = interiorOrientation.getPrincipleDistance().getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_c); 
			A.set(1, column, collinearityEquation.par_ys_c); 
		}

		
		// exterior orientation of the image
		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_X0);
			A.set(1, column, collinearityEquation.par_ys_X0);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_Y0);
			A.set(1, column, collinearityEquation.par_ys_Y0);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_Z0);
			A.set(1, column, collinearityEquation.par_ys_Z0);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column); 
			A.set(0, column, collinearityEquation.par_xs_omega);
			A.set(1, column, collinearityEquation.par_ys_omega);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_phi);
			A.set(1, column, collinearityEquation.par_ys_phi);
		}

		column = exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn();
		if (column >= 0 && column != Integer.MAX_VALUE) {
			columns.add(column);
			A.set(0, column, collinearityEquation.par_xs_kappa);
			A.set(1, column, collinearityEquation.par_ys_kappa);
		}
		
		// apply distortion model modifications
		Collection<DistortionModel> distortionModels = camera.getDistortionModels();
		for (DistortionModel distortionModel : distortionModels) {
			switch (distortionModel.getType()) {
			case AFFINITY_AND_SHEAR:
				AffinityShearDistortionModelFactory.apply((AffinityShearDistortionModel)distortionModel, collinearityEquation, columns, A, w);
				break;
			case DISTANCE_DISTORTION:
				RadialDistanceDistortionModelFactory.apply((RadialDistanceDistortionModel)distortionModel, collinearityEquation, columns, A, w);
				break;
			case RADIAL_DISTORTION:
				RadiallySymmetricDistortionModelFactory.apply((RadiallySymmetricDistortionModel)distortionModel, collinearityEquation, columns, A, w);
				break;
			case TANGENTIAL_DISTORTION:
				TangentialDistortionModelFactory.apply((TangentialDistortionModel)distortionModel, collinearityEquation, columns, A, w);
				break;
			case ZERNIKE_POLYNOMIAL:
				ZernikeDistortionModelFactory.apply((ZernikeDistortionModel)distortionModel, collinearityEquation, columns, A, w);
				break;
			}
		}
		
		return stackNormalEquationSystem(NEQ, neq, A, P, w, new ArrayList<Integer>(columns), diagonalWeighting);
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
