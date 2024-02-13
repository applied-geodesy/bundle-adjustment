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

package org.applied_geodesy.adjustment.bundle.dlt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ImageCoordinate;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.Constant;
import org.applied_geodesy.adjustment.DefaultValue;
import org.applied_geodesy.adjustment.MathExtension;
import org.applied_geodesy.adjustment.NormalEquationSystem;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.MatrixSingularException;
import no.uib.cipr.matrix.UpperSymmBandMatrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;
import no.uib.cipr.matrix.Vector;

public class DirectLinearTransformation {
	public enum RestrictionType {
		IDENTICAL_PRINCIPLE_DISTANCE,
		ROTATION_WITHOUT_SHEAR,
		FIXED_PRINCIPLE_DISTANCE_X,
		FIXED_PRINCIPLE_DISTANCE_Y,
		FIXED_PRINCIPAL_POINT_X,
		FIXED_PRINCIPAL_POINT_Y
	}

	private static double SQRT_EPS = Math.sqrt(Constant.EPS);
	
	private static int maximalNumberOfIterations = DefaultValue.getMaximalNumberOfIterations(),
			           numberOfUnknownParameters = 11;

	private DirectLinearTransformation() {}
	
	public static boolean adjust(DLTCoefficients coefficients, Map<String, ObjectCoordinate> objectCoordinates, RestrictionType ...restrictions) {
		restrictions = validateRestrictions(restrictions);
		prepareUnknwonParameters(coefficients);
		
		Image image = coefficients.getReference();
		
		double sumSqrDistanceWorld = 0; 
		double sumSqrDistanceImage = 0;
				
		List<ImageCoordinate> homologousImageCoordinates = new ArrayList<ImageCoordinate>();

		for (ImageCoordinate imageCoordinate : image) {
			String name = imageCoordinate.getObjectCoordinate().getName();
			if (objectCoordinates.containsKey(name)) {
				ObjectCoordinate objectCoordinate = objectCoordinates.get(name);
				homologousImageCoordinates.add(imageCoordinate);
				
				double x = imageCoordinate.getX().getValue();
				double y = imageCoordinate.getY().getValue();
				
				double X = objectCoordinate.getX().getValue();
				double Y = objectCoordinate.getY().getValue();
				double Z = objectCoordinate.getZ().getValue();
				
				sumSqrDistanceWorld += X*X + Y*Y + Z*Z; 
				sumSqrDistanceImage += x*x + y*y;
			}
		}
		
		if (homologousImageCoordinates.size() < 6) {
			try {
				throw new MatrixSingularException("Error, insufficient number of homologous points (" + homologousImageCoordinates.size() + " vs. 6) in image #" + image.getId());
			}
			catch(Exception e) {
				e.printStackTrace();
				return false;
			}
		}

		double scale = sumSqrDistanceImage > 0 ? Math.sqrt(sumSqrDistanceWorld / sumSqrDistanceImage) : 1.0;
		double maxAbsDx = 0.0;

		int runs = maximalNumberOfIterations - 1;
		boolean isEstimated = false, estimateCompleteModel = false, isConverge = true;

		if (maximalNumberOfIterations == 0) 
			estimateCompleteModel = isEstimated = true;

		try {
			boolean includeRestrictions = false;
			do {
				maxAbsDx = 0.0;
				
				// erzeuge Normalgleichung
				NormalEquationSystem neq = createNormalEquation(coefficients, homologousImageCoordinates, objectCoordinates, scale, includeRestrictions ? restrictions : new RestrictionType[0]);

				// Nutze Vorkonditionierung
				NormalEquationSystem.applyPrecondition(neq);

				DenseVector n = neq.getVector();
				UpperSymmPackMatrix N = neq.getMatrix();
				Vector dx = n;

				estimateCompleteModel = isEstimated || restrictions.length == 0;
				try {
					// Solve Nx = n in-place, i.e., n <-- dx 
					MathExtension.solve(N, n, !includeRestrictions ? numberOfUnknownParameters : n.size(), false);
					NormalEquationSystem.applyPrecondition(neq.getPreconditioner(), null, dx);

					n = null;
					N = null;
				}
				catch (Exception e) {
					e.printStackTrace();
					return false;
				}

				maxAbsDx = updateUnknownParameters(coefficients, dx);
				dx = null;
				includeRestrictions = true;

				if (Double.isInfinite(maxAbsDx) || Double.isNaN(maxAbsDx)) {
					return false;
				}
				else if (maxAbsDx <= SQRT_EPS && runs > 0 && (includeRestrictions || restrictions.length == 0)) {
					isEstimated = true;
				}
				else if (runs-- <= 1) {
					if (estimateCompleteModel)
						isConverge = false;
					isEstimated = true;
				}
			}
			while (!estimateCompleteModel);
		}
		catch (OutOfMemoryError e) {
			e.printStackTrace();
			return false;
		}
		
		expandUnknownParameters(coefficients, scale);
		return isConverge;
	}
	
	private static double updateUnknownParameters(DLTCoefficients coefficients, Vector dx) {
		double maxAbsDx = 0;
		for (UnknownParameter<?> unknownParameter : coefficients) {
			int column = unknownParameter.getColumn();
			if (column >= 0 && column < Integer.MAX_VALUE) {
				double value  = unknownParameter.getValue();
				double dvalue = dx.get(column);
				maxAbsDx = Math.max(Math.abs(dvalue), maxAbsDx);
				unknownParameter.setValue(value + dvalue);
			}
		}
		return maxAbsDx;
	}
	
	private static void expandUnknownParameters(DLTCoefficients coefficients, double scale) {
		for (UnknownParameter<?> unknownParameter : coefficients) {
			if (unknownParameter.getParameterType() == ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14 || unknownParameter.getParameterType() == ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24)
				continue;
			
			double value = unknownParameter.getValue();
			unknownParameter.setValue(value / scale);
		}
		
		double b11 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11).getValue();
		double b12 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12).getValue();
		double b13 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13).getValue();
		double b14 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14).getValue();
		
		double b21 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21).getValue();
		double b22 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22).getValue();
		double b23 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23).getValue();
		double b24 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24).getValue();
		
		double b31 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31).getValue();
		double b32 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32).getValue();
		double b33 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33).getValue();
		
		double b = b31*b31 + b32*b32 + b33*b33;
		
		double x0 = (b11*b31 + b12*b32 + b13*b33)/b;
		double y0 = (b21*b31 + b22*b32 + b23*b33)/b;

		double cx = Math.sqrt((b11*b11 + b12*b12 + b13*b13)/b - x0*x0);
		double cy = Math.sqrt((b21*b21 + b22*b22 + b23*b23)/b - y0*y0);

		double r11 = -(x0*b31 - b11)/Math.sqrt(b)/cx;
		double r12 = -(y0*b31 - b21)/Math.sqrt(b)/cy;
		double r13 = -b31/Math.sqrt(b);

		double r21 = -(x0*b32 - b12)/Math.sqrt(b)/cx;
		double r22 = -(y0*b32 - b22)/Math.sqrt(b)/cy;
		double r23 = -b32/Math.sqrt(b);

		double r31 = -(x0*b33 - b13)/Math.sqrt(b)/cx;
		double r32 = -(y0*b33 - b23)/Math.sqrt(b)/cy;
		double r33 = -b33/Math.sqrt(b);
		
		double detR = r11*r22*r33 + r12*r23*r31 + r13*r21*r32 - r13*r22*r31 - r11*r23*r32 - r12*r21*r33;

		if (detR < 0) {
			r11 = -r11;
			r12 = -r12;
			r13 = -r13;
			
			r21 = -r21;
			r22 = -r22;
			r23 = -r23;
		    
			r31 = -r31;
			r32 = -r32;
			r33 = -r33;
		}
		
		double omega = Math.atan2(-r23, r33);
		double phi   = Math.asin(r13);
		double kappa = Math.atan2(-r12, r11);
		
		DenseMatrix F = new DenseMatrix(new double[][] {{b11, b12, b13}, {b21, b22, b23}, {b31, b32, b33}});
		DenseVector f = new DenseVector(new double[] {-b14, -b24, -1.0});
		DenseVector t = new DenseVector(f.size());
		F.solve(f, t);
		
		if (coefficients.get(ParameterType.PRINCIPAL_DISTANCE).getColumn() != Integer.MAX_VALUE)
			coefficients.get(ParameterType.PRINCIPAL_DISTANCE).setValue(0.5 * (cx + cy));
		if (coefficients.get(ParameterType.PRINCIPAL_POINT_X).getColumn() != Integer.MAX_VALUE)
			coefficients.get(ParameterType.PRINCIPAL_POINT_X).setValue(x0);
		if (coefficients.get(ParameterType.PRINCIPAL_POINT_Y).getColumn() != Integer.MAX_VALUE)
			coefficients.get(ParameterType.PRINCIPAL_POINT_Y).setValue(y0);

		coefficients.get(ParameterType.CAMERA_COORDINATE_X).setValue(t.get(0));
		coefficients.get(ParameterType.CAMERA_COORDINATE_Y).setValue(t.get(1));
		coefficients.get(ParameterType.CAMERA_COORDINATE_Z).setValue(t.get(2));
		
		coefficients.get(ParameterType.CAMERA_OMEGA).setValue(omega);
		coefficients.get(ParameterType.CAMERA_PHI).setValue(phi);
		coefficients.get(ParameterType.CAMERA_KAPPA).setValue(kappa);
	}
	
	private static RestrictionType[] validateRestrictions(RestrictionType ...restrictions) {
		// avoid over-constrained system
		Set<RestrictionType> restrictionSet = new LinkedHashSet<RestrictionType>(Arrays.asList(restrictions));
		if (restrictionSet.contains(RestrictionType.FIXED_PRINCIPLE_DISTANCE_X) && 
				restrictionSet.contains(RestrictionType.FIXED_PRINCIPLE_DISTANCE_Y) &&
				restrictionSet.contains(RestrictionType.IDENTICAL_PRINCIPLE_DISTANCE)) {
			restrictionSet.remove(RestrictionType.IDENTICAL_PRINCIPLE_DISTANCE);
		}
		return restrictionSet.toArray(new RestrictionType[restrictionSet.size()]);
	}
	
	private static void prepareUnknwonParameters(DLTCoefficients coefficients) {
		int column = 0;
		for (UnknownParameter<?> unknownParameter : coefficients) {
			unknownParameter.setValue(0);
			
			switch(unknownParameter.getParameterType()) {
			case DIRECT_LINEAR_TRANSFORMATION_B11:
			case DIRECT_LINEAR_TRANSFORMATION_B12:
			case DIRECT_LINEAR_TRANSFORMATION_B13:
			case DIRECT_LINEAR_TRANSFORMATION_B14:
			case DIRECT_LINEAR_TRANSFORMATION_B21:
			case DIRECT_LINEAR_TRANSFORMATION_B22:
			case DIRECT_LINEAR_TRANSFORMATION_B23:
			case DIRECT_LINEAR_TRANSFORMATION_B24:
			case DIRECT_LINEAR_TRANSFORMATION_B31:
			case DIRECT_LINEAR_TRANSFORMATION_B32:
			case DIRECT_LINEAR_TRANSFORMATION_B33:
				unknownParameter.setColumn(column++);
				break;
			default:
				unknownParameter.setColumn(unknownParameter.getColumn() == Integer.MAX_VALUE ? Integer.MAX_VALUE : -1);
				break;
			}
		}
		Image image = coefficients.getReference();
		Camera camera = image.getReference();
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		UnknownParameter<?> interiorOrientationParam[] = new UnknownParameter<?>[] {
			interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE),
			interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X),
			interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y)
		};
		
		for (UnknownParameter<?> unknownParameter : interiorOrientationParam) {
			coefficients.get(unknownParameter.getParameterType()).setValue(unknownParameter.getValue());
			if (unknownParameter.getColumn() == Integer.MAX_VALUE)
				coefficients.get(unknownParameter.getParameterType()).setColumn(Integer.MAX_VALUE);
		}
	}

	private static NormalEquationSystem createNormalEquation(DLTCoefficients coefficients, List<ImageCoordinate> homologousImageCoordinates, Map<String, ObjectCoordinate> objectCoordinates, double scale, RestrictionType ...restrictions) {
		int numberOfRestrictions = restrictions.length;
		
		UpperSymmPackMatrix N = new UpperSymmPackMatrix(numberOfUnknownParameters + numberOfRestrictions);
		UpperSymmBandMatrix V = new UpperSymmBandMatrix( N.numRows(), 0 );
		DenseVector n         = new DenseVector(N.numRows());
		
		for (int i = 0; i < homologousImageCoordinates.size(); i++) {
			ImageCoordinate imageCoordinate = homologousImageCoordinates.get(i);
			String name = imageCoordinate.getObjectCoordinate().getName();
			ObjectCoordinate objectCoordinate = objectCoordinates.get(name);

			double x = imageCoordinate.getX().getValue();
			double y = imageCoordinate.getY().getValue();

			double X = objectCoordinate.getX().getValue();
			double Y = objectCoordinate.getY().getValue();
			double Z = objectCoordinate.getZ().getValue();

			DLTPartialDerivativeFactory.addPartialNormalEquationOfDLTParameters(N, n, coefficients, x, y, X/scale, Y/scale, Z/scale);
		}
		
		if (numberOfRestrictions > 0)
			DLTPartialDerivativeFactory.setParameterRestrictions(N, n, numberOfUnknownParameters, coefficients, restrictions);
		
		for (int column = 0; column < N.numColumns(); column++) {
			double value = N.get(column, column);
			V.set(column, column, value > Constant.EPS ? 1.0 / Math.sqrt(value) : 1.0);
		}
		
		return new NormalEquationSystem(N, n, V);
	}
}
