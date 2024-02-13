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
import java.util.Collections;
import java.util.List;

import org.applied_geodesy.adjustment.bundle.ImageCoordinate;
import org.applied_geodesy.adjustment.bundle.dlt.DirectLinearTransformation.RestrictionType;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.UpperSymmPackMatrix;

class DLTPartialDerivativeFactory {
	private DLTPartialDerivativeFactory() {}
	
	static double getMisclosure(ParameterType observationParamType, DLTCoefficients coefficients, double x, double y, double X, double Y, double Z) {
		if (observationParamType == ParameterType.IMAGE_COORDINATE_X)
			return x - getCollinearityEquationValue(observationParamType, coefficients, x, y, X, Y, Z);
		else if (observationParamType == ParameterType.IMAGE_COORDINATE_Y)
			return y - getCollinearityEquationValue(observationParamType, coefficients, x, y, X, Y, Z);
		return 0;
	}
	
	public static double getCollinearityEquationValue(ParameterType observationParamType, DLTCoefficients coefficients, double x, double y, double X, double Y, double Z) {
		UnknownParameter<DLTCoefficients> B11 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11);
		UnknownParameter<DLTCoefficients> B12 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12);
		UnknownParameter<DLTCoefficients> B13 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13);
		UnknownParameter<DLTCoefficients> B14 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14);
		
		UnknownParameter<DLTCoefficients> B21 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21);
		UnknownParameter<DLTCoefficients> B22 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22);
		UnknownParameter<DLTCoefficients> B23 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23);
		UnknownParameter<DLTCoefficients> B24 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24);
		
		UnknownParameter<DLTCoefficients> B31 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31);
		UnknownParameter<DLTCoefficients> B32 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32);
		UnknownParameter<DLTCoefficients> B33 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33);
		
		
		double b11 = B11.getValue();
		double b12 = B12.getValue();
		double b13 = B13.getValue();
		double b14 = B14.getValue();
		
		double b21 = B21.getValue();
		double b22 = B22.getValue();
		double b23 = B23.getValue();
		double b24 = B24.getValue();
		
		double b31 = B31.getValue();
		double b32 = B32.getValue();
		double b33 = B33.getValue();
        
        if (observationParamType == ParameterType.IMAGE_COORDINATE_X)
			return X*b11 + Y*b12 + Z*b13 + b14 - x*X*b31 - x*Y*b32 - x*Z*b33;
		else if (observationParamType == ParameterType.IMAGE_COORDINATE_Y)
			return X*b21 + Y*b22 + Z*b23 + b24 - y*X*b31 - y*Y*b32 - y*Z*b33;
		return 0;
	}
			

	static void setParameterRestrictions(UpperSymmPackMatrix NEQ, DenseVector neq, int rowIndex, DLTCoefficients coefficients, RestrictionType ...restrictions) {
		UnknownParameter<DLTCoefficients> B11 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11);
		UnknownParameter<DLTCoefficients> B12 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12);
		UnknownParameter<DLTCoefficients> B13 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13);
		UnknownParameter<DLTCoefficients> B14 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14);
		
		UnknownParameter<DLTCoefficients> B21 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21);
		UnknownParameter<DLTCoefficients> B22 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22);
		UnknownParameter<DLTCoefficients> B23 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23);
		UnknownParameter<DLTCoefficients> B24 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24);
		
		UnknownParameter<DLTCoefficients> B31 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31);
		UnknownParameter<DLTCoefficients> B32 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32);
		UnknownParameter<DLTCoefficients> B33 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33);
		
		double c  = coefficients.get(ParameterType.PRINCIPAL_DISTANCE).getValue();
		double x0 = coefficients.get(ParameterType.PRINCIPAL_POINT_X).getValue();
		double y0 = coefficients.get(ParameterType.PRINCIPAL_POINT_Y).getValue();
		
		
		double b11 = B11.getValue();
		double b12 = B12.getValue();
		double b13 = B13.getValue();
		// double b14 = B14.getValue();
		double b21 = B21.getValue();
		double b22 = B22.getValue();
		double b23 = B23.getValue();
		// double b24 = B24.getValue();
		double b31 = B31.getValue();
		double b32 = B32.getValue();
		double b33 = B33.getValue();
		
		double sb11 = b11*b11;
		double sb12 = b12*b12;
		double sb13 = b13*b13;
		
		double sb21 = b21*b21;
		double sb22 = b22*b22;
		double sb23 = b23*b23;
		
		double sb31 = b31*b31;
		double sb32 = b32*b32;
		double sb33 = b33*b33;
		
		double b1 = sb11 + sb12 + sb13;
		double b2 = sb21 + sb22 + sb23;
		double b3 = sb31 + sb32 + sb33;
		double bx = b11*b31 + b12*b32 + b13*b33;
		double by = b21*b31 + b22*b32 + b23*b33;
		
		for (RestrictionType restriction : restrictions) {
			switch (restriction) {
			case FIXED_PRINCIPAL_POINT_X:
				NEQ.set(B11.getColumn(), rowIndex, b31/b3);
				NEQ.set(B12.getColumn(), rowIndex, b32/b3);
				NEQ.set(B13.getColumn(), rowIndex, b33/b3);
				NEQ.set(B14.getColumn(), rowIndex, 0);
				NEQ.set(B21.getColumn(), rowIndex, 0);
				NEQ.set(B22.getColumn(), rowIndex, 0);
				NEQ.set(B23.getColumn(), rowIndex, 0);
				NEQ.set(B24.getColumn(), rowIndex, 0);
				NEQ.set(B31.getColumn(), rowIndex, -(2.0*b31*(b12*b32 + b13*b33) + b11*(sb31 - sb32 - sb33))/b3/b3);
				NEQ.set(B32.getColumn(), rowIndex, -(2.0*b32*(b11*b31 + b13*b33) + b12*(sb32 - sb31 - sb33))/b3/b3);
				NEQ.set(B33.getColumn(), rowIndex, -(2.0*b33*(b11*b31 + b12*b32) + b13*(sb33 - sb32 - sb31))/b3/b3);
				
				neq.set(rowIndex++, x0 - bx/b3);
				
				break;
			case FIXED_PRINCIPAL_POINT_Y:
				NEQ.set(B11.getColumn(), rowIndex, 0);
				NEQ.set(B12.getColumn(), rowIndex, 0);
				NEQ.set(B13.getColumn(), rowIndex, 0);
				NEQ.set(B14.getColumn(), rowIndex, 0);
				NEQ.set(B21.getColumn(), rowIndex, b31/b3);
				NEQ.set(B22.getColumn(), rowIndex, b32/b3);
				NEQ.set(B23.getColumn(), rowIndex, b33/b3);
				NEQ.set(B24.getColumn(), rowIndex, 0);
				NEQ.set(B31.getColumn(), rowIndex, -(2.0*b31*(b22*b32 + b23*b33) + b21*(sb31 - sb32 - sb33))/b3/b3);
				NEQ.set(B32.getColumn(), rowIndex, -(2.0*b32*(b21*b31 + b23*b33) + b22*(sb32 - sb31 - sb33))/b3/b3);
				NEQ.set(B33.getColumn(), rowIndex, -(2.0*b33*(b21*b31 + b22*b32) + b23*(sb33 - sb32 - sb31))/b3/b3);
				
				neq.set(rowIndex++, y0 - by/b3);
				
				break;
			case FIXED_PRINCIPLE_DISTANCE_X:
				NEQ.set(B11.getColumn(), rowIndex, 2.0*(b11*(sb32 + sb33) - b31*(b12*b32 + b13*b33))/b3/b3);
				NEQ.set(B12.getColumn(), rowIndex, 2.0*(b12*(sb31 + sb33) - b32*(b11*b31 + b13*b33))/b3/b3);
				NEQ.set(B13.getColumn(), rowIndex, 2.0*(b13*(sb31 + sb32) - b33*(b11*b31 + b12*b32))/b3/b3);
				NEQ.set(B14.getColumn(), rowIndex, 0);
				NEQ.set(B21.getColumn(), rowIndex, 0);
				NEQ.set(B22.getColumn(), rowIndex, 0);
				NEQ.set(B23.getColumn(), rowIndex, 0);
				NEQ.set(B24.getColumn(), rowIndex, 0);
				NEQ.set(B31.getColumn(), rowIndex, 4.0*(b31*bx*bx - 0.5*b3*(b31*b1 + bx*b11))/(b3*b3*b3));
				NEQ.set(B32.getColumn(), rowIndex, 4.0*(b32*bx*bx - 0.5*b3*(b32*b1 + bx*b12))/(b3*b3*b3));
				NEQ.set(B33.getColumn(), rowIndex, 4.0*(b33*bx*bx - 0.5*b3*(b33*b1 + bx*b13))/(b3*b3*b3));

				neq.set(rowIndex++, c*c - b1/b3 + bx*bx/b3/b3);
				
				break;
			case FIXED_PRINCIPLE_DISTANCE_Y:
				NEQ.set(B11.getColumn(), rowIndex, 0);
				NEQ.set(B12.getColumn(), rowIndex, 0);
				NEQ.set(B13.getColumn(), rowIndex, 0);
				NEQ.set(B14.getColumn(), rowIndex, 0);
				NEQ.set(B21.getColumn(), rowIndex, 2.0*(b21*(sb32 + sb33) - b31*(b22*b32 + b23*b33))/b3/b3);
				NEQ.set(B22.getColumn(), rowIndex, 2.0*(b22*(sb31 + sb33) - b32*(b21*b31 + b23*b33))/b3/b3);
				NEQ.set(B23.getColumn(), rowIndex, 2.0*(b23*(sb31 + sb32) - b33*(b21*b31 + b22*b32))/b3/b3);
				NEQ.set(B24.getColumn(), rowIndex, 0);
				NEQ.set(B31.getColumn(), rowIndex, 4.0*(b31*by*by - 0.5*b3*(b31*b2 + by*b21))/(b3*b3*b3));
				NEQ.set(B32.getColumn(), rowIndex, 4.0*(b32*by*by - 0.5*b3*(b32*b2 + by*b22))/(b3*b3*b3));
				NEQ.set(B33.getColumn(), rowIndex, 4.0*(b33*by*by - 0.5*b3*(b33*b2 + by*b23))/(b3*b3*b3));

				neq.set(rowIndex++, c*c - b2/b3 + by*by/b3/b3);
				
				break;
			case IDENTICAL_PRINCIPLE_DISTANCE:
				NEQ.set(B11.getColumn(), rowIndex,  2.0*(b11*sb32 - b12*b31*b32 + b11*sb33 - b13*b31*b33));
				NEQ.set(B12.getColumn(), rowIndex,  2.0*(b12*sb31 - b11*b32*b31 + b12*sb33 - b13*b32*b33));
				NEQ.set(B13.getColumn(), rowIndex,  2.0*(b13*sb31 - b11*b33*b31 + b13*sb32 - b12*b33*b32));
				NEQ.set(B14.getColumn(), rowIndex,  0.0);
				NEQ.set(B21.getColumn(), rowIndex, -2.0*(b21*sb32 - b22*b31*b32 + b21*sb33 - b23*b31*b33));
				NEQ.set(B22.getColumn(), rowIndex, -2.0*(b22*sb31 - b21*b32*b31 + b22*sb33 - b23*b32*b33));
				NEQ.set(B23.getColumn(), rowIndex, -2.0*(b23*sb31 - b21*b33*b31 + b23*sb32 - b22*b33*b32));
				NEQ.set(B24.getColumn(), rowIndex,  0.0);
				NEQ.set(B31.getColumn(), rowIndex,  2.0*(b31*sb12 - b11*b32*b12 + b31*sb13 - b11*b33*b13 - b31*sb22 + b21*b32*b22 - b31*sb23 + b21*b33*b23));
				NEQ.set(B32.getColumn(), rowIndex,  2.0*(b32*sb11 - b12*b31*b11 + b32*sb13 - b12*b33*b13 - b32*sb21 + b22*b31*b21 - b32*sb23 + b22*b33*b23));
				NEQ.set(B33.getColumn(), rowIndex,  2.0*(b33*sb11 - b13*b31*b11 + b33*sb12 - b13*b32*b12 - b33*sb21 + b23*b31*b21 - b33*sb22 + b23*b32*b22));

				neq.set(rowIndex++, -b3 * (b1 - b2) + bx*bx - by*by);
				
				break;
			case ROTATION_WITHOUT_SHEAR:
				NEQ.set(B11.getColumn(), rowIndex, -b21*sb32 + b22*b31*b32 - b21*sb33 + b23*b31*b33);
				NEQ.set(B12.getColumn(), rowIndex, -b22*sb31 + b21*b32*b31 - b22*sb33 + b23*b32*b33);
				NEQ.set(B13.getColumn(), rowIndex, -b23*sb31 + b21*b33*b31 - b23*sb32 + b22*b33*b32);
				NEQ.set(B14.getColumn(), rowIndex, 0.0);
				NEQ.set(B21.getColumn(), rowIndex, -b11*sb32 + b12*b31*b32 - b11*sb33 + b13*b31*b33);
				NEQ.set(B22.getColumn(), rowIndex, -b12*sb31 + b11*b32*b31 - b12*sb33 + b13*b32*b33);
				NEQ.set(B23.getColumn(), rowIndex, -b13*sb31 + b11*b33*b31 - b13*sb32 + b12*b33*b32);
				NEQ.set(B24.getColumn(), rowIndex, 0.0);
				NEQ.set(B31.getColumn(), rowIndex, b11*b22*b32 + b12*b21*b32 - 2.0*b12*b22*b31 + b11*b23*b33 + b13*b21*b33 - 2.0*b13*b23*b31);
				NEQ.set(B32.getColumn(), rowIndex, b11*b22*b31 - 2.0*b11*b21*b32 + b12*b21*b31 + b12*b23*b33 + b13*b22*b33 - 2.0*b13*b23*b32);
				NEQ.set(B33.getColumn(), rowIndex, b11*b23*b31 - 2.0*b11*b21*b33 + b13*b21*b31 - 2.0*b12*b22*b33 + b12*b23*b32 + b13*b22*b32);

				neq.set(rowIndex++, b3 * (b11*b21 + b12*b22 + b13*b23) - bx * by);
				
				break;			
			}
		}
	}
	
	static void addPartialNormalEquationOfDLTParameters(UpperSymmPackMatrix NEQ, DenseVector neq, DLTCoefficients coefficients, double x, double y, double X, double Y, double Z) {
		UnknownParameter<DLTCoefficients> B11 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11);
		UnknownParameter<DLTCoefficients> B12 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12);
		UnknownParameter<DLTCoefficients> B13 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13);
		UnknownParameter<DLTCoefficients> B14 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14);
		
		UnknownParameter<DLTCoefficients> B21 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21);
		UnknownParameter<DLTCoefficients> B22 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22);
		UnknownParameter<DLTCoefficients> B23 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23);
		UnknownParameter<DLTCoefficients> B24 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24);
		
		UnknownParameter<DLTCoefficients> B31 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31);
		UnknownParameter<DLTCoefficients> B32 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32);
		UnknownParameter<DLTCoefficients> B33 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33);
		
		
		double b11 = B11.getValue();
		double b12 = B12.getValue();
		double b13 = B13.getValue();
		double b14 = B14.getValue();
		
		double b21 = B21.getValue();
		double b22 = B22.getValue();
		double b23 = B23.getValue();
		double b24 = B24.getValue();
		
		double b31 = B31.getValue();
		double b32 = B32.getValue();
		double b33 = B33.getValue();
	
		DenseVector w = new DenseVector(2);
		DenseMatrix A = new DenseMatrix(w.size(), neq.size());
		
		List<Integer> columns = new ArrayList<Integer>(11);
		int column = -1;

		// DLT B1X
		column = B11.getColumn();
		columns.add(column);
		A.set(0, column, X);

		column = B12.getColumn();
		columns.add(column);
		A.set(0, column, Y);

		column = B13.getColumn();
		columns.add(column);
		A.set(0, column, Z);

		column = B14.getColumn();
		columns.add(column);
		A.set(0, column, 1.0);

		// DLT B2X
		column = B21.getColumn();
		columns.add(column);
		A.set(1, column, X);

		column = B22.getColumn();
		columns.add(column);
		A.set(1, column, Y);

		column = B23.getColumn();
		columns.add(column);
		A.set(1, column, Z);

		column = B24.getColumn();
		columns.add(column);
		A.set(1, column, 1.0);

		// DLT B3X
		column = B31.getColumn();
		columns.add(column);
		A.set(0, column, -x*X);
		A.set(1, column, -y*X);

		column = B32.getColumn();
		columns.add(column);
		A.set(0, column, -x*Y);
		A.set(1, column, -y*Y);

		column = B33.getColumn();
		columns.add(column);
		A.set(0, column, -x*Z);
		A.set(1, column, -y*Z);


		w.set(0, x - (X*b11 + Y*b12 + Z*b13 + b14 - x*X*b31 - x*Y*b32 - x*Z*b33));
        w.set(1, y - (X*b21 + Y*b22 + Z*b23 + b24 - y*X*b31 - y*Y*b32 - y*Z*b33));
        
        Collections.sort(columns);
        for (int row = 0; row < A.numRows(); row++) {
        	for (int columnATIdx = 0; columnATIdx < columns.size(); columnATIdx++) {
        		int colAT = columns.get(columnATIdx);
        		double aT = A.get(row, colAT);
        		neq.add(colAT, aT * w.get(row));

        		if (NEQ != null) {
        			for (int columnAIdx = columnATIdx; columnAIdx < columns.size(); columnAIdx++) {
        				int colA = columns.get(columnAIdx);	
        				double a = A.get(row, colA);
        				NEQ.add(colAT, colA, aT * a);
        			}
        		}
        	}	
		}
	}
	
	public static void addPartialNormalEquationOfUnknownPosition(UpperSymmPackMatrix NEQ, DenseVector neq, ImageCoordinate imageCoordinate, DLTCoefficients coefficients) {
		double x = imageCoordinate.getX().getValue();
		double y = imageCoordinate.getY().getValue();
		
		UnknownParameter<DLTCoefficients> B11 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B11);
		UnknownParameter<DLTCoefficients> B12 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B12);
		UnknownParameter<DLTCoefficients> B13 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B13);
		UnknownParameter<DLTCoefficients> B14 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B14);
		
		UnknownParameter<DLTCoefficients> B21 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B21);
		UnknownParameter<DLTCoefficients> B22 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B22);
		UnknownParameter<DLTCoefficients> B23 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B23);
		UnknownParameter<DLTCoefficients> B24 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B24);
		
		UnknownParameter<DLTCoefficients> B31 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B31);
		UnknownParameter<DLTCoefficients> B32 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B32);
		UnknownParameter<DLTCoefficients> B33 = coefficients.get(ParameterType.DIRECT_LINEAR_TRANSFORMATION_B33);

		double b11 = B11.getValue();
		double b12 = B12.getValue();
		double b13 = B13.getValue();
		double b14 = B14.getValue();
		
		double b21 = B21.getValue();
		double b22 = B22.getValue();
		double b23 = B23.getValue();
		double b24 = B24.getValue();
		
		double b31 = B31.getValue();
		double b32 = B32.getValue();
		double b33 = B33.getValue();
	
		DenseVector w = new DenseVector(2);
		DenseMatrix A = new DenseMatrix(w.size(), 3);
		
		A.set(0, 0, b11 - x * b31);
		A.set(0, 1, b12 - x * b32);
		A.set(0, 2, b13 - x * b33);
		
		A.set(1, 0, b21 - y * b31);
		A.set(1, 1, b22 - y * b32);
		A.set(1, 2, b23 - y * b33);
		
		w.set(0, -b14 + x);
		w.set(1, -b24 + y);
		
		for (int row = 0; row < A.numRows(); row++) {
        	for (int colAT = 0; colAT < A.numColumns(); colAT++) {
        		double aT = A.get(row, colAT);
        		neq.add(colAT, aT * w.get(row));

        		if (NEQ != null) {
        			for (int colA = colAT; colA < A.numColumns(); colA++) {
        				double a = A.get(row, colA);
        				NEQ.add(colAT, colA, aT * a);
        			}
        		}
        	}	
		}
	}
}
