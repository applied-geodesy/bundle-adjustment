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

package org.applied_geodesy.util.io.writer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.UpperSymmPackMatrix;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.MatlabType;
import us.hebi.matlab.mat.types.Matrix;
import us.hebi.matlab.mat.types.Struct;

public class MatlabResultWriter extends BundleAdjustmentResultWriter {

	public MatlabResultWriter(String exportPathAndFileBaseName) {
		super(exportPathAndFileBaseName);
	}

	@Override
	public void export(BundleAdjustment bundleAdjustment) throws NullPointerException, IllegalArgumentException, IOException {
		if (bundleAdjustment == null)
			throw new NullPointerException("Error, bundle adjustment object cannot be null!");
		
		String exportPathAndFileBaseName = this.getExportPathAndFileBaseName();
		
		if (exportPathAndFileBaseName == null)
			throw new NullPointerException("Error, export path cannot be null!");
		
		File binFile = new File(exportPathAndFileBaseName + ".mat");
		
		int numberOfObservations      = bundleAdjustment.getNumberOfObservations();
		int numberOfDatumConditions   = bundleAdjustment.getNumberOfDatumConditions();
		int numberOfUnknownParameters = bundleAdjustment.getNumberOfUnknownParameters();
		int degreeOfFreedom           = bundleAdjustment.getDegreeOfFreedom();
		
		Collection<ObjectCoordinate> objectCoordinates = bundleAdjustment.getObjectCoordinates();
		Collection<Camera> cameras = bundleAdjustment.getCameras();
		UpperSymmPackMatrix cofactor = bundleAdjustment.getCofactorMatrix();
		
		boolean exportDispersionMatrix = !(cofactor == null || cofactor.numRows() < (numberOfUnknownParameters + numberOfDatumConditions));
		
		double sigma2aprio = bundleAdjustment.getVarianceFactorApriori();
		double sigma2apost = bundleAdjustment.getVarianceFactorAposteriori();

		Struct coordinates          = Mat5.newStruct(1, objectCoordinates.size());
		Struct interiorOrientations = Mat5.newStruct(1, cameras.size() * 13);
		
		List<Integer> indices = exportDispersionMatrix ? new ArrayList<Integer>(objectCoordinates.size() * 3 + cameras.size() * 13) : null;

		int structIndex = 0;
		int columnIndex = 1;
		for (ObjectCoordinate objectCoordinate : objectCoordinates) {
			UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
			UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
			UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();
			String name = objectCoordinate.getName();
			
			int columnX = X.getColumn();
			int columnY = Y.getColumn();
			int columnZ = Z.getColumn();

			coordinates.set("name", structIndex, Mat5.newString(name));
			
			coordinates.set("X", structIndex, newDouble(X.getValue()));
			coordinates.set("Y", structIndex, newDouble(Y.getValue()));
			coordinates.set("Z", structIndex, newDouble(Z.getValue()));
			
			if (columnX >= 0 && columnX < Integer.MAX_VALUE) {
				indices.add(columnX);
				columnX = columnIndex++;
			}
			else
				columnX = -1;
			
			if (columnY >= 0 && columnY < Integer.MAX_VALUE) {
				indices.add(columnY);
				columnY = columnIndex++;
			}
			else
				columnY = -1;
			
			if (columnZ >= 0 && columnZ < Integer.MAX_VALUE) {
				indices.add(columnZ);
				columnZ = columnIndex++;
			}
			else
				columnZ = -1;
			
			coordinates.set("covx", structIndex, newInteger(columnX));
			coordinates.set("covy", structIndex, newInteger(columnY));
			coordinates.set("covz", structIndex, newInteger(columnZ));
			
			structIndex++;
		}
		
		structIndex = 0;
		for (Camera camera : cameras) {
			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			for (UnknownParameter<InteriorOrientation> param : interiorOrientation) {
				interiorOrientations.set("cam_id", structIndex, newLong(camera.getId()));
				interiorOrientations.set("name",   structIndex, Mat5.newString(param.getParameterType().name().toLowerCase()));
				interiorOrientations.set("value",  structIndex, newDouble(param.getValue()));
				if (exportDispersionMatrix) {
					int column = param.getColumn();
					
					if (column >= 0 && column < cofactor.numColumns()) {
						indices.add(column);
						column = columnIndex++;
					}
					else
						column = -1;
					
					interiorOrientations.set("cov", structIndex, newInteger(column));
				}
				structIndex++;
			}
		}
	
		
		MatFile matFile = Mat5.newMatFile();
		
		matFile.addArray("variance_of_unit_weight_prio", newDouble(sigma2aprio));
		matFile.addArray("variance_of_unit_weight_post", newDouble(sigma2apost));
		matFile.addArray("degree_of_freedom",            newInteger(degreeOfFreedom));
		matFile.addArray("number_of_observations",       newInteger(numberOfObservations));
		matFile.addArray("number_of_unknowns",           newInteger(numberOfUnknownParameters));
		
		matFile.addArray("coordinates",           coordinates);
		matFile.addArray("interior_orientations", interiorOrientations);
		
		if (exportDispersionMatrix) {
			Matrix dispersion = Mat5.newMatrix(indices.size(), indices.size(), MatlabType.Double);
			for (int rowIdx = 0; rowIdx < indices.size(); rowIdx++) {
				int row = indices.get(rowIdx);
				double var = cofactor.get(row, row);
				dispersion.setDouble(rowIdx, rowIdx, var);
				for (int columnIdx = rowIdx + 1; columnIdx < indices.size(); columnIdx++) {
					int column = indices.get(columnIdx);
					double covar = cofactor.get(row, column);
					dispersion.setDouble(rowIdx, columnIdx, covar);
					dispersion.setDouble(columnIdx, rowIdx, covar);
				}
			}
			matFile.addArray("dispersion", dispersion);
		}

		Mat5.writeToFile(matFile, binFile);		
	}

	private static Matrix newInteger(int value) {
		Matrix matrix = Mat5.newMatrix(1, 1, MatlabType.Int32);
		matrix.setInt(0, 0, value);
		return matrix;
	}
	
	private static Matrix newLong(long value) {
		Matrix matrix = Mat5.newMatrix(1, 1, MatlabType.Int64);
		matrix.setLong(0, 0, value);
		return matrix;
	}
	
	private static Matrix newDouble(double value) {
		Matrix matrix = Mat5.newMatrix(1, 1, MatlabType.Double);
		matrix.setDouble(0, 0, value);
		return matrix;
	}
}
