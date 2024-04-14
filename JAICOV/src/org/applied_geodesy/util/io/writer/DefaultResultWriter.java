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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;

import no.uib.cipr.matrix.UpperSymmPackMatrix;

public class DefaultResultWriter extends BundleAdjustmentResultWriter {

	public DefaultResultWriter(String exportPathAndFileBaseName) {
		super(exportPathAndFileBaseName);
	}

	@Override
	public void export(BundleAdjustment bundleAdjustment) throws NullPointerException, IllegalArgumentException, IOException {
		if (bundleAdjustment == null)
			throw new NullPointerException("Error, bundle adjustment object cannot be null!");
		
		String exportPathAndFileBaseName = this.getExportPathAndFileBaseName();
		
		if (exportPathAndFileBaseName == null)
			throw new NullPointerException("Error, export path cannot be null!");

		List<Integer> indices = this.exportCovarianceInformation(bundleAdjustment, new File(exportPathAndFileBaseName + ".info"));
		this.exportCovarianceMatrix(bundleAdjustment, indices, new File(exportPathAndFileBaseName + ".cxx"));
	}

	/**
	 * Exports coordinates of object points and related index in dispersion matrix
	 * @param exportPathAndFileBaseName
	 * @return isWritten
	 * @throws IOException 
	 */
	private List<Integer> exportCovarianceInformation(BundleAdjustment bundleAdjustment, File file) throws IOException {
		Collection<ObjectCoordinate> objectCoordinates = bundleAdjustment.getObjectCoordinates();
		
		List<Integer> indices = new ArrayList<Integer>(objectCoordinates.size() * 3);
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			//Name, Type(XYZ), Coordinate, Row/Column in NES
			String format = "%25s\t%5s\t%35.15f\t%10d%n";
			int columnIndex = 0;
			for (ObjectCoordinate objectCoordinate : objectCoordinates) {
				UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
				UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
				UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();
				String name = objectCoordinate.getName();
				
				int columnX = X.getColumn();
				int columnY = Y.getColumn();
				int columnZ = Z.getColumn();
				
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

				pw.printf(Locale.ENGLISH, format, name, 'X', X.getValue(), columnX);
				pw.printf(Locale.ENGLISH, format, name, 'Y', Y.getValue(), columnY);
				pw.printf(Locale.ENGLISH, format, name, 'Z', Z.getValue(), columnZ);
			}
		}
		finally {
			if (pw != null) {
				pw.close();
			}
		}
		
		return indices;
	}

	/**
	 * Exports dispersion matrix
	 * @param exportPathAndFileBaseName
	 * @return isWritten
	 */
	private void exportCovarianceMatrix(BundleAdjustment bundleAdjustment, List<Integer> indices, File file) throws IOException {
		int numberOfDatumConditions   = bundleAdjustment.getNumberOfDatumConditions();
		int numberOfUnknownParameters = bundleAdjustment.getNumberOfUnknownParameters();
		UpperSymmPackMatrix cofactor  = bundleAdjustment.getCofactorMatrix();
		boolean exportDispersionMatrix = !(cofactor == null || cofactor.numRows() < (numberOfUnknownParameters + numberOfDatumConditions));
		
		if (!exportDispersionMatrix)
			return;

		PrintWriter pw = null;
		double sigma2apost = bundleAdjustment.getVarianceFactorAposteriori();

		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			for (int rowIdx = 0; rowIdx < indices.size(); rowIdx++) {
				int row = indices.get(rowIdx);
				for (int columnIdx = 0; columnIdx < indices.size(); columnIdx++) {
					int column = indices.get(columnIdx);
					double value = cofactor.get(row, column);
					pw.printf(Locale.ENGLISH, "%+35.15f  ", sigma2apost * value);
				}
				pw.println();
			}
		} 
		finally {
			if (pw != null) {
				pw.close();
			}
		}
	}
}
