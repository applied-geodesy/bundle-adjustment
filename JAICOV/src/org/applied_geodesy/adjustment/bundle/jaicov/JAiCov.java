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

package org.applied_geodesy.adjustment.bundle.jaicov;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Collection;
import java.util.Locale;
import java.util.Set;

import org.applied_geodesy.adjustment.EstimationStateType;
import org.applied_geodesy.adjustment.EstimationType;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment.MatrixInversion;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ImageCoordinate;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.util.io.AICONReportFileReader;

import no.uib.cipr.matrix.Matrix;

public class JAiCov implements PropertyChangeListener {

	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		System.out.println("Info: " + evt.getPropertyName() + " " +evt.getOldValue() + " --> " + evt.getNewValue());
	}

	public static void main(String[] args) throws Exception {
		// Prevent warning of (missing) native implementation
		System.setProperty("com.github.fommil.netlib.BLAS",   "com.github.fommil.netlib.F2jBLAS");
		System.setProperty("com.github.fommil.netlib.LAPACK", "com.github.fommil.netlib.F2jLAPACK");
		System.setProperty("com.github.fommil.netlib.ARPACK", "com.github.fommil.netlib.F2jARPACK");

		long t = System.currentTimeMillis();
		MatrixInversion estimateDispersionMatrix = MatrixInversion.FULL;
		
		// Read adjustment report from Aicon Studio 3D
		AICONReportFileReader reader = new AICONReportFileReader("example/example.htm");

		// Create an adjustment object using a specific file reader
		BundleAdjustment adjustment = reader.readAndImport();
		
		// Get all camera objects from reader
		Collection<Camera> cameras = reader.getCameras();

		// Specify the points, which are used to define the datum of the frame
		// !!!this step is not mandatory!!!
		for (Camera camera : cameras) {
			// Call the images taken from this camera object
			for (Image image : camera) {
				// Extract the image coordinates
				for (ImageCoordinate imageCoordinate : image) {
					// Select corresponding object coordinates
					ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
					if (objectCoordinate.getName().length() > 3)
						objectCoordinate.setDatum(Boolean.FALSE);
				}
			}
		}

		// Add a listener to get notifications of the adjustment process
		adjustment.addPropertyChangeListener(new JAiCov());
		// Select the estimation type, e.g., L2Norm or Simulation
		adjustment.setEstimationType(EstimationType.L2NORM);
		// Use NONE, if the dispersion matrix is not required
		adjustment.setInvertNormalEquation(estimateDispersionMatrix);

		// Call estimate to start the bundle adjustment
		EstimationStateType estimationStateType = adjustment.estimateModel();

		// Check the result
		if (estimationStateType != EstimationStateType.ERROR_FREE_ESTIMATION) {
			System.err.println("Error, bundle adjustment fails...");
		}
		else {
			System.out.println("Bundle adjustment finished successfully...");

			// derive dispersion of parameters
			Matrix D = adjustment.getCofactorMatrix().scale(adjustment.getVarianceFactorAposteriori());
			String template = "%10s\t%+16.5f\t%+16.5f\t%+16.5f\t%+8.5f\t%+8.5f\t%+8.5f\t%1s";

			// print coordinates of object points and related uncertainties
			Set<ObjectCoordinate> objectCoordinates = adjustment.getObjectCoordinates();
			for (ObjectCoordinate objectCoordinate : objectCoordinates) {
				UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
				UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
				UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();

				// indicate datum (d) and object points (o)
				char datum = objectCoordinate.isDatum() ? 'd' : 'o';

				double x = X.getValue();
				double y = Y.getValue();
				double z = Z.getValue();

				double ux = 0, uy = 0, uz = 0;

				if (estimateDispersionMatrix != MatrixInversion.NONE && X.getColumn() >= 0 && Y.getColumn() >= 0 && Z.getColumn() >= 0 && X.getColumn() != Integer.MAX_VALUE && Y.getColumn() != Integer.MAX_VALUE && Z.getColumn() != Integer.MAX_VALUE) {
					ux = Math.sqrt(Math.abs(D.get(X.getColumn(), X.getColumn())));
					uy = Math.sqrt(Math.abs(D.get(Y.getColumn(), Y.getColumn())));
					uz = Math.sqrt(Math.abs(D.get(Z.getColumn(), Z.getColumn())));
				}
				String formattedStr = String.format(Locale.ENGLISH, template, objectCoordinate.getName(), x, y, z, ux, uy, uz, datum);
				System.out.println(formattedStr);
			}
			
			for (Camera camera : cameras) {
				InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
				for (UnknownParameter<InteriorOrientation> interiorOrientationParameter : interiorOrientation) {
					System.out.println(String.format(Locale.ENGLISH, "%-25s = %+12.8f %s", interiorOrientationParameter.getParameterType().name(), interiorOrientationParameter.getValue(), interiorOrientationParameter.getColumn() == Integer.MAX_VALUE  ? "fixed" : ""));
				}
				System.out.println();
			}

			// print some statistical parameters
			System.out.println("Number of observations:           " + adjustment.getNumberOfObservations());
			System.out.println("Number of unknown parameters:     " + adjustment.getNumberOfUnknownParameters());
			System.out.println("Degree of freedom:                " + adjustment.getDegreeOfFreedom());
			System.out.println("Variances of unit weight:         1.0 : " + adjustment.getVarianceFactorAposteriori() / adjustment.getVarianceFactorApriori());
			System.out.println("Variances of unit weight (ratio): " + adjustment.getVarianceFactorApriori() + " : " + adjustment.getVarianceFactorAposteriori());
			System.out.println("Estimation time:                  " + ((System.currentTimeMillis() - t)/1000.0) + " sec");	
		}
	}
}
