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

package org.applied_geodesy.adjustment.bundle.example;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

import org.applied_geodesy.adjustment.EstimationStateType;
import org.applied_geodesy.adjustment.EstimationType;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ImageCoordinate;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.ScaleBar;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment.MatrixInversion;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.DirectlyObservedParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.util.io.reader.aicon.EORFileReader;
import org.applied_geodesy.util.io.reader.aicon.IORFileReader;
import org.applied_geodesy.util.io.reader.aicon.OBCFileReader;
import org.applied_geodesy.util.io.reader.aicon.PHCFileReader;
import org.applied_geodesy.util.io.reader.aicon.ScaleFileReader;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.UpperSPDPackMatrix;

public class ExampleFlatFiles implements PropertyChangeListener {

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
		
		final String basepath = "example/example";
		
		// Read coordinates of object points
		OBCFileReader obcReader = new OBCFileReader(basepath + ".obc");
		Map<String, ObjectCoordinate> coordinates = obcReader.readAndImport();
		
		// Read scalebar definition
		ScaleFileReader scaleReader = new ScaleFileReader(basepath + ".scale", coordinates);
		Collection<ScaleBar> scaleBars = scaleReader.readAndImport();
		
		// Read interior orientation and create camera object
		IORFileReader iorReader = new IORFileReader(basepath + ".ior");
		Camera camera = iorReader.readAndImport();
		
		// Fixing some interior orientation parameters
		camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A3).setColumn(Integer.MAX_VALUE);
		
		camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C1).setColumn(Integer.MAX_VALUE);
		camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C2).setColumn(Integer.MAX_VALUE);
		
		camera.getInteriorOrientation().get(ParameterType.DISTANCE_DISTORTION_D1).setColumn(Integer.MAX_VALUE);
		camera.getInteriorOrientation().get(ParameterType.DISTANCE_DISTORTION_D2).setColumn(Integer.MAX_VALUE);
		camera.getInteriorOrientation().get(ParameterType.DISTANCE_DISTORTION_D3).setColumn(Integer.MAX_VALUE);
		
		// Read taken images and exterior orientation parameters of images
		EORFileReader eorReader = new EORFileReader(basepath + ".eor", camera);
		eorReader.readAndImport();

		// Read image coordinates and add to camera
		PHCFileReader phcReader = new PHCFileReader(basepath + ".phc", camera, coordinates);
		phcReader.readAndImport();
		
		// Specify the points, which are used to define the frame datum using a random stochastic model
		// !!! This step is not mandatory, just for demonstrations !!!
		ArrayList<ObservationParameter<? extends UnknownParameter<?>>> observedCoordinates = new ArrayList<ObservationParameter<? extends UnknownParameter<?>>>();
		// Random stochastic model
		Random random = new Random();
		double sigma0 = 0.001;
		for (Image image : camera) {
			// Call the images taken from this camera object
			// and extract the image coordinates
			for (ImageCoordinate imageCoordinate : image) {
				// Select corresponding object coordinates
				ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
				if (objectCoordinate.getName().length() > 3)
					objectCoordinate.setDatum(Boolean.FALSE);
				
				if (objectCoordinate.isDatum()) {
					objectCoordinate.setDatum(Boolean.FALSE);
					
					ObservationParameter<UnknownParameter<ObjectCoordinate>> obsX = new ObservationParameter<UnknownParameter<ObjectCoordinate>>(objectCoordinate.getX());
					ObservationParameter<UnknownParameter<ObjectCoordinate>> obsY = new ObservationParameter<UnknownParameter<ObjectCoordinate>>(objectCoordinate.getY());
					ObservationParameter<UnknownParameter<ObjectCoordinate>> obsZ = new ObservationParameter<UnknownParameter<ObjectCoordinate>>(objectCoordinate.getZ());
					
					// If only variances (or standard deviations) are known, simply add the corresponding values directly to the observations to save memory
					// obsX.setVariance(Math.pow(random.nextGaussian(0, sigma0), 2));
					// obsY.setVariance(Math.pow(random.nextGaussian(0, sigma0), 2));
					// obsZ.setVariance(Math.pow(random.nextGaussian(0, sigma0), 2));
					
					observedCoordinates.add(obsX);
					observedCoordinates.add(obsY);
					observedCoordinates.add(obsZ);
				}
			}
		}

		// Create a fully populated dispersion matrix D = U'*U
		int size = observedCoordinates.size();
		DenseMatrix U = new DenseMatrix(size, size);
		for (MatrixEntry element : U) 
			element.set(random.nextGaussian(0, sigma0));

		DenseMatrix V = new DenseMatrix(U, Boolean.TRUE);
		UpperSPDPackMatrix dispersion = new UpperSPDPackMatrix(size);
		U.transAmult(V, dispersion);
		DirectlyObservedParameterGroup observedObjectCoordinate = new DirectlyObservedParameterGroup(dispersion, observedCoordinates);	
		
		
		// Create an adjustment object
		BundleAdjustment adjustment = new BundleAdjustment();
		// Add camera(s) and images to adjustment object
		adjustment.add(camera);
		// Add scale bar(s) to adjustment object
		for (ScaleBar scaleBar : scaleBars)
			adjustment.add(scaleBar);
		// Add observed coordinates to adjustment object
		adjustment.add(observedObjectCoordinate);
		
		
		// Add a listener to get notifications of the adjustment process (not required)
		adjustment.addPropertyChangeListener(new ExampleFlatFiles());
		// Select the estimation type, e.g., L2Norm or Simulation
		adjustment.setEstimationType(EstimationType.L2NORM);
		// Use NONE, if the dispersion matrix is not required
		adjustment.setInvertNormalEquation(MatrixInversion.REDUCED);
		// Call estimate to start the bundle adjustment
		EstimationStateType estimationStateType = adjustment.estimateModel();

		// Check the result
		if (estimationStateType != EstimationStateType.ERROR_FREE_ESTIMATION) {
			System.err.println("Error, bundle adjustment fails...");
		}
		else {
			System.out.println("Bundle adjustment finished successfully...");

			// derive dispersion of parameters
			Matrix D = adjustment.getCofactorMatrix();
			if (D != null)
				D.scale(adjustment.getVarianceFactorAposteriori());
			
			String template = "%10s\t%+16.5f\t%+16.5f\t%+16.5f\t%+12.5f\t%+12.5f\t%+12.5f\t%1s";

			// print coordinates of object points and related uncertainties
			Collection<ObjectCoordinate> objectCoordinates = adjustment.getObjectCoordinates();
			for (ObjectCoordinate objectCoordinate : objectCoordinates) {
				UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
				UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
				UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();

				// indicate observed point (o) and new object points (n)
				char datum = (objectCoordinate.getName().length() > 3) ? 'n' : 'o';

				double x = X.getValue();
				double y = Y.getValue();
				double z = Z.getValue();

				double ux = 0, uy = 0, uz = 0;

				if (D != null && X.getColumn() >= 0 && Y.getColumn() >= 0 && Z.getColumn() >= 0 && X.getColumn() != Integer.MAX_VALUE && Y.getColumn() != Integer.MAX_VALUE && Z.getColumn() != Integer.MAX_VALUE) {
					ux = Math.sqrt(Math.abs(D.get(X.getColumn(), X.getColumn())));
					uy = Math.sqrt(Math.abs(D.get(Y.getColumn(), Y.getColumn())));
					uz = Math.sqrt(Math.abs(D.get(Z.getColumn(), Z.getColumn())));
				}
				System.out.println(String.format(Locale.ENGLISH, template, objectCoordinate.getName(), x, y, z, ux, uy, uz, datum));
			}
			System.out.println();
			
			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			for (UnknownParameter<InteriorOrientation> interiorOrientationParameter : interiorOrientation) {
				System.out.println(String.format(Locale.ENGLISH, "%-25s = %+12.10f %s", interiorOrientationParameter.getParameterType().name(), interiorOrientationParameter.getValue(), interiorOrientationParameter.getColumn() == Integer.MAX_VALUE  ? "fixed" : ""));
			}

			System.out.println();
			
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
