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
import java.util.Collection;
import java.util.Locale;
import java.util.Map;

import org.applied_geodesy.adjustment.EstimationStateType;
import org.applied_geodesy.adjustment.EstimationType;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.ScaleBar;
import org.applied_geodesy.adjustment.bundle.BundleAdjustment.MatrixInversion;
import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadiallySymmetricDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.ZernikeDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.util.io.reader.aicon.EORFileReader;
import org.applied_geodesy.util.io.reader.aicon.IORFileReader;
import org.applied_geodesy.util.io.reader.aicon.OBCFileReader;
import org.applied_geodesy.util.io.reader.aicon.PHCFileReader;
import org.applied_geodesy.util.io.reader.aicon.ScaleFileReader;

import no.uib.cipr.matrix.Matrix;

public class ExampleDistortionModel implements PropertyChangeListener {

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
		IORFileReader iorReader = new IORFileReader(basepath + ".ior", DistortionModel.Type.ZERNIKE_GRADIENT, DistortionModel.Type.ZERNIKE_X, DistortionModel.Type.ZERNIKE_Y);
		Camera camera = iorReader.readAndImport();
		
		// fixing interior camera parameters such as principle distance because of correlations with specific Zernike polynomial e.g. c and Z(4)
		camera.getInteriorOrientation().getPrincipleDistance().setValue(28);
		camera.getInteriorOrientation().getPrincipleDistance().setColumn(Integer.MAX_VALUE);	

		// fixing parameters of standard distortion models to apply Zernike approach
		RadiallySymmetricDistortionModel radiallySymmetricDistortionModel = (RadiallySymmetricDistortionModel)camera.getDistortionModel(DistortionModel.Type.RADIAL_DISTORTION);
		for (UnknownParameter<?> distortionParameter : radiallySymmetricDistortionModel) {
			distortionParameter.setValue(0);
			distortionParameter.setColumn(Integer.MAX_VALUE);
		}
		
		// specify Zernike polynomials
		ZernikeDistortionModel.Gradient zernikeDistortion = (ZernikeDistortionModel.Gradient)camera.getDistortionModel(DistortionModel.Type.ZERNIKE_GRADIENT);
		
		// apply radially symmetric distortions only
		for (int i = 1, order = 0; i < 6; i++) {
			order += i*4;
			zernikeDistortion.add(order);
		}	
		
		// Read taken images and exterior orientation parameters of images
		EORFileReader eorReader = new EORFileReader(basepath + ".eor", camera);
		eorReader.readAndImport();

		// Read image coordinates and add to camera
		PHCFileReader phcReader = new PHCFileReader(basepath + ".phc", camera, coordinates);
		phcReader.readAndImport();
				
		// Create an adjustment object
		BundleAdjustment adjustment = new BundleAdjustment();
		// Add camera(s) and images to adjustment object
		adjustment.add(camera);
		// Add scale bar(s) to adjustment object
		for (ScaleBar scaleBar : scaleBars)
			adjustment.add(scaleBar);
		
		// Add a listener to get notifications of the adjustment process (not required)
		adjustment.addPropertyChangeListener(new ExampleDistortionModel());
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
			for (UnknownParameter<?> unknownParameter : interiorOrientation) {
				System.out.println(String.format(Locale.ENGLISH, "%-27s = %+15.10f %s", unknownParameter.getParameterType().name(), unknownParameter.getValue(), unknownParameter.getColumn() == Integer.MAX_VALUE  ? "fixed" : ""));
			}
			Collection<DistortionModel> distortionModels = camera.getDistortionModels();
			for (DistortionModel model : distortionModels) {
				for (UnknownParameter<?> unknownParameter : model) {
					int order = (unknownParameter instanceof PolynomialCoefficient) ? ((PolynomialCoefficient<?>)unknownParameter).getOrder() : -1;
					System.out.println(String.format(Locale.ENGLISH, "%-27s = %+15.10f %s", unknownParameter.getParameterType().name() + (order < 0 ? "" : "(" + order + ")"), unknownParameter.getValue(), unknownParameter.getColumn() == Integer.MAX_VALUE  ? "fixed" : ""));
				}
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
