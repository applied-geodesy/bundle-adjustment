# Bundle-Adjustment

## Photogrammetry

Photogrammetry is a technique to obtain spatial coordinates of observed objects by analysing planar coordinates from taken images. 
The collinearity equations yield the functional relation between the planar image coordinates, i.e., the observations, and the spatial 
coordinates, i.e., the parameters to be estimated. 
If the number of observations exceeds the number of required observations to estimate the unknown parameters uniquely, 
usually a least-squares adjustment is performed. 

## Bundle Adjustment
In most photogrammetric applications, a bundle adjustment solves the 
overdetermined system of equations, and yields the spatial coordinates in a consistent frame.
Beside the spatial coordinates, the bundle adjustment also provides auxiliary parameters like specific camera parameters and isometric parameters, 
also known as interior and exterior orientation. In most software packages, the fully populated dispersion matrix of the estimated parameters 
is not derived, but only the standard deviations are provided.
Thus, there is a lack of statistical parameters to evaluate the reliability of the adjustment. Certainly, there are many applications 
that do not need such statistical parameters. However, a small number of applications exists that calls for a rigorous adjustment 
and the fully populated dispersion matrix - especially, if the results of the bundle adjustment are treated as incoming data in 
further analysis steps. For instance, in the framework of deformation analysis, the dispersion of the parameters is required to 
obtain reliable test statistics.

## Metrology and Close-Range Photogrammetry
In the field of metrology and close-range photogrammetry, usually marked points or special targets are used for precise measurements. 
Moreover, the number of points to be estimated is discrete, and so is the dimension of the dispersion matrix. The computational 
effort is justifiable, and the benefit of a small but high accurate set of points preponderate against a dense but less precise 
point cloud. For such specific applications, this tiny bundle adjustment library is developed.

## Example

The following code snippet demonstrates the use case of the library.

```java
// Create an adjustment object using your specific source file reader
BundleAdjustment adjustment = reader.readAndImport();
		
// call estimate to start the bundle adjustment
EstimationStateType estimationStateType = adjustment.estimateModel();
		
// Check the result
if (estimationStateType != EstimationStateType.ERROR_FREE_ESTIMATION) {
	System.err.println("Error, bundle adjustment fails...");
}
else {
	System.out.println("Bundle adjustment finished successfully...");
			
	// derive dispersion of parameters
	Matrix D = adjustment.getCofactorMatrix().scale(adjustment.getVarianceFactorAposteriori());
	String template = "%10s\t%+16.5f\t%+16.5f\t%+16.5f\t%+8.5f\t%+8.5f\t%+8.5f";

	// print coordinates of object points and related uncertainties
	for (ObjectCoordinate objectCoordinate : objectCoordinates) {
		UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
		UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
		UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();

		double x = X.getValue();
		double y = Y.getValue();
		double z = Z.getValue();
		double ux = 0, uy = 0, uz = 0;

		if (X.getColumn() >= 0 && Y.getColumn() >= 0 && Z.getColumn() >= 0) {
			ux = Math.sqrt(Math.abs(D.get(X.getColumn(), X.getColumn())));
			uy = Math.sqrt(Math.abs(D.get(Y.getColumn(), Y.getColumn())));
			uz = Math.sqrt(Math.abs(D.get(Z.getColumn(), Z.getColumn())));
		}
		String formattedStr = String.format(Locale.ENGLISH, template, objectCoordinate.getName(), x, y, z, ux, uy, uz);
		System.out.println(formattedStr);
	}

	// print some statistical parameters
	System.out.println("Degree of freedom:          " + adjustment.getDegreeOfFreedom());
	System.out.println("Variance of unit weight:    " + adjustment.getVarianceFactorApriori() + " : " + adjustment.getVarianceFactorAposteriori());
}
```
