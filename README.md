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

## Field of application
The library was developed within the international project [GeoMetre](https://www.ptb.de/empir2018/geometre/home/), a joint research project within the European Metrology Research Programme EMPIR (Grant Number: 18SIB01, Funding: [10.13039/100014132](https://doi.org/10.13039/100014132)). The bundle adjustment was applied in the framework of reference point determination of laser telescopes used for satellite laser ranging (SLR). In order to detect smallest deformations of the main reflector as well as the subreflector of radio telescopes used for very long baseline interferometry (VLBI), this library was used for the rigorous adjustment of photogrammetric measurements.

## References
- Lösler, M., Eschelbach, C., Klügel, T., Riepl, S.: *ILRS Reference Point Determination using Close Range Photogrammetry.* Applied Sciences, 11(6), 2785, 2021. DOI: [10.3390/app11062785](https://doi.org/10.3390/app11062785)
- Lösler, M., Eschelbach, C., Klügel, T.: *Close Range Photogrammetry for High-Precision Reference Point Determination: A Proof of Concept at Satellite Observing System Wettzell.* In: Freymueller, J. T., Sánchez, L. (eds.): Geodesy for a Sustainable Earth, Scientific Assembly of the International Association of Geodesy (IAG), Springer, Berlin, 2022. DOI: [10.1007/1345_2022_141](https://doi.org/10.1007/1345_2022_141)
- Eschelbach, C., Lösler, M.: *A Feasibility Study for Accelerated Reference Point Determination Using Close Range Photogrammetry.* 5th Joint International Symposium on Deformation Monitoring (JISDM), 20-22 June 2022, Polytechnic University of Valencia (UPV), Valencia, Spain, 2022. DOI: [10.4995/JISDM2022.2022.13417](https://doi.org/10.4995/JISDM2022.2022.13417)
- Lösler, M., Eschelbach, C., Greiwe, A., Brechtken, R., Plötz, C., Kronschnabl, G., Neidhardt, A.: *Ray Tracing-Based Delay Model for Compensating Gravitational Deformations of VLBI Radio Telescopes.* Journal of Geodetic Science, 12(1), pp. 165-184, 2022. DOI: [10.515/jogs-2022-0141](https://doi.org/10.515/jogs-2022-0141)

## System requirements
The bundle adjustment is written in the platform-independent programming language Java and, therefore, the software is runnable at each platform and operation system that provides a Java Runtime Environment (JRE). The JRE can be found for several platforms at Oracles [download page](https://java.oracle.com) or at the [OpenJDK](https://openjdk.java.net)-project pages.

## Example
The following [code snippet](https://github.com/applied-geodesy/bundle-adjustment/blob/main/JAICOV/src/org/applied_geodesy/adjustment/bundle/jaicov/JAiCov.java) demonstrates the use case of the library.

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
