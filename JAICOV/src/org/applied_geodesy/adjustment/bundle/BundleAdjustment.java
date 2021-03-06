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

package org.applied_geodesy.adjustment.bundle;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import org.applied_geodesy.adjustment.Constant;
import org.applied_geodesy.adjustment.DefaultValue;
import org.applied_geodesy.adjustment.EstimationStateType;
import org.applied_geodesy.adjustment.EstimationType;
import org.applied_geodesy.adjustment.MathExtension;
import org.applied_geodesy.adjustment.NormalEquationSystem;
import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameter;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.adjustment.defect.DefectType;
import org.applied_geodesy.adjustment.defect.RankDefect;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.MatrixNotSPDException;
import no.uib.cipr.matrix.MatrixSingularException;
import no.uib.cipr.matrix.UpperSymmBandMatrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;
import no.uib.cipr.matrix.Vector;

public class BundleAdjustment {
	public enum MatrixInversion {
		NONE,
		FULL,
		REDUCED
	};
	
	private class DispersionMatrixExportProperties {
		private String exportPathAndFileName;
		private boolean includeDatumConditions = Boolean.FALSE;
		
		DispersionMatrixExportProperties(String exportPathAndFileName, boolean includeDatumConditions) {
			this.exportPathAndFileName  = exportPathAndFileName;
			this.includeDatumConditions = includeDatumConditions;
		}
	}
	
	private final PropertyChangeSupport change = new PropertyChangeSupport(this);
	private EstimationStateType currentEstimationStatus = EstimationStateType.BUSY;
	private EstimationType estimationType = EstimationType.L2NORM;
	private DispersionMatrixExportProperties dispersionMatrixExportProperties = null;
	
	private static double SQRT_EPS = Math.sqrt(Constant.EPS);
	private int maximalNumberOfIterations = DefaultValue.getMaximalNumberOfIterations(),
			iterationStep                 = 0,
			numberOfUnknownParameters     = 0,
			numberOfObservations          = 0,
			numberOfInteriorOrientationParameters = 0;

	private boolean interrupt                    = false,
			applyAposterioriVarianceOfUnitWeight = true,
			useCentroidedCoordinates             = true,
			calculateStochasticParameters        = false;
	
	private MatrixInversion invertNormalEquationMatrix = MatrixInversion.FULL;

	private double maxDx    = Double.MIN_VALUE,
			omega           = 0.0,
			currentMaxAbsDx = this.maxDx,
			sigma2apriori   = 1.0;
	
	private ObjectCoordinate centroid = new ObjectCoordinate(this.getClass().getSimpleName(), 0, 0, 0);

	private UpperSymmPackMatrix Qxx = null;

	private RankDefect rankDefect = new RankDefect();

	private List<Camera> cameras = new ArrayList<Camera>();
	private List<ObservationParameter<?>> observations = new ArrayList<ObservationParameter<?>>();
	private List<UnknownParameter<?>> unknownParameters = new ArrayList<UnknownParameter<?>>();
	private Map<UnknownParameter<?>, LinkedHashSet<ObservationParameter<?>>> observationsOfUnknownParameters = new LinkedHashMap<UnknownParameter<?>, LinkedHashSet<ObservationParameter<?>>>();
	private Set<ObjectCoordinate> objectCoordinates = new LinkedHashSet<ObjectCoordinate>();
	private Set<ScaleBar> scaleBars = new LinkedHashSet<ScaleBar>();

	public BundleAdjustment() {}
	
	private void centroidCoordinates(boolean invert) throws UnsupportedOperationException {
		if (!invert) {
			double x0 = 0, y0 = 0, z0 = 0;
			int cntX = 0, cntY = 0, cntZ = 0;
			for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
				ParameterType paramType = unknownParameter.getParameterType();
				switch (paramType) {
				case CAMERA_COORDINATE_X:
				case OBJECT_COORDINATE_X:
					x0 += unknownParameter.getValue();
					cntX++;
					break;
				case CAMERA_COORDINATE_Y:
				case OBJECT_COORDINATE_Y:
					y0 += unknownParameter.getValue();
					cntY++;
					break;
				case CAMERA_COORDINATE_Z:
				case OBJECT_COORDINATE_Z:
					z0 += unknownParameter.getValue();
					cntZ++;
					break;
				default:
					break;
				}
			}

			if (cntX == cntY && cntX == cntZ && cntY == cntZ && cntX > 0) {
				x0 /= cntX;
				y0 /= cntY;
				z0 /= cntZ;
				this.centroid.getX().setValue(x0);
				this.centroid.getY().setValue(y0);
				this.centroid.getZ().setValue(z0);
			}
			else 
				throw new UnsupportedOperationException(this.getClass().getSimpleName() +" Error, the numbers of coordinate components are un-equal or zero [" + cntX + ", " + cntY + ", " + cntZ + "]");
		}
		
		double sign = invert ? 1.0 : -1.0;
		
		double x0 = sign * this.centroid.getX().getValue();
		double y0 = sign * this.centroid.getY().getValue();
		double z0 = sign * this.centroid.getZ().getValue();
		
		for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
			ParameterType paramType = unknownParameter.getParameterType();
			switch (paramType) {
			case CAMERA_COORDINATE_X:
			case OBJECT_COORDINATE_X:
				unknownParameter.setValue(unknownParameter.getValue() + x0);
				break;
			case CAMERA_COORDINATE_Y:
			case OBJECT_COORDINATE_Y:
				unknownParameter.setValue(unknownParameter.getValue() + y0);
				break;
			case CAMERA_COORDINATE_Z:
			case OBJECT_COORDINATE_Z:
				unknownParameter.setValue(unknownParameter.getValue() + z0);
				break;
			default:
				break;
			}
		}
	}

	public EstimationStateType estimateModel() {
		this.currentEstimationStatus = EstimationStateType.BUSY;
		this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);

		this.maxDx = Double.MIN_VALUE;
		this.currentMaxAbsDx = this.maxDx;

		int runs = this.maximalNumberOfIterations-1;
		boolean isEstimated = false, estimateCompleteModel = false, isConverge = true;

		if (this.maximalNumberOfIterations == 0)
			estimateCompleteModel = isEstimated = true;

		this.sigma2apriori = this.sigma2apriori > 0 ? this.sigma2apriori : 1.0;

		this.prepareUnknwonParameters();
		if (this.useCentroidedCoordinates)
			this.centroidCoordinates(false);
		
		try {
			// Reset aller Iterationseinstellungen
			this.maxDx = Double.MIN_VALUE;
			this.currentMaxAbsDx = this.maxDx;
			runs = this.maximalNumberOfIterations - 1;
			isEstimated = false;
			estimateCompleteModel = false;
			isConverge = true;

			do {
				this.maxDx = Double.MIN_VALUE;
				this.iterationStep = this.maximalNumberOfIterations-runs;
				this.currentEstimationStatus = EstimationStateType.ITERATE;
				this.change.firePropertyChange(this.currentEstimationStatus.name(), this.maximalNumberOfIterations, this.iterationStep);

				// erzeuge Normalgleichung
				NormalEquationSystem neq = this.createNormalEquation();

				// Nutze Vorkonditionierung
				this.applyPrecondition(neq);

				if (this.interrupt || neq == null) {
					this.currentEstimationStatus = EstimationStateType.INTERRUPT;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					this.interrupt = false;
					return this.currentEstimationStatus;
				}
				DenseVector n = neq.getVector();
				UpperSymmPackMatrix N = neq.getMatrix();
				Vector dx = n;

				estimateCompleteModel = isEstimated;				
				try {
					if (estimateCompleteModel) {
						this.calculateStochasticParameters = estimateCompleteModel;

						if (this.invertNormalEquationMatrix != MatrixInversion.NONE) {
							this.currentEstimationStatus = EstimationStateType.INVERT_NORMAL_EQUATION_MATRIX;
							this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
						}
						
						// In-Situ Invertierung der NGL: N <-- Qxx, n <-- dx 
						if (this.invertNormalEquationMatrix == MatrixInversion.REDUCED) {
							int numRows = this.numberOfInteriorOrientationParameters + this.objectCoordinates.size() * 3 + this.rankDefect.getDefect();
							// reduziere Gleichungssystem
							this.reduceNormalEquationMatrix(neq);
							// In-Situ Invertierung der NGL: N <-- Qxx, n <-- dx 
							MathExtension.solve(N, n, numRows, Boolean.TRUE);
						}
						else {
							MathExtension.solve(N, n, this.invertNormalEquationMatrix == MatrixInversion.FULL);
						}

						this.applyPrecondition(neq.getPreconditioner(), N);
						this.Qxx = N;


						if (this.calculateStochasticParameters) {
							this.currentEstimationStatus = EstimationStateType.ESTIAMTE_STOCHASTIC_PARAMETERS;
							this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
						}
					}
					else
						// Loese Nx=n und ueberschreibe n durch die Loesung x
						MathExtension.solve(N, n, false);
					this.applyPrecondition(neq.getPreconditioner(), dx);

					n = null;
					N = null;
				}
				catch (MatrixSingularException | MatrixNotSPDException | IllegalArgumentException | ArrayIndexOutOfBoundsException e) {
					e.printStackTrace();
					this.currentEstimationStatus = EstimationStateType.SINGULAR_MATRIX;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					return this.currentEstimationStatus;
				}
				catch (Exception e) {
					e.printStackTrace();
					this.currentEstimationStatus = EstimationStateType.INTERRUPT;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					return this.currentEstimationStatus;
				}

				this.updateModel(dx, estimateCompleteModel);
				dx = null;

				if (this.interrupt) {
					this.currentEstimationStatus = EstimationStateType.INTERRUPT;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					this.interrupt = false;
					return this.currentEstimationStatus;
				}

				if (Double.isInfinite(this.maxDx) || Double.isNaN(this.maxDx)) {
					this.currentEstimationStatus = EstimationStateType.SINGULAR_MATRIX;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					return this.currentEstimationStatus;
				}
				else if (this.maxDx <= SQRT_EPS && runs > 0) {
					isEstimated = true;
					this.currentEstimationStatus = EstimationStateType.CONVERGENCE;
					//if (this.estimationType != EstimationType.UNSCENTED_TRANSFORMATION)
					this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxDx);
				}
				else if (runs-- <= 1) {
					if (estimateCompleteModel) {
						this.currentEstimationStatus = EstimationStateType.NO_CONVERGENCE;
						this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxDx);
						isConverge = false;
					}
					isEstimated = true;
				}
				else {
					this.currentEstimationStatus = EstimationStateType.CONVERGENCE;
					//if (this.estimationType != EstimationType.UNSCENTED_TRANSFORMATION)
					this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxDx);
				}

				// Sollten nur stochastische Punkte enthalten sein, ist maxDx MIN_VALUE
				if (this.maxDx > Double.MIN_VALUE)
					this.currentMaxAbsDx = this.maxDx;
				else
					this.maxDx = this.currentMaxAbsDx;
			}
			while (!estimateCompleteModel);
			
			if (this.useCentroidedCoordinates)
				this.centroidCoordinates(true);
			
			// export dispersion matrix to file
			if (!this.exportCovarianceMatrix()) 
				System.err.println("Fehler, Varianz-Kovarianz-Matrix konnte nicht exportiert werden.");

		}
		catch (OutOfMemoryError e) {
			e.printStackTrace();
			this.currentEstimationStatus = EstimationStateType.OUT_OF_MEMORY;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
			return this.currentEstimationStatus;
		}

		if (!isConverge) {
			this.currentEstimationStatus = EstimationStateType.NO_CONVERGENCE;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxDx);
		}
		else if (this.currentEstimationStatus.getId() == EstimationStateType.BUSY.getId() || this.calculateStochasticParameters) {
			this.currentEstimationStatus = EstimationStateType.ERROR_FREE_ESTIMATION;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxDx);
		}

		return this.currentEstimationStatus;
	}
	
	private void updateModel(Vector dx, boolean updateCompleteModel) {
		this.omega = 0;
		for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
			int column = unknownParameter.getColumn();
			this.maxDx = Math.max(this.maxDx, Math.abs(dx.get(column))); 
			unknownParameter.setValue(unknownParameter.getValue() + dx.get(column));
		}
		if (updateCompleteModel) {
			for (ObservationParameter<?> observation : this.observations) {
				double v = this.estimationType == EstimationType.SIMULATION ? 0.0 :  PartialDerivativeFactory.getMisclosure(observation);
				this.omega += v * v * this.sigma2apriori / observation.getVariance();
			}
		}
	}

//	private void addDatumConditionColumns(UpperSymmPackMatrix N) {
//		// center of mass
//		double x0 = 0, y0 = 0, z0 = 0;
//		int count = 0;
//		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
//			if (objectCoordinate.isDatum()) {
//				x0 += objectCoordinate.getX().getValue();
//				y0 += objectCoordinate.getY().getValue();
//				z0 += objectCoordinate.getZ().getValue();
//				count++;
//			}
//		}
//
//		if (count < 3)
//			throw new IllegalArgumentException(this.getClass() + " Error, not enought object points to realise the frame datum!");
//
//		x0 = x0 / (double)count;
//		y0 = y0 / (double)count;
//		z0 = z0 / (double)count;
//
//		int defectSize = this.rankDefect.getDefect();
//		int row = N.numRows() - defectSize;
//
//		// Position in condition matrix
//		int defectRow = 0;
//		int tx   = this.rankDefect.estimateTranslationX() ? defectRow++ : -1;
//		int ty   = this.rankDefect.estimateTranslationY() ? defectRow++ : -1;
//		int tz   = this.rankDefect.estimateTranslationZ() ? defectRow++ : -1;
//		int rx   = this.rankDefect.estimateRotationX()    ? defectRow++ : -1;
//		int ry   = this.rankDefect.estimateRotationY()    ? defectRow++ : -1;
//		int rz   = this.rankDefect.estimateRotationZ()    ? defectRow++ : -1;
//		int mxyz = this.rankDefect.estimateScale()        ? defectRow++ : -1;
//
//		// Sum of squared row
//		double normColumn[] = new double[defectSize];
//
//		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
//			if (!objectCoordinate.isDatum())
//				continue;
//
//			int columnX = objectCoordinate.getX().getColumn();
//			int columnY = objectCoordinate.getY().getColumn();
//			int columnZ = objectCoordinate.getZ().getColumn();
//
//			double x = objectCoordinate.getX().getValue() - x0;
//			double y = objectCoordinate.getY().getValue() - y0;
//			double z = objectCoordinate.getZ().getValue() - z0;
//
//			if (tx >= 0) {
//				N.set(columnX, row + tx, 1.0);
//
//				normColumn[tx] += 1.0; 
//			}
//
//			if (ty >= 0) {
//				N.set(columnY, row + ty, 1.0);
//
//				normColumn[ty] += 1.0; 
//			}
//
//			if (tz >= 0) {
//				N.set(columnZ, row + tz, 1.0);
//
//				normColumn[tz] += 1.0; 
//			}
//
//			if (rx >= 0) {
//				N.set(columnY, row + rx,  z );
//				N.set(columnZ, row + rx, -y );
//
//				normColumn[rx] += z*z + y*y;
//			}
//
//			if (ry >= 0) {
//				N.set(columnX, row + ry, -z );
//				N.set(columnZ, row + ry,  x );
//
//				normColumn[ry] += z*z + x*x;
//			}
//
//			if (rz >= 0) {
//				N.set(columnX, row + rz,  y );
//				N.set(columnY, row + rz, -x );	
//
//				normColumn[rz] += x*x + y*y;
//			}
//
//			if (mxyz >= 0) {
//				N.set(columnX, row+mxyz, x );
//				N.set(columnY, row+mxyz, y );
//				N.set(columnZ, row+mxyz, z );
//
//				normColumn[mxyz] += x*x + y*y + z*z;
//			}
//		}
//
//		// Normieren der Spalten
//		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
//			if (!objectCoordinate.isDatum())
//				continue;
//
//			int columnX = objectCoordinate.getX().getColumn();
//			int columnY = objectCoordinate.getY().getColumn();
//			int columnZ = objectCoordinate.getZ().getColumn();
//
//			if (tx >= 0)
//				N.set(columnX, row + tx, N.get(columnX, row + tx) / Math.sqrt(normColumn[tx]));
//
//			if (ty >= 0)
//				N.set(columnY, row + ty, N.get(columnY, row + ty) / Math.sqrt(normColumn[ty]));
//
//			if (tz >= 0)
//				N.set(columnZ, row + tz, N.get(columnZ, row + tz) / Math.sqrt(normColumn[tz]));
//
//			if (rx >= 0) {
//				N.set(columnY, row + rx, N.get(columnY, row + rx) / Math.sqrt(normColumn[rx]) );
//				N.set(columnZ, row + rx, N.get(columnZ, row + rx) / Math.sqrt(normColumn[rx]) );
//			}
//
//			if (ry >= 0) {
//				N.set(columnX, row + ry, N.get(columnX, row + ry) / Math.sqrt(normColumn[ry]) );
//				N.set(columnZ, row + ry, N.get(columnZ, row + ry) / Math.sqrt(normColumn[ry]) );
//			}
//
//			if (rz >= 0) {
//				N.set(columnX, row + rz, N.get(columnX, row + rz) / Math.sqrt(normColumn[rz]) );
//				N.set(columnY, row + rz, N.get(columnY, row + rz) / Math.sqrt(normColumn[rz]) );	
//			}
//
//			if (mxyz >= 0) {
//				N.set(columnX, row + mxyz, N.get(columnX, row + mxyz) / Math.sqrt(normColumn[mxyz]) );
//				N.set(columnY, row + mxyz, N.get(columnY, row + mxyz) / Math.sqrt(normColumn[mxyz]) );
//				N.set(columnZ, row + mxyz, N.get(columnZ, row + mxyz) / Math.sqrt(normColumn[mxyz]) );
//			}
//		}
//	}

	private void addDatumConditionRows(UpperSymmPackMatrix N) {
		// center of mass
		double x0 = 0, y0 = 0, z0 = 0;
		int count = 0;
		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
			if (objectCoordinate.isDatum()) {
				x0 += objectCoordinate.getX().getValue();
				y0 += objectCoordinate.getY().getValue();
				z0 += objectCoordinate.getZ().getValue();
				count++;
			}
		}

		if (count < 3)
			throw new IllegalArgumentException(this.getClass() + " Error, not enought object points to realise the frame datum!");

		x0 = x0 / (double)count;
		y0 = y0 / (double)count;
		z0 = z0 / (double)count;

		int defectSize = this.rankDefect.getDefect();
		int row = 0;

		// Position in condition matrix
		int defectRow = 0;
		int tx   = this.rankDefect.estimateTranslationX() ? defectRow++ : -1;
		int ty   = this.rankDefect.estimateTranslationY() ? defectRow++ : -1;
		int tz   = this.rankDefect.estimateTranslationZ() ? defectRow++ : -1;
		int rx   = this.rankDefect.estimateRotationX()    ? defectRow++ : -1;
		int ry   = this.rankDefect.estimateRotationY()    ? defectRow++ : -1;
		int rz   = this.rankDefect.estimateRotationZ()    ? defectRow++ : -1;
		int mxyz = this.rankDefect.estimateScale()        ? defectRow++ : -1;

		// Sum of squared row
		double normColumn[] = new double[defectSize];

		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
			if (!objectCoordinate.isDatum())
				continue;

			int columnX = objectCoordinate.getX().getColumn();
			int columnY = objectCoordinate.getY().getColumn();
			int columnZ = objectCoordinate.getZ().getColumn();

			double x = objectCoordinate.getX().getValue() - x0;
			double y = objectCoordinate.getY().getValue() - y0;
			double z = objectCoordinate.getZ().getValue() - z0;

			if (tx >= 0) {
				N.set(row + tx, columnX, 1.0);

				normColumn[tx] += 1.0; 
			}

			if (ty >= 0) {
				N.set(row + ty, columnY, 1.0);

				normColumn[ty] += 1.0; 
			}

			if (tz >= 0) {
				N.set(row + tz, columnZ, 1.0);

				normColumn[tz] += 1.0; 
			}

			if (rx >= 0) {
				N.set(row + rx, columnY,  z );
				N.set(row + rx, columnZ, -y );

				normColumn[rx] += z*z + y*y;
			}

			if (ry >= 0) {
				N.set(row + ry, columnX, -z );
				N.set(row + ry, columnZ,  x );

				normColumn[ry] += z*z + x*x;
			}

			if (rz >= 0) {
				N.set(row + rz, columnX,  y );
				N.set(row + rz, columnY, -x );	

				normColumn[rz] += x*x + y*y;
			}

			if (mxyz >= 0) {
				N.set(row+mxyz, columnX, x );
				N.set(row+mxyz, columnY, y );
				N.set(row+mxyz, columnZ, z );

				normColumn[mxyz] += x*x + y*y + z*z;
			}
		}

		// Normieren der Spalten
		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
			if (!objectCoordinate.isDatum())
				continue;

			int columnX = objectCoordinate.getX().getColumn();
			int columnY = objectCoordinate.getY().getColumn();
			int columnZ = objectCoordinate.getZ().getColumn();

			if (tx >= 0)
				N.set(row + tx, columnX, N.get(row + tx, columnX) / Math.sqrt(normColumn[tx]));

			if (ty >= 0)
				N.set(row + ty, columnY, N.get(row + ty, columnY) / Math.sqrt(normColumn[ty]));

			if (tz >= 0)
				N.set(row + tz, columnZ, N.get(row + tz, columnZ) / Math.sqrt(normColumn[tz]));

			if (rx >= 0) {
				N.set(row + rx, columnY, N.get(row + rx, columnY) / Math.sqrt(normColumn[rx]) );
				N.set(row + rx, columnZ, N.get(row + rx, columnZ) / Math.sqrt(normColumn[rx]) );
			}

			if (ry >= 0) {
				N.set(row + ry, columnX, N.get(row + ry, columnX) / Math.sqrt(normColumn[ry]) );
				N.set(row + ry, columnZ, N.get(row + ry, columnZ) / Math.sqrt(normColumn[ry]) );
			}

			if (rz >= 0) {
				N.set(row + rz, columnX, N.get(row + rz, columnX) / Math.sqrt(normColumn[rz]) );
				N.set(row + rz, columnY, N.get(row + rz, columnY) / Math.sqrt(normColumn[rz]) );	
			}

			if (mxyz >= 0) {
				N.set(row + mxyz, columnX, N.get(row + mxyz, columnX) / Math.sqrt(normColumn[mxyz]) );
				N.set(row + mxyz, columnY, N.get(row + mxyz, columnY) / Math.sqrt(normColumn[mxyz]) );
				N.set(row + mxyz, columnZ, N.get(row + mxyz, columnZ) / Math.sqrt(normColumn[mxyz]) );
			}
		}
	}

	private void addObservation(UnknownParameter<?> unknownParameter, ObservationParameter<?> observation) {
		if (!this.observationsOfUnknownParameters.containsKey(unknownParameter))
			this.observationsOfUnknownParameters.put(unknownParameter, new LinkedHashSet<ObservationParameter<?>>());
		this.observationsOfUnknownParameters.get(unknownParameter).add(observation);

		if (!this.observations.contains(observation)) {
			this.observations.add(observation);
			this.sigma2apriori = Math.min(this.sigma2apriori, observation.getVariance());
		}
	}

	private void addUnknownParameter(UnknownParameter<?> unknownParameter) {
		if (unknownParameter.getColumn() == -1 && !this.unknownParameters.contains(unknownParameter)) {
			unknownParameter.setColumn( this.numberOfUnknownParameters++ );
			this.unknownParameters.add( unknownParameter );
		}
	}
	
	/**
	 * @deprecated - Debugging only
	 * @param coordinates
	 */
	public void add(List<ObjectCoordinate> coordinates) {
		for (ObjectCoordinate coordinate : coordinates) {
			this.addUnknownParameter(coordinate.getX());
			this.addUnknownParameter(coordinate.getY());
			this.addUnknownParameter(coordinate.getZ());
		}
	}
	
	public void add(Camera camera) {
		this.cameras.add(camera);
	}

	public void add(ScaleBar scaleBar) {
		this.scaleBars.add(scaleBar);
	}
	
	private void prepareUnknwonParameters() {
		for (Camera camera : this.cameras) {
			// add parameters of interior orientation
			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			
			for (Image image : camera) {
				// add parameters of exterior orientation
				ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();

				// add pixel coordinatens as observations
				for (ImageCoordinate imageCoordinate : image) {
					imageCoordinate.getX().setRow( this.numberOfObservations++ );
					imageCoordinate.getY().setRow( this.numberOfObservations++ );

					ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
					this.objectCoordinates.add(objectCoordinate);

					this.addUnknownParameter(objectCoordinate.getX());
					this.addUnknownParameter(objectCoordinate.getY());
					this.addUnknownParameter(objectCoordinate.getZ());

					this.addObservation(objectCoordinate.getX(), imageCoordinate.getX());
					this.addObservation(objectCoordinate.getX(), imageCoordinate.getY());

					this.addObservation(objectCoordinate.getY(), imageCoordinate.getX());
					this.addObservation(objectCoordinate.getY(), imageCoordinate.getY());

					this.addObservation(objectCoordinate.getZ(), imageCoordinate.getX());
					this.addObservation(objectCoordinate.getZ(), imageCoordinate.getY());

					for (UnknownParameter<InteriorOrientation> unknownParameter : interiorOrientation) {
						this.addObservation(unknownParameter, imageCoordinate.getX());
						this.addObservation(unknownParameter, imageCoordinate.getY());
					}

					for (UnknownParameter<ExteriorOrientation> unknownParameter : exteriorOrientation) {
						this.addObservation(unknownParameter, imageCoordinate.getX());
						this.addObservation(unknownParameter, imageCoordinate.getY());
					}

				}
			}
		}
		
		for (Camera camera : this.cameras) {
			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			for (UnknownParameter<InteriorOrientation> unknownParameter : interiorOrientation) {
				if (unknownParameter.getColumn() == -1) {
					this.addUnknownParameter(unknownParameter);
					if (unknownParameter.getColumn() != -1)	
						this.numberOfInteriorOrientationParameters++;
				}
			}
		}
		
		for (Camera camera : this.cameras) {
			for (Image image : camera) {
				// add parameters of exterior orientation
				ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
				for (UnknownParameter<ExteriorOrientation> unknownParameter : exteriorOrientation)
					this.addUnknownParameter(unknownParameter);
			}
		}
		
		for (ScaleBar scaleBar : this.scaleBars) {
			scaleBar.getLength().setRow( this.numberOfObservations++ );

			ObjectCoordinate objectCoordinateA = scaleBar.getObjectCoordinateA();
			ObjectCoordinate objectCoordinateB = scaleBar.getObjectCoordinateB();

			this.addUnknownParameter(objectCoordinateA.getX());
			this.addUnknownParameter(objectCoordinateA.getY());
			this.addUnknownParameter(objectCoordinateA.getZ());

			this.addUnknownParameter(objectCoordinateB.getX());
			this.addUnknownParameter(objectCoordinateB.getY());
			this.addUnknownParameter(objectCoordinateB.getZ());

			this.addObservation(objectCoordinateA.getX(), scaleBar.getLength());
			this.addObservation(objectCoordinateA.getY(), scaleBar.getLength());
			this.addObservation(objectCoordinateA.getZ(), scaleBar.getLength());

			this.addObservation(objectCoordinateB.getX(), scaleBar.getLength());
			this.addObservation(objectCoordinateB.getY(), scaleBar.getLength());
			this.addObservation(objectCoordinateB.getZ(), scaleBar.getLength());
		}
		
		this.detectRankDefect();

		// re-number 
		int numberOfDatumConditions = this.rankDefect.getDefect();
		if (numberOfDatumConditions > 0) {
			for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
				unknownParameter.setColumn(unknownParameter.getColumn() + numberOfDatumConditions);
			}
		}
	}

	/**
	 * @return - Normalgleichungssystem
	 */
	public NormalEquationSystem createNormalEquation() {

		UpperSymmPackMatrix N = new UpperSymmPackMatrix( this.numberOfUnknownParameters + this.rankDefect.getDefect());
		UpperSymmBandMatrix V = new UpperSymmBandMatrix( N.numRows(), 0 );
		DenseVector n = new DenseVector( N.numRows() );
		
		for (ObservationParameter<?> observation : this.observations) {
			PartialDerivativeFactory.getPartialDerivative(this.sigma2apriori, N, n, observation);
		}

		this.addDatumConditionRows(N);
		
//		// Vorkonditionierer == Wurzel der Hauptdiagonale
		for (int column = 0; column < N.numColumns(); column++) {
			double value = N.get(column, column);
			V.set(column, column, value > Constant.EPS ? 1.0 / Math.sqrt(value) : 1.0);
		}
		
		
		if (this.estimationType == EstimationType.SIMULATION)
			n.zero();
		return new NormalEquationSystem(N, n, V);
	}

	private void detectRankDefect() {
		this.rankDefect.reset();

		this.rankDefect.setTranslationX(DefectType.FREE);
		this.rankDefect.setTranslationY(DefectType.FREE);
		this.rankDefect.setTranslationZ(DefectType.FREE);

		this.rankDefect.setRotationX(DefectType.FREE);
		this.rankDefect.setRotationY(DefectType.FREE);
		this.rankDefect.setRotationZ(DefectType.FREE);

		this.rankDefect.setScale(this.scaleBars.isEmpty() ? DefectType.FREE : DefectType.FIXED);
	}
	
	private void applyPrecondition(NormalEquationSystem neq) {
		if (neq.getPreconditioner() != null) {
			this.applyPrecondition(neq.getPreconditioner(), neq.getMatrix());
			this.applyPrecondition(neq.getPreconditioner(), neq.getVector());
		}
	}

	private void applyPrecondition(UpperSymmBandMatrix V, UpperSymmPackMatrix M) {
		for (int row = 0; row < M.numRows(); row++) {
			for (int column = row; column < M.numColumns(); column++) {
				M.set(row, column, V.get(column, column) * M.get(row, column) * V.get(row, row));
			}
		}
		
//		for (int i = 0; i < this.unknownParameters.size(); i++) {
//			UnknownParameter<?> unknownParameteri = this.unknownParameters.get(i);
//			int row = unknownParameteri.getColumn();
//			for (int j = i; j < this.unknownParameters.size(); j++) {
//				UnknownParameter<?> unknownParameterj = this.unknownParameters.get(j);
//				int column = unknownParameterj.getColumn();
//				M.set(row, column, V.get(column, column) * M.get(row, column) * V.get(row, row));
//			}
//		}
	}
	
	private void applyPrecondition(UpperSymmBandMatrix V, Vector m) {
		for (int row = 0; row < m.size(); row++) {
			m.set(row, V.get(row, row) * m.get(row));
		}
		
//		for (int i = 0; i < this.unknownParameters.size(); i++) {
//			UnknownParameter<?> unknownParameteri = this.unknownParameters.get(i);
//			int row = unknownParameteri.getColumn();
//			m.set(row, V.get(row, row) * m.get(row));
//		}
	}

	public Set<ObjectCoordinate> getObjectCoordinates() {
		return this.objectCoordinates;
	}
	
	public int getDegreeOfFreedom() {
		return this.numberOfObservations - this.numberOfUnknownParameters + this.rankDefect.getDefect();
	}
	
	public double getVarianceFactorAposteriori() {
		int degreeOfFreedom = this.getDegreeOfFreedom();
		return degreeOfFreedom > 0 && this.omega > 0 && this.estimationType != EstimationType.SIMULATION && this.applyAposterioriVarianceOfUnitWeight ? Math.abs(this.omega/(double)degreeOfFreedom) : 1.0;
	}
	
	public double getVarianceFactorApriori() {
		return this.sigma2apriori;
	}
	
	public void setCovarianceExportPathAndBaseName(String path, boolean exportDatumConditions) {
		this.dispersionMatrixExportProperties = path == null ? null : new DispersionMatrixExportProperties(path, exportDatumConditions);
	}
	
	public void setEstimationType(EstimationType estimationType) throws UnsupportedOperationException {
		if (estimationType == EstimationType.L2NORM || estimationType == EstimationType.SIMULATION)
			this.estimationType = estimationType;
		else
			throw new UnsupportedOperationException(this.getClass().getSimpleName() + " Error, this estimation type is not supported! " + estimationType);
	}
	
	/**
	 * Legt fest, ob die Normalgleichung invertiert werden soll
	 * um die Kofaktormatrix zu bestimmen. Sind nur die Parameter
	 * von Interesse, dann kann hierdurch z.T. erhebliche Rechenzeit 
	 * gesparrt werden.
	 * 
	 * @param invert
	 */
	public void setInvertNormalEquation(MatrixInversion invert) {
		this.invertNormalEquationMatrix = invert;
	}
	
	private boolean exportCovarianceMatrix() {
		if (this.dispersionMatrixExportProperties == null)
			return true;

		String exportPathAndFileName   = this.dispersionMatrixExportProperties.exportPathAndFileName;
		boolean includeDatumConditions = this.dispersionMatrixExportProperties.includeDatumConditions;

		File coVarMatrixFile = new File(exportPathAndFileName + ".cxx");
		File coVarInfoFile   = new File(exportPathAndFileName + ".info");

		return this.exportCovarianceMatrixInfoToFile(coVarInfoFile, includeDatumConditions) && this.exportCovarianceMatrixToFile(coVarMatrixFile, includeDatumConditions);
	}
	
	public UpperSymmPackMatrix getCofactorMatrix() {
		return this.Qxx;
	}

	/**
	 * Schreibt punktbezogene Informationen zur CoVar raus
	 * @param f
	 * @return isWritten
	 */
	private boolean exportCovarianceMatrixInfoToFile(File f, boolean includeDatumConditions) {
		// noch keine Loesung vorhanden
		if (f == null) //  || this.Qxx == null
			return false;

		this.currentEstimationStatus = EstimationStateType.EXPORT_COVARIANCE_INFORMATION;
		this.change.firePropertyChange(this.currentEstimationStatus.name(), null, f.toString());

		int indexDatum = includeDatumConditions ? 0 : this.rankDefect.getDefect();
		boolean isComplete = false;
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter( f )));
			//Pkt,Type(XYZ),Coord,Row in NGL
			String format = "%25s\t%5s\t%35.15f\t%10d%n";
			for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
				UnknownParameter<ObjectCoordinate> X = objectCoordinate.getX();
				UnknownParameter<ObjectCoordinate> Y = objectCoordinate.getY();
				UnknownParameter<ObjectCoordinate> Z = objectCoordinate.getZ();
				String name = objectCoordinate.getName();
				
				pw.printf(Locale.ENGLISH, format, name, 'X', X.getValue(), X.getColumn() - indexDatum);
				pw.printf(Locale.ENGLISH, format, name, 'Y', Y.getValue(), Y.getColumn() - indexDatum);
				pw.printf(Locale.ENGLISH, format, name, 'Z', Z.getValue(), Z.getColumn() - indexDatum);

			}
			isComplete = true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		finally {
			if (pw != null) {
				pw.close();
			}
		}

		return isComplete;
	}

	/**
	 * Schreibt die CoVar raus
	 * @param f
	 * @return isWritten
	 */
	private boolean exportCovarianceMatrixToFile(File f, boolean includeDatumConditions) {
		// noch keine Loesung vorhanden
		if (f == null || this.Qxx == null || this.Qxx.numRows() < this.numberOfUnknownParameters)
			return false;

		this.currentEstimationStatus = EstimationStateType.EXPORT_COVARIANCE_MATRIX;
		this.change.firePropertyChange(this.currentEstimationStatus.name(), null, f.toString());

		boolean isComplete = false;
		PrintWriter pw = null;
		double sigma2apost = this.getVarianceFactorAposteriori();

		int numberOfDatumConditions = this.rankDefect.getDefect();
		int size = this.invertNormalEquationMatrix == MatrixInversion.REDUCED ? this.numberOfInteriorOrientationParameters + this.objectCoordinates.size() * 3 : this.Qxx.numRows() - numberOfDatumConditions;
		size = includeDatumConditions ? size + numberOfDatumConditions : size;
		
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter( f )));
			for (int i = includeDatumConditions ? 0 : numberOfDatumConditions; i < size; i++) {
				for (int j = includeDatumConditions ? 0 : numberOfDatumConditions; j < size; j++) {
					pw.printf(Locale.ENGLISH, "%+35.15f  ", sigma2apost*this.Qxx.get(i, j));
				}
				pw.println();	
			}

			isComplete = true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		finally {
			if (pw != null) {
				pw.close();
			}
		}
		return isComplete;
	}
	
	public void useCentroidedCoordinates(boolean useCentroidedCoordinates) {
		this.useCentroidedCoordinates = useCentroidedCoordinates;
	}
	
	public void applyAposterioriVarianceOfUnitWeight(boolean applyAposterioriVarianceOfUnitWeight) {
		this.applyAposterioriVarianceOfUnitWeight = applyAposterioriVarianceOfUnitWeight;
	}
	
	private void reduceNormalEquationMatrix(NormalEquationSystem neq) throws IllegalArgumentException, MatrixSingularException {
		for (Camera camera : this.cameras)
			this.reduceNormalEquationMatrix(neq, camera);
	}
	
	private void reduceNormalEquationMatrix(NormalEquationSystem neq, Camera camera) throws IllegalArgumentException, MatrixSingularException {
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		
		List<UnknownParameter<InteriorOrientation>> unknownInteriorOrientation = new ArrayList<UnknownParameter<InteriorOrientation>>(10);
		for (UnknownParameter<InteriorOrientation> unknownParameter : interiorOrientation) {
			if (unknownParameter.getColumn() < 0 || unknownParameter.getColumn() == Integer.MAX_VALUE)
				continue;
			unknownInteriorOrientation.add(unknownParameter);
		}

		for (Image image : camera)
			this.reduceNormalEquationMatrix(neq, image, unknownInteriorOrientation);
	}
	
	private void reduceNormalEquationMatrix(NormalEquationSystem neq, Image image, List<UnknownParameter<InteriorOrientation>> unknownInteriorOrientation) throws IllegalArgumentException, MatrixSingularException {
		UpperSymmPackMatrix N = neq.getMatrix();
		DenseVector n = neq.getVector();
		
		ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
		List<UnknownParameter<ExteriorOrientation>> unknownExteriorOrientation = new ArrayList<UnknownParameter<ExteriorOrientation>>(10);
		
		for (UnknownParameter<ExteriorOrientation> unknownParameter : exteriorOrientation) {
			if (unknownParameter.getColumn() < 0 || unknownParameter.getColumn() == Integer.MAX_VALUE)
				continue;
			unknownExteriorOrientation.add(unknownParameter);
		}
		
		UpperSymmPackMatrix N22 = new UpperSymmPackMatrix(unknownExteriorOrientation.size());
		DenseVector n2 = new DenseVector(N22.numRows());
		for (int rowN22 = 0; rowN22 < N22.numRows(); rowN22++) {
			UnknownParameter<ExteriorOrientation> unknownParameterRow = unknownExteriorOrientation.get(rowN22);
			int rowN = unknownParameterRow.getColumn();
			
			for (int columnN22 = rowN22; columnN22 < N22.numColumns(); columnN22++) {
				UnknownParameter<ExteriorOrientation> unknownParameterColumn = unknownExteriorOrientation.get(columnN22);
				int columnN = unknownParameterColumn.getColumn();
				
				N22.set(rowN22, columnN22, N.get(rowN, columnN));
			}
			n2.set(rowN22, n.get(rowN));
		}
		
		MathExtension.inv(N22);

		List<UnknownParameter<?>> unknownParameters = new ArrayList<UnknownParameter<?>>();
		unknownParameters.addAll(unknownInteriorOrientation);
		
		for (ImageCoordinate imageCoordinate : image) {
			ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
			unknownParameters.add(objectCoordinate.getX());
			unknownParameters.add(objectCoordinate.getY());
			unknownParameters.add(objectCoordinate.getZ());
		}
		
		for (UnknownParameter<?> unknownParameterRow : unknownParameters) {
			int rowN = unknownParameterRow.getColumn();
			DenseVector n12 = new DenseVector(N22.numColumns());
			for (int columnN22 = 0; columnN22 < N22.numColumns(); columnN22++) {
				double dot = 0;
				for (int rowN22 = 0; rowN22 < N22.numRows(); rowN22++) {
					int columnN12 = unknownExteriorOrientation.get(rowN22).getColumn();
					dot += N.get(rowN, columnN12) * N22.get(rowN22, columnN22);
				}
				n12.set(columnN22, dot);
			}
			
			n.add(rowN, -n12.dot(n2));

			for (UnknownParameter<?> unknownParameterColumn : unknownParameters) {
				int columnN = unknownParameterColumn.getColumn();
				double dot = 0;
				for (int rowN22 = 0; rowN22 < N22.numRows(); rowN22++) {
					int rowN12 = unknownExteriorOrientation.get(rowN22).getColumn();
					dot += n12.get(rowN22) * N.get(columnN, rowN12);
				}
				N.add(rowN, columnN, -dot);
			}
		}
	}
	
	public void addPropertyChangeListener(PropertyChangeListener listener) {
		this.change.addPropertyChangeListener(listener);
	}

	public void removePropertyChangeListener(PropertyChangeListener listener) {
		this.change.removePropertyChangeListener(listener);
	}
}
