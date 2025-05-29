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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
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
import org.applied_geodesy.adjustment.bundle.parameter.ObservationParameterGroup;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.adjustment.bundle.parameter.UnknownParameter;
import org.applied_geodesy.adjustment.defect.DefectType;
import org.applied_geodesy.adjustment.defect.RankDefect;
import org.applied_geodesy.util.io.writer.AdjustmentResultWritable ;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
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
	
	private final PropertyChangeSupport change = new PropertyChangeSupport(this);
	private EstimationStateType currentEstimationStatus = EstimationStateType.BUSY;
	private EstimationType estimationType = EstimationType.L2NORM;
	private AdjustmentResultWritable  adjustmentResultWriter = null;
	
	private static double SQRT_EPS = Math.sqrt(Constant.EPS);
	private int maximalNumberOfIterations = DefaultValue.getMaximalNumberOfIterations(),
			iterationStep                 = 0,
			numberOfUnknownParameters     = 0,
			numberOfObservations          = 0,
			numberOfInteriorOrientationParameters = 0;

	private boolean interrupt                    = false,
			applyAposterioriVarianceOfUnitWeight = true,
			useCentroidedCoordinates             = true,
			calculateStochasticParameters        = false,
			deriveFirstAdaptedDampingValue       = false;
	
	private MatrixInversion invertNormalEquationMatrix = MatrixInversion.FULL;

	private double maxAbsDx     = Double.MIN_VALUE,
			omega               = 0.0,
			lastValidmaxAbsDx   = 0.0,
			dampingValue        = 0.0,
			adaptedDampingValue = 0.0,
			sigma2apriori       = 1.0;
	
	private ObjectCoordinate centroid = new ObjectCoordinate(this.getClass().getSimpleName(), 0, 0, 0);

	private UpperSymmPackMatrix Qxx = null;

	private RankDefect rankDefect = new RankDefect();

	private List<Camera> cameras = new ArrayList<Camera>();
	private Set<ObservationParameterGroup<?>> observationGroups = new LinkedHashSet<ObservationParameterGroup<?>>();
	private Set<UnknownParameter<?>> unknownParameters = new LinkedHashSet<UnknownParameter<?>>();
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
		
		this.deriveFirstAdaptedDampingValue = this.dampingValue > 0;
		this.adaptedDampingValue = 0;

		this.maxAbsDx          = 0.0;
		this.lastValidmaxAbsDx = 0.0;

		int runs = this.maximalNumberOfIterations-1;
		boolean isEstimated = false, estimateCompleteModel = false, isConverge = true;

		if (this.maximalNumberOfIterations == 0) {
			estimateCompleteModel = isEstimated = true;
			this.adaptedDampingValue = 0;
		}

		this.sigma2apriori = this.sigma2apriori > 0 ? this.sigma2apriori : 1.0;

		this.prepareUnknownParameters();
		if (this.useCentroidedCoordinates)
			this.centroidCoordinates(false);
		
		try {
			do {
				this.maxAbsDx      = 0.0;
				this.iterationStep = this.maximalNumberOfIterations-runs;
				this.currentEstimationStatus = EstimationStateType.ITERATE;
				this.change.firePropertyChange(this.currentEstimationStatus.name(), this.maximalNumberOfIterations, this.iterationStep);

				// create normal equation system
				NormalEquationSystem neq = this.createNormalEquation();

				// apply pre-condition to normal equation system
				NormalEquationSystem.applyPrecondition(neq);

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
						
						// in-place solution of NES, i.e., N <-- Qxx, n <-- dx 
						if (this.invertNormalEquationMatrix == MatrixInversion.REDUCED) {
							int numRows = this.numberOfInteriorOrientationParameters + this.objectCoordinates.size() * 3 + this.rankDefect.getDefect();
							// reduced system of equations
							this.reduceNormalEquationMatrix(neq);
							// in-place solution of NES, i.e., N <-- Qxx, n <-- dx 
							MathExtension.solve(N, n, numRows, Boolean.TRUE);
						}
						else {
							MathExtension.solve(N, n, this.invertNormalEquationMatrix == MatrixInversion.FULL);
						}

						NormalEquationSystem.applyPrecondition(neq.getPreconditioner(), N, dx);
						this.Qxx = N;


						if (this.calculateStochasticParameters) {
							this.currentEstimationStatus = EstimationStateType.ESTIAMTE_STOCHASTIC_PARAMETERS;
							this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
						}
					}
					else {
						// Solve Nx = n in-place, i.e., n <-- dx 
						MathExtension.solve(N, n, false);
						NormalEquationSystem.applyPrecondition(neq.getPreconditioner(), null, dx);
					}

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

				if (Double.isInfinite(this.maxAbsDx) || Double.isNaN(this.maxAbsDx)) {
					this.currentEstimationStatus = EstimationStateType.SINGULAR_MATRIX;
					this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
					return this.currentEstimationStatus;
				}
				else if (this.maxAbsDx <= SQRT_EPS && runs > 0 && this.adaptedDampingValue == 0) {
					isEstimated = true;
					this.currentEstimationStatus = EstimationStateType.CONVERGENCE;
					//if (this.estimationType != EstimationType.UNSCENTED_TRANSFORMATION)
					this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxAbsDx);
				}
				else if (runs-- <= 1) {
					if (estimateCompleteModel) {
						this.currentEstimationStatus = EstimationStateType.NO_CONVERGENCE;
						this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxAbsDx);
						isConverge = false;
					}
					isEstimated = true;
				}
				else {
					this.currentEstimationStatus = EstimationStateType.CONVERGENCE;
					//if (this.estimationType != EstimationType.UNSCENTED_TRANSFORMATION)
					this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxAbsDx);
				}
				
				if (isEstimated || this.adaptedDampingValue <= SQRT_EPS || runs < this.maximalNumberOfIterations * 0.5 + 1)
					this.adaptedDampingValue = 0.0;
			}
			while (!estimateCompleteModel);
			
			if (this.useCentroidedCoordinates)
				this.centroidCoordinates(true);
			
			// export adjustment results to file
			try {
				this.exportAdjustmentResults();
			} catch (NullPointerException | IllegalArgumentException | IOException e) {
				e.printStackTrace();
				this.currentEstimationStatus = EstimationStateType.EXPORT_ADJUSTMENT_RESULTS_FAILED;
				this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
				return this.currentEstimationStatus;
			} 
		}
		catch (OutOfMemoryError e) {
			e.printStackTrace();
			this.currentEstimationStatus = EstimationStateType.OUT_OF_MEMORY;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), false, true);
			return this.currentEstimationStatus;
		}

		if (!isConverge) {
			this.currentEstimationStatus = EstimationStateType.NO_CONVERGENCE;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxAbsDx);
		}
		else if (this.currentEstimationStatus.getId() == EstimationStateType.BUSY.getId() || this.calculateStochasticParameters) {
			this.currentEstimationStatus = EstimationStateType.ERROR_FREE_ESTIMATION;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), SQRT_EPS, this.maxAbsDx);
		}

		return this.currentEstimationStatus;
	}
	
	private void updateModel(Vector dx, boolean updateCompleteModel) {
		if (this.adaptedDampingValue > 0) {
			// reduce the step length
			double alpha = 0.25 * Math.pow(this.adaptedDampingValue, -0.05);
			alpha = Math.min(alpha, 0.75);
			dx.scale(alpha);
			
			double prevOmega = this.omega;
			double curOmega  = this.getOmega(dx);
			prevOmega = prevOmega <= 0 ? Double.MAX_VALUE : prevOmega;
			// Check if the current solution converts - if not, reject
			boolean lmaConverge = prevOmega >= curOmega;
			this.omega = curOmega;
			double lastAdaptedDampingValue = this.adaptedDampingValue;
			if (lmaConverge) {
				this.adaptedDampingValue *= 0.2;
			}
			else {
				this.adaptedDampingValue *= 5.0;
				
				// To avoid infinity
				if (this.adaptedDampingValue > 1.0/SQRT_EPS) {
					this.adaptedDampingValue = 1.0/SQRT_EPS;
					// force an update within the next iteration 
					this.omega = 0.0;
				}
			}
			
			this.currentEstimationStatus = EstimationStateType.LEVENBERG_MARQUARDT_STEP;
			this.change.firePropertyChange(this.currentEstimationStatus.name(), lastAdaptedDampingValue, this.adaptedDampingValue);
			
			if (!lmaConverge) {
				// Current solution is NOT an improvement --> cancel procedure
				this.maxAbsDx = this.lastValidmaxAbsDx;
				dx.zero();
				return;
			}
		}
		
		// Estimate of sum of squared residuals from last adjustment step before parameter update
		if (updateCompleteModel) 
			this.omega = this.estimationType == EstimationType.SIMULATION ? 0.0 : this.getOmega(dx);
		// updating model parameters --> estimated maxAbsDx
		this.maxAbsDx = this.updateUnknownParameters(dx);
		this.lastValidmaxAbsDx = this.maxAbsDx;
		
//		if (updateCompleteModel) {
//			this.omega = 0.0;
//			if (this.estimationType != EstimationType.SIMULATION) {
//				for (ObservationParameterGroup<?> observations : this.observationGroups)
//					this.omega += PartialDerivativeFactory.getWeightedSumOfSquaredResiduals(this.sigma2apriori, observations);
//			}
//		}
	}
	
	/**
	 * Perform update of the model parameters
	 * X = X0 + dx
	 * @param dx
	 * @return max(|dx|)
	 */
	private double updateUnknownParameters(Vector dx) {
		double maxAbsDx = 0;
		for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
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
	
	/**
	 * Estimates the residuals of the misclosures v = A*dx - w
	 * and returns omega = v'*P*v
	 * 
	 * 
	 * @param dx
	 * @return omega = v'*P*v
	 */
	private double getOmega(Vector dx) {
		double omega = 0.0;
		
		for (ObservationParameterGroup<?> observations : this.observationGroups) {
			GaussMarkovEquations APw = PartialDerivativeFactory.getPartialDerivative(this.sigma2apriori, null, new DenseVector(dx.size()), observations);
			if (APw == null)
				continue;
			
			Matrix A = APw.getJacobian();
			Matrix P = APw.getWeights();
			Vector v = APw.getgetMisclosure();
			// multAdd: y = alpha*A*x + y
			A.multAdd(-1, dx, v);
			Vector Pv = new DenseVector(v.size());
			P.mult(v, Pv);

			omega += v.dot(Pv);
		}
		return omega;
	}
	
	private void addDatumConditionRows(UpperSymmPackMatrix N) {
		// center of mass
		double x0 = 0, y0 = 0, z0 = 0;
		int count = 0;
		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
			int columnX = objectCoordinate.getX().getColumn();
			int columnY = objectCoordinate.getY().getColumn();
			int columnZ = objectCoordinate.getZ().getColumn();
			
			if (!objectCoordinate.isDatum() || columnX == Integer.MAX_VALUE || columnY == Integer.MAX_VALUE || columnZ == Integer.MAX_VALUE)
				continue;
			
			x0 += objectCoordinate.getX().getValue();
			y0 += objectCoordinate.getY().getValue();
			z0 += objectCoordinate.getZ().getValue();
			count++;

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
			int columnX = objectCoordinate.getX().getColumn();
			int columnY = objectCoordinate.getY().getColumn();
			int columnZ = objectCoordinate.getZ().getColumn();
			
			if (!objectCoordinate.isDatum() || columnX == Integer.MAX_VALUE || columnY == Integer.MAX_VALUE || columnZ == Integer.MAX_VALUE)
				continue;

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

		// Normalize the columns
		for (ObjectCoordinate objectCoordinate : this.objectCoordinates) {
			int columnX = objectCoordinate.getX().getColumn();
			int columnY = objectCoordinate.getY().getColumn();
			int columnZ = objectCoordinate.getZ().getColumn();
			
			if (!objectCoordinate.isDatum() || columnX == Integer.MAX_VALUE || columnY == Integer.MAX_VALUE || columnZ == Integer.MAX_VALUE)
				continue;

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
	
	private void addObservationGroup(ObservationParameterGroup<?> observations) {
		if (!this.observationGroups.contains(observations)) {
			this.observationGroups.add(observations);
			for (ObservationParameter<?> observation : observations)
				this.sigma2apriori = Math.min(this.sigma2apriori, observation.getVariance());
		}
	}

	private void addUnknownParameter(UnknownParameter<?> unknownParameter) {
		if (unknownParameter.getColumn() == -1 && !this.unknownParameters.contains(unknownParameter)) {
			unknownParameter.setColumn( this.numberOfUnknownParameters++ );
			this.unknownParameters.add( unknownParameter );
		}
	}
	
	public void add(Camera camera) {
		this.cameras.add(camera);
	}

	public void add(ScaleBar scaleBar) {
		this.scaleBars.add(scaleBar);
	}
	
	private void prepareUnknownParameters() {
		for (Camera camera : this.cameras) {		
			for (Image image : camera) {
				for (ImageCoordinate imageCoordinate : image) {
					imageCoordinate.getX().setRow( this.numberOfObservations++ );
					imageCoordinate.getY().setRow( this.numberOfObservations++ );

					ObjectCoordinate objectCoordinate = imageCoordinate.getObjectCoordinate();
					this.objectCoordinates.add(objectCoordinate);

					// add object coordinates as unknowns
					this.addUnknownParameter(objectCoordinate.getX());
					this.addUnknownParameter(objectCoordinate.getY());
					this.addUnknownParameter(objectCoordinate.getZ());

					// add pixel coordinates as observations
					this.addObservationGroup(imageCoordinate);
				}
			}
		}
		
		for (Camera camera : this.cameras) {
			// add parameters of interior orientation
			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			for (UnknownParameter<InteriorOrientation> unknownParameter : interiorOrientation) {
				if (unknownParameter.getColumn() == -1)
					this.numberOfInteriorOrientationParameters++;
				this.addUnknownParameter(unknownParameter);
			}
		}
		
		for (Camera camera : this.cameras) {
			for (Image image : camera) {
				// add parameters of exterior orientation per image
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

			this.addObservationGroup(scaleBar);
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
	 * Create system of normal equations Nx = n, were N contains datum conditions. Moreover
	 * the preconditioner, i.e., the root of the main diagonal of N, is derived
	 * @return Normal equation system
	 */
	public NormalEquationSystem createNormalEquation() {

		UpperSymmPackMatrix N = new UpperSymmPackMatrix( this.numberOfUnknownParameters + this.rankDefect.getDefect() );
		UpperSymmBandMatrix V = new UpperSymmBandMatrix( N.numRows(), 0 );
		DenseVector n = new DenseVector( N.numRows() );
		
		for (ObservationParameterGroup<?> observations : this.observationGroups) {
			PartialDerivativeFactory.getPartialDerivative(this.sigma2apriori, N, n, observations);
		}

		this.addDatumConditionRows(N);
		
		if (this.deriveFirstAdaptedDampingValue) {
//			double maxElement = 0;
//			for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
//				int column = unknownParameter.getColumn();
//				if (column < 0)
//					continue;
//				maxElement = Math.max(maxElement, N.get(column, column));
//			}
			// derive first damping value for LMA
			this.adaptedDampingValue = this.dampingValue;// * maxElement;
			this.deriveFirstAdaptedDampingValue = false;
		}
		
		if (this.adaptedDampingValue > 0) {
			for (UnknownParameter<?> unknownParameter : this.unknownParameters) {
				int column = unknownParameter.getColumn();
				if (column < 0)
					continue;
				N.add(column, column, this.adaptedDampingValue * N.get(column, column));
				//N.add(column, column, this.adaptedDampingValue);
			}
		}
		
		// Preconditioner == Root of the main diagonal
		for (int column = 0; column < N.numColumns(); column++) {
			double value = N.get(column, column);
			V.set(column, column, value > Constant.EPS ? 1.0 / Math.sqrt(value) : 1.0);
		}
		
		if (this.estimationType == EstimationType.SIMULATION)
			n.zero();
		
		return new NormalEquationSystem(N, n, V);
	}

	private void detectRankDefect() {
		boolean hasScaleBars = !this.scaleBars.isEmpty();
		this.rankDefect.reset();
		
		this.rankDefect.setTranslationX(DefectType.FREE);
		this.rankDefect.setTranslationY(DefectType.FREE);
		this.rankDefect.setTranslationZ(DefectType.FREE);
		
		this.rankDefect.setRotationX(DefectType.FREE);
		this.rankDefect.setRotationY(DefectType.FREE);
		this.rankDefect.setRotationZ(DefectType.FREE);
		
		this.rankDefect.setScale(hasScaleBars ? DefectType.FIXED : DefectType.FREE);
		
		int countFixedX = 0;
		int countFixedY = 0;
		int countFixedZ = 0;
		
		for (ObjectCoordinate coordinate : this.objectCoordinates) {
			countFixedX += coordinate.getX().getColumn() == Integer.MAX_VALUE ? 1 : 0;
			countFixedY += coordinate.getY().getColumn() == Integer.MAX_VALUE ? 1 : 0;
			countFixedZ += coordinate.getZ().getColumn() == Integer.MAX_VALUE ? 1 : 0;
			
			if (this.rankDefect.estimateTranslationX() && countFixedX > 0)
				this.rankDefect.setTranslationX(DefectType.FIXED);
			if (this.rankDefect.estimateTranslationY() && countFixedY > 0)
				this.rankDefect.setTranslationY(DefectType.FIXED);
			if (this.rankDefect.estimateTranslationZ() && countFixedZ > 0)
				this.rankDefect.setTranslationZ(DefectType.FIXED);
			
			if (!hasScaleBars && (countFixedX >= 2 || countFixedY >= 2 || countFixedZ >= 2))
				this.rankDefect.setScale(DefectType.FIXED);
			
			//if (countFixedX >= 3 && countFixedY >= 3 && countFixedZ >= 3)
			if (countFixedX > 0 && countFixedY > 0 && countFixedZ > 0 && 
					(hasScaleBars && countFixedX + countFixedY + countFixedZ >= 6 || !hasScaleBars && countFixedX + countFixedY + countFixedZ >= 7)) {
				this.rankDefect.setRotationX(DefectType.FIXED);
				this.rankDefect.setRotationY(DefectType.FIXED);
				this.rankDefect.setRotationZ(DefectType.FIXED);
				
				break;
			}
		}
		
		if (this.rankDefect.estimateTranslationX() || this.rankDefect.estimateTranslationY() || this.rankDefect.estimateTranslationZ() ||
				this.rankDefect.estimateRotationX() || this.rankDefect.estimateRotationY() || this.rankDefect.estimateRotationZ() ||
				this.rankDefect.estimateScale()) {

			for (Camera camera : this.cameras) {
				for (Image image : camera) {
					// add parameters of exterior orientation
					ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
					if (this.rankDefect.estimateRotationX() && exteriorOrientation.get(ParameterType.CAMERA_OMEGA).getColumn() == Integer.MAX_VALUE)
						this.rankDefect.setRotationX(DefectType.FIXED);
					if (this.rankDefect.estimateRotationY() && exteriorOrientation.get(ParameterType.CAMERA_PHI).getColumn() == Integer.MAX_VALUE)
						this.rankDefect.setRotationY(DefectType.FIXED);
					if (this.rankDefect.estimateRotationZ() && exteriorOrientation.get(ParameterType.CAMERA_KAPPA).getColumn() == Integer.MAX_VALUE)
						this.rankDefect.setRotationZ(DefectType.FIXED);

					countFixedX += exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).getColumn() == Integer.MAX_VALUE ? 1 : 0;
					countFixedY += exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).getColumn() == Integer.MAX_VALUE ? 1 : 0;
					countFixedZ += exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).getColumn() == Integer.MAX_VALUE ? 1 : 0;

					if (this.rankDefect.estimateTranslationX() && countFixedX > 0)
						this.rankDefect.setTranslationX(DefectType.FIXED);
					if (this.rankDefect.estimateTranslationY() && countFixedY > 0)
						this.rankDefect.setTranslationY(DefectType.FIXED);
					if (this.rankDefect.estimateTranslationZ() && countFixedZ > 0)
						this.rankDefect.setTranslationZ(DefectType.FIXED);
					
					if (!hasScaleBars && (countFixedX >= 2 || countFixedY >= 2 || countFixedZ >= 2))
						this.rankDefect.setScale(DefectType.FIXED);
					
					if (!this.rankDefect.estimateTranslationX() && !this.rankDefect.estimateTranslationY() && !this.rankDefect.estimateTranslationZ() &&
						!this.rankDefect.estimateRotationX() && !this.rankDefect.estimateRotationY() && !this.rankDefect.estimateRotationZ() &&
						!this.rankDefect.estimateScale())
						break;

					//if (countFixedX >= 3 && countFixedY >= 3 && countFixedZ >= 3)
					if (countFixedX > 0 && countFixedY > 0 && countFixedZ > 0 && 
							(hasScaleBars && countFixedX + countFixedY + countFixedZ >= 6 || !hasScaleBars && countFixedX + countFixedY + countFixedZ >= 7)) {
						this.rankDefect.setRotationX(DefectType.FIXED);
						this.rankDefect.setRotationY(DefectType.FIXED);
						this.rankDefect.setRotationZ(DefectType.FIXED);
						
						break;
					}
				}
				if (!this.rankDefect.estimateTranslationX() && !this.rankDefect.estimateTranslationY() && !this.rankDefect.estimateTranslationZ() &&
						!this.rankDefect.estimateRotationX() && !this.rankDefect.estimateRotationY() && !this.rankDefect.estimateRotationZ() &&
						!this.rankDefect.estimateScale())
						break;
			}
		}
	}
	
	/**
	 * Returns a collection containing the adjusted object coordinates
	 * @return coordinates
	 */
	public Collection<ObjectCoordinate> getObjectCoordinates() {
		return this.objectCoordinates;
	}
	
	/**
	 * Returns the number of observations
	 * @return number of observations
	 */
	public int getNumberOfObservations() {
		return this.numberOfObservations;
	}
	
	/**
	 * Returns the number of adjusted parameters 
	 * @return number of parameters  
	 */
	public int getNumberOfUnknownParameters() {
		return this.numberOfUnknownParameters;
	}
	
	/**
	 * Returns the defect of the normal equation system, i.e. nullity(N)
	 * @return defect of normal equations
	 */
	public int getNumberOfDatumConditions() {
		return this.rankDefect.getDefect();
	}
	
	/**
	 * Returns the overall degree of freedom of the adjustment, i.e., dof = n - u + d
	 * @return degree of freedom
	 */
	public int getDegreeOfFreedom() {
		return this.numberOfObservations - this.numberOfUnknownParameters + this.rankDefect.getDefect();
	}
	
	/**
	 * Returns the a-posteriori variance of the unit weight, if the option applyAposterioriVarianceOfUnitWeight 
	 * is set to true and the estimation type is set to L2NORM, in any other case, the a-priori variance of the 
	 * unit weight is returned 
	 * @return variance of unit weight
	 */
	public double getVarianceFactorAposteriori() {
		int degreeOfFreedom = this.getDegreeOfFreedom();
		return degreeOfFreedom > 0 && this.omega > 0 && this.estimationType != EstimationType.SIMULATION && this.applyAposterioriVarianceOfUnitWeight ? Math.abs(this.omega/(double)degreeOfFreedom) : this.sigma2apriori;
	}
	
	/**
	 * Returns the a-priori variance of the unit weight
	 * @return variance of unit weight
	 */
	public double getVarianceFactorApriori() {
		return this.sigma2apriori;
	}
	
	/**
	 * Returns a collection of the cameras used within the adjustment
	 * @return cameras
	 */
	public Collection<Camera> getCameras() {
		return this.cameras;
	}
	
	/**
	 * Returns a collection of the scale bars
	 * @return scale bars
	 */
	public Collection<ScaleBar> getScaleBars() {
		return this.scaleBars;
	}
	
	/**
	 * Defines the result writer to export specific adjustment results
	 * @param adjustmentResultWriter
	 */
	public void setAdjustmentResultWriter(AdjustmentResultWritable adjustmentResultWriter) {
		this.adjustmentResultWriter = adjustmentResultWriter;
	}
	
	/**
	 * Defines the estimation type, i.e. L2NORM or SIMULATION
	 * @param estimationType
	 * @throws UnsupportedOperationException
	 */
	public void setEstimationType(EstimationType estimationType) throws UnsupportedOperationException {
		if (estimationType == EstimationType.L2NORM || estimationType == EstimationType.SIMULATION)
			this.estimationType = estimationType;
		else
			throw new UnsupportedOperationException(this.getClass().getSimpleName() + " Error, this estimation type is not supported! " + estimationType);
	}
	
	/**
	 * Defines whether the normal equation should be inverted to determine the co-factor matrix. 
	 * If only the adjusted parameters are of interest, this option can save considerable 
	 * computing time. 
	 * 
	 * @param invert
	 */
	public void setInvertNormalEquation(MatrixInversion invert) {
		this.invertNormalEquationMatrix = invert;
	}
	
	/**
	 * Returns whether the normal equation should be inverted or even not
	 * @return invert
	 */
	public MatrixInversion getInvertNormalEquation() {
		return this.invertNormalEquationMatrix;
	}
	
	/**
	 * Writes adjustment results
	 * @throws IOException 
	 * @throws IllegalArgumentException 
	 * @throws NullPointerException 
	 */
	private void exportAdjustmentResults() throws NullPointerException, IllegalArgumentException, IOException {
		if (this.adjustmentResultWriter == null)
			return;

		this.currentEstimationStatus = EstimationStateType.EXPORT_ADJUSTMENT_RESULTS;
		this.change.firePropertyChange(this.currentEstimationStatus.name(), null, this.adjustmentResultWriter.toString());
		this.adjustmentResultWriter.export(this);
	}

	/**
	 * Returns the co-factor matrix, if inversion mode is not set to NONE
	 * @return co-factor
	 */
	public UpperSymmPackMatrix getCofactorMatrix() {
		return (this.invertNormalEquationMatrix == MatrixInversion.NONE) ? null : this.Qxx;
	}

	public void useCentroidedCoordinates(boolean useCentroidedCoordinates) {
		this.useCentroidedCoordinates = useCentroidedCoordinates;
	}
	
	public void applyAposterioriVarianceOfUnitWeight(boolean applyAposterioriVarianceOfUnitWeight) {
		this.applyAposterioriVarianceOfUnitWeight = applyAposterioriVarianceOfUnitWeight;
	}
	
	public void setLevenbergMarquardtDampingValue(double lambda) {
		this.dampingValue = Math.abs(lambda);
	}
	
	public double getLevenbergMarquardtDampingValue() {
		return this.dampingValue;
	}
	
	private void reduceNormalEquationMatrix(NormalEquationSystem neq) throws IllegalArgumentException, MatrixSingularException {
		for (Camera camera : this.cameras)
			this.reduceNormalEquationMatrix(neq, camera);
	}
	
	private void reduceNormalEquationMatrix(NormalEquationSystem neq, Camera camera) throws IllegalArgumentException, MatrixSingularException {
		InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
		
		List<UnknownParameter<InteriorOrientation>> unknownInteriorOrientation = new ArrayList<UnknownParameter<InteriorOrientation>>(13);
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
		List<UnknownParameter<ExteriorOrientation>> unknownExteriorOrientation = new ArrayList<UnknownParameter<ExteriorOrientation>>(6);
		
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
	
	public void interrupt() {
		this.interrupt = true;
	}
	
	public void addPropertyChangeListener(PropertyChangeListener listener) {
		this.change.addPropertyChangeListener(listener);
	}

	public void removePropertyChangeListener(PropertyChangeListener listener) {
		this.change.removePropertyChangeListener(listener);
	}
}
