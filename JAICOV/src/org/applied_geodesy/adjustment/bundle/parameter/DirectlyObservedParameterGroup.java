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

package org.applied_geodesy.adjustment.bundle.parameter;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.applied_geodesy.adjustment.MathExtension;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixNotSPDException;
import no.uib.cipr.matrix.UpperSPDPackMatrix;
import no.uib.cipr.matrix.UpperSymmBandMatrix;

public class DirectlyObservedParameterGroup implements ObservationParameterGroup<UnknownParameter<?>> {

	private Set<ObservationParameter<? extends UnknownParameter<?>>> observedParameters;
	private UpperSPDPackMatrix weightMatrix = null;
	private double sigma2apriori = -1;
	
	public DirectlyObservedParameterGroup(List<ObservationParameter<? extends UnknownParameter<?>>> observedParameters) throws IllegalArgumentException {
		this.observedParameters = new LinkedHashSet<ObservationParameter<? extends UnknownParameter<?>>>(observedParameters);
		
		if (observedParameters.size() != this.getNumberOfParameters())
			throw new IllegalArgumentException("Error, array contains duplicate observations.");
	}
	
	public DirectlyObservedParameterGroup(UpperSPDPackMatrix dispersionMatrix, List<ObservationParameter<? extends UnknownParameter<?>>> observedParameters) throws IllegalArgumentException {	
		this(observedParameters);
		
		if (dispersionMatrix.numColumns() != this.getNumberOfParameters())
			throw new IllegalArgumentException("Error, number of observations and number of rows/columns in dispersion matrix are unequal: " + this.getNumberOfParameters() + " vs. " + dispersionMatrix.numColumns());
		
		int row = 0;
		for (ObservationParameter<? extends UnknownParameter<?>> observedParameter : this.observedParameters)
			observedParameter.setVariance(dispersionMatrix.get(row, row++)); 
		
		// Do we need a deep copy (matrix will be overwritten during the adjustment process)?
		this.weightMatrix = new UpperSPDPackMatrix(dispersionMatrix, Boolean.FALSE);
	}

	public boolean hasFullyPopulatedWeightMatrix() {
		return this.weightMatrix != null;
	}
	
	public Matrix getWeightMatrix(double sigma2apriori) throws MatrixNotSPDException, IllegalArgumentException {
		if (sigma2apriori <= 0)
			throw new IllegalArgumentException("Error, variance of unit weight must be positive: " + sigma2apriori);
		
		if (!this.hasFullyPopulatedWeightMatrix()) {
			int nop = this.getNumberOfParameters();
			Matrix matrix = new UpperSymmBandMatrix(nop, 0);
			int row = 0;
			for (ObservationParameter<? extends UnknownParameter<?>> observedParameter : this.observedParameters)
				matrix.set(row, row++, sigma2apriori/observedParameter.getVariance());
			
			return matrix;
		}				
		else {
			// Matrix is already inverted (i.e. a weight matrix) but variance of unit weight has changed during iteration
			if (this.sigma2apriori > 0 && this.sigma2apriori != sigma2apriori)
				this.weightMatrix.scale(this.sigma2apriori/sigma2apriori);

			if (this.sigma2apriori < 0)
				MathExtension.inv((UpperSPDPackMatrix)this.weightMatrix.scale(1.0/sigma2apriori));

			this.sigma2apriori = sigma2apriori;
			return this.weightMatrix;
		}
	}
	
	@Override
	public int getNumberOfParameters() {
		return this.observedParameters.size();
	}

	@Override
	public Iterator<ObservationParameter<? extends UnknownParameter<?>>> iterator() {
		return observedParameters.iterator();
	}
}
