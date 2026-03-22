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

package org.applied_geodesy.util.io.reader.aicon;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;

import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.camera.distortion.AffinityShearDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.DistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.RadiallySymmetricDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.distortion.TangentialDistortionModel;
import org.applied_geodesy.adjustment.bundle.camera.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.PolynomialCoefficient;
import org.applied_geodesy.util.io.reader.SourceFileReader;

public class IORFileReader extends SourceFileReader<Camera> {
	private Camera camera = null;
	private InteriorOrientation interiorOrientation;
	private int lineCounter = 0;
	private final int[] LINE_LENGHTS = new int[] {8, 1, 2, 2, 4};
	
	public IORFileReader(String fileName) {
		this(new File(fileName).toPath());
	}

	public IORFileReader(File sf) {
		this(sf.toPath());
	}
	
	public IORFileReader(Path path) {
		super(path);
		this.reset();
	}
	
	@Override
	public void reset() {
		this.lineCounter = 0;
	}

	@Override
	public Camera readAndImport() throws IOException, SQLException {
		this.reset();
		this.ignoreLinesWhichStartWith("#");
		
		super.read();
		return this.camera;
	}

	@Override
	public void parse(String line) {
		line = line.trim();
		try {
			// 1. Zeile 
			//    1. Spalte Kameranummer
			//    2. Spalte Interner Parameter
			//    3. Spalte Kamerakonstante ck
			//    4. Spalte Hauptpunktlage xh
			//    5. Spalte Hauptpunktlage yh
			//    6. Spalte radialsymmetrische Verzeichnung A1
			//    7. Spalte radialsymmetrische Verzeichnung A2
			//    8. Spalte zweiter Nulldurchgang der Verzeichnungskurve R0
			// 2. Zeile 
			//    1. Spalte radialsymmetrische Verzeichnung A3
			// 3. Zeile 
			//    1. Spalte radialasym. und tang. Verzeichnung B1
			//    2. Spalte radialasym. und tang. Verzeichnung B2
			// 4. Zeile 
			//    1. Spalte Scherung und Affinität C1
			//    2. Spalte Scherung und Affinität C2
			// 5. Zeile 
			//    1. Spalte Breite des Sensors [mm]
			//    2. Spalte Höhe des Sensors [mm]
			//    3. Spalte Anzahl Pixel in x -Richtung (Breite)
			//    4. Spalte Anzahl Pixel in y -Richtung (Höhe)
			
			String columns[] = line.split("\\s+");
			if (this.lineCounter >= LINE_LENGHTS.length || columns.length < LINE_LENGHTS[this.lineCounter] || this.lineCounter > 0 && this.camera == null)
				return;
			
			DistortionModel distortionModel = null;
			PolynomialCoefficient<?> polynomialCoefficient = null;
			
			switch(this.lineCounter++) {
			case 0: // 1. Zeile
				long camid = Long.parseLong(columns[0].trim());

				double c  = Double.parseDouble(columns[2].trim());
				double x0 = Double.parseDouble(columns[3].trim());
				double y0 = Double.parseDouble(columns[4].trim());

				double A1 = Double.parseDouble(columns[5].trim());
				double A2 = Double.parseDouble(columns[6].trim());
				
				double r0 = Double.parseDouble(columns[7].trim());
				
				this.camera = new Camera(camid, r0, DistortionModel.Type.RADIAL_DISTORTION, DistortionModel.Type.TANGENTIAL_DISTORTION, DistortionModel.Type.AFFINITY_AND_SHEAR);
				this.interiorOrientation = camera.getInteriorOrientation();
				
				this.interiorOrientation.getPrincipleDistance().setValue(-c);
				this.interiorOrientation.getPrincipleDistance().setColumn(-1);

				this.interiorOrientation.getPrinciplePointX().setValue(x0);
				this.interiorOrientation.getPrinciplePointX().setColumn(-1);

				this.interiorOrientation.getPrinciplePointY().setValue(y0);
				this.interiorOrientation.getPrinciplePointY().setColumn(-1);
				
				distortionModel = this.camera.getDistortionModel(DistortionModel.Type.RADIAL_DISTORTION);
				
				polynomialCoefficient = ((RadiallySymmetricDistortionModel)distortionModel).add(1);
				polynomialCoefficient.setValue(A1);
				polynomialCoefficient.setColumn(-1);
				
				polynomialCoefficient = ((RadiallySymmetricDistortionModel)distortionModel).add(2);
				polynomialCoefficient.setValue(A2);
				polynomialCoefficient.setColumn(-1);
				
				break;
			case 1: // 2. Zeile
				double A3 = Double.parseDouble(columns[0].trim());
				
				distortionModel = this.camera.getDistortionModel(DistortionModel.Type.RADIAL_DISTORTION);
				
				polynomialCoefficient = ((RadiallySymmetricDistortionModel)distortionModel).add(3);
				polynomialCoefficient.setValue(A3);
				polynomialCoefficient.setColumn(-1);

				break;
			case 2: // 3. Zeile
				double Bx = Double.parseDouble(columns[0].trim());
				double By = Double.parseDouble(columns[1].trim());
				
				distortionModel = this.camera.getDistortionModel(DistortionModel.Type.TANGENTIAL_DISTORTION);

				((TangentialDistortionModel)distortionModel).getBx().setValue(Bx);
				((TangentialDistortionModel)distortionModel).getBx().setColumn(-1);
				
				((TangentialDistortionModel)distortionModel).getBy().setValue(By);
				((TangentialDistortionModel)distortionModel).getBy().setColumn(-1);
				
				break;
			case 3: // 4. Zeile
				double Cx = Double.parseDouble(columns[0].trim());
				double Cy = Double.parseDouble(columns[1].trim());
				
				distortionModel = this.camera.getDistortionModel(DistortionModel.Type.AFFINITY_AND_SHEAR);
				
				((AffinityShearDistortionModel)distortionModel).getCx().setValue(Cx);
				((AffinityShearDistortionModel)distortionModel).getCx().setColumn(-1);
				
				((AffinityShearDistortionModel)distortionModel).getCy().setValue(Cy);
				((AffinityShearDistortionModel)distortionModel).getCy().setColumn(-1);

				break;
			}
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}
	}
}