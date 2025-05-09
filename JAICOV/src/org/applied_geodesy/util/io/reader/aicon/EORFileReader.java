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

import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;
import org.applied_geodesy.util.io.reader.SourceFileReader;

public class EORFileReader extends SourceFileReader<Camera> {
	private final Camera camera;

	public EORFileReader(Camera camera) {
		this.camera = camera;
		this.reset();
	}
	
	public EORFileReader(String fileName, Camera camera) {
		this(new File(fileName).toPath(), camera);
	}

	public EORFileReader(File sf, Camera camera) {
		this(sf.toPath(), camera);
	}
	
	public EORFileReader(Path path, Camera camera) {
		super(path);
		this.camera = camera;
		this.reset();
	}
	
	@Override
	public void reset() {	}

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
			//  1. Spalte Bildnummer
			//  2. Spalte Kameranummer
			//  3. Spalte X- Koordinate
			//  4. Spalte Y- Koordinate
			//  5. Spalte Z- Koordinate
			//  6. Spalte Drehwinkel Omega
			//  7. Spalte Drehwinkel Phi
			//  8. Spalte Drehwinkel Kappa
			//  9. Spalte Drehreihenfolge 0 = CAP-Rotation
			// 10. Spalte Status des Bildes 
			//       0 = nicht aktiv, > 0 = aktiv
			// 11. Spalte Status der Bildorientierung
			//       1 = Bild nicht orientiert
			//       2 = Kamerastandort aus Vororientierung 
			//       3 = Kamerastandort aus Bündelausgleichung
			//
			
			String columns[] = line.split("\\s+");
			if (columns.length < 11)
				return;
			
			long camid = Long.parseLong(columns[1].trim());
			
			boolean capRotation =  columns[8].trim().equalsIgnoreCase("0");
			boolean enable      = !columns[9].trim().equalsIgnoreCase("0");
			boolean orient      = !columns[10].trim().equalsIgnoreCase("1");

			if (!enable || !capRotation || !orient || camid != this.camera.getId())
				return;
			
			long imgid = Long.parseLong(columns[0].trim());

			double X0 = Double.parseDouble(columns[2].trim());
			double Y0 = Double.parseDouble(columns[3].trim());
			double Z0 = Double.parseDouble(columns[4].trim());
			
			double OMEGA = Double.parseDouble(columns[5].trim());
			double PHI   = Double.parseDouble(columns[6].trim());
			double KAPPA = Double.parseDouble(columns[7].trim());

			Image image = this.camera.add(imgid);
			
			ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).setValue(X0);
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).setValue(Y0);
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).setValue(Z0);

			exteriorOrientation.get(ParameterType.CAMERA_OMEGA).setValue(OMEGA);
			exteriorOrientation.get(ParameterType.CAMERA_PHI).setValue(PHI);
			exteriorOrientation.get(ParameterType.CAMERA_KAPPA).setValue(KAPPA);
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}	
	}
}