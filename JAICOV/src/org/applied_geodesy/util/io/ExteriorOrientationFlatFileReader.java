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

package org.applied_geodesy.util.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;

import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

public class ExteriorOrientationFlatFileReader extends SourceFileReader<Camera> {
	private final Camera camera;
	private int imgid = 1;
	public ExteriorOrientationFlatFileReader(Camera camera) {
		this.camera = camera;
		this.reset();
	}
	
	public ExteriorOrientationFlatFileReader(String fileName, Camera camera) {
		this(new File(fileName).toPath(), camera);
	}

	public ExteriorOrientationFlatFileReader(File sf, Camera camera) {
		this(sf.toPath(), camera);
	}
	
	public ExteriorOrientationFlatFileReader(Path path, Camera camera) {
		super(path);
		this.camera = camera;
		this.reset();
	}
	
	@Override
	public void reset() {
		this.imgid = 1;
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
			
			
//		     6.9023  -1976.8731  -1394.6327     2.15936712 -0.02770565 -0.03094882
//		    61.7070  -1987.8665  -1376.9822     2.13014527 -0.06227612 -1.57543791
//		     3.9796  -1983.6757  -1410.3555     2.08371542 -0.04995296  1.52097805
//		    75.4013  -1802.6736  -1501.9709     2.13684490 -0.09142014 -3.11116284
			
			String columns[] = line.split("\\s+");
			if (columns.length < 6)
				return;
			
			double X0 = Double.parseDouble(columns[0].trim());
			double Y0 = Double.parseDouble(columns[1].trim());
			double Z0 = Double.parseDouble(columns[2].trim());
			
			double OMEGA = Double.parseDouble(columns[3].trim());
			double PHI   = Double.parseDouble(columns[4].trim());
			double KAPPA = Double.parseDouble(columns[5].trim());
			
			Image image = this.camera.add(this.imgid);
	
			ExteriorOrientation exteriorOrientation = image.getExteriorOrientation();
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_X).setValue(X0);
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Y).setValue(Y0);
			exteriorOrientation.get(ParameterType.CAMERA_COORDINATE_Z).setValue(Z0);
			
			exteriorOrientation.get(ParameterType.CAMERA_OMEGA).setValue(OMEGA);
			exteriorOrientation.get(ParameterType.CAMERA_PHI).setValue(PHI);
			exteriorOrientation.get(ParameterType.CAMERA_KAPPA).setValue(KAPPA);
			
			this.imgid++;
		}
		catch (Exception err) {
			err.printStackTrace();
			// nichts, Beobachtung unbrauchbar...
			return;
		}	
	}
}