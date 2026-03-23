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

package org.applied_geodesy.util.io.reader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;

import org.applied_geodesy.adjustment.bundle.camera.Camera;
import org.applied_geodesy.adjustment.bundle.camera.orientation.InteriorOrientation;

public class InteriorOrientationFlatFileReader extends SourceFileReader<Camera> {
	private final Camera camera;
	public InteriorOrientationFlatFileReader(Camera camera) {
		this.camera = camera;
		this.reset();
	}
	
	public InteriorOrientationFlatFileReader(String fileName, Camera camera) {
		this(new File(fileName).toPath(), camera);
	}

	public InteriorOrientationFlatFileReader(File sf, Camera camera) {
		this(sf.toPath(), camera);
	}
	
	public InteriorOrientationFlatFileReader(Path path, Camera camera) {
		super(path);
		this.camera = camera;
		this.reset();
	}
	
	@Override
	public void reset() {}

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
		    // 1 7.244192e-003	1.164305e-001	2.868147e+001
			
			String columns[] = line.split("\\s+");
			if (columns.length < 4)
				return;
			
			long camid = Long.parseLong(columns[0].trim());

			if (camid != this.camera.getId())
				throw new IllegalArgumentException("Error, camera-id mismatch: " + this.camera.getId() + " vs. " + camid + "!");

			int col = 1;
			double x0 = Double.parseDouble(columns[col++].trim());
			double y0 = Double.parseDouble(columns[col++].trim());
			double c  = Double.parseDouble(columns[col++].trim());
			
			InteriorOrientation interiorOrientation = this.camera.getInteriorOrientation();
			interiorOrientation.getPrinciplePointX().setValue(x0);
			interiorOrientation.getPrinciplePointY().setValue(y0);
			interiorOrientation.getPrincipleDistance().setValue(c);
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}	
	}
}
