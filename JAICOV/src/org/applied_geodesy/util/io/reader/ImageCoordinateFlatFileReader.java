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
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;

public class ImageCoordinateFlatFileReader extends SourceFileReader<Camera> {
	private final Map<String, ObjectCoordinate> objectCoordinates;
	private final Camera camera;
	
	public ImageCoordinateFlatFileReader(Camera camera, Map<String, ObjectCoordinate> objectCoordinates) {
		this.objectCoordinates = objectCoordinates;
		this.camera = camera;
		this.reset();
	}
	
	public ImageCoordinateFlatFileReader(String fileName, Camera camera, Map<String, ObjectCoordinate> objectCoordinates) {
		this(new File(fileName).toPath(), camera, objectCoordinates);
	}

	public ImageCoordinateFlatFileReader(File sf, Camera camera, Map<String, ObjectCoordinate> objectCoordinates) {
		this(sf.toPath(), camera, objectCoordinates);
	}
	
	public ImageCoordinateFlatFileReader(Path path, Camera camera, Map<String, ObjectCoordinate> objectCoordinates) {
		super(path);
		this.objectCoordinates = objectCoordinates;
		this.camera = camera;
		this.reset();
	}
	
	@Override
	public void reset() {
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
			
			
			// 1  16      116         -9.457563      2.724932      0.000500    0.000500
			// 1  16      117         -9.919849     -0.616718      0.000500    0.000500
			// 1  16      120         -5.538143      3.031979      0.000500    0.000500
			
			String columns[] = line.split("\\s+");
			if (columns.length < 7)
				return;

			long camid = Long.parseLong(columns[0].trim());
			long imgid = Long.parseLong(columns[1].trim());
			
			if (camid != this.camera.getId())
				return;
			
			String name = columns[2].trim();
			double xp = Double.parseDouble(columns[3].trim());
			double yp = Double.parseDouble(columns[4].trim());
			double sx = Double.parseDouble(columns[5].trim());
			double sy = Double.parseDouble(columns[6].trim());
			double rho = columns.length > 6 ? Double.parseDouble(columns[7].trim()) : 0.0;
			
			Image image = this.camera.add(imgid);				
			if (this.objectCoordinates.containsKey(name)) {
				image.add(this.objectCoordinates.get(name), xp, yp, sx, sy, rho);
			}
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}	
	}
}
