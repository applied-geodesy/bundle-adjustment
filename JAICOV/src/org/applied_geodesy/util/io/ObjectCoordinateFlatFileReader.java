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
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;

public class ObjectCoordinateFlatFileReader extends SourceFileReader<Map<String, ObjectCoordinate>> {
	private Map<String, ObjectCoordinate> coordinates = new LinkedHashMap<String, ObjectCoordinate>();
	
	public ObjectCoordinateFlatFileReader() {
		this.reset();
	}
	
	public ObjectCoordinateFlatFileReader(String fileName) {
		this(new File(fileName).toPath());
	}

	public ObjectCoordinateFlatFileReader(File sf) {
		this(sf.toPath());
	}
	
	public ObjectCoordinateFlatFileReader(Path path) {
		super(path);
		this.reset();
	}
	
	@Override
	public void reset() {
		if (this.coordinates == null) 
			this.coordinates = new LinkedHashMap<String, ObjectCoordinate>();
		
		this.coordinates.clear();
	}

	@Override
	public Map<String, ObjectCoordinate> readAndImport() throws IOException, SQLException {
		this.reset();
		this.ignoreLinesWhichStartWith("#");
		
		super.read();
		return this.coordinates;
	}

	@Override
	public void parse(String line) {
		line = line.trim();
		
		ObjectCoordinate coordinate = null;
		try {
			
			String columns[] = line.split("\\s+");
			if (columns.length < 4)
				return;
			
			String name = columns[0].trim(); 
			double x = Double.parseDouble(columns[1].trim());
			double y = Double.parseDouble(columns[2].trim());
			double z = Double.parseDouble(columns[3].trim());
			
			boolean datum = columns.length > 4 && columns[4].trim().equalsIgnoreCase("1");
			
			coordinate = new ObjectCoordinate(name, x, y, z);
			coordinate.setDatum(datum);
			this.coordinates.put(name, coordinate);
		}
		catch (Exception err) {
			err.printStackTrace();
			// nichts, Beobachtung unbrauchbar...
			return;
		}	
	}
}
