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
import java.util.LinkedHashMap;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.util.io.reader.SourceFileReader;

public class OBCFileReader extends SourceFileReader<Map<String, ObjectCoordinate>> {
	private Map<String, ObjectCoordinate> coordinates = new LinkedHashMap<String, ObjectCoordinate>();
	
	public OBCFileReader() {
		this.reset();
	}
	
	public OBCFileReader(String fileName) {
		this(new File(fileName).toPath());
	}

	public OBCFileReader(File sf) {
		this(sf.toPath());
	}
	
	public OBCFileReader(Path path) {
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
			//  1. Spalte Punktnummer
			//  2. Spalte Objektkoordinate X
			//  3. Spalte Objektkoordinate Y
			//  4. Spalte Objektkoordinate Z
			//  5. Spalte Standardabweichung in X
			//  6. Spalte Standardabweichung in Y
			//  7. Spalte Standardabweichung in Z
			//  8. Spalte Strahlen = Anzahl der Bilder, in denen der Punkt gemessen wurde.
			//  9. Spalte Status des Bildes 1 = aktiv; 0 = inaktiv
			// 10. Spalte Parameter NEU; Bündelausgleichung: Neupunkt wenn ungleich 0
			// 11. Spalte Datumspunkt, wenn ungleich 0

			String columns[] = line.split("\\s+");
			if (columns.length < 4)
				return;
			
			boolean enable = columns.length < 11 || !columns[8].trim().equalsIgnoreCase("0");
			
			if (!enable)
				return;
			
			String name = columns[0].trim(); 
			double x = Double.parseDouble(columns[1].trim());
			double y = Double.parseDouble(columns[2].trim());
			double z = Double.parseDouble(columns[3].trim());
			
			coordinate = new ObjectCoordinate(name, x, y, z);
			this.coordinates.put(name, coordinate);
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}	
	}
}