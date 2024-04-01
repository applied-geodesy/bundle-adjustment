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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.ScaleBar;

public class ScaleBarFlatFileReader extends SourceFileReader<List<ScaleBar>> {
	private final Map<String, ObjectCoordinate> objectCoordinates;
	private List<ScaleBar> scaleBars = new ArrayList<ScaleBar>();
	
	public ScaleBarFlatFileReader(Map<String, ObjectCoordinate> objectCoordinates) {
		this.objectCoordinates = objectCoordinates;
		this.reset();
	}
	
	public ScaleBarFlatFileReader(String fileName, Map<String, ObjectCoordinate> objectCoordinates) {
		this(new File(fileName).toPath(), objectCoordinates);
	}

	public ScaleBarFlatFileReader(File sf, Map<String, ObjectCoordinate> objectCoordinates) {
		this(sf.toPath(), objectCoordinates);
	}
	
	public ScaleBarFlatFileReader(Path path, Map<String, ObjectCoordinate> objectCoordinates) {
		super(path);
		this.objectCoordinates = objectCoordinates;
		this.reset();
	}
	
	@Override
	public void reset() {
		if (this.scaleBars == null)
			this.scaleBars = new ArrayList<ScaleBar>();
		
		this.scaleBars.clear();
	}

	@Override
	public List<ScaleBar> readAndImport() throws IOException, SQLException {
		this.reset();
		this.ignoreLinesWhichStartWith("#");
		
		super.read();
		return this.scaleBars;
	}

	@Override
	public void parse(String line) {
		line = line.trim();
		
		ScaleBar scaleBar = null;
		try {
			
			
//			506     507        789.8480    0.0100
//		    510     511        789.7940    0.0100
			
			String columns[] = line.split("\\s+");
			if (columns.length < 4)
				return;
			
			String nameA = columns[0].trim();
			String nameB = columns[1].trim();
	
			if (!this.objectCoordinates.containsKey(nameA) || !this.objectCoordinates.containsKey(nameB))
				return;
			
			double length = Double.parseDouble(columns[2].trim());
			double sigma  = Double.parseDouble(columns[3].trim());
			
			scaleBar = new ScaleBar(this.objectCoordinates.get(nameA), this.objectCoordinates.get(nameB),length,sigma);
			
			this.scaleBars.add(scaleBar);
		}
		catch (Exception err) {
			err.printStackTrace();
			// nichts, Beobachtung unbrauchbar...
			return;
		}	
	}
}