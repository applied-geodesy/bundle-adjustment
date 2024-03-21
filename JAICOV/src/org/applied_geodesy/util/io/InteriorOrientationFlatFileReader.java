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
import org.applied_geodesy.adjustment.bundle.orientation.InteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

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
		    // 1 7.244192e-003	1.164305e-001	2.868147e+001	-1.097712e-004	1.535086e-007	0.000000e+000	8.273977e-006	-1.054657e-005	-7.008010e-005	-3.126270e-005	0.000000e+000	0.000000e+000	0.000000e+000
			
			String columns[] = line.split("\\s+");
			if (columns.length < 14)
				return;
			
			long camid = Long.parseLong(columns[0].trim());

			if (camid != this.camera.getId())
				throw new IllegalArgumentException("Error, camera-id mismatch: " + this.camera.getId() + " vs. " + camid + "!");

			int col = 1;
			double x0 = Double.parseDouble(columns[col++].trim());
			double y0 = Double.parseDouble(columns[col++].trim());
			double c  = Double.parseDouble(columns[col++].trim());

			double A1 = Double.parseDouble(columns[col++].trim());
			double A2 = Double.parseDouble(columns[col++].trim());
			double A3 = Double.parseDouble(columns[col++].trim());

			double B1 = Double.parseDouble(columns[col++].trim());
			double B2 = Double.parseDouble(columns[col++].trim());

			double C1 = Double.parseDouble(columns[col++].trim());
			double C2 = Double.parseDouble(columns[col++].trim());
			
			double D1 = Double.parseDouble(columns[col++].trim());
			double D2 = Double.parseDouble(columns[col++].trim());
			double D3 = Double.parseDouble(columns[col++].trim());

			InteriorOrientation interiorOrientation = camera.getInteriorOrientation();
			interiorOrientation.get(ParameterType.PRINCIPAL_POINT_X).setValue(x0);
			interiorOrientation.get(ParameterType.PRINCIPAL_POINT_Y).setValue(y0);
			interiorOrientation.get(ParameterType.PRINCIPAL_DISTANCE).setValue(c);

			interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A1).setValue(A1);
			interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A2).setValue(A2);
			interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).setValue(A3);
			interiorOrientation.get(ParameterType.RADIAL_DISTORTION_A3).setColumn(Integer.MAX_VALUE);

			interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B1).setValue(B1);
			interiorOrientation.get(ParameterType.TANGENTIAL_DISTORTION_B2).setValue(B2);

			interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).setValue(C1);
			interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C1).setColumn(Integer.MAX_VALUE);
			interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).setValue(C2);
			interiorOrientation.get(ParameterType.AFFINITY_AND_SHEAR_C2).setColumn(Integer.MAX_VALUE);
			
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).setValue(D1);
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D1).setColumn(Integer.MAX_VALUE);
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).setValue(D2);
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D2).setColumn(Integer.MAX_VALUE);
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).setValue(D3);
			interiorOrientation.get(ParameterType.DISTANCE_DISTORTION_D3).setColumn(Integer.MAX_VALUE);
		}
		catch (Exception err) {
			err.printStackTrace();
			return;
		}	
	}
}
