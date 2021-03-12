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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.applied_geodesy.adjustment.bundle.BundleAdjustment;
import org.applied_geodesy.adjustment.bundle.Camera;
import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.ScaleBar;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

public class AICONReportFileReader extends SourceFileReader<BundleAdjustment> {
	
	private enum ContentType {
		INTERIOR_ORIENTATION,
		EXTERIOR_ORIENTATION,
		OBJECT_COORDINATES,
		IMAGE_COORDINATES,
		SCALE_BARS,
		UNDEFINED;
	}

	private ContentType currentContentType = ContentType.UNDEFINED;
	private Camera camera = null;
	private Image image = null;
	private Map<Integer, Camera> cameras = new LinkedHashMap<Integer, Camera>();
	private Map<Integer, Image> images = new LinkedHashMap<Integer, Image>();
	private List<ScaleBar> scaleBars = new ArrayList<ScaleBar>();
	private Map<String, ObjectCoordinate> objectCoordinates = new LinkedHashMap<String, ObjectCoordinate>();
	private Map<String, ObjectCoordinate> datumCoordinates = new HashMap<String, ObjectCoordinate>();
	public AICONReportFileReader() {
		this.reset();
	}
	
	public AICONReportFileReader(String fileName) {
		this(new File(fileName).toPath(), null);
	}

	public AICONReportFileReader(File sf) {
		this(sf.toPath(), null);
	}
	
	public AICONReportFileReader(Path path) {
		this(path, null);
	}
	
	public AICONReportFileReader(String fileName, Map<String, ObjectCoordinate> datumCoordinates) {
		this(new File(fileName).toPath(), datumCoordinates);
	}

	public AICONReportFileReader(File sf, Map<String, ObjectCoordinate> datumCoordinates) {
		this(sf.toPath(), datumCoordinates);
	}
	
	public AICONReportFileReader(Path path, Map<String, ObjectCoordinate> datumCoordinates) {
		super(path);
		this.reset();
		if (datumCoordinates != null)
			this.datumCoordinates = datumCoordinates;
	}
	
	@Override
	public void reset() {
		if (this.cameras == null)
			this.cameras = new LinkedHashMap<Integer, Camera>();
		if (this.objectCoordinates == null)
			this.objectCoordinates = new LinkedHashMap<String, ObjectCoordinate>();
		if (this.images == null)
			this.images = new LinkedHashMap<Integer, Image>();
		if (this.scaleBars == null)
			this.scaleBars = new ArrayList<ScaleBar>();

		this.cameras.clear();
		this.objectCoordinates.clear();
		this.images.clear();
		this.scaleBars.clear();
		this.currentContentType = ContentType.UNDEFINED;
	}

	@Override
	public BundleAdjustment readAndImport() throws IOException, SQLException {
		this.reset();
		super.read();
		
		BundleAdjustment adjustment = new BundleAdjustment();
		for (Camera camera : this.cameras.values())
			adjustment.add(camera);
		for (ScaleBar scaleBar : this.scaleBars)
			adjustment.add(scaleBar);

		return adjustment;
	}

	@Override
	public void parse(String line) {
		line = line.trim();
		try {
			if (line.contains("#Start") || line.contains("zum Anfang"))
				this.currentContentType = ContentType.UNDEFINED;
			
			if (line.contains("name=\"interior_orientations\"") || line.contains("*** Innere Orientierungen ***"))
				this.currentContentType = ContentType.INTERIOR_ORIENTATION;
			
			if (line.contains("name=\"exterior_orientations\"") || line.contains("ussere Orientierungen ***"))
				this.currentContentType = ContentType.EXTERIOR_ORIENTATION;
			
			if (line.contains("name=\"object_points\"") || line.contains("*** Objektpunkte ***"))
				this.currentContentType = ContentType.OBJECT_COORDINATES;
			
			if (line.contains("name=\"image_coordinates\"") || line.contains("*** Bildkoordinaten ***"))
				this.currentContentType = ContentType.IMAGE_COORDINATES;
			
			if (line.contains("name=\"distances\"") || line.contains("*** Strecken ***"))
				this.currentContentType = ContentType.SCALE_BARS;
			
			
			
			switch(this.currentContentType) {
			case INTERIOR_ORIENTATION:
				this.readInteriorOrientationParameter(line);
				break;
			case EXTERIOR_ORIENTATION:
				this.readExteriorOrientationParameter(line);
				break;
			case OBJECT_COORDINATES:
				this.readObjectCoordinates(line);
				break;
			case IMAGE_COORDINATES:
				this.readImageCoordinates(line);
				break;
			case SCALE_BARS:
				this.readScaleBars(line);
				break;
			default:
				break;
			}
		}
		catch (Exception err) {
			err.printStackTrace();
			// nichts, Beobachtung unbrauchbar...
			return;
		}	
	}
	
	private void readScaleBars(String line) {
	    //506     507        789.8480    -0.0117     0.1771     0.0100       3.71  ---
		//510     511        789.7940     0.0117     0.1771     0.0100       3.71  ---
		
		String regexp = new String("^\\w+\\s+\\w+\\s+[\\d\\.+-]+.+");
		if (!line.matches(regexp))
			return;
		
		String columns[] = line.split("\\s+");

		if (columns.length < 7)
			return;
		
		String nameA = columns[0].trim();
		String nameB = columns[1].trim();
		
		if (!this.objectCoordinates.containsKey(nameA) || !this.objectCoordinates.containsKey(nameB))
			return;
		
		double value = Double.parseDouble(columns[2]);
		double sigma = Double.parseDouble(columns[5]);
		
		ScaleBar scaleBar = new ScaleBar(this.objectCoordinates.get(nameA), this.objectCoordinates.get(nameB), value, sigma);
		this.scaleBars.add(scaleBar);
	}
	
	private void readImageCoordinates(String line) {
		   //116      2         -0.555783     -8.932040   -0.001288    0.007369    0.000500    0.000500   1.00   1.00  3.44 19.67  ***
		   //117      2          3.852806     -8.677704   -0.021392    0.013159    0.000500    0.000500   1.00   1.00 57.10 35.13  ***
		   //115      1        -11.228082      1.363158   -0.000100   -0.000288    0.000500    0.000500   0.91   0.94  0.28  0.79 
		   //116      1         -7.280413     -0.104988    0.000015   -0.000155    0.000500    0.000500   0.93   0.94  0.04  0.43 
		if (line.endsWith("***"))
			return;
		
		String regexp = new String("^\\w+\\s+\\d+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+");
		if (!line.matches(regexp))
			return;
		
		String columns[] = line.split("\\s+");
		if (columns.length != 12)
			return;
		
		String name = columns[0].trim();
		int imgId = Integer.parseInt(columns[1]);
		
		if (!this.objectCoordinates.containsKey(name) || !this.images.containsKey(imgId))
			return;
		
		double xp = Double.parseDouble(columns[2]);
		double yp = Double.parseDouble(columns[3]);
		double sx = Double.parseDouble(columns[6]);
		double sy = Double.parseDouble(columns[7]);
		
		Image image = this.images.get(imgId);
		image.add(this.objectCoordinates.get(name), xp, yp, sx, sy);		
	}
	
	private void readObjectCoordinates(String line) {
	     //115        -876.3998       3.2146      90.0435    0.0164    0.0105    0.0100       22        1
	     //116        -523.4996       4.8431     -71.4074    0.0129    0.0090    0.0085       28        3
	     //117        -505.8577     197.4650    -371.5378    0.0128    0.0088    0.0084       30        3
		
		String regexp = new String("^\\w+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+\\d+\\s+\\d+");
		if (!line.matches(regexp))
			return;
		
		String columns[] = line.split("\\s+");
		if (columns.length != 9)
			return;
		
		String name = columns[0].trim();
		double x = Double.parseDouble(columns[1].trim());
		double y = Double.parseDouble(columns[2].trim());
		double z = Double.parseDouble(columns[3].trim());
		
		ObjectCoordinate coordinate = new ObjectCoordinate(name, x, y, z);
		coordinate.setDatum(this.datumCoordinates == null || this.datumCoordinates.isEmpty());

		if (this.datumCoordinates != null && !this.datumCoordinates.isEmpty() && this.datumCoordinates.containsKey(name))
			coordinate = this.datumCoordinates.get(name);
		
		this.objectCoordinates.put(name, coordinate);
	}
	
	private void readExteriorOrientationParameter(String line) {
		//1            1      6.9023  -1976.8731  -1394.6327     0.0313     0.0413     0.0356        313
		String regexpXYZ   = new String("^\\d+\\s+\\d+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+\\d+");
		//air  rad    2.15936712 -0.02770565 -0.03094882   0.000023   0.000023   0.000035   0.000320   0.000297
		String regexpAngle = new String("^air\\s+rad\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.+-]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+[\\d\\.]+");
		
		if (line.matches(regexpXYZ)) {
			String columns[] = line.split("\\s+");
			Camera camera = this.cameras.get(Integer.parseInt(columns[1]));
			if (camera == null)
				return;

			int imgId = Integer.parseInt(columns[0]);
			this.image = camera.add(imgId);
			this.image.getExteriorOrientation().get(ParameterType.CAMERA_COORDINATE_X).setValue(Double.parseDouble(columns[2]));
			this.image.getExteriorOrientation().get(ParameterType.CAMERA_COORDINATE_Y).setValue(Double.parseDouble(columns[3]));
			this.image.getExteriorOrientation().get(ParameterType.CAMERA_COORDINATE_Z).setValue(Double.parseDouble(columns[4]));
			this.images.put(imgId, this.image);
		}
		else if (this.image != null && line.matches(regexpAngle)) {
			String columns[] = line.split("\\s+");

			this.image.getExteriorOrientation().get(ParameterType.CAMERA_OMEGA).setValue(Double.parseDouble(columns[2]));
			this.image.getExteriorOrientation().get(ParameterType.CAMERA_PHI).setValue(Double.parseDouble(columns[3]));
			this.image.getExteriorOrientation().get(ParameterType.CAMERA_KAPPA).setValue(Double.parseDouble(columns[4]));
		}
	}
	
	private void readInteriorOrientationParameter(String line) {
		if (!line.contains(":"))
			return;
		
		String columns[] = line.split("[:\\s]+");
		if (columns.length != 3)
			return;
		
		String type = columns[0].trim();
		if (type.equalsIgnoreCase("Kamera/R0")) {
			int camId = Integer.parseInt(columns[1].trim());
			double r0 = Double.parseDouble(columns[2].trim());
			this.camera = new Camera(camId, r0);
			this.cameras.put(camId, this.camera);
		}
		if (this.camera == null)
			return;
		
		double value  = Double.parseDouble(columns[1].trim());
		boolean fixed = columns[2].equalsIgnoreCase("fest");
		
		switch(type) {
			case "Ck":
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_DISTANCE).setValue(-value);
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_DISTANCE).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "Xh":
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_POINT_X).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_POINT_X).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "Yh":
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_POINT_Y).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.PRINCIPAL_POINT_Y).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "A1":
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A1).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A1).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "A2":
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A2).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A2).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "A3":
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A3).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.RADIAL_DISTORTION_A3).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "B1":
				this.camera.getInteriorOrientation().get(ParameterType.TANGENTIAL_DISTORTION_B1).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.TANGENTIAL_DISTORTION_B1).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "B2":
				this.camera.getInteriorOrientation().get(ParameterType.TANGENTIAL_DISTORTION_B2).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.TANGENTIAL_DISTORTION_B2).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "C1":
				this.camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C1).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C1).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
			case "C2":
				this.camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C2).setValue(value);
				this.camera.getInteriorOrientation().get(ParameterType.AFFINITY_AND_SHEAR_C2).setColumn(fixed ? Integer.MAX_VALUE : -1);
				break;
		}
	}
}
