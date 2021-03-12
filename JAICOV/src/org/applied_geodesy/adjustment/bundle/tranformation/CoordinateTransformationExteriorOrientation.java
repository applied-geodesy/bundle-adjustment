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

package org.applied_geodesy.adjustment.bundle.tranformation;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.applied_geodesy.adjustment.bundle.Image;
import org.applied_geodesy.adjustment.bundle.ObjectCoordinate;
import org.applied_geodesy.adjustment.bundle.orientation.ExteriorOrientation;
import org.applied_geodesy.adjustment.bundle.parameter.ParameterType;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;

public class CoordinateTransformationExteriorOrientation {
	private static CoordinateTransformationExteriorOrientation trans = new CoordinateTransformationExteriorOrientation();
	private UpperSymmPackMatrix covariance;
	private List<ObjectCoordinate> transformedCoordinates;
	
	private CoordinateTransformationExteriorOrientation() {}

	public static CoordinateTransformationExteriorOrientation getInstance() {
		return trans;
	}

	public void transform(Set<ObjectCoordinate> objectCoordinatesToTransform, Map<Image, ArrayList<Image>> imagesToAlign, double sigma2, UpperSymmPackMatrix CoVar) {
		//int numberOfObjectCoordinates = 3 * objectCoordinatesToTransform.size();
		int rows = 0;
		int columns = 3 * objectCoordinatesToTransform.size();
		for (Map.Entry<Image, ArrayList<Image>> entry : imagesToAlign.entrySet()) {
			//Image referenceImage = entry.getKey();
			columns += 6;
			ArrayList<Image> images = entry.getValue();
			columns += 6 * images.size();
			for (Image image : images) {
				for (ObjectCoordinate objectCoordinateSrc : objectCoordinatesToTransform) {
					// skip points not visible in current image
					if (image.get(objectCoordinateSrc) == null) 
						continue;
					rows += 3;
				}
			}
			//rows += numberOfObjectCoordinates * images.size();
		}
		columns = CoVar.numColumns();
		LinkedSparseMatrix J = new LinkedSparseMatrix(rows, columns);
//		DenseMatrix J = new DenseMatrix(rows, columns);
		this.transformedCoordinates = new ArrayList<ObjectCoordinate>(rows/3);
		int row = 0;
		for (Map.Entry<Image, ArrayList<Image>> entry : imagesToAlign.entrySet()) {
			Image referenceImage = entry.getKey();
			ArrayList<Image> images = entry.getValue();
			ExteriorOrientation exteriorOrientationTrg = referenceImage.getExteriorOrientation();

			for (Image image : images) {
				ExteriorOrientation exteriorOrientationSrc = image.getExteriorOrientation();
				//System.out.println("Ref: " + referenceImage.getId()+"; Img: " + image.getId());
				for (ObjectCoordinate objectCoordinateSrc : objectCoordinatesToTransform) {
					// skip points not visible in current image
					if (image.get(objectCoordinateSrc) == null) {
						//System.err.println("Note: point " + objectCoordinateSrc.getName() + " not visible in image " + image.getId() +"!");
						continue;
					}
					
					int rowsInJacobian[] = new int[] {
							row++,
							row++,
							row++,
					};
					String name = objectCoordinateSrc.getName() + " " + image.getId() + " " + referenceImage.getId();
					ObjectCoordinate objectCoordinateTrg = setPartialDerivations(name, rowsInJacobian, J, objectCoordinateSrc, exteriorOrientationTrg, exteriorOrientationSrc);
					this.transformedCoordinates.add(objectCoordinateTrg);
					//System.out.println(objectCoordinateTrg.getName()+"  "+objectCoordinateTrg.getX().getValue()+"  "+objectCoordinateTrg.getY().getValue()+"  "+objectCoordinateTrg.getZ().getValue());
				}
				
			}
		}
//		System.out.println("J " + J.numRows() + "x" + J.numColumns());
//		LinkedSparseMatrix JC = new LinkedSparseMatrix(J.numRows(), CoVar.numRows());
//		J.mult(CoVar, JC);
//		System.out.println("JCJT");		
//		this.covariance = new UpperSymmPackMatrix(J.numRows());
//		JC.transBmult(sigma2, J, this.covariance);

		//Matrix CJT = new DenseMatrix(CoVar.numRows(), J.numRows());
		LinkedSparseMatrix CJT = new LinkedSparseMatrix(CoVar.numRows(), J.numRows());
		CoVar.transBmult(J, CJT);

		this.covariance = new UpperSymmPackMatrix(J.numRows());
		J.mult(sigma2, CJT, this.covariance);
//		
//		try {
//			org.applied_geodesy.util.io.MatrixWriter.write(new java.io.File("j.txt"), J);
//		} catch (java.io.IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}

	public Matrix getCovarianceMatrix() {
		return this.covariance;
	}
	
	public List<ObjectCoordinate> getTransformedCoordinates() {
		return this.transformedCoordinates;
	}
	
	private static ObjectCoordinate setPartialDerivations(String name, int rows[], Matrix J, ObjectCoordinate objectCoordinateSrc, ExteriorOrientation exteriorOrientationTrg, ExteriorOrientation exteriorOrientationSrc) {
		boolean isReferenceImage = exteriorOrientationTrg == exteriorOrientationSrc;
		
		int rowX = rows[0];
		int rowY = rows[1];
		int rowZ = rows[2];
		
		double XiSrc = objectCoordinateSrc.getX().getValue();
		double YiSrc = objectCoordinateSrc.getY().getValue();
		double ZiSrc = objectCoordinateSrc.getZ().getValue();
		
		int columnXiSrc = objectCoordinateSrc.getX().getColumn();
		int columnYiSrc = objectCoordinateSrc.getY().getColumn();
		int columnZiSrc = objectCoordinateSrc.getZ().getColumn();
		
		ObjectCoordinate objectCoordinateTrg = null;
		name = name == null || name.isBlank() ? objectCoordinateSrc.getName() : name;
		
		if (isReferenceImage) {
			J.set(rowX, columnXiSrc, 1.0);
			J.set(rowY, columnYiSrc, 1.0);
			J.set(rowZ, columnZiSrc, 1.0);

			objectCoordinateTrg = new ObjectCoordinate(name, XiSrc, YiSrc, ZiSrc);
			objectCoordinateTrg.getX().setColumn(rowX);
			objectCoordinateTrg.getY().setColumn(rowY);
			objectCoordinateTrg.getZ().setColumn(rowZ);
		}
		else {

			double X0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_X).getValue();
			double Y0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_Y).getValue();
			double Z0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_Z).getValue();

			double omegaTrg = exteriorOrientationTrg.get(ParameterType.CAMERA_OMEGA).getValue();
			double phiTrg   = exteriorOrientationTrg.get(ParameterType.CAMERA_PHI).getValue();
			double kappaTrg = exteriorOrientationTrg.get(ParameterType.CAMERA_KAPPA).getValue();

			int columnX0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
			int columnY0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
			int columnZ0Trg = exteriorOrientationTrg.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();

			int columnOmegaTrg = exteriorOrientationTrg.get(ParameterType.CAMERA_OMEGA).getColumn();
			int columnPhiTrg   = exteriorOrientationTrg.get(ParameterType.CAMERA_PHI).getColumn();
			int columnKappaTrg = exteriorOrientationTrg.get(ParameterType.CAMERA_KAPPA).getColumn();

			double cosOmegaTrg = Math.cos(omegaTrg);
			double sinOmegaTrg = Math.sin(omegaTrg);

			double cosPhiTrg = Math.cos(phiTrg);
			double sinPhiTrg = Math.sin(phiTrg);

			double cosKappaTrg = Math.cos(kappaTrg);
			double sinKappaTrg = Math.sin(kappaTrg);

			double r11Trg =  cosPhiTrg * cosKappaTrg;
			double r12Trg = -cosPhiTrg * sinKappaTrg;
			double r13Trg =  sinPhiTrg;

			double r21Trg =  cosOmegaTrg*sinKappaTrg + sinOmegaTrg*sinPhiTrg*cosKappaTrg;
			double r22Trg =  cosOmegaTrg*cosKappaTrg - sinOmegaTrg*sinPhiTrg*sinKappaTrg;
			double r23Trg = -sinOmegaTrg*cosPhiTrg;

			double r31Trg = sinOmegaTrg*sinKappaTrg - cosOmegaTrg*sinPhiTrg*cosKappaTrg;
			double r32Trg = sinOmegaTrg*cosKappaTrg + cosOmegaTrg*sinPhiTrg*sinKappaTrg;
			double r33Trg = cosOmegaTrg*cosPhiTrg;



			double X0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_X).getValue();
			double Y0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_Y).getValue();
			double Z0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_Z).getValue();

			double omegaSrc = exteriorOrientationSrc.get(ParameterType.CAMERA_OMEGA).getValue();
			double phiSrc   = exteriorOrientationSrc.get(ParameterType.CAMERA_PHI).getValue();
			double kappaSrc = exteriorOrientationSrc.get(ParameterType.CAMERA_KAPPA).getValue();

			int columnX0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_X).getColumn();
			int columnY0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_Y).getColumn();
			int columnZ0Src = exteriorOrientationSrc.get(ParameterType.CAMERA_COORDINATE_Z).getColumn();

			int columnOmegaSrc = exteriorOrientationSrc.get(ParameterType.CAMERA_OMEGA).getColumn();
			int columnPhiSrc   = exteriorOrientationSrc.get(ParameterType.CAMERA_PHI).getColumn();
			int columnKappaSrc = exteriorOrientationSrc.get(ParameterType.CAMERA_KAPPA).getColumn();

			double cosOmegaSrc = Math.cos(omegaSrc);
			double sinOmegaSrc = Math.sin(omegaSrc);

			double cosPhiSrc = Math.cos(phiSrc);
			double sinPhiSrc = Math.sin(phiSrc);

			double cosKappaSrc = Math.cos(kappaSrc);
			double sinKappaSrc = Math.sin(kappaSrc);

			double r11Src =  cosPhiSrc * cosKappaSrc;
			double r12Src = -cosPhiSrc * sinKappaSrc;
			double r13Src =  sinPhiSrc;

			double r21Src =  cosOmegaSrc*sinKappaSrc + sinOmegaSrc*sinPhiSrc*cosKappaSrc;
			double r22Src =  cosOmegaSrc*cosKappaSrc - sinOmegaSrc*sinPhiSrc*sinKappaSrc;
			double r23Src = -sinOmegaSrc*cosPhiSrc;

			double r31Src = sinOmegaSrc*sinKappaSrc - cosOmegaSrc*sinPhiSrc*cosKappaSrc;
			double r32Src = sinOmegaSrc*cosKappaSrc + cosOmegaSrc*sinPhiSrc*sinKappaSrc;
			double r33Src = cosOmegaSrc*cosPhiSrc;




			double dXi = XiSrc - X0Src;
			double dYi = YiSrc - Y0Src;
			double dZi = ZiSrc - Z0Src;

			double dXiTrg = r11Src * dXi + r21Src * dYi + r31Src * dZi;
			double dYiTrg = r12Src * dXi + r22Src * dYi + r32Src * dZi;
			double dZiTrg = r13Src * dXi + r23Src * dYi + r33Src * dZi;

			double XiTrg = X0Trg + r11Trg * dXiTrg + r12Trg * dYiTrg + r13Trg * dZiTrg;
			double YiTrg = Y0Trg + r21Trg * dXiTrg + r22Trg * dYiTrg + r23Trg * dZiTrg;
			double ZiTrg = Z0Trg + r31Trg * dXiTrg + r32Trg * dYiTrg + r33Trg * dZiTrg;


			objectCoordinateTrg = new ObjectCoordinate(name, XiTrg, YiTrg, ZiTrg);
			objectCoordinateTrg.getX().setColumn(rowX);
			objectCoordinateTrg.getY().setColumn(rowY);
			objectCoordinateTrg.getZ().setColumn(rowZ);

			// params = [X0Trg Y0Trg Z0Trg  omegaTrg phiTrg kappaTrg  X0Src Y0Src Z0Src omegaSrc phiSrc kappaSrc  XiSrc YiSrc ZiSrc];
			/** X-Komponente **/
			// Referenzbild (Zielsytstem)
			J.set(rowX, columnX0Trg,     1.0);
			J.set(rowX, columnY0Trg,     0.0);
			J.set(rowX, columnZ0Trg,     0.0);
			J.set(rowX, columnOmegaTrg,  0.0);
			J.set(rowX, columnPhiTrg,    cosPhiTrg*(dXi*r13Src+dYi*r23Src+dZi*r33Src)-cosKappaTrg*r13Trg*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)+r13Trg*sinKappaTrg*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src));
			J.set(rowX, columnKappaTrg, -r11Trg*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)-cosPhiTrg*sinKappaTrg*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src));
			// Referenzbild (Startsystem)
			J.set(rowX, columnX0Src,    -r11Src*r11Trg-r13Src*r13Trg+cosPhiTrg*r12Src*sinKappaTrg);
			J.set(rowX, columnY0Src,    -r11Trg*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+cosPhiTrg*sinKappaTrg*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+cosPhiSrc*r13Trg*sinOmegaSrc);
			J.set(rowX, columnZ0Src,    -r11Trg*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)-r13Trg*r33Src+cosPhiTrg*sinKappaTrg*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc));
			J.set(rowX, columnOmegaSrc, -r13Trg*(dYi*r33Src+cosPhiSrc*dZi*sinOmegaSrc)-r11Trg*(dYi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)-dZi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc))-r12Trg*(dYi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)-dZi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)));
			J.set(rowX, columnPhiSrc,   -r11Trg*(cosKappaSrc*dXi*r13Src+cosOmegaSrc*dZi*r11Src-dYi*r11Src*sinOmegaSrc)+r13Trg*(cosPhiSrc*dXi-cosOmegaSrc*dZi*r13Src+dYi*r13Src*sinOmegaSrc)-cosPhiTrg*sinKappaTrg*(dXi*r13Src*sinKappaSrc+dYi*r12Src*sinOmegaSrc+dZi*r33Src*sinKappaSrc));
			J.set(rowX, columnKappaSrc,  r11Trg*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)-r12Trg*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src));
			// Objektkoordinate
			J.set(rowX, columnXiSrc,     r11Src*r11Trg+r13Src*r13Trg+cosPhiSrc*cosPhiTrg*sinKappaSrc*sinKappaTrg);
			J.set(rowX, columnYiSrc,     r11Trg*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+r12Trg*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+r13Trg*r23Src);
			J.set(rowX, columnZiSrc,     r11Trg*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+r12Trg*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+r13Trg*r33Src);

			/** Y-Komponente **/
			// Referenzbild (Zielsytstem)
			J.set(rowY, columnX0Trg,     0.0);
			J.set(rowY, columnY0Trg,     1.0);
			J.set(rowY, columnZ0Trg,     0.0);
			J.set(rowY, columnOmegaTrg, -(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)-(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)-r33Trg*(dXi*r13Src+dYi*r23Src+dZi*r33Src));
			J.set(rowY, columnPhiTrg,    r13Trg*sinOmegaTrg*(dXi*r13Src+dYi*r23Src+dZi*r33Src)+r11Trg*sinOmegaTrg*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)-cosPhiTrg*sinOmegaTrg*sinKappaTrg*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src));
			J.set(rowY, columnKappaTrg, -(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)+(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src));
			// Referenzbild (Startsystem)
			J.set(rowY, columnX0Src,    -r11Src*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)+cosPhiSrc*sinKappaSrc*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)+cosPhiTrg*r13Src*sinOmegaTrg);
			J.set(rowY, columnY0Src,    -(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)-(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)+cosPhiTrg*r23Src*sinOmegaTrg);
			J.set(rowY, columnZ0Src,    -(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)-(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)+cosPhiTrg*r33Src*sinOmegaTrg);
			J.set(rowY, columnOmegaSrc, -(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)*(dYi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)-dZi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc))-(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)*(dYi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)-dZi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc))-r23Trg*(dYi*r33Src+cosPhiSrc*dZi*sinOmegaSrc));
			J.set(rowY, columnPhiSrc,   -(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)*(cosKappaSrc*dXi*r13Src+cosOmegaSrc*dZi*r11Src-dYi*r11Src*sinOmegaSrc)+(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)*(dXi*r13Src*sinKappaSrc+dYi*r12Src*sinOmegaSrc+dZi*r33Src*sinKappaSrc)-cosPhiTrg*sinOmegaTrg*(cosPhiSrc*dXi-cosOmegaSrc*dZi*r13Src+dYi*r13Src*sinOmegaSrc));
			J.set(rowY, columnKappaSrc,  (cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)-(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src));
			// Objektkoordinate
			J.set(rowY, columnXiSrc,     r11Src*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)+r12Src*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)+r13Src*r23Trg);
			J.set(rowY, columnYiSrc,     (cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)+(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)+cosPhiSrc*cosPhiTrg*sinOmegaSrc*sinOmegaTrg);
			J.set(rowY, columnZiSrc,     r23Trg*r33Src+(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)*(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)+(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)*(cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg));

			/** Z-Komponente **/
			// Referenzbild (Zielsytstem)
			J.set(rowZ, columnX0Trg,     0.0);
			J.set(rowZ, columnY0Trg,     0.0);
			J.set(rowZ, columnZ0Trg,     1.0);
			J.set(rowZ, columnOmegaTrg,  (cosOmegaTrg*sinKappaTrg+cosKappaTrg*r13Trg*sinOmegaTrg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)+(cosOmegaTrg*cosKappaTrg-r13Trg*sinOmegaTrg*sinKappaTrg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src)-cosPhiTrg*sinOmegaTrg*(dXi*r13Src+dYi*r23Src+dZi*r33Src));
			J.set(rowZ, columnPhiTrg,   -cosOmegaTrg*r13Trg*(dXi*r13Src+dYi*r23Src+dZi*r33Src)-cosOmegaTrg*r11Trg*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)-cosOmegaTrg*r12Trg*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src));
			J.set(rowZ, columnKappaTrg,  (cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)-(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src));
			// Referenzbild (Startsystem)
			J.set(rowZ, columnX0Src,    -r11Src*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)-r13Src*r33Trg+cosPhiSrc*sinKappaSrc*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg));
			J.set(rowZ, columnY0Src,    -(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)-(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)+cosPhiSrc*r33Trg*sinOmegaSrc);
			J.set(rowZ, columnZ0Src,    -r33Src*r33Trg-(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)-(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg));
			J.set(rowZ, columnOmegaSrc, -(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)*(dYi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)-dZi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc))-(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)*(dYi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)-dZi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc))-r33Trg*(dYi*r33Src+cosPhiSrc*dZi*sinOmegaSrc));
			J.set(rowZ, columnPhiSrc,    r33Trg*(cosPhiSrc*dXi-cosOmegaSrc*dZi*r13Src+dYi*r13Src*sinOmegaSrc)-(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)*(cosKappaSrc*dXi*r13Src+cosOmegaSrc*dZi*r11Src-dYi*r11Src*sinOmegaSrc)+(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)*(dXi*r13Src*sinKappaSrc+dYi*r12Src*sinOmegaSrc+dZi*r33Src*sinKappaSrc));
			J.set(rowZ, columnKappaSrc, -(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)*(dYi*(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)+dZi*(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)+dXi*r11Src)+(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)*(dYi*(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)+dZi*(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)+dXi*r12Src));
			// Objektkoordinate
			J.set(rowZ, columnXiSrc,     r11Src*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)+r12Src*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)+r13Src*r33Trg);
			J.set(rowZ, columnYiSrc,     r23Src*r33Trg+(cosOmegaSrc*sinKappaSrc+cosKappaSrc*r13Src*sinOmegaSrc)*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg)+(cosOmegaSrc*cosKappaSrc-r13Src*sinOmegaSrc*sinKappaSrc)*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg));
			J.set(rowZ, columnZiSrc,     r33Src*r33Trg+(cosKappaSrc*sinOmegaSrc+cosOmegaSrc*r13Src*sinKappaSrc)*(cosKappaTrg*sinOmegaTrg+cosOmegaTrg*r13Trg*sinKappaTrg)+(sinOmegaSrc*sinKappaSrc-cosOmegaSrc*cosKappaSrc*r13Src)*(sinOmegaTrg*sinKappaTrg-cosOmegaTrg*cosKappaTrg*r13Trg));
		}
		return objectCoordinateTrg;
	}
	
}
