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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;

import no.uib.cipr.matrix.Matrix;

public class MatrixWriter {
	
	private MatrixWriter() {}

	public static void write(File f, Matrix M) throws IOException {
		write(f, 1.0, M);
	}
	
	public static void write(File f, double scale, Matrix M) throws IOException {
		int columns = M.numColumns();
		int rows = M.numRows();
		
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter( f )));

			for (int i=0; i<rows; i++) {
				for (int j=0; j<columns; j++) {
					pw.printf(Locale.ENGLISH, "%+35.15f  ", scale * M.get(i,j));
				}
				pw.println();	
			}
		} 
		finally {
			if (pw != null) {
				pw.close();
			}
		}
	}
}
