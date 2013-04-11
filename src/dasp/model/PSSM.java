/*
 File: PSSM.java

 Copyright (c) 2009, Wake Forest University

 This library is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published
 by the Free Software Foundation; either version 3.0 of the License, or
 any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 documentation provided hereunder is on an "as is" basis, and 
 neither Wake Forest University nor the University of Califonia
 have any obligations to provide maintenance, support,
 updates, enhancements or modifications.  In no event shall
 Wake Forest University or the University of Califonia
 be liable to any party for direct, indirect, special,
 incidental or consequential damages, including lost profits, arising
 out of the use of this software and its documentation, even if the
 Wake Forest University and/or the University of Califonia
 have been advised of the possibility of such damage.  See
 the GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

package dasp.model;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

/**
 * A PSSM records the frequency of amino acids at a particular point
 * in an alignment, given the normal frequency of amino acids in the PDB
 * This is taken from Ryan Huff's theses (see page 20):
 *
 * foreach column j in the alignment:
 *   record the observed counts of each residue, <i>i</i>: <i>q<sub>i</sub></i>
 *   record the pseudocounts of each residue, i: p<sub>i</sub> = b X f<sub>i</sub>,
 *          where f<sub>i</sub> is the overall frequency of residue i in the database,
 *          and b is the pseudocount weight (currently set to .1)
 *   foreach residue i:
 *      compute the score: S<sub>i</sub> = (q<sub>i</sub> + p<sub>i</sub>) / (n + b),
 *             where n is the number of sequences in the alignment
 *   enter the logarithm (base 2) of S<sub>i</sub> in cell [i, j] of the PSSM
 *
 * See also: Bioinformatics: Sequence and Genome Analysis. David W. Mount, 2001
 */
public class PSSM implements Comparable {
	private double PSSMmatrix[][] = null;
	private Alignment alignment;
	private static final double PSEUDOCOUNT_WEIGHT = 0.1;
	private static final String aaList = "ACDEFGHIKLMNPQRSTVWY";
	
	/**
	 * Steps to getting the PSSM object loaded
	 * 1) load the amino acid frequencies from the PDB database
	 *	Hashtable aaFreqTable = loadAAfreq(aaFreqFile, codebase);
	 * 2) create the PSSMs for each alignment object containing a profile fragment/motif
	 * LinkedList pssmList = createPSSMs(msaList.size(), id, numPDBs, aaFreqTable);
	 * 3) sort the PSSMs by length -- have to search database using longest motif to shortest
	 * sort(pssmList);
	 */
	public PSSM (Alignment profileFragAlignment) {
		alignment = profileFragAlignment;
		PSSMmatrix = new double [20][alignment.getAlignmentWidth()];
		createPSSM();
	}

	public double[][] getPSSM() {
		return PSSMmatrix;
	}

	public int getWidth() {
		if (PSSMmatrix == null) return -1;

		return PSSMmatrix[0].length;
	}

	public String toString() {
		String pssm = "";
		for(int aa = 0; aa < 20; aa++) {
			pssm += numToAA(aa)+"\t";
			for(int position = 0; position < getWidth(); position++){
				pssm += ""+PSSMmatrix[aa][position]+"\t";
			}
			pssm += "\n";
		}
		return pssm;
	}
	
	//create the PSSMs for each motif
	/* Modified by Amy Olex 6/24/09
	 * This method supposedly calculates the PSSM for each motif of the ASP.
	 * I am modifying this method so that the PSSM for each motif is printed to a file.
	 * Each motif has been aligned by Clustal and is in a .aln file.
	 **/
	/**
	 * Calculates the PSSM for each profile fragment alignment
	 * NEEDS LOTS OF WORK...Just plopped in code.
	 *
	 */
	private void createPSSM() {
		int alignmentLength = alignment.getAlignmentWidth();
		int numSeqs = alignment.getAlignment().size();

		//initilize the PSSM matrix to zero. k = rows (AAs), m = columns(motif positions)
		for(int aa = 0; aa < 20; aa++)
			for(int position = 0; position < alignmentLength; position++)
				PSSMmatrix[aa][position] = 0;

		//ok, we are looping through each row of each column, so the first motif position is loop through all the way down followed by the second etc. Instead of looping through each column of each row.
		for (String alignedString: alignment.getAlignment()) {
			String upperString = alignedString.toUpperCase();
			for(int position = 0; position < alignmentLength; position++) {
				char aa = upperString.charAt(position);
				//if the cur char is not a -, Z or X then increment the corresponding PSSM location in that column
				if (isAmino(aa))
					PSSMmatrix[getAANum(aa)][position]++;
			}
		}

		//back to k = AAs and m = motif positions
		for(int aa = 0; aa < 20; aa++) {
			for(int position = 0; position < alignmentLength; position++){
				NumberFormat nf = NumberFormat.getNumberInstance();
				nf.setMaximumFractionDigits(0);
				//b is a global variable and is equal to 0.1
				PSSMmatrix[aa][position] += PSEUDOCOUNT_WEIGHT*aaFreq[aa];
				PSSMmatrix[aa][position] /= ((double)(PSEUDOCOUNT_WEIGHT + numSeqs));
				double temp = Math.log(PSSMmatrix[aa][position]);
				PSSMmatrix[aa][position] = Double.parseDouble(nf.format(temp));
			}
		}
	}

	public static boolean isAmino(char aa) {
		if (aa == '-' || aa == 'Z' || aa == 'X')
			return false;

		return true;
	}

	public static int getAANum(char aa) {
		return aaList.indexOf(aa);
	}

	public static char numToAA(int anum) {
		return aaList.charAt(anum);
	}

    public String getProfileFragAlignment(){
        return alignment.getAlignmentAsString();
    }


	private static final double[] aaFreq = {
				0.07843351606835683,	// ALA
				0.015193800942076749,	// CYS
				0.0582383486327026,		// ASP
				0.06565599729916563,	// GLU
				0.04031959873318007,	// PHE
				0.07247078115203447,	// GLY
				0.023592109705319678,	// HIS
				0.05868044949600502,	// ILE
				0.0605951481439802,		// LYS
				0.08932848897963121,	// LEU
				0.022564828062955164,	// MET
				0.044698808739128335,	// ASN
				0.04393357233573943,	// PRO
				0.03837596257415237,	// GLN
				0.04933041814703471,	// ARG
				0.05985724161214089,	// SER
				0.057519733774898316,	// THR
				0.07181807951384982,	// VAL
				0.013960741443338746,	// TRP
				0.035430767004806844	// TYR
	};

    public int compareTo(Object pssm) {
        int thisLength = this.getWidth();
        PSSM cmp = (PSSM)pssm;
        int cmpLength = cmp.getWidth();

        if(thisLength > cmpLength)
            return 1;
        else if(thisLength < cmpLength)
            return -1;

        return 0;
    }
		
}
