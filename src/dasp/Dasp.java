/*
 File: Dasp.java

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

package dasp;

import java.lang.Double;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import dasp.algorithms.FastAFileSearcher;
import dasp.algorithms.DBSearch;
import dasp.algorithms.ClustalAlign;
import dasp.algorithms.RyansPSSMSearch;
import dasp.model.ActiveSiteProfile;
import dasp.model.ActiveSiteSignature;
import dasp.model.Alignment;
import dasp.model.DBSearchResult;
import dasp.model.PSSM;
import dasp.model.SearchResult;
import java.util.Collection;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The DASP mainline.  Dasp expects as input a file which contains a series of PDB identifiers
 * and the "key" residues, a cutoff for the radius around each key residue to search,
 * and a database where the PDB files are stored.
 *
 *		Basic flow of the DASP algorithm:
 *		1. Read our input file and create the Active Site Signatures
 *		2. Align the Active Site Signatures to create our Active Site Profile
 *		3. Split the Active Site Profile into the individual fragments
 *		4. Realign the fragments
 *		5. Create the Position-Specific Scoring Matrices
 *		6. Search the sequence database using the PSSMs
 *		7. Combine the results from the searches
 *
 * Where the input file to DASP to dasp is of the form
 * 	pdbID#chain: residueList
 * where residueList is a comma-speparated list of residues
 */
public class Dasp {
	
	// Parameterize these
	// private static String database = "nr.gz";
	private static String databaseDir = "/databases/mol/blast/db/";
	private static String database = "nr";
	private static String pdbDatabase = null;
	private static String inputFile = null;
	private static double cutoff = 1e-50;
	private static double radius = 10;
	private static String outputFile = null;
	private static boolean vFlag = false;
	private static File profilePath = null;

	public Dasp () {
	}
	
	/**
	 * This is the method that will be called when Dasp is run from
	 * the command line.  The command-line arguments are:
	 *
	 * <b>-i</b> <i>filename</i>	The name of the input file
	 * <b>-r</b> <i>radius</i>	The radius for inclusion into the active site signature
	 * <b>-c</b> <i>cutoff</i>	The cutoff value for searching the sequence database
	 * <b>-p</b> <i>pdbDatabase</i>	The path to the pdb database
	 * <b>-o</b> <i>filename</i>	The file to wite the active stie profile into
	 * <b>-d</b> <i>database</i>	The database to use for the sequence search
	 * <b>-h</b> the help text
	 */
	public static void main(String[] args) throws InterruptedException{
		File pdbPath = null;
	
		//Before starting algorithms, check for args
		//If there are args then there have to be on two or three args
		//First is always a file name, second is always the number of lines
		//and the third is always a new profile radius.

		GetOpt opts = new GetOpt(args, "i:c:p:o:d:r:P:hv");

		int result;
		while ((result = opts.getopt()) >= 0) {
			switch((char)result) {
			case 'i':
				inputFile = opts.optArg;
				break;

			case 'r':
				try {
					radius = new Double(opts.optArg);
				} catch (Exception e) {
					System.err.println("Radius argument must be a float");
					System.exit(1);
				}
				break;

			case 'c':
				try {
					cutoff = new Double(opts.optArg);
				} catch (Exception e) {
					System.err.println("Cutoff argument must be a float");
					System.exit(1);
				}
				break;

			case 'o':
				outputFile = opts.optArg;
				break;

			case 'd':
				database = opts.optArg;
				try {
					File f = new File(database);
					if (!f.isAbsolute()) {
						database = databaseDir+database;
						f = new File(database);
					}
				} catch (Exception e) {
					System.err.println("Unable to open database file '"+database+"': "+e.getMessage());
					System.exit(1);
				}
				break;

			case 'p':
				pdbDatabase = opts.optArg;
				try {
					pdbPath = new File(pdbDatabase);
					if (!pdbPath.isDirectory())
						throw new Exception("PDB argument must point to a directory");
				} catch (Exception e) {
					System.err.println("Unable to open directory '"+pdbDatabase+"': "+e.getMessage());
					System.exit(1);
				}
			case 'P':
				String proFile = opts.optArg;
				try {
					profilePath = new File(proFile);
				} catch (Exception e) {
					System.err.println("Unable to open profile '"+proFile+"': "+e.getMessage());
					System.exit(1);
				}

			case 'v':
				vFlag = true;
				break;

			case 'h':
				usage();
				System.exit(0);

			default:
				usage();
				System.exit(0);
			}
		}

		if (inputFile == null && profilePath == null) {
			System.err.println("Must provide an input file!");
			usage();
			System.exit(1);
		}

		PrintStream outputStream = System.out;
		if (outputFile != null) {
			try {
				outputStream = new PrintStream(outputFile);
			} catch (Exception e) {
				System.err.println("Unable to output file '"+outputFile+"': "+e.getMessage());
				System.exit(1);
			}
		}

		File dbFile = new File(database);

		List<ActiveSiteSignature> asSigList = new ArrayList();
		ActiveSiteProfile	profile;

		if (profilePath == null) {
			// 	1. Read our input file and create the Active Site Signatures
			try {
				BufferedReader reader = new BufferedReader(new FileReader(inputFile));
				String line = null;
				while ((line = reader.readLine()) != null) {
					ActiveSiteSignature sig = new ActiveSiteSignature(line, pdbPath);
					asSigList.add(sig);
					if (vFlag) {
						System.out.println("Active Site Signature for: "+line.trim()+" is "+sig.getSignature(radius));
					}
				}
			} catch (Exception e) {
				System.err.println("Unable to read input file "+inputFile+": "+e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}

			// 	2. Align the Active Site Signatures to create our Active Site Profile
			//  After this call profile will contain the ASP score, alignment and identity score
			profile = new ActiveSiteProfile(asSigList, new ClustalAlign(), radius);
			if (vFlag) {
				System.out.println("ActiveSiteProfile: \n"+profile.getAlignmentAsString());
			}
		} else {
			profile = new ActiveSiteProfile(profilePath);
			if (vFlag) {
				System.out.println("ActiveSiteProfile: \n"+profile.getAlignmentAsString());
			}
		}

		// 	3. Split the Active Site Profile into the individual fragments
		List<Alignment>fragmentList = profile.findProfileFragments(new ClustalAlign());

		// 	4. Realign the fragments --if all of the structures from step 3 are set up correctly then this should work just fine.
		for (Alignment a: fragmentList) { 
			a.doAlign(); 
			if (vFlag) {
				//System.out.println("Fragment alignment: \n"+a.getAlignmentAsString());
			}
		}

		// 	5. Create the Position-Specific Scoring Matrices
		List<PSSM>pssmList = new ArrayList();
		for (Alignment a: fragmentList) { 
			PSSM pssm = new PSSM(a);
			pssmList.add(pssm); 
		}

		Collections.sort(pssmList); //sorts shortest to longest on pssm width
		Collections.reverse(pssmList);

		// 	6. Search the sequence database using the PSSMs
		DBSearch searcher = new FastAFileSearcher();
		try {
			List<DBSearchResult> searchResults = searcher.search(dbFile, pssmList, new RyansPSSMSearch(), cutoff);

			Collections.sort(searchResults);  //search results are sorted on pvalue

			outputStream.println("\n\n"+searchResults.size()+" Total Database Search Results\n");
			outputStream.print("Pvalue\tSeqName\tPseudosig\t");
			for(int i=0; i<pssmList.size(); i++){
				outputStream.print("Index"+i+"\tMatchingSubSeq"+i+"\t");
			}
			for(DBSearchResult r: searchResults){
				outputStream.print("\n"+r.getPval()+"\t"+r.getName()+"\t"+r.getPseudoSig()+"\t");
				SearchResult[] matches = r.getPssmMatches();
				for(SearchResult s: matches){
					outputStream.print(s.getIndex()+"\t"+r.getFullSeq().substring(s.getIndex(), s.getIndex()+s.getPssmLength())+"\t");
				}
			}

			//print out motifs in order
			int idx=0;
			for(PSSM p: pssmList){
				outputStream.println("\nProfileFragment "+idx+":\n"+p.getProfileFragAlignment());
				idx++;
			}

		} catch (Exception ex) {
			Logger.getLogger(Dasp.class.getName()).log(Level.SEVERE, null, ex);
			ex.printStackTrace();
		}

		// 	7. Combine the results from the searches
		// 	TODO: search.combineResults();

		// TODO: output final results
	}

	private static void usage() {
		System.out.println("Usage: dasp -i filename [-r n] [-c n.nn] [-p dir] [-o file] [-d db] [-h][-v]");
		System.out.println("arguments: ");
	 	System.out.println("    -i filename	The name of the input file");
	 	System.out.println("    -r radius	The radius for inclusion into the active site signature");
	 	System.out.println("    -c cutoff	The cutoff value for searching the sequence database");
	 	System.out.println("    -p pdbDatabase	The path to the pdb database");
	 	System.out.println("    -o filename	The file to wite the active stie profile into");
	 	System.out.println("    -d database	The database to use for the sequence search");
	 	System.out.println("    -h the help text");
	 	System.out.println("    -v print verbose output");
	}
}
