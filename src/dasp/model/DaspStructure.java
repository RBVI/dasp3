/*
 File: DaspStructure.java

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

import java.io.File;
import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

/**
 * DaspStructure is the main class that contains information about a PDB
 * structure.
 */
public class DaspStructure {
	String databasePrefix = null; // A path to the PDB file hierarchy
	List<Residue> residues = null; // the residues for this structure
	Map<String,Residue> residueMap = null; // the residues for this structure
	String chainID = null; // the chain ID
	String pdbFile = null;
    int nextResIdx = -1; //the enumerator for the nextResidue Method

	/**
	 * Create a new structure from a pdb file.  If there is a file separator
	 * in the PDB file, we assume that the file is a full path to an independent
	 * PDB file.  If there is no slash, we will use the standard PDB director
	 * structure to find the identifier.
	 *
	 * @param pdbDatabase the path to the PDB database
	 * @param pdbFile the pdb identifier or path to the pdb file
	 */
	public DaspStructure(File pdbDatabase, String pdbFile) throws IOException {
		this.chainID = null;
		this.pdbFile = pdbFile;
		getPDB(pdbFile, pdbDatabase);
	}

	/**
	 * Create a new structure from a pdb file.  If there is a file separator
	 * in the PDB file, we assume that the file is a full path to an independent
	 * PDB file.  If there is no slash, we will use the standard PDB director
	 * structure to find the identifier.  This form of the constructor takes
	 * a chain ID to pull a specific chain from the file
	 *
	 * @param pdbDatabase the path to the PDB database
	 * @param pdbFile the pdb identifier or path to the pdb file
	 * @param chain the chain ID
	 */
	public DaspStructure(File pdbDatabase, String pdbFile, String chain) throws IOException {
		this.chainID = chain;
		this.pdbFile = pdbFile;
		// System.out.println("Getting pdb file: "+pdbFile+"#"+chain+" from database: "+pdbDatabase);
		getPDB(pdbFile, pdbDatabase);
	}

	/**
	 * Return the amino acid sequence for this PDB file.
	 *
	 * @return the sequence
	 */
	public String getSequence() {
		String sequence = "";
		for (Residue r: residues) {
			sequence += r.getAA();
		}
		return sequence;
	}
    
    /**
     * Return the pbd id and chain as a string for this sequence
     * Does not account ofr a possible 'pdb' in front of the pdb id.
     * 
     * @return a String containing the PDB ID and chain
     */
    public String getPdbId() {
        String pdbID = "";
        String chain;

        if(chainID != null)
            chain = chainID;
        else
            chain = "";

        if(pdbFile == null)
        {
            //output some sort of error message?
        }
        if (pdbFile.indexOf("/") > 0) {
            if(pdbFile.indexOf(".") > 0){
                pdbID = pdbFile.substring(pdbFile.lastIndexOf("/"),pdbFile.indexOf("."))+chain;  
            }
            else
                pdbID = pdbFile.substring(pdbFile.lastIndexOf("/"))+chain;
		}
        else
            pdbID = pdbFile+chain;
        
        
        return pdbID;
    }

    /**
     * Returns the number of residues in the sequence.
     *
     * @return the length of the sequence
     */
    public int getSeqLength()
    {
        return residues.size();
    }

		/**
		 * Return the list of residues
		 *
		 * @return the list of residues
		 */
		public List<Residue> getResidues() { return residues; }
    
    /**
     * Enumerates through sequence of residues and returns the next Residue object.
     *
     * Is there a better way to implement this?  Using something already avaliable?
     * I didn't want the other method to have direct access to the residue seqeunce list.
     *
     * @return the next residue object in the sequence.
     */
    public Residue nextResidue()
    {
        nextResIdx++;
        return residues.get(nextResIdx);
    }

    /**
     * Tests to see if there is at least one more Residue object in the sequence
     *
     * @return true if more Residue objects exist, false otherwise.
     */
    public boolean hasMoreResidues()
    {
        if(nextResIdx >= residues.size())
            return false;
        else
            return true;
    }

    /**
     * Resets the enumerator variable to -1
     */
    public void resetEnum()
    {
        nextResIdx = -1;
    }

	/**
	 * Return a Residue object for a residue within this PDB file.
	 *
	 * @param residueID the ID for the residue to get from the structure.  The
	 * residue ID can be the residue # as an int, or an amino acid code followed
	 * by an integer value (e.g. A123)
	 * @return a residue object corresponding to that residue ID
	 */
	public Residue getResidue(String residueID) {
		// System.out.println("Searching for residue: '"+residueID+"'");
		String residueStr = null;
		char aaCode = 'X';
		int residueIndex = 0;
		residueID = residueID.trim();
		if (residueID.matches("[a-zA-Z][0-9]+[a-zA-Z]*")) {
			// System.out.println("Found aa code: "+residueID.charAt(0));
			residueStr = residueID.substring(1);
			aaCode = residueID.charAt(0);
		// } else if (residueID.matches("[0-9]+[a-zA-Z]*")) {
		// 	System.out.println("Found trailing characters: "+residueID.charAt(0));
		// 	residueStr = residueID.substring(1);
		} else {
			residueStr = residueID;
		}

		// System.out.println("Residue string: "+residueStr+" aa code = "+aaCode);
		// Get the residue from our list
		Residue res = residueMap.get(residueStr);
		if (aaCode == 'X' || res.getAA() == Character.toUpperCase(aaCode))
			return res;

		// AA code didn't match
		return null;
	}

	public Residue getResidue(int index) {
		if (residueMap.containsKey(String.valueOf(index)))
			return residueMap.get(String.valueOf(index));
		return null;
	}

	/**
	 * This method reads and parses the pdb file and creates the internal
	 * structures to represent the protein structure.  Note that as a by-product
	 * of this method, the individual Residue objects are created and the 
	 * corresponding geometric centers are calculated.
	 *
	 * @param fileName the pdb file name
	 * @param pdbDatabase the pdb database directory
	 * @throws IOException on either read or parsing errors
	 */
	private void getPDB(String fileName, File pdbDatabase) throws IOException {
		File pdbFile = new File(fileName);
		String filePath = null;
		if (!pdbFile.isAbsolute()) {
			// Need to figure out the prefix and find the file
			filePath = findPDBPath(fileName, pdbDatabase);
		} else {
			filePath = pdbFile.getPath();;
		}
		// System.out.println("getPDB: "+filePath);

		PDBFileReader pdbreader = new PDBFileReader();
		try {
			// System.out.println("Getting structure: "+filePath);
			Structure struct = pdbreader.getStructure(filePath);
			Chain chain = null;

			if (chainID != null && struct.hasChain(chainID)) {
				// System.out.println("Getting chain: "+chainID);
				chain = struct.getChainByPDB(chainID);
			} else {
				chain = struct.getChain(0);
			}

			residues = new ArrayList();
			residueMap = new HashMap();
			for (Group g: chain.getAtomGroups()) {
				// System.out.println("Got group "+g);
				try {
					Residue r = new Residue(g);
					residues.add(r);
					r.setIndex(residues.size()-1);
					// System.out.println("Read residue: "+r+" position: "+(residues.size()-1));
					residueMap.put(g.getPDBCode(),r);
				} catch (Exception e) {
					// Not an amino acid
				}
			}
		} catch (Exception e) {
			throw new IOException("Can't read PDB file "+filePath+": "+e.getMessage());
		}
	}

	private String findPDBPath(String fileName, File pdbDatabase) {
		// A PDB file name could be "pdbnmmp.ent" or "pdbnmmp" or "nmmp" 
		// where n is a number, mm is a two letter index, and p is either a
		// letter or a number.
		String filePath = pdbDatabase.getPath();
		// Make the filename lower case (convention)
		fileName = fileName.toLowerCase();

		int indexOffset = 1;
		if (fileName.startsWith("pdb"))
			indexOffset = 4;

		filePath = filePath + "/" + fileName.substring(indexOffset, indexOffset+2) + "/";

        if (indexOffset == 1)
			filePath = filePath + "pdb";

		filePath = filePath + fileName;

		if (!fileName.endsWith(".ent"))
			filePath = filePath + ".ent";

		return filePath;
	}
}
