/*
 File: ActiveSiteSignature.java

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
import java.util.List;

/**
 * An ActiveSiteSignature is a sequence that is formed by concatenating
 * all of the fragments in a structure that are within a defined distance
 * of one or more key residues.
 */
public class ActiveSiteSignature {
	private DaspStructure struct = null;
	private List<Residue> keyResidues = null;
	private List<SequenceFragment> fragments = null;
	private String signature = null;
	private double radius = 10.0;

	public ActiveSiteSignature (String inputString, File database) throws IOException {
		// Parse input string
		parseInput(inputString, database);
	}

	public ActiveSiteSignature (DaspStructure struct, List<Residue> residues) {
		this.keyResidues = residues;
		this.struct = struct;
	}

	public String getSignature(double profileRadius) {
		if (signature == null || radius != profileRadius) {
			// calculate the ASSig
			//Should this be performed in the constructor?
			this.radius = profileRadius;
			this.fragments = calculateASSig(profileRadius);

			this.signature = alternate(fragments);
		}

		return signature;
	}

	public String getPdbId(){
		return struct.getPdbId();
	}

	public String getFragmentsInContext() {
		String pdbSequence = struct.getSequence();
		for (SequenceFragment frag: fragments) {
			frag.makeUpper(pdbSequence);
		}
		return pdbSequence;
	}

	private String alternate(List<SequenceFragment>frags) {
		boolean toUpper = true;
		String result = "";
		for (SequenceFragment frag: frags) {
			result += frag.getSequence(toUpper);
			if (toUpper) 
				toUpper = false;
			else
				toUpper = true;
		}
		return result;
	}

	/**
	 * parseInput will read the input line and create an appropriate
	 * DaspStructure object with the key residues.
	 *
	 * Input format:
	 *	 structure#chain: residueList
	 *  where
	 *		structure is a structural identifier (a pdbID or a file path)
	 *		residueList is a comma-separated list of residues.  A residue is either
	 *			a simple integer or a residue name followed by an integer.
	 *
	 * @param input the input string
	 * @throws IOException if we get a read error
	 */
	private void parseInput(String input, File database) throws IOException {
		String pdbFile = null;
		String chain = null;
        
		int split1 = input.indexOf(':');

		if (split1 < 0)
			throw new IOException("No ':' in input string: "+input);

		int chainStart = input.substring(0,split1).indexOf('#');
		if (chainStart > 0) {
			chain = input.substring(chainStart+1, split1);
			pdbFile = input.substring(0,chainStart);
		} else {
			pdbFile = input.substring(0, split1);
		}

		String[] residues = input.substring(split1+1).split(",");
		if (residues.length == 0)
			throw new IOException("No residues in input string: "+input);

		struct = new DaspStructure(database, pdbFile, chain);

		keyResidues = new ArrayList();
		for (int residue = 0; residue < residues.length; residue++) {
			Residue r = struct.getResidue(residues[residue]);
			if (r != null) {
				keyResidues.add(r);
			} else {
				throw new IOException("PDB file "+pdbFile+" doesn't have residue "+residues[residue]);
			}
		}
	}

	/**
   * Extracts the fragments that are part of the active site signature
   * 
   * @param radius  The profile radius around the key residues
   * @return  Returns the list of fragments that are within the defined radius of at least 1 key residue.
   */
	private List<SequenceFragment> calculateASSig(double radius) {
			List<SequenceFragment> returnFrags = new ArrayList();
	
			//not sure this is the best way to start this off.
			SequenceFragment tempFrag = new SequenceFragment();  //starts the first fragment object
			int prevResId = 0;  //initializes the previous residue ID
			for (Residue curRes: struct.getResidues())
			{ 
				for(Residue keyRes: keyResidues)
				{
					// System.out.println("Comparing: "+curRes.getAA()+" to "+keyRes.getAA());
					//get centers
					//compare centers
					double EDdist = dist(keyRes, curRes);
					
					//test for within radius
					if(EDdist <= radius && EDdist >= -radius)
					{
							if(prevResId > 0 && curRes.getIndex() != prevResId+1)
							{
								returnFrags.add(tempFrag);
								tempFrag = new SequenceFragment();
							}

							tempFrag.addFragRes(curRes); //if true, add residue to fragments list and break
							prevResId = curRes.getIndex();
							break;
					} 
					//if false move onto next keyRes
				}

			}
			returnFrags.add(tempFrag); // Add the last fragment

			//By this point we should have all our fragments that are within the radius
			//of the active site keyResidues, they should be in sequencial order,
			//and each of the consecutive fragments is stored in a different fragment object.
	
			return returnFrags;
		}
	
	/**
	 * Calculates the Euclidean distance between two residue centers
	 * 
	 * @param r1 First residue
	 * @param r2 Second residue
	 * @return the Euclidean distance between the 2 residue centers as a double value.
	 */
	public double dist(Residue r1, Residue r2){
		Location loc1 = r1.getCenter();
		Location loc2 = r2.getCenter();
		
		return Math.sqrt(Math.pow((loc2.getX() - loc1.getX()), 2) + 
				Math.pow((loc2.getY() - loc1.getY()), 2) +
				Math.pow((loc2.getZ() - loc1.getZ()), 2));
	} //end dist

}
