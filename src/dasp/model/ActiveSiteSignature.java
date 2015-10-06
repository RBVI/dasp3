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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * An ActiveSiteSignature is a sequence that is formed by concatenating
 * all of the fragments in a structure that are within a defined distance
 * of one or more key residues.
 */
public class ActiveSiteSignature {
	private DaspStructure struct = null;
	private List<Residue> keyResidues = null;
	private List<SequenceFragment> fragments = null;
	private List<SequenceFragment> residueFrags = null;
	private List<SequenceFragment> otherFrags = null;
	private Map<Integer, ActiveSiteSignature> fragSigCache = null;
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

	public ActiveSiteSignature (DaspStructure struct, List<Residue> residues, SequenceFragment frag) {
		this.keyResidues = residues;
		this.struct = struct;
		fragments = new ArrayList<SequenceFragment>();
		fragments.add(frag);
		signature = frag.getSequence(false);
	}

	public String getSignature(double profileRadius) {
		if (signature != null && keyResidues == null)
			return signature;

		if (signature == null || radius != profileRadius) {
			// calculate the ASSig
			//Should this be performed in the constructor?
			this.radius = profileRadius;
			this.fragments = calculateASSig(profileRadius);

			this.signature = alternate(fragments);
		}

		// System.out.println("signature = "+signature);

		return signature;
	}

	public String getPdbId(){
		return struct.getPdbId();
	}

	public List<SequenceFragment> getFragments() {
		return fragments;
	}

	public String getFragmentsInContext() {
		String pdbSequence = struct.getSequence();
		for (SequenceFragment frag: fragments) {
			frag.makeUpper(pdbSequence);
		}
		return pdbSequence;
	}

	public ActiveSiteSignature getFragmentAsSig(int frag) {
		if (fragSigCache == null)
			fragSigCache = new HashMap<Integer, ActiveSiteSignature>();

		if (fragSigCache.containsKey(frag))
			return fragSigCache.get(frag);

		SequenceFragment newFrag = null;
		if (frag >= fragments.size()) {
			newFrag = new SequenceFragment();
			fragments.add(newFrag);
			otherFrags.add(newFrag);
		} else if (frag >= residueFrags.size()) {
			newFrag = otherFrags.get(frag-residueFrags.size());
		} else {
			newFrag = residueFrags.get(frag);
		} 

		ActiveSiteSignature aSig = new ActiveSiteSignature(struct, keyResidues, newFrag);
		fragSigCache.put(frag, aSig);
		return aSig;
	}

	public int keyFragCount() {
		return residueFrags.size();
	}

	public int fragCount() {
		return residueFrags.size()+otherFrags.size();
	}

	public List<Residue> getKeyResidues() {
		return keyResidues;
	}

	public void removeFragment(int fragment) {
		if (fragment >= fragments.size())
			return;
		if (fragment < residueFrags.size()) {
			residueFrags.remove(fragment);
		} else if (fragment >= residueFrags.size()) {
			otherFrags.remove(fragment-residueFrags.size());
		}
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
        
		String[] tokens = input.trim().split(":");

		if (tokens.length < 2)
			throw new IOException("No ':' in input string: "+input);

		int chainStart = tokens[0].indexOf('#');
		if (chainStart > 0) {
			chain = tokens[0].substring(chainStart+1);
			pdbFile = tokens[0].substring(0,chainStart);
		} else {
			pdbFile = tokens[0];
		}

		if (tokens.length == 2) {
			String[] residues = tokens[1].split(",");
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
		} else if (tokens.length == 3) {
			// We're getting the signature -- just read it and
			// skip all the work calculating it
			signature = tokens[2].trim();
			struct = new DaspStructure(database, pdbFile, chain);
			// TODO: Create SequenceFragments!!!!!
			// Read through the signature a character at a time
			// to create the fragments
		}
	}

	/**
   * Extracts the fragments that are part of the active site signature
   * 
   * @param radius  The profile radius around the key residues
   * @return  Returns the list of fragments that are within the defined radius of at least 1 key residue.
   */
	private List<SequenceFragment> calculateASSig(double radius) {
		List<SequenceFragment> returnFrags = new ArrayList<SequenceFragment>();
		residueFrags = new ArrayList<SequenceFragment>();
		otherFrags = new ArrayList<SequenceFragment>();

		//not sure this is the best way to start this off.
		SequenceFragment tempFrag = new SequenceFragment();  //starts the first fragment object
		int prevResId = 0;  //initializes the previous residue ID
		boolean haveKey = false;
		for (Residue curRes: struct.getResidues())
		{ 
			for (Residue keyRes: keyResidues)
			{
				// System.out.println("Comparing: "+curRes.getAA()+" to "+keyRes.getAA());
				//get centers
				//compare centers
				double EDdist = keyRes.dist(curRes);
				
				//test for within radius
				if(EDdist <= radius && EDdist >= -radius)
				{
						if(prevResId > 0 && curRes.getIndex() != prevResId+1)
						{
							if (tempFrag.getResidueCount() >= 4) {
								returnFrags.add(tempFrag);
								if (haveKey)
									residueFrags.add(tempFrag);
								else
									otherFrags.add(tempFrag);
							}
							tempFrag = new SequenceFragment();
							haveKey = false;
						}

						int keyIndex = keyResidues.indexOf(curRes);
						tempFrag.addFragRes(curRes, keyRes, keyIndex); //if true, add residue to fragments list and break
						if (keyIndex >= 0) haveKey = true;

						prevResId = curRes.getIndex();
						break;
				} 
				//if false move onto next keyRes
			}

		}

		if (tempFrag.getResidueCount() >= 4) {
			if (haveKey)
				residueFrags.add(tempFrag);
			else
				otherFrags.add(tempFrag);

			returnFrags.add(tempFrag); // Add the last fragment
		}

		//By this point we should have all our fragments that are within the radius
		//of the active site keyResidues, they should be in sequencial order,
		//and each of the consecutive fragments is stored in a different fragment object.

		return returnFrags;
	}
	
	/**
 	 * Extend an alignment string.  NOTE: this method only operates on the
 	 * first fragment.  It assumes that this is being used in fragment-oriented
 	 * mode.
 	 *
 	 * @param sequence the alignment string (with dashes)
 	 * @return extended sequence
 	 */
	public String extendAlignment(String alignmentString) {
		if (fragments.size() > 1) return null;
		// System.out.println("Extending "+alignmentString+" from "+struct.getPdbId());
		SequenceFragment frag = fragments.get(0);
		String newAlignment = "";

		int leadingDashes = countLeadingDashes(alignmentString);
		if (leadingDashes == alignmentString.length())
			return alignmentString;

		int trailingDashes = countTrailingDashes(alignmentString);

		Residue leadingResidue = frag.getFragRes(0);
		Residue trailingResidue = frag.getFragRes(frag.getResidueCount()-1);

		int leadingIndex = leadingResidue.getIndex();
		int trailingIndex = trailingResidue.getIndex();
		// System.out.println("leadingIndex = "+leadingIndex+" leading dashes = "+leadingDashes);
		// System.out.println("trailingIndex = "+trailingIndex+" trailing dashes = "+trailingDashes);
		for (int r = leadingIndex-leadingDashes; r<leadingIndex; r++) {
			if (r >=0 && struct.getResidues().get(r) != null) {
				newAlignment += struct.getResidues().get(r).getAA();
			} else {
				newAlignment += "-";
			}
		}

		newAlignment += alignmentString.substring(leadingDashes, alignmentString.length()-trailingDashes);

		if (trailingDashes > 0) {
			int endIndex = trailingIndex+trailingDashes+1;
			for (int r = trailingIndex+1; r < endIndex; r++) {
				// System.out.println("Looking for residue "+r+": "+struct.getResidues().get(r).getAA());
				// System.out.println("   Residue count = "+frag.getResidueCount());
				if (r < struct.getResidues().size() && struct.getResidues().get(r) != null) {
					// System.out.println("Appending....");
					newAlignment += struct.getResidues().get(r).getAA();
				} else
					newAlignment += "-";
			}
		}
		// System.out.println("New alignment = "+newAlignment);
		return newAlignment;
	}

	private int countLeadingDashes(String alignmentString) {
		int count = 0;
		for (int i = 0; i < alignmentString.length(); i++) {
			if (alignmentString.charAt(i) != '-')
				break;
			count++;
		}
		return count;
	}

	private int countTrailingDashes(String alignmentString) {
		int count = 0;
		for (int i = alignmentString.length()-1; i >= 0; i--) {
			if (alignmentString.charAt(i) != '-')
				break;
			count++;
		}
		return count;
	}
}
