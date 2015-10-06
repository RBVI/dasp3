/**
 * Residue objects contain information about the amino acids that make up
 * protein structures.  This information includes the amino acid type, and
 * the location (as an average of the locations of its member atoms).
 */

package dasp.model;

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

public class Residue {
	private static Map<String,Character>aaMap = null;
	private char aminoAcid = 'X';
  private String PDBcode = null; //residue number+insertion code - these can be
                           //non-sequential, so just use as a unique identifier.
                           //This needs looking into.  What I really need is just the residue ID as a number.
                           //in the code that uses this I pretend it is a number.
	private Group thisGroup = null;
	private int index = 0;
	Location center = null;
	
	public Residue(Group g) throws Exception {
		if (Residue.aaMap == null)
			Residue.initializeAAMap();

		thisGroup = g;

		// Get the location
		double xSum = 0.0;
		double ySum = 0.0;
		double zSum = 0.0;

		// Need to get the 3 character name and convert it ourselves to handle
		// non-standard amino acids
		// System.out.println("g.getPDBName = "+g.getPDBName());
		Character oneLetterCode = Residue.getAA(g.getPDBName());
		if (oneLetterCode == null) 
			throw new Exception(g.getPDBName()+" is not an amino acid");

		aminoAcid = Residue.getAA(g.getPDBName()).charValue();
		// System.out.println("Amino acid = "+aminoAcid);
		PDBcode = g.getPDBCode();

		// Just to make sure...
		aminoAcid = Character.toUpperCase(aminoAcid);

		for (Atom a: g.getAtoms()) {
			xSum += a.getX();
			ySum += a.getY();
			zSum += a.getZ();
		}
		center = new Location(xSum/g.size(), ySum/g.size(), zSum/g.size());
	}

	public void setIndex(int index) { this.index = index; }
	public int getIndex() { return index; }

	public char getAA() {
		return aminoAcid;
	}

	public Location getCenter()
	{
		return center;
	}

	/**
	 * Calculates the Euclidean distance between this residue center and
	 * a second residue center
	 * 
	 * @param r Second residue
	 * @return the Euclidean distance between the 2 residue centers as a double value.
	 */
	public double dist(Residue r){
		Location loc1 = getCenter();
		Location loc2 = r.getCenter();
		
		return Math.sqrt(Math.pow((loc2.getX() - loc1.getX()), 2) + 
				Math.pow((loc2.getY() - loc1.getY()), 2) +
				Math.pow((loc2.getZ() - loc1.getZ()), 2));
	} //end dist

	public String getPDBCode() {
		return this.PDBcode;
	}

	private static void initializeAAMap() {
		aaMap = new HashMap();
		aaMap.put("ALA",'A');
		aaMap.put("ARG",'R');
		aaMap.put("ASN",'N');
		aaMap.put("ASP",'D');
		aaMap.put("CYS",'C');
		aaMap.put("GLU",'E');
		aaMap.put("GLN",'Q');
		aaMap.put("GLY",'G');
		aaMap.put("HIS",'H');
		aaMap.put("ILE",'I');
		aaMap.put("LEU",'L');
		aaMap.put("LYS",'K');
		aaMap.put("MET",'M');
		aaMap.put("PHE",'F');
		aaMap.put("PRO",'P');
		aaMap.put("SER",'S');
		aaMap.put("THR",'T');
		aaMap.put("TRP",'W');
		aaMap.put("TYR",'Y');
		aaMap.put("VAL",'V');

		// Non-standard amino acids
		aaMap.put("ASH",'D');
		aaMap.put("CYX",'C');
		aaMap.put("GLH",'E');
		aaMap.put("HID",'H');
		aaMap.put("HIE",'H');
		aaMap.put("HIP",'H');
		aaMap.put("HYP",'P');
		aaMap.put("MSE",'M');
		aaMap.put("SNC",'C');
		aaMap.put("+N",'N');
	}

	private static Character getAA(String threeLetter) {
		if (aaMap.containsKey(threeLetter))
			return aaMap.get(threeLetter);
		return null;
	}
}
