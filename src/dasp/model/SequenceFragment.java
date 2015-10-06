
package dasp.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class SequenceFragment {

	List<Residue> fragResidues = null;
	Set<Residue> keySet = null;
	Map<Integer, Residue> keyResidues = null;

	public SequenceFragment() {
		fragResidues = new ArrayList<Residue>();
		keySet = new HashSet<Residue>();
		keyResidues = new HashMap<Integer, Residue>();
	}

	public void addFragRes(Residue r, Residue key, int keyNumber)
	{
		fragResidues.add(r);
		keySet.add(key);
		if (keyNumber >=0 ) keyResidues.put(keyNumber, r);
	}
	
	public Residue getFragRes(int resIdx)
	{
		return fragResidues.get(resIdx);
	}

	public int getResidueCount() {
		return fragResidues.size();
	}

	public String makeUpper(String seq) {
		return null;
	}

	public String getSequence(boolean toUpper) {
		if (fragResidues.size() == 0)
			return "-";
		String seq = "";
		for (Residue r: fragResidues) {
			seq+=r.getAA();
		}
		if (toUpper)
			return seq.toUpperCase();

		return seq.toLowerCase();
	}

	public int getClosestKey() {
		// Calculate the centroid of all of our residues
		double xCent = 0.0;
		double yCent = 0.0;
		double zCent = 0.0;
		for (Residue r: fragResidues) {
			xCent += r.getCenter().getX();
			yCent += r.getCenter().getY();
			zCent += r.getCenter().getZ();
		}
		xCent /= fragResidues.size();
		yCent /= fragResidues.size();
		zCent /= fragResidues.size();

		// Now, for each key residue, get the distance between
		// that residue and our centroid.  Return the shortest
		double closest = Double.MAX_VALUE;
		Integer closestKey = null;
		for (Integer keyIndex: keyResidues.keySet()) {
			Residue k = keyResidues.get(keyIndex);
			Location loc2 = k.getCenter();
			double dist = Math.sqrt(Math.pow((loc2.getX() - xCent), 2) + 
				Math.pow((loc2.getY() - yCent), 2) +
				Math.pow((loc2.getZ() - zCent), 2));
			if (dist < closest) {
				closest = dist;
				closestKey = keyIndex;
			}
		}
		return closestKey.intValue();
	}
}
