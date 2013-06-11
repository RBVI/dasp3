
package dasp.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SequenceFragment {

	List<Residue> fragResidues = null;
	Map<Integer, Residue> keyResidues = null;

	public SequenceFragment() {
		fragResidues = new ArrayList<Residue>();
		keyResidues = new HashMap<Integer, Residue>();
	}

	public void addFragRes(Residue r, int keyNumber)
	{
		fragResidues.add(r);
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
}
