
package dasp.model;

import java.util.List;
import java.util.ArrayList;
import dasp.algorithms.Align;

public class ProfileFragment extends Alignment {
	List<Residue> fragResidues = null;

	public ProfileFragment (Align alnMethod) {
		fragResidues = new ArrayList();
		method = alnMethod;
	}

	public void addFragRes(Residue r)
	{
		fragResidues.add(r);
	}
	
	public Residue getFragRes(int resIdx)
	{
		return fragResidues.get(resIdx);
	}

	public String makeUpper(String seq) {
		return null;
	}

	public String getSequence(boolean toUpper) {
		String seq = "";
		for (Residue r: fragResidues) {
			seq+=r.getAA();
		}
		if (toUpper)
			return seq.toUpperCase();

		return seq.toLowerCase();
	}
}
