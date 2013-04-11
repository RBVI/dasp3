
package dasp.model;

import java.util.List;
import java.util.ArrayList;

public class Fragment {

	List<Residue> fragResidues = null;

	public Fragment() {
		fragResidues = new ArrayList();
	}

	public void addFragRes(Residue r)
	{
		fragResidues.add(r);
	}
	
	public Residue getFragRes(int resIdx)
	{
		return fragResidues.get(resIdx);
	}

    //need to implement this!!!!
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
