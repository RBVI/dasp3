/**
 * The DBSearch interface provides a general method to 
 * use position-specific scoring matrices to search sequence
 * databases.
 *
 */

package dasp.algorithms;

import java.io.File;
import java.io.IOException;
import java.util.List;

import dasp.model.PSSM;
import dasp.model.DBSearchResult;

public interface DBSearch {
	public List<DBSearchResult> search(File database, List<PSSM>pssmList, 
	                                   PSSMSearch searchAlg, double threashold) throws Exception;
}
