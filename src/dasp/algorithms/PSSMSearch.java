/**
 * The PSSMSearch interface provides a general method to 
 * use position-specific scoring matrices to search sequence
 * databases.
 *
 */

package dasp.algorithms;

import java.util.List;

import dasp.model.PSSM;
import dasp.model.SearchResult;

interface PSSMSearch {
	SearchResult[] search(List<PSSM>pssmList, String sequence, int xCount, boolean debug);
}
