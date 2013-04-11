
package dasp.model;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class DBSearchResult implements Comparable{

    private double pvalueScore = -1; //final pvalue score
    private String pseudoSig = ""; //constructed pseudoSignature based on sorted array of matches
    //ISSUE: to get the pseudo signature I have to order this array of matches based on
    //pssm length.  However, to stay in sync with the PSSM list in main it needs to be unordered.
    private SearchResult[] pssmMatches = null; //this needs to remain unordered to stay in sync with the PSSM list
    private String seqName = "";
    private String fullSeq = "";
    private String formattedFullSeq = ""; //the full sequence in all lowercase except for the pssm matches.
    
    //the unsorted array SearchResult objects. It is left unsorted so that it remains in the same order as the PSSM list
    //should we have the PSSM list saved in here?  I am thinking not as we can use the one in the main class
    
	public DBSearchResult(double pvalueScore, SearchResult[] pssmMatches, String seqName, String fullSeq) {
        this.pvalueScore = pvalueScore;
        this.pssmMatches = pssmMatches;
        this.seqName = seqName;
        this.fullSeq = fullSeq;
        
        generatePseudoSig();
        formatFullSeq();
	}

    /**
     * Gets the array of PSSM matches for this result.
     * 
     * @return The SearchResult array of pssm matches.
     */
    public SearchResult[] getPssmMatches() { return pssmMatches; }

    /**
     * Gets the pvalue score for this result.
     *
     * @return Returns the pvalue of this result as a double.
     */
    public double getPval(){ return pvalueScore; }

    /**
     * Gets the name of the sequence for this result.
     *
     * @return Returns the name of the sequence as a String.
     */
    public String getName(){ return seqName; }

    /**
     * Gets the full sequence for this search result.
     *
     * @return Returns the full sequence as a String.
     */
    public String getFullSeq(){ return fullSeq; }

    /**
     * Gets the pseudo signature which is fragment formatted.
     * 
     * @return Returns the pseudo signature as a String.
     */
    public String getPseudoSig(){ return pseudoSig; }

    /**
     * Gets the full sequence formatted so that the fragments of the pseudo signature
     * are in context.
     *
     * @return Returns the full sequence formatted with pseudo sig fragments in context.
     */
    public String getSigInContext(){ return formattedFullSeq; }

    /**
     * Used to sort a list of DBSearchResults by increasing pvalue
     *
     * @param result Another DBSearchResult object
     * @return Returns 1, 0 or -1 if this pvalue is greater, equal or less than
     * the pvalue of the passed result, respectively.
     */
    public int compareTo(Object result) {
        DBSearchResult cmp = (DBSearchResult)result;
        double thisPval = this.pvalueScore;

        double cmpPval = cmp.getPval();

        if(thisPval > cmpPval)
            return 1;
        else if(thisPval < cmpPval)
            return -1;

        return 0;
    }

    /**
     * Puts the pseudo signature fragments in context of the full sequence
     */
    private void formatFullSeq() {
        //convert fullSeq to all lowercase
        formattedFullSeq = fullSeq.toLowerCase();
        //loop through each match and convert to uppercase

    }

    /**
     * Generates the fragment formatted pseudo signature for this significant
     * search result.
     */
    private void generatePseudoSig() {

        SearchResult[] tempResults = new SearchResult[pssmMatches.length];
        for(int i=0;i<pssmMatches.length;i++){
            tempResults[i] = pssmMatches[i];
        }

        List<SearchResult> list = Arrays.asList(tempResults);
        Collections.sort(list);

		boolean toUpper = true;
		for (SearchResult s: list) {
			if (toUpper){
                pseudoSig += fullSeq.substring(s.getIndex(), s.getIndex()+s.getPssmLength()).toUpperCase();
                toUpper = false;
            } 
			else{
                pseudoSig += fullSeq.substring(s.getIndex(), s.getIndex()+s.getPssmLength()).toLowerCase();
                toUpper = true;
            }
		}
    }


}
