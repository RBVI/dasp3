
package dasp.model;

public class SearchResult implements Comparable{
	int index;
    int pssmLength;
	double pvalue;

	public SearchResult (int index, double pvalue, int pssmLength) {
		this.index = index;
		this.pvalue = pvalue;
        this.pssmLength = pssmLength;
	}

	public void setIndex(int index) { this.index = index; }
	public void setPvalue(double pvalue) { this.pvalue = pvalue; }
	public int getIndex() { return index; }
	public double getPvalue() { return pvalue; }
	public void setPssmLength(int length) { this.pssmLength = length; }
	public int getPssmLength() { return pssmLength; };

	public int compareTo(Object result) {
		SearchResult cmp = (SearchResult)result;
		int thisIndex = this.getIndex();
		int cmpIndex = cmp.getIndex();

		if(thisIndex > cmpIndex)
			return 1;
		else if(thisIndex < cmpIndex)
			return -1;
		else if(thisIndex==-1 && cmpIndex==-1){
			System.out.println("ERROR: cannot sort pseudoSig fragments because of non-existent index.");
		}
		else if(thisIndex==cmpIndex){
			System.out.println("ERROR: cannot sort pseudoSig fragments have overlapping start index.");
		}

		return 0;
	}

}
