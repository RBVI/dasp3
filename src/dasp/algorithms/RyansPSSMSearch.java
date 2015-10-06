
package dasp.algorithms;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;

import dasp.model.PSSM;
import dasp.model.SearchResult;
import java.util.ArrayList;
import java.util.Collections;

public class RyansPSSMSearch implements PSSMSearch {
	private boolean includeX = false;

	public RyansPSSMSearch (boolean includeX) {
		this.includeX = includeX;
	}

	//match the PSSMs to a database sequence
	//Note PSSM list should already be in sequential order from largest PSSM to smallest
	public SearchResult[] search(List<PSSM> list, String seq, int xCount, boolean debug){
		SearchResult[] pvalues = new SearchResult[list.size()];

		char[] seqArray = seq.toCharArray();
		double aaIncrement = 1.0/(double)(seqArray.length-xCount);

		double[] seqFreqs;
		if (includeX)
			seqFreqs = new double[21];
		else
			seqFreqs = new double[20];

		//count the number of occurences of each residue in the sequence
		for(int position = 0; position < seq.length(); position++){
			int aa = list.get(0).getAANum(seqArray[position]);
			if (aa == -1) continue;
			seqFreqs[aa] += aaIncrement;
		}

		int pssmNumber = 0;
		//match each PSSM to the sequence...
		for (PSSM pssm: list) {
			double[][] matrix = pssm.getPSSM();

			double max = -1*Double.MAX_VALUE;
			int max_index = -1;
			//find the position where the current PSSM matches the sequence
			for(int position = 0; position < seqArray.length - matrix[0].length + 1; position++) {
				double score = -1*Double.MAX_VALUE;
				boolean badRes = false;
				for(int pssmColumn = 0; pssmColumn < matrix[0].length; pssmColumn++){
					if (includeX || !(seqArray[position+pssmColumn] == 'X')) {
						int index = pssm.getAANum(seqArray[position+pssmColumn]); 
						if(pssmColumn == 0)
							score = matrix[index][pssmColumn];
						else
							score += matrix[index][pssmColumn];
					} else {
						badRes = true;
						if (debug)
							System.out.println("Found bad residue at: "+(position+pssmColumn));
					}
				}

				if((score > max) && notTaken(pvalues, position, pssm, list) && !badRes){
					if (debug)
						System.out.println("Found match at "+position+" score = "+score);
					max_index = position;
					max = score;
				}
			}

			//if we matched the PSSM to a valid position...
			if(max_index != -1){
				//calculate the p-value of the match
				double pvalue = calcPvalue(matrix, max, seqFreqs);
				if (debug)
					System.out.println("Pvalue for match at "+max_index+" is = "+pvalue);
				//normalize the p-value, but we need very high precision to do it
				BigDecimal bigPvalue = new BigDecimal(pvalue, MathContext.DECIMAL128);
				BigDecimal minusOnePvalue = BigDecimal.ONE.subtract(bigPvalue, MathContext.DECIMAL128);
				BigDecimal bigPpvalue = 
					BigDecimal.ONE.subtract(
							minusOnePvalue.pow(seqArray.length 
					                       - xCount 
					                       - matrix[0].length + 1, 
					                       MathContext.DECIMAL128), 
              MathContext.DECIMAL128);
				double ppvalue = bigPpvalue.doubleValue();
				if (debug) {
					System.out.println("Normalized pvalue for match at "+max_index+" is = "+ppvalue);
					if (ppvalue == 0.0) {
						System.out.println("Normalization power: "+(seqArray.length - xCount - matrix[0].length + 1));
						System.out.println("Normalization result: "+Math.pow(1.0-pvalue, (double)(seqArray.length - xCount - matrix[0].length + 1)));
					}
				}

				pvalues[pssmNumber] = new SearchResult(max_index, ppvalue, matrix[0].length);
			}
			else{
				pvalues[pssmNumber] = new SearchResult(max_index, -1*Double.MAX_VALUE, matrix[0].length);
			}
			pssmNumber++;
		}
		boolean valid = true;
		for(int i = 0; i < pvalues.length; i++){
			int key = pvalues[i].getIndex();
			if(key == -1)
				valid = false;
		}
		if(valid)
			return pvalues;
		else
			return null;
	}


	/**
	 * Returns whether or not a PSSM has already matched at a position.
	 * I don't think this algorithm is correct, but I am leaving it as is so it
	 * works the same as the original DASP.  We need to come back to this algortihm
	 * and make sure it is doing what we want it to.
	 *
	 * @param pvalues   The array of SearchResults already found
	 * @param start     The position we are at in the sequence iteration
	 * @param curPSSM   The current PSSM object we are trying to match
	 * @param pssmList  The full list of PSSM objects
	 * @return          True if this position has not been taken, false otherwise
	 */
	private boolean notTaken(SearchResult[] pvalues, int start, PSSM curPSSM, List<PSSM> pssmList) {
		boolean flag = true;
		int pssmNumber = 0;
		for (PSSM matrix: pssmList) {
			if (matrix == curPSSM)
                break;
            
			int s = pvalues[pssmNumber].getIndex();
			if(s != -1){
				if(start <= s){
					if((start + curPSSM.getWidth()) >= s){
						flag = false;
					}
				}
				else{
					if((s + matrix.getWidth()) >= start){
						flag = false;
					}
				}
			}
            pssmNumber++;
		}
		return flag;
	}

	//get the p-value for the score of 'max'
	private double calcPvalue(double[][] matrix, double max, double[] seqFreqs){
		int numCols = matrix[0].length;
		LinkedList<Map<Double,Double>> PGFs = new LinkedList();

		//create the probability generating functions for each column of the PSSM
		//loop through each column of the pssm
		for(int pssmColumn = 0; pssmColumn < numCols; pssmColumn++){
			Map<Double,Double> colPGF = Collections.synchronizedMap(new HashMap());
			//for each amino acid (pssm rows), so going down a column and
			//looking at all amino acid frequencies for each pssm position.
			for(int pssmRow = 0; pssmRow < matrix.length; pssmRow++){
				//get the pssm score for this amino acid
				double score = matrix[pssmRow][pssmColumn];
				//turn it into a Double object
				Double sObj = new Double(score);
				//if the PGF for this column is empty or it does not contain a
				//key for the pssm score then create an entry for this pssm score
				if(colPGF.isEmpty() || !colPGF.containsKey(sObj)){
					//get this amino acids frequency in the sequence
					Double freq = new Double(seqFreqs[pssmRow]);
					//colPGF contains the pssm score for this amino acid in thie position
					//and the overall frequency of this amino acid in the sequence
					colPGF.put(sObj, freq);
				}
				else{
					double fncn = ((Double)colPGF.get(sObj)).doubleValue();
					//basically, if this pssm score already exists then it adds the previous amino acid seq frequency to this frequency
					colPGF.put(sObj, new Double(fncn + seqFreqs[pssmRow]));
				}
			}
			PGFs.add(colPGF);
		}
        
		//create the final probabilty generating function for the entire PSSM
		while(PGFs.size() > 1){
			Map<Double,Double> colA = PGFs.get(0);
			Map<Double,Double> colB = PGFs.get(1);
			Map<Double,Double> C = Collections.synchronizedMap(new HashMap());

			for (Double scoreAObj: colA.keySet()) {
				double scoreA = scoreAObj.doubleValue();
				double fncnA = colA.get(scoreAObj).doubleValue();
				for (Double scoreBObj: colB.keySet()) {
					double scoreB = scoreBObj.doubleValue();
					double sum = scoreA + scoreB;
					Double sObj = new Double(sum);
					if(C.isEmpty() || !C.containsKey(sObj))
						C.put(sObj, new Double(
							fncnA*colB.get(scoreBObj).doubleValue()));
					else{
						double fncn = C.get(sObj).doubleValue();
						fncn += (fncnA * colB.get(scoreBObj).doubleValue());
						C.put(sObj, new Double(fncn));
					}
				}

				/**
				 * This next commented out line was causing a ConcurrentModificationException
				 * We were deleting a key while looping through them.
				 * Ryan does this in his code too, but he is looping through an enum object
				 * of the keys but removing this key from the hashtable, so 2 different objects.
				 * I'm not sure why this is done in this code here or even in his code
				 * as the colA objects are all reset anyway after we are done looping.
				 * I currently have it commented out as I'm not sure what to do with it.
				 * I don't think we need this line as I don't think it does anything at all.
				 *
				 */
				//colA.remove(scoreAObj);
                
			}
			PGFs.remove(colA);
			PGFs.remove(colB);
			PGFs.addFirst(C);
		}
		Map<Double,Double> finalPGF = PGFs.getLast();
		//reach into the PGF and grab the p-value for our score
		return finalPGF.get(new Double(max)).doubleValue();
	}
}
