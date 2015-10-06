/**
 * The PSSMSearch interface provides a general method to 
 * use position-specific scoring matrices to search sequence
 * databases.
 *
 */

package dasp.algorithms;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import dasp.model.PSSM;
import dasp.model.DBSearchResult;
import dasp.model.SearchResult;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;

import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;

public class FastAFileSearcher implements DBSearch {
	public List<DBSearchResult> search(File database, List<PSSM>pssmList, 
									   PSSMSearch searchAlg, double threshold, int numThreads) throws Exception {
		//
		// Open database - done
		// Find a sequence - done
		// 	 Create a search result for that sequence
		//	 Massage the sequence
		//	 Search using searchAlg
		//	 Apply results to result
		//	 If search significant, add to list
		
		//open the database file for searching
		//get the fasta file iterator
		SequenceIterator fastaIter = parseDBfile(database, "fasta");

		//Create a database search results list
		//only significant sequence search results are added to this object
		List<DBSearchResult> DBresults = new ArrayList();
		int sequenceNumber = 0;

		ExecutorService[] threadPools = new ExecutorService[numThreads];
		for (int pool = 0; pool < threadPools.length; pool++)
			threadPools[pool] = Executors.newFixedThreadPool(1);

		//iterate through each sequence in the file
		int maxQueueDepth = 1000;
		while(fastaIter.hasNext()){
			for (int queueCount = 0; queueCount < maxQueueDepth; queueCount++) {
				for (int pool = 0; pool < threadPools.length; pool++) {
					if (fastaIter.hasNext()) {
						Sequence curSeq = fastaIter.nextSequence();
						Runnable r = new ParallelSearcher(curSeq, pssmList, searchAlg, threshold, DBresults, sequenceNumber++);
						threadPools[pool].submit(r);
					}
				}
			}

			// OK, now wait for all the threads to complete
			for (int pool = 0; pool < threadPools.length; pool++) {
				threadPools[pool].shutdown();
				try {
					boolean result = threadPools[pool].awaitTermination(7, TimeUnit.DAYS);
				} catch (Exception e) {}
				threadPools[pool] = Executors.newFixedThreadPool(1);
			}
		}
		
		return DBresults;
	}

  /**
   * Get a Sequence Iterator over all the sequences in the file.
   * SeqIOTools.fileToBiojava() returns an Object. If the file read
   * is an alignment format like MSF and Alignment object is returned
   * otherwise a SequenceIterator is returned.
	 *
	 * @param inFile the file to read from
	 * @param fileFormat a string that represents the file format to read (e.g. "fasta").
	 * @return an iterator over all of the sequences in the file
   */
	private SequenceIterator parseDBfile(File inFile, String fileFormat) throws IOException {

		String alphabet = "protein";
		SequenceIterator iter = null;
		try {
		  iter = (SequenceIterator)SeqIOTools.fileToBiojava(fileFormat,alphabet, new BufferedReader(new FileReader(inFile)));
		}
		catch (BioException ex) {
			throw new IOException(ex.getMessage());
		}
		catch (FileNotFoundException ex) {
			throw new IOException(ex.getMessage());
		}
		return iter;
	}

	//QFAST algorithm: for finding the p-value of a product of p-values
	private double QFAST(int n, double p){
		double x = 0;
		if(p == 0)
			return 0;
		if(n > 1)
			x = -1*Math.log(p);
		double t = p;
		double q = p;
		for(int i = 1; i < n; i++){
			t *= (x/(double)i);
			q += t;
		}
		return q;
	}

	class ParallelSearcher implements Runnable {
		final Sequence curSeq;
		final List<PSSM> pssmList;
		final PSSMSearch searchAlg;
		final double threshold;
		final List<DBSearchResult> DBresults;
		int xCount = 0;
		int sequenceNumber;

		public ParallelSearcher(Sequence curSeq, List<PSSM>pssmList, 
		                        PSSMSearch searchAlg, double threshold, List<DBSearchResult> DBresults,
		                        int sequenceNumber) {
			this.curSeq = curSeq;
			this.pssmList = pssmList;
			this.searchAlg = searchAlg;
			this.threshold = threshold;
			this.DBresults = DBresults;
			this.sequenceNumber = sequenceNumber;
		}

		public void run() {
			/**
			 * ISSUE: getName() seems to truncate at the first space in the
			 * fasta name line.  Have not looked at how to get around that.
			 * Currently I have pre-inserted 4 $ chars for each space in the
			 * test DB fasta file, and then I parse those out and replace with spaces before
			 * I print out the name.  I don't like to do that, but we need
			 * all the information on each line.  Will look into fixing this issue later.
			 */
			boolean debug = false;
			String seqName = curSeq.getName().replaceAll("\\${4}", " ");
			// System.out.println("\n\nSeqName:"+seqName+"\nFullSeq:"+curSeq.seqString());

			//We have a seqeunce, now munch it and get the xCount
			this.xCount = 0;
			String munchSeq = munchSequence(curSeq.seqString());
			if (debug)
				System.out.println("\n\nSeqName:"+seqName+"\nmunchSeq:"+munchSeq);

			//munchSeq will be null if it is a DNA sequence, so we want to skip that.
			if(munchSeq!=null){
				// XXX Should this be in a separate thread?
				//now we can call the search routine
				SearchResult[] seqResults = searchAlg.search(pssmList, munchSeq, xCount, debug);

				//The search alg will return null for the seqResults if one or
					//more of the profile fragments did not match to the seqeunce.
				if(seqResults!=null){
					//Need to get the product of the pvalues to pass to QFAST
					double product=1.0;
					for(SearchResult r: seqResults){
						product *= r.getPvalue();
					}
					double finalPval = QFAST(pssmList.size(), product);

					//print out results we have for this sequence.
					// System.out.println("\nSearch Results for sequence:"+seqName+"\nFinal Pvalue="+finalPval);
					//int sIndex = 0;
					//for(SearchResult s: seqResults){
					//	System.out.println("Match "+sIndex+":Index "+s.getIndex()+":Length "+s.getPssmLength()+":Pvalue "+s.getPvalue()+":SubStr "+munchSeq.substring(s.getIndex(), s.getIndex()+s.getPssmLength()));
					//	sIndex++;
					//}


					//if the pval is significant then create and DBSearchResults object and add to list.
					if(finalPval > 0.0 && finalPval < threshold){
						System.out.println("Seq "+seqName+" PASSED.  Final pValue = "+finalPval);
						DBSearchResult result = new DBSearchResult(finalPval, seqResults, seqName, curSeq.seqString());
						synchronized (DBresults) {
							DBresults.add(result);
						}
					}
				}
			}
			// System.out.print("Finished "+sequenceNumber+"\r");
		}

		private String munchSequence(String sequence) {
			char[] inputSequence = sequence.toCharArray();
			char[] outputSequence = new char[inputSequence.length+1];
			int outputIndex = 0;
			boolean notDNA = false;
			for (int inputIndex = 0; inputIndex < inputSequence.length; inputIndex++) {
				char seq = inputSequence[inputIndex];
				switch (seq) {
				case 'Z':
				case 'B':
				case 'J':
				case 'O':
				case 'U':
				case 'X':
					outputSequence[outputIndex++] = 'X';
					this.xCount++;
					break;
	
				case 'D':
				case 'E':
				case 'F':
				case 'H':
				case 'I':
				case 'K':
				case 'L':
				case 'M':
				case 'N':
				case 'P':
				case 'Q':
				case 'R':
				case 'S':
				case 'V':
				case 'W':
				case 'Y':
					notDNA = true;
				case 'A':
				case 'C':
				case 'T':
				case 'G':
					outputSequence[outputIndex++] = seq;
					break;
	
				}
			}
	
			outputSequence[outputIndex] = 0;
	
			if (notDNA)
				return new String(outputSequence).trim();
	
			return null;
		}
	}
}
