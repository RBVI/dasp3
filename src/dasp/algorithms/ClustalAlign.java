/**
 * This is the implementation of an alignment using ClustalW.
 */

package dasp.algorithms;

import dasp.model.Alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.PrintStream;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ClustalAlign implements Align {
	private double alignmentScore = 0.0;
	private Alignment alignment = null;
	private double percentID = 0.0;
	private Map<String, String> alignMap;
	private int match=0,strong=0,weak=0,totalGaps=0,alnLength=0;

	/**
	 * Will need to store the alignment or files the alignemnt is saved to?
	 */
	public ClustalAlign() {
		alignMap = new HashMap();
	}

	/**
	 * align actually performs the alignment.  This will almost certainly throw some sort
	 * of exception, depending on the alignment algorithm or program
	 * 
	 * @param model the data model we're using
	 * @return the alignment score
	 */
	public double align(List<String> sequences) {

		/**
		 * Can possibly use the ClustalW java client, but I'm not sure how to get that working.
		 * See the EBI CLustalW site for more information.
		 */
	
		/**
		 * First need to save all seqs to a FASTA file for use with Clustal.
		 */
		//Need to find a way to get the original file name for use with this, so that each Clustal run can have a unique name.
		File tmpFastaFile = null;
		File alignmentFile = null;
		try {
			tmpFastaFile = File.createTempFile("tmpFastaFile",".fasta");
			alignmentFile = File.createTempFile("tmpFastaFile",".aln");
		} catch (IOException e) {
				System.err.println("Error creating temp files.");
				e.printStackTrace();
				System.exit(-1);
		}

		try{
			FileOutputStream FOS = new FileOutputStream(tmpFastaFile);
			PrintStream PS = new PrintStream(FOS);
	
			for (String seq: sequences){
				String sp[] = seq.split("\t");
				PS.println(">" + sp[0]);
				PS.println(sp[1]);
			}
			FOS.close();
		} catch(IOException e) {
				System.err.println("Error writing fasta signatures to file.");
				e.printStackTrace();
				System.exit(-1);
		}
        
 		try {
             
			String s = null;
			//Need to make sure the file exists first!
			// XXX The location of the clustalw binary should be a property
			String cmd[] = {"/usr/local/bin/clustalw2", "-INFILE="+tmpFastaFile.getPath(), "-OUTFILE="+alignmentFile.getPath()};
			Process p = Runtime.getRuntime().exec(cmd);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
			FileOutputStream FOS = new FileOutputStream("ClustalScreenOutput.cw");
			PrintStream PS = new PrintStream(FOS);

			//prints Clustal screen output to a file
			while((s = stdInput.readLine()) != null){
				//System.out.println(s);
				PS.println(s);
			}
			//prints clustal error messages to a file AND to the screen
			while((s = stdError.readLine()) != null){
				System.err.println(s);
				PS.println(s);
			}
			//Need to close process streams
			p.getErrorStream().close();
			p.getInputStream().close();
			p.getOutputStream().close();

		} catch(IOException e){
			System.err.println("Exception happened running CLUSTALW: ");
			e.printStackTrace();
			System.exit(-1);
		}
         
		/**
		 * At this point we should have the following files saved in the outDir:
		 * tmpFastaFile.aln
		 * tmpFastaFile.dnd
		 * tmpFastaFile.cw
		 *
		 * The next step DASP takes is to calculate the ASP score for the alignment
		 * and appends information to the .aln file.
		 * I think if this appended file is to be the important one, we should
		 * have it named something unique so it doesn't get overridden.  Not sure
		 * how to do that at the moment.
		 *
		 * I think the ASP score is calculated from the .aln file.  Should we calculate
		 * that here or in the ActiveSiteProfile class?  I think it would make sense to calculate it here
		 * but I'm not sure all alignemnt methods will give output compatable with the algorithm.
		 * I will put in the method here, then we can change it later if we need to.
		 * I guess the score this method returns should be related to the specific alignment algorithm
		 * However, clustal's score is always just the percent ID I think.
		 * Have to think about this more, but either way I need a method to claculate the ASP score.
		 * Maybe that should be part of the Align interface.
		 *
		 */

		parseAlnFile(alignmentFile, sequences);

		formatAlignment(sequences);

		return getASPscore(alignmentFile, sequences.size());
	}

   /**
    * Return the ASP score as defined in Cammer et al for this alignment.
    * method basically parses the .aln file and processes it by blocks.
    * Is there an easier way to do this?  Try to find out later as this should work for now.
    * Need to look closer at this as it came from the identityTableMaker class which is pairwise ASP scores.
    * It looks like it will work for the full alignment though, so not sure why the DASP code does it differently.
    * 
    *
    * @return the ASP score
    */
  public double getASPscore(File alnFile, int numSequences){

		double Si = 1.0;
		double Ss = 0.2;
		double Sw = 0.1;
		double Sg = -0.5;
		return (Si*match + Ss*strong + Sw*weak + Sg*totalGaps)/alnLength;
  }
  
  /**
	 * Return the % identity for all of the sequences in the alignment.
	 *
	 * @return the % identity
	 */
  public double getIdentity() {
    return 0.0;
  }


  /**
	 * Return the alignment itself with each sequence as a single string and with dashes inserted
	 * to account for insertions.
     * Need to include the PDB id or structure id at the start of each line followed by some
     * delimiter like spaces.  findProfileFragments() currently assumes spaces.
     * Basically reads in the .aln file and formats it appropriatly.  Should we do the fragment thing here or elsewhere?
     * We have the fragments, so why not do them at this step?
	 *
	 * @return the alignment
	 */
  public Map<String, String> getAlignment() { return alignMap; }


	/**
 	 * Return the alignment itself as an HTML-formatted string
 	 * to account for insertions.
 	 *
 	 * @param alignment the string alignment
 	 * @return the HTML formatted alignment
 	 */
	public String getHTMLAlignment(Map<String, String> alignment) { return null; }

	public Align NewAlignFactory() {
		return new ClustalAlign();
	}

	/**
	 * Parse the ALN file
	 */
	private void parseAlnFile(File alnFile, List<String>sequences) {
		for (String seq: sequences) {
			String sp[] = seq.split("\t");
			alignMap.put(sp[0],"");
		}
		alignMap.put(CONSERVATION, "");

		int numSequences = sequences.size();

		String lastLine = null, inLine = null, fourthLine = null;
		int index; //needs to be more descriptive, what is last?

		try{
			FileReader fr = new FileReader(alnFile);
			LineNumberReader lnr = new LineNumberReader(fr);

			while((inLine = lnr.readLine()) != null){
				if(!inLine.equals("") && lnr.getLineNumber() != 1) {
					// System.out.println(lnr.getLineNumber() + "\t" + inLine);

					String key = addLine(inLine, alignMap);

					// put the block data in cache
					for(int i=0; i< numSequences; i++){
						inLine = lnr.readLine();
						key = addLine(inLine, alignMap);
					}

					// manipualte block data
					// get identity, strong, weak
					lastLine = alignMap.get(CONSERVATION);
					for(int i=0; i<lastLine.length(); i++){
						switch(lastLine.charAt(i)){
							case '*': {
								match++;
								break;
							}
							case ':': {
								strong++;
								break;
							}
							case '.': {
								weak++;
								break;
							}
						}
					}

					// get gaps
					int[] gaps = null;
					for(String line: alignMap.values()) {
						if (gaps == null)  {
							alnLength = line.length();
							gaps = new int [alnLength];
							for (int i = 0; i < gaps.length; i++) { gaps[i] = 0; }
						}
						char[] temp = line.toCharArray();
						for(int j=0; j<temp.length; j++){
							if(temp[j] == '-'){
								gaps[j] = 1;
							}
						}
					}
					for (int i=0; i<gaps.length; i++){
						if(gaps[i] == 1){
							totalGaps++;
						}
					}

					// get N
				}
			}
			lnr.close();
			fr.close();
		} catch (IOException x) {
			System.out.println("ASP Scoring Error");
		}


	}

	private String addLine(String input, Map<String,String>alignMap) {
		// See if we have a sequence line or a total line
		String key = CONSERVATION;
		// System.out.println(input);
		if (!input.startsWith(" ")) {
			int start = input.indexOf(' ');
			key = input.substring(0,start);
			input = input.substring(start);
		}
		input = input.trim();

		String newLine = alignMap.get(key)+input;
		alignMap.put(key, newLine);
		return key;
	}

	private void formatAlignment(List<String> sequences) {
		for (String seq: sequences) {
			String sp[] = seq.split("\t");
			String alignedSeq = alignMap.get(sp[0]);
			alignMap.put(sp[0], formatSequence(sp[1], alignedSeq));
		}
	}

	private String formatSequence(String formatted, String aligned) {
		char seq[] = new char[aligned.length()];
		int formattedPointer = 0;
		for (int i = 0; i < aligned.length(); i++) {
			if (aligned.charAt(i) == '-')
				seq[i] = '-';
			else
				seq[i] = formatted.charAt(formattedPointer++);
		}
		return new String(seq);
	}

	private void printAlignment() {
		for (String str: alignMap.keySet()) { 
			if (!str.equals(CONSERVATION))
				System.out.println(str+": "+alignMap.get(str)); 
		}
		
		System.out.println(CONSERVATION+": "+alignMap.get(CONSERVATION)); 
		System.out.print("identity: " + match);
		System.out.print("  strong: " + strong);
		System.out.print("  weak: " + weak);
		System.out.print("  gap: " + totalGaps);
		System.out.print("  total: " + alnLength + "\n");
	}
}
