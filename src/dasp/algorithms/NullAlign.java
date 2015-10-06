/**
 * This is an implementation of an alignment that actually does nothing.
 * This assumes that the input is already aligned.
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

public class NullAlign implements Align {
	private double alignmentScore = 0.0;
	private Alignment alignment = null;
	private double percentID = 0.0;
	private Map<String, String> alignMap;
	private int match=0,strong=0,weak=0,totalGaps=0,alnLength=0;

	/**
	 * Will need to store the alignment or files the alignemnt is saved to?
	 */
	public NullAlign() {
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
		for (String seq: sequences) {
			String[] line = seq.split("\t");
			alignMap.put(line[0], line[1]);
		}

		return 1.0;
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
		return new NullAlign();
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
