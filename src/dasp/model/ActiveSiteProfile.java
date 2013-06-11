/*
 File: Dasp.java

 Copyright (c) 2009, Wake Forest University

 This library is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published
 by the Free Software Foundation; either version 3.0 of the License, or
 any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 documentation provided hereunder is on an "as is" basis, and 
 neither Wake Forest University nor the University of Califonia
 have any obligations to provide maintenance, support,
 updates, enhancements or modifications.  In no event shall
 Wake Forest University or the University of Califonia
 be liable to any party for direct, indirect, special,
 incidental or consequential damages, including lost profits, arising
 out of the use of this software and its documentation, even if the
 Wake Forest University and/or the University of Califonia
 have been advised of the possibility of such damage.  See
 the GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

package dasp.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.LinkedList;

import dasp.algorithms.Align;

/**
 * An ActiveSiteProfile is an alignment that is formed by aligning
 * the ActiveSiteSignatures of a group of structures.
 */
public class ActiveSiteProfile extends Alignment {

	// NOTE: this assumes that each profile contains a single column (one fragment)
	static public ActiveSiteProfile concatenate(List<ActiveSiteProfile> profileList) {
		int sequenceCount = profileList.get(0).getSequences().size();
		List<String> newAlignment = new ArrayList<String>();
		for (int row = 0; row < sequenceCount; row++) {
			boolean upper = true;
			String aln = "";
			for (ActiveSiteProfile profile: profileList) {
				String seq = profile.getAlignmentString(row);
				if (upper) {
					aln += seq.toUpperCase();
					upper = false;
				} else {
					aln += seq.toLowerCase();
					upper = true;
				}
			}
			newAlignment.add(aln);
		}

		return new ActiveSiteProfile(profileList.get(0), newAlignment);
	}

	public ActiveSiteProfile(ActiveSiteProfile template, List<String>alignment) {
		super();
		List<String> sequences = template.getSequences();
		for (int strIndex = 0; strIndex < sequences.size(); strIndex++) {
			String seq = sequences.get(strIndex);
			String[] splitStr = seq.split("\t");
			String pdbId = splitStr[0].trim();
			addSequence(pdbId+"\t"+alignment.get(strIndex));
			alignmentMap.put(pdbId, alignment.get(strIndex));
			// NOTE: not updating name map at this point
		}
	}

	public ActiveSiteProfile(List<ActiveSiteSignature>signatureList, Align alignmentMethod, double radius) {
		super();
		method = alignmentMethod;

		for (ActiveSiteSignature sig: signatureList) {
		  //sig.getSignature returns a String with fragments in upper and lower case.
		  //for running through Clustal it needs to be in fasta format.
		  //possibly implement that following so that each Sequence entry is already in fasta format?
		  //Actually, the Alignment class assumes the sequence string is in fasta format.
		  //possibly should add a getFastSignature() method to ActiveSiteSignature() class, and have the signature calculated in advance.
		  // System.out.println("ASP adding sequence: "+sig.getPdbId()+"\t"+sig.getSignature(radius));
		  addSequence(sig.getPdbId()+"\t"+sig.getSignature(radius));
			// System.out.println("Adding sequence: "+sig.getPdbId()+": "+sig.getSignature(radius));
		}

		doAlign();
		// System.out.println(getAlignmentAsString());

		for (ActiveSiteSignature sig: signatureList) {
			if (sig.getFragments().size() > 1) break;

			String alignment = alignmentMap.get(sig.getPdbId());
			String newAlignment = sig.extendAlignment(alignment);
			alignmentMap.put(sig.getPdbId(), newAlignment);
		}
		// System.out.println(getAlignmentAsString());
	}

	public ActiveSiteProfile(File profilePath) {
		super();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(profilePath));
			String line = null;
			int lineNumber = 0;
			while ((line = reader.readLine()) != null) {
				lineNumber++;
				if (line.startsWith("#")) continue;
				line = line.trim();
				if (line.length() == 0) continue;
				System.out.println("Line: "+line);
				String[] split = line.split("\t");
				if (split.length != 2) {
					System.err.println("Input error in "+profilePath+" line #"+lineNumber+": input not tab-delimited"); 
					continue;
				}
				addSequence(split[0]+"\t"+split[1]);
			}
		} catch (Exception e) {
			System.err.println("Unable to read input file "+profilePath+": "+e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	/**
	 * Identifies the motifs (profile fragments) of the ASP
	 *
	 * Take the active site profile that has sequence fragments identified and find the profile 
	 * fragments or motifs within the alignment.
	 * Steps:
	 * 1) We start with the fragment formatted alignment stored in 'alignmentMap'
	 * 2) We need to apply Ryan's algorithm for finding motifs to the list of aligned strings.
	 * 3) once a motif is identified we need to put that in a new Alignment object.
	 * 4) each fragment needs to then be realigned.
	 *
	 * See Amy's powerpoint on the motif identification algorithm for a chart of how this works in general.
	 *
	 * Currently this code has basically been pasted into here with the initial
	 * section modified to get it to work with the current alignmentMap structure set up.
	 * TO-DO: Need to get it to return a List of alignments for use in the PSSM searching.
	 *
	 * @return a list of Alignment objects where each contains the fragments in the motif.
	 * Note the profile fragments returned are not yet aligned.  The aligning is done in Main.
	 */
	public List<Alignment> findProfileFragments(Align alnMethod) {

		List<Alignment>FragList = null;
		int numPDBs = seqs.size();

		List<Alignment> list = new ArrayList(); // this is a list of arrays that hold each fragment of the profile.

		//Since the rest of the code is based off the Linked List and String array structures
		//for simplicity sake I will just convert the 'alignmentMap' into those structures.
		//we can change it later into something better when we have time.

		String sequences[] = new String[numPDBs];
		List<String> names = new ArrayList();

		int m = 0;
		for(String key: alignmentMap.keySet()){
			if(key.equals(Align.CONSERVATION))
				continue;
			names.add(key);
			sequences[m] = alignmentMap.get(key);
			m++;
		}

		int end = sequences[0].length();  //the length of the alignment?
		int column = 0;

		// for each column {
		//   boxSize = (shortestFrag > 3)
		//   for each row {
		//      frag = getFrags(row, column, boxSize)
		//   }
    // }

		//loops through each column i of the alignment
		while(column < end){
			Alignment alignment = new ProfileFragment(alnMethod.NewAlignFactory());
			int minBox = end;
			int minIndex = -1; //this is the row/sequence with the current minimum length fragment starting at position i
			//loops through each sequence in the alignemnt (each row k)
			for(int row = 0; row < sequences.length; row++){
				int fragmentLength = pieceOfFrag(sequences[row], column); //temp must be the length of the fragment.
				if(fragmentLength >= 4){
					if(fragmentLength < minBox){
						minIndex = row;
						minBox = fragmentLength;
					}
				}
			}
			//if there is a valid fragment for at least one of the signatures then do this procedure.
			if(minIndex != -1) {
				int boxEnd = getEnd(sequences[minIndex], column);  //is this index position of the motif bounds?
				int counter = 1;
				//loops through each sequence again for position i
				for(int row = 0; row < sequences.length; row++) {
					if(row != minIndex) {
						int boxColumn = column;
						int tempBound = -1;
						//this right here is why the end fragment is used instead of the longest on within the bounds.
						//
						while(tempBound == -1 && boxColumn <= boxEnd) {
							tempBound = pieceOfFrag(sequences[row], boxColumn);
							if(tempBound >= 3) {
								tempBound = getEnd(sequences[row], boxColumn);
	
								//HAHAHA!  This is the problem right here!  So say the bounding box
								//is only 3 long with positions 1, 2 and 3.  Say with sequence k at position 1
								//you only have a fragment that is 2 long, so it is skipped.  Same at position 2 the
								//fragment is now only 1 long.  BUT at position 3 the full fragment is 5 long,
								//however, the last 4 letters are outside the bounding box.  So, this method here
								//just chops off that last 4 letters and only includes the first letter in the final motif.
								//With LONG motifs this may be ok, but with short motifs like 3 or 4 this sounds
								//like a really bad way to do it!  The motif should be at least 3 long INSIDE
								//the bounding box to be accepted.  Otherwise you can end up with motifs like I have
								//with the kinases where maybe 4 sequences have a length 3 motif, but the rest only
								//only have 1 or 2 letters.  I would think this would be insignificant, especially
								//for alignment by Clustal.
								if(tempBound > boxEnd){
									tempBound = boxEnd;
								}
							}
							else
								tempBound = -1;
							boxColumn++;
						}
						if(tempBound != -1) {//if tempBound is not -1 then we are adding it to the motif
							counter++;
							alignment.addSequence(names.get(row), unAlignString(sequences[row], boxColumn-1, tempBound));
						}
						//
					} else {
						alignment.addSequence(names.get(minIndex), unAlignString(sequences[minIndex], column, boxEnd));
					}
				}
				//if at least 50% of the sequences are included in the motif then consider it a motif for searching
				if(counter > Math.ceil((double)sequences.length/(double)2)) {
					list.add(alignment);
				}
				column = (boxEnd + 1); //increment i to be just after the last motif box searched
				//I actually think this is another bad way to do it because if the motif was not added then we
				//should just increment i by 1 in case another motif is found at that position.
			} else
				column++;
		}  //end while i < end

		return list;

		//This list is a list of arrays containing the first index and last index of each motif block
		//in the alignment.  The names in h don't seem to be used, so I'm wondering why they are there.
		//After this the extractResidues() method is called which basically take the stored indexes
		//and goes through each alignment sequence using the original fragTable hashtable for the full
		//aligned signatures and pdb ids.  All the letters within the specified index range for a motif
		//are written directly to a fasta file.
		//I think we should be able to do this on the fly in this method and save this data to an object
		//instead of writing it directly to a file.
		//The extractResidues() method is pasted below and commented out for reference.
	}

	//check if a piece of a fragment exists at the current position in the profile
	private int pieceOfFrag(String s, int index){
		boolean u = Character.isUpperCase(s.charAt(index));
		int gapCounter = 0;
		int curPos = index + 1;
		int length = 1;
		if(!Character.isLetter(s.charAt(index)))
			return 0;
		if(u){
			while(curPos < s.length() && (!Character.isLetter(s.charAt(curPos)) ||
							Character.isUpperCase(s.charAt(curPos)))){
				if(!Character.isLetter(s.charAt(curPos)))
					gapCounter++;
				else{
					gapCounter = 0;
					length++;
				}
				curPos++;
			}
		}
		else{
			while(curPos < s.length() && (!Character.isLetter(s.charAt(curPos)) ||
							Character.isLowerCase(s.charAt(curPos)))){
				if(!Character.isLetter(s.charAt(curPos)))
					gapCounter++;
				else{
					gapCounter = 0;
					length++;
				}
				curPos++;
			}
		}
		return length;
	}

	//check where the current fragment ends in the profile
	private int getEnd(String s, int index){
		boolean u = Character.isUpperCase(s.charAt(index));
		int gapCounter = 0;
		int curPos = index + 1;
		if(u){
			while(curPos < s.length() && (!Character.isLetter(s.charAt(curPos)) ||
							Character.isUpperCase(s.charAt(curPos)))){
				if(!Character.isLetter(s.charAt(curPos)))
					gapCounter++;
				else
					gapCounter = 0;
				curPos++;
			}
		}
		else{
			while(curPos < s.length() && (!Character.isLetter(s.charAt(curPos)) ||
							Character.isLowerCase(s.charAt(curPos)))){
				if(!Character.isLetter(s.charAt(curPos)))
					gapCounter++;
				else
					gapCounter = 0;
				curPos++;
			}
		}
		if(gapCounter == 0)
			return (curPos - 1);
		else
			return (curPos - 1 - gapCounter);
	}

	private String unAlignString (String alignedString, int start, int end) {
		String str = alignedString.substring(start, end+1).replace("-","");
		return str;
	}
}
