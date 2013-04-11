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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dasp.algorithms.Align;

public class Alignment {
	protected Align method = null;
	protected List<String>seqs = null;
	protected Map<String,String>nameMap = null;
	protected Map<String,String>alignmentMap = null;
	private double percentIdentity = 0.0;
	private double alignmentScore = 0.0;

	private int match = 0;
	private int strong = 0;
	private int weak = 0;
	private int gaps = 0;
	private int length = 0;

	public Alignment () {
		seqs = new ArrayList();
		this.alignmentMap = new HashMap();
		this.nameMap = new HashMap();
	}

	public void doAlign() {
		alignmentScore = method.align(seqs);
		alignmentMap = method.getAlignment();
		percentIdentity = method.getIdentity();
	}

	/**
	 * addSequence can be used to add an addtional sequence to the alignment.  Note that
	 * the alignment is *not* recalculated when you add a new sequence.
	 *
	 * @param seq the sequence to add as a FASTA-formatted string
	 */
	public void addSequence (String seq) {
		seqs.add(seq);
		updateNameMap(nameMap, seq);
	}

	/**
	 * addSequence can be used to add an addtional sequence to the alignment.  Note that
	 * the alignment is *not* recalculated when you add a new sequence.
	 *
	 * @param name the name (identifier) of the sequence
	 * @param rawSeq the sequence to add as an unformatted string
	 */
	public void addSequence (String name, String rawSeq) {
		addSequence(name+"\t"+rawSeq);
	}


	/**
	 * deleteSequence can be used to remove a sequence from the alignment.  Note that
	 * the alignment is *not* recalculated when you remove a new sequence.
	 *
	 * @param seq the sequence to remove as a FASTA-formatted string
	 */
	public void deleteSequence(String seq) { }

	/**
	 * addSequences can be used to add a list of addtional sequences to the alignment.  Note that
	 * the alignment is *not* recalculated when you add new sequences.
	 *
	 * @param seqList the sequences to add as a list of FASTA-formatted strings
	 */
	public void addSequences (List<String> seqList) {
		seqs.addAll(seqList);
		updateNameMap(nameMap, seqList);
	}

	/**
 	 * Return all of the sequences as FASTA-formatted strings
 	 *
 	 * @return FASTA-formatted sequences as a list
 	 */
	public List<String> getSequences() {return seqs;}

	/**
 	 * Return the % identity for all of the sequences in the alignment.
 	 *
 	 * @return the % identity
 	 */
	public double getIdentity() {
		return percentIdentity;
	}

	/**
 	 * Return the alignment score that results from the alignment.  Usually this will be the
 	 * expectation value (e-value) or some other measure of the probably that this alignment
 	 * could result from a random assembly of sequences.
 	 *
 	 * @return the score
 	 */
	public double getScore() {
		return alignmentScore;
	}

	/**
 	 * Return the alignment itself as a string with newlines between sequences and dashes inserted
 	 * to account for insertions.
 	 *
 	 * @return the alignment
 	 */
	public String getAlignmentAsString() { 
		String result = "";
		for (String key: alignmentMap.keySet()) {
			if (key.equals(Align.CONSERVATION))
				continue;
			result += key+": "+alignmentMap.get(key)+"\n";
		}
		return result;
	}

	/**
 	 * Return the alignment itself as a list of strings with dashes inserted
 	 * to account for insertions.
 	 *
 	 * @return the alignment as a list of strings
 	 */
	public List<String> getAlignment() { 
		List<String>alignment = new ArrayList();
		for (String key: alignmentMap.keySet())
			if (!key.equals(Align.CONSERVATION)) alignment.add(alignmentMap.get(key));
		return alignment;
	}

	/**
 	 * Return the alignment string for just one of the proteins in the alignment.  This assumes that
 	 * each of the entries in the fasta file were named.
 	 *
 	 * @param name the name of this protein (usually from the fasta file)
 	 * @return the alignment string for name
 	 */
	public String getAlignmentByName(String name) { 
		if (alignmentMap.containsKey(name))
			return alignmentMap.get(name); 
		return null;
	}

	/**
	 * Return the entire alignment map, including the conservation line.
	 *
	 * @return the alignment map
	 */
	public Map<String,String> getAlignmentMap() {
		return alignmentMap;
	}

	/**
 	 * Return the alignment itself as an HTML-formatted string
 	 * to account for insertions.
 	 *
 	 * @return the HTML formatted alignment
 	 */
	public String getHTMLAlignment() { 
		if (alignmentMap != null)
			return method.getHTMLAlignment(getAlignmentMap());
		return null; 
	}

	/**
	 * Return the width of the alignment
	 *
	 * @return the alignment width
	 */
	public int getAlignmentWidth() {
		if (alignmentMap != null && alignmentMap.size() > 0) {
			for (String key: alignmentMap.keySet()) {
				if (!key.equals(Align.CONSERVATION)) {
					return alignmentMap.get(key).length();
				}
			}
		}
		return -1;
	}

	/**
 	 * This method updates a name map from a fasta sequence
 	 *
 	 * @param map the name map to update
 	 * @param seq the sequence that needs to be updated
 	 * @return the name map
 	 */
	private void updateNameMap(Map<String,String>map, String seq) {
		String name = fastaName(seq);
		String rawSeq = fastaSeq(seq);
		if (name != null)
			map.put(name, rawSeq);
	}

	/**
 	 * This method updates a name map from a list of fasta sequences
 	 *
 	 * @param seqs list of sequences
 	 * @return the name map
 	 */
	private void updateNameMap(Map<String,String>map, List<String>seqs) {
		for (String seq: seqs) {
			String name = fastaName(seq);
			String rawSeq = fastaSeq(seq);
			if (name != null)
				map.put(name, rawSeq);
		}
	}

	private String fastaName(String fastaString) {
		int nameIndex = fastaString.indexOf("> ");
		if (nameIndex == -1) return null;
		int nameEnd = fastaString.indexOf("\n", nameIndex);
		return fastaString.substring(nameIndex, nameEnd);
	}

	private String fastaSeq(String fastaString) {
		String rawSequence = fastaString;
		int nameIndex = fastaString.indexOf("> ");
		if (nameIndex >= 0){
			int nameEnd = fastaString.indexOf("\n", nameIndex)+1;
			rawSequence = fastaString.substring(nameEnd);
		}

		// Now strip all spaces and newlines
		rawSequence = rawSequence.replaceAll(" ", "");
		rawSequence = rawSequence.replaceAll("\n", "");
		return rawSequence;
	}

	private String seq2FastA(String name, String rawSeq) {
		return ("> "+name+"\n"+rawSeq);
	}
}
