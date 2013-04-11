/**
 * The Align interface provides a general method to call various sequence alignment algorithms
 *
 */

package dasp.algorithms;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import dasp.model.Alignment;

public interface Align {
	public static final String CONSERVATION = "Conservation";
	/**
	 * align actually performs the alignment.  This will almost certainly throw some sort
	 * of exception, depending on the alignment algorithm or program
	 *
	 * @param sequences list of sequences in FASTA format
	 * @return the alignment score
	 */
	public double align(List<String>sequences);

  /**
	 * Return the % identity for all of the sequences in the alignment.
	 *
	 * @return the % identity
	 */
  public double getIdentity();

  /**
   * Return the ASP score as defined in Cammer et al for this alignment.
   *
   * @return the ASP score
   */
  public double getASPscore(File alnFile, int numSequences);

	/**
 	 * Return the alignment itself as a list of aligned sequences with dashes
 	 * to account for insertions.
 	 *
 	 * @param alignedSequences list of aligned sequences
 	 * @return list of aligned sequences
 	 */
	public Map<String, String> getAlignment();

	/**
 	 * Return the alignment itself as an HTML-formatted string
 	 * to account for insertions.
 	 *
 	 * @param alignedSequences list of aligned sequences
 	 * @return the HTML formatted alignment
 	 */
	public String getHTMLAlignment(Map<String, String>alignedSequences);

	/**
	 * This method is useful when you need to get a fresh alignment algorithm
	 * object from an existing object (i.e. from inside of a loop).
	 *
	 * @return new Align object
	 */
	public Align NewAlignFactory();
}
