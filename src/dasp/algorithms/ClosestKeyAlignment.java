/**
 * Align fragments based on the closest key residue.  Essentially, align them such
 * that those closest to residue 1 are aligned first, etc.
 *
 */

package dasp.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import dasp.model.ActiveSiteProfile;
import dasp.model.ActiveSiteSignature;

public class ClosestKeyAlignment {
	private Map<ActiveSiteProfile, List<ActiveSiteSignature>> profileMap = null;
	private Map<ActiveSiteSignature, List<ActiveSiteProfile>> signatureMap = null;
	private SortedMap<Double, List<ActiveSiteProfile>> orderedMap = null;
	private int fragCount;
	private double radius;

	public ClosestKeyAlignment (int fragCount, double radius, List<ActiveSiteSignature> signatures) {
		this.fragCount = fragCount;
		this.radius = radius;
		profileMap = new HashMap<ActiveSiteProfile, List<ActiveSiteSignature>>();
		signatureMap = new HashMap<ActiveSiteSignature, List<ActiveSiteProfile>>();
		orderedMap = new TreeMap<Double, List<ActiveSiteProfile>>();
		getClosestKeyAlignments(signatures);
	}

	public void getClosestKeyAlignments(List<ActiveSiteSignature> sigList) {
		// For each key residue
		//    For each signature
		//       Get the closest fragment & add it to the list
		//    Create a profile & update
	}

/*
	public void getAllAlignments(List<ActiveSiteSignature> fragList, 
	                             List<ActiveSiteSignature> sigList) {

		ActiveSiteSignature sig = sigList.get(0);
		for (int fragIndex = sig.keyFragCount(); fragIndex < fragCount; fragIndex++) {
			List<ActiveSiteSignature> sList = new ArrayList<ActiveSiteSignature>(fragList);
			sList.add(sig.getFragmentAsSig(fragIndex));
			if (sigList.size() > 1) {
				getAllAlignments(sList, sigList.subList(1, sigList.size()));
			} else {
				ActiveSiteProfile prof = new ActiveSiteProfile(sList, new ClustalAlign(), radius);
				prof.updateAlignmentScore();
				updateMaps(prof, sList);
			}
		}
		return;
	}
*/

	public List<ActiveSiteProfile>getBestProfiles() {
		return bestProfiles();
	}

	private void updateMaps(ActiveSiteProfile prof, List<ActiveSiteSignature>tryList) {
		List<ActiveSiteProfile> profList = null;
/*
		System.out.println("\nAdding profile fragment: ");
		System.out.println(prof.getAlignmentAsString());
		System.out.println("  Score: "+prof.getScore());
		System.out.println("Using signatures: ");
		for (ActiveSiteSignature s: tryList) {
			System.out.println("   "+s.getSignature(radius));
		}
*/
		if (!orderedMap.containsKey(prof.getScore()))
			profList = new ArrayList<ActiveSiteProfile>();
		else
			profList = orderedMap.get(prof.getScore());

		profList.add(prof);
		orderedMap.put(prof.getScore(), profList);

		profileMap.put(prof, tryList);

		for (ActiveSiteSignature sig: tryList) {
			List<ActiveSiteProfile> pList = null;
			if (!signatureMap.containsKey(sig))
				pList = new ArrayList<ActiveSiteProfile>();
			else
				pList = signatureMap.get(sig);
			pList.add(prof);
			signatureMap.put(sig, pList);
		}
		
	}

	private List<ActiveSiteProfile> bestProfiles() {
		List<ActiveSiteProfile> profileList = new ArrayList<ActiveSiteProfile>();
		while (orderedMap.size() > 0) {
			// System.out.println("Map size = "+orderedMap.size());
			double bestScore = orderedMap.lastKey();
			ActiveSiteProfile bestProfile = orderedMap.get(bestScore).get(0);

			// System.out.println("Best profile is: ");
			// System.out.println(bestProfile.getAlignmentAsString());

			removeOverlaps(bestProfile, bestScore);
			profileList.add(bestProfile);
			// System.out.println("Map size = "+orderedMap.size());
		}
		return profileList;
	}

	private void removeOverlaps(ActiveSiteProfile profile, double bestScore) {
		List<ActiveSiteSignature> conflictList = profileMap.get(profile);

		for (ActiveSiteSignature sig: conflictList) {
			List<ActiveSiteProfile> conflictProf = signatureMap.get(sig);
			for (ActiveSiteProfile prof: conflictProf) {
				removeProfile(prof.getScore(), prof);
			}
			signatureMap.remove(sig);
		}
		removeProfile(bestScore, profile);
	}

	private void removeProfile(double score, ActiveSiteProfile profile) {
		// System.out.println("Removing profile: ");
		// System.out.println(profile.getAlignmentAsString());
		// System.out.println("Score = "+score);
		// Remove from orderedMap
		if (orderedMap.containsKey(score)) {
			List<ActiveSiteProfile> profList = orderedMap.get(score);
			if (profList.size() == 1 && profList.get(0) == profile) {
				// System.out.println("Removing last profile: ");
				// System.out.println(profList.get(0).getAlignmentAsString());
				// System.out.println("Score = "+profList.get(0).getScore());
				orderedMap.remove(score);
			} else {
				profList.remove(profile);
				orderedMap.put(score, profList);
			}
		}

		// Remove from profileMap
		if (profileMap.containsKey(profile))
			profileMap.remove(profile);
	}
}
