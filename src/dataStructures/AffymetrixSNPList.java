package dataStructures;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * List of {@link AffymetrixSNPList} organized in a map indexed by chromosomes 
 * @author Julien Lajugie
 */
public class AffymetrixSNPList {

	private final Map<String, List<AffymetrixSNP>> affySNPMap;


	/**
	 * Creates an instance of {@link AffymetrixSNPList}
	 */
	public AffymetrixSNPList() {
		this.affySNPMap = new HashMap<>();
	}


	/**
	 * Load the list from an affymetrix file
	 * @param affyFile
	 * @throws IOException
	 */
	public void loadAffymetrixFile(File affyFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(affyFile));
			String line = null;
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					AffymetrixSNP snpToAdd = AffymetrixSNP.createFromAffymetrixLine(line);
					if (snpToAdd != null) {
						addSNP(snpToAdd);
					}
				}
			}
			sort();
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}

	
	/**
	 * Adds a specified {@link AffymetrixSNP} to the map
	 * @param snpToAdd {@link AffymetrixSNP}
	 */
	public void addSNP(AffymetrixSNP snpToAdd) {
		// if the list doesn't contain the chromosome we add it
		if (!affySNPMap.containsKey(snpToAdd.getChromosome())) {
			List<AffymetrixSNP> listToAdd = new ArrayList<>();
			listToAdd.add(snpToAdd);
			affySNPMap.put(snpToAdd.getChromosome(), listToAdd);						
		} else {
			affySNPMap.get(snpToAdd.getChromosome()).add(snpToAdd);
		}
	}
	
	
	/**
	 * @param chromosome
	 * @param position
	 * @return the {@link AffymetrixSNP} on the specified chromosome, at the specified location if it exists.  Return null otherwise
	 */
	public AffymetrixSNP get(String chromosome, int position) {
		if (!affySNPMap.containsKey(chromosome)) {
			return null;
		}
		List<AffymetrixSNP> snpList = affySNPMap.get(chromosome);
		AffymetrixSNP snpToSearch = new AffymetrixSNP(chromosome, position, "AA");
		int index = Collections.binarySearch(snpList, snpToSearch);
		if (index >= 0) {
			return snpList.get(index);
		} else {
			return null;
		}		
	}
	
	
	/**
	 * Sort the lists of {@link AffymetrixSNP}
	 */
	public void sort() {
		for (List<AffymetrixSNP> list: affySNPMap.values()) {
			Collections.sort(list);
		}
	}
	
	
	/**
	 * @return the total number of SNPs
	 */
	public int SNPCount() {
		int count = 0;
		for (List<AffymetrixSNP> list: affySNPMap.values()) {
			count += list.size();
		}
		return count;
	}
}
