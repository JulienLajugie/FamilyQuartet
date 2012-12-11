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
 * This class represents lists of {@link SegmentalDuplication} organised by chromosome 
 * @author Julien Lajugie
 */
public class SegmentalDuplicationList {

	private final Map<String, List<SegmentalDuplication>> segDupListMap; // lists of segmental duplications organised by chromosome


	/**
	 * Creates an instance of {@link SegmentalDuplicationList} 
	 */
	public SegmentalDuplicationList() {
		this.segDupListMap = new HashMap<String, List<SegmentalDuplication>>();
	}


	/**
	 * @param file bed or bgr file
	 * @throws IOException
	 */
	public void loadBedOrBgr(File file) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't care about the comment lines 
				if ((line.trim().charAt(0) != '#') && 
						((line.length() < 5) || (!line.substring(0, 5).equalsIgnoreCase("track"))) && 
						((line.length() < 7) || (!line.substring(0, 7).equalsIgnoreCase("browser")))) {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int startPosition = Integer.parseInt(splitLine[1].trim());
					int stopPosition = Integer.parseInt(splitLine[2].trim());
					SegmentalDuplication segDupToAdd = new SegmentalDuplication(startPosition, stopPosition);
					addDuplication(chromosome, segDupToAdd);
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
	 * @param file bed or bgr file
	 * @throws IOException
	 */
	public void loadBedOrBgrWithScore(File file) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't care about the comment lines 
				if ((line.trim().charAt(0) != '#') && 
						((line.length() < 5) || (!line.substring(0, 5).equalsIgnoreCase("track"))) && 
						((line.length() < 7) || (!line.substring(0, 7).equalsIgnoreCase("browser")))) {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int startPosition = Integer.parseInt(splitLine[1].trim());
					int stopPosition = Integer.parseInt(splitLine[2].trim());
					double score = Double.parseDouble(splitLine[3].trim());
					SegmentalDuplication segDupToAdd = new SegmentalDuplication(startPosition, stopPosition, score);
					addDuplication(chromosome, segDupToAdd);
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
	 * Adds a segmental duplication to the map
	 * @param chromosome chromosome of the duplication
	 * @param duplicationToAdd
	 */
	public void addDuplication(String chromosome, SegmentalDuplication duplicationToAdd) {
		// if the list doesn't contain the chromosome we add it
		if (!segDupListMap.containsKey(chromosome)) {
			List<SegmentalDuplication> listToAdd = new ArrayList<>();
			listToAdd.add(duplicationToAdd);
			segDupListMap.put(chromosome, listToAdd);						
		} else {
			segDupListMap.get(chromosome).add(duplicationToAdd);
		}
	}


	/**
	 * Sorts the lists of block in position order
	 */
	public void sort() {
		for (List<SegmentalDuplication> currentList: segDupListMap.values()) {
			Collections.sort(currentList);
		}
	}


	/**
	 * @param variant a {@link Variant}
	 * @return true if the variant is located inside a segmental duplication. False otherwise
	 */
	public boolean isInSegmentalDuplication(Variant variant) {
		return (getBlock(variant.getChromosome(), variant.getPosition()) != null);
	}


	/**
	 * @param chromosome a chromosome
	 * @param position a position
	 * @return the segmental duplication that contains the specified chromosome and position
	 */
	public SegmentalDuplication getBlock(String chromosome, int position) {
		if (segDupListMap.containsKey(chromosome)) {
			List<SegmentalDuplication> chromosomeBlockList = segDupListMap.get(chromosome);
			for (SegmentalDuplication currentDuplication: chromosomeBlockList) {
				if ((position >= currentDuplication.getStartPosition()) 
						&& (position <= currentDuplication.getStopPosition())) {
					return currentDuplication;
				}
			}
		}
		return null;
	}


	/**
	 * @param chromosome a chromosome
	 * @param block a {@link SegmentalDuplication}
	 * @return true if the list contains a bock that overlap with the specified {@link SegmentalDuplication} on the specified chromosome
	 */
	public boolean containsBlockOverlapping(String chromosome, SegmentalDuplication block) {
		return getBlockOverlapping(chromosome, block) != null;
	}

	
	/**
	 * @param chromosome a chromosome
	 * @param block a {@link SegmentalDuplication}
	 * @return a block that overlap with the specified {@link SegmentalDuplication} on the specified chromosome. Null if there is none
	 */
	public SegmentalDuplication getBlockOverlapping(String chromosome, SegmentalDuplication block) {
		if (segDupListMap.containsKey(chromosome)) {
			List<SegmentalDuplication> chromosomeBlockList = segDupListMap.get(chromosome);
			for (SegmentalDuplication currentDuplication: chromosomeBlockList) {
				if (((block.getStartPosition() >= currentDuplication.getStartPosition()) && (block.getStartPosition() <= currentDuplication.getStopPosition()))
						|| ((block.getStopPosition() >= currentDuplication.getStartPosition()) && (block.getStopPosition() <= currentDuplication.getStopPosition()))
						|| ((block.getStartPosition() <= currentDuplication.getStartPosition()) && (block.getStopPosition() >= currentDuplication.getStopPosition()))) {
					return currentDuplication;
				}
			}
		}
		return null;
	}	
	
	

	/**
	 * @return the segmental duplication blocks organized in a map sorted per chromosome
	 */
	public final Map<String, List<SegmentalDuplication>> getBlocks() {
		return segDupListMap;
	}
}
