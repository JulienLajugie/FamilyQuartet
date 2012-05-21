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
	 * Creates an instance of {@link SegmentalDuplicationList} from the specified file
	 * @param file a bed or a begraph file
	 * @throws IOException 
	 */
	public SegmentalDuplicationList(File file) throws IOException {
		this.segDupListMap = createSegmentalDuplicationMap(file);
	}


	/**
	 * @param file bed or bgr file
	 * @return lists of segmental duplication organized by chromosome from the specified file
	 * @throws IOException
	 */
	private Map<String, List<SegmentalDuplication>> createSegmentalDuplicationMap(File file) throws IOException {
		BufferedReader reader = null;
		Map<String, List<SegmentalDuplication>> segDupListMap = new HashMap<String, List<SegmentalDuplication>>();
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
					// if the list doesn't contain the chromosome we add it
					if (!segDupListMap.containsKey(chromosome)) {
						List<SegmentalDuplication> listToAdd = new ArrayList<>();
						listToAdd.add(segDupToAdd);
						segDupListMap.put(chromosome, listToAdd);						
					} else {
						segDupListMap.get(chromosome).add(segDupToAdd);
					}
				}
			}
			for (List<SegmentalDuplication> currentList: segDupListMap.values()) {
				Collections.sort(currentList);
			}
			return segDupListMap;
		} finally {
			if (reader != null) {
				reader.close();
			}
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
						&& (position < currentDuplication.getStopPosition())) {
					return currentDuplication;
				}
			}
		}
		return null;
	}
}
