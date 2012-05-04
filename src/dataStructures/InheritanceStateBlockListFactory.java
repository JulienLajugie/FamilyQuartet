package dataStructures;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Factory that creates {@link InheritanceStateBlockList} instances
 * @author Julien Lajugie
 */
public class InheritanceStateBlockListFactory {
	
	
	/**
	 * Creates an {@link InheritanceStateBlockList} from a smoothed_blocks_with_intercalated_partial_blocks file
	 * This file is an output file from the software ISCA (0.1.6) - Jared Roach <jroach@systemsbiology.org>
	 * @param iscaBlockFile output file from the software ISCA
	 * @return an {@link InheritanceStateBlockList} loaded from the ISCA file
	 * @throws IOException
	 */
	public static InheritanceStateBlockList<QuartetInheritanceState> createFromISCAFile(File iscaBlockFile) throws IOException {
		BufferedReader reader = null;
		Map<String, List<InheritanceStateBlock<QuartetInheritanceState>>> ISBlockMap = new HashMap<String, List<InheritanceStateBlock<QuartetInheritanceState>>>();
		try {
			reader = new BufferedReader(new FileReader(iscaBlockFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't care about the comment lines 
				if (line.trim().charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					String binaryState = splitLine[3].trim();
					int startPosition = Integer.parseInt(splitLine[5].trim());
					int stopPosition = Integer.parseInt(splitLine[6].trim());
					QuartetInheritanceState blockState = QuartetInheritanceState.valueOfBinaryState(binaryState);
					InheritanceStateBlock<QuartetInheritanceState> blockToAdd = new QuartetInheritanceStateBlock(blockState, chromosome, startPosition, stopPosition);
					// if the list doesn't contain the chromosome we add it
					if (!ISBlockMap.containsKey(chromosome)) {
						List<InheritanceStateBlock<QuartetInheritanceState>> listToAdd = new ArrayList<>();
						listToAdd.add(blockToAdd);
						ISBlockMap.put(chromosome, listToAdd);
					} else {
						ISBlockMap.get(chromosome).add(blockToAdd);
					}
				}
			}
			return new InheritanceStateBlockList<>(ISBlockMap);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Creates an {@link InheritanceStateBlockList} from a bedgraph file
	 * @param bgrFile bedgraph file
	 * @return an {@link InheritanceStateBlockList} loaded from a quartet bgr file
	 * @throws IOException
	 */
	public static InheritanceStateBlockList<QuartetInheritanceState> createFromQuartetBgrFile(File bgrFile) throws IOException {
		BufferedReader reader = null;
		Map<String, List<InheritanceStateBlock<QuartetInheritanceState>>> ISBlockMap = new HashMap<String, List<InheritanceStateBlock<QuartetInheritanceState>>>();
		try {
			reader = new BufferedReader(new FileReader(bgrFile));
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
					int stateScore = (int) Double.parseDouble(splitLine[3].trim());
					QuartetInheritanceState blockState = QuartetInheritanceState.valueOf(stateScore);
					InheritanceStateBlock<QuartetInheritanceState> blockToAdd = new QuartetInheritanceStateBlock(blockState, chromosome, startPosition, stopPosition);
					// if the list doesn't contain the chromosome we add it
					if (!ISBlockMap.containsKey(chromosome)) {
						List<InheritanceStateBlock<QuartetInheritanceState>> listToAdd = new ArrayList<>();
						listToAdd.add(blockToAdd);
						ISBlockMap.put(chromosome, listToAdd);						
					} else {
						ISBlockMap.get(chromosome).add(blockToAdd);
					}
				}
			}	
			return new InheritanceStateBlockList<>(ISBlockMap);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Creates an {@link InheritanceStateBlockList} from a bedgraph file
	 * @param bgrFile bedgraph file
	 * @return an {@link InheritanceStateBlockList} loaded from a cross trio bgr file
	 * @throws IOException
	 */
	public static InheritanceStateBlockList<CrossTriosInheritanceState> createFromCrossTriosBgrFile(File bgrFile) throws IOException {
		BufferedReader reader = null;
		Map<String, List<InheritanceStateBlock<CrossTriosInheritanceState>>> ISBlockMap = new HashMap<String, List<InheritanceStateBlock<CrossTriosInheritanceState>>>();
		try {
			reader = new BufferedReader(new FileReader(bgrFile));
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
					int stateScore = (int) Double.parseDouble(splitLine[3].trim());
					CrossTriosInheritanceState blockState = CrossTriosInheritanceState.valueOf(stateScore);
					InheritanceStateBlock<CrossTriosInheritanceState> blockToAdd = new CrossTriosInheritanceStateBlock(blockState, chromosome, startPosition, stopPosition);
					// if the list doesn't contain the chromosome we add it
					if (!ISBlockMap.containsKey(chromosome)) {
						List<InheritanceStateBlock<CrossTriosInheritanceState>> listToAdd = new ArrayList<>();
						listToAdd.add(blockToAdd);
						ISBlockMap.put(chromosome, listToAdd);						
					} else {
						ISBlockMap.get(chromosome).add(blockToAdd);
					}
				}
			}	
			return new InheritanceStateBlockList<>(ISBlockMap);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
