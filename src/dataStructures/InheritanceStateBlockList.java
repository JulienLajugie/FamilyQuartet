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
 * This class represent a list of {@link InheritanceStateBlock}
 * @author Julien Lajugie
 *
 */
public class InheritanceStateBlockList {


	private final Map<String, List<InheritanceStateBlock>> ISBlockMap; // list of blocks organised by chromosome


	/**
	 * Creates an instance of {@link InheritanceStateBlockList}
	 */
	public InheritanceStateBlockList() {
		ISBlockMap = new HashMap<String, List<InheritanceStateBlock>>();
	}


	/**
	 * Creates an {@link InheritanceStateBlockList} from a smoothed_blocks_with_intercalated_partial_blocks file
	 * This file is an output file from the software ISCA (0.1.6) - Jared Roach <jroach@systemsbiology.org>
	 * @param iscaBlockFile output file from the software ISCA
	 * @throws IOException
	 */
	public void loadFromISCAFile(File iscaBlockFile) throws IOException {
		BufferedReader reader = null;
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
					InheritanceState blockState = binaryStateToBlockState(binaryState);
					InheritanceStateBlock blockToAdd = new InheritanceStateBlock(blockState, chromosome, startPosition, stopPosition);

					// if the list doesn't contain the chromosome we add it
					if (!ISBlockMap.containsKey(chromosome)) {
						List<InheritanceStateBlock> listToAdd = new ArrayList<InheritanceStateBlock>();
						listToAdd.add(blockToAdd);
						ISBlockMap.put(chromosome, listToAdd);						
					} else {
						ISBlockMap.get(chromosome).add(blockToAdd);
					}
				}
			}			
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Creates an {@link InheritanceStateBlockList} from a bedgraph file
	 * @param bgrFile bedgraph file
	 * @throws IOException
	 */
	public void loadFromBgrFile(File bgrFile) throws IOException {
		BufferedReader reader = null;
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
					int binaryState = (int) Double.parseDouble(splitLine[3].trim());
					InheritanceState blockState = intStateToBlockState(binaryState);
					InheritanceStateBlock blockToAdd = new InheritanceStateBlock(blockState, chromosome, startPosition, stopPosition);
					// if the list doesn't contain the chromosome we add it
					if (!ISBlockMap.containsKey(chromosome)) {
						List<InheritanceStateBlock> listToAdd = new ArrayList<InheritanceStateBlock>();
						listToAdd.add(blockToAdd);
						ISBlockMap.put(chromosome, listToAdd);						
					} else {
						ISBlockMap.get(chromosome).add(blockToAdd);
					}
				}
			}			
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Analyzes the current variant to generate stats on the number of MIE / SCE / Not informative variant per block
	 * @param variant variant to analyze
	 */
	public void analyzeVariant(Variant variant) {
		InheritanceStateBlock variantBlock = getBlock(variant);
		if (variantBlock != null) {
			variantBlock.analyzeVariant(variant);
		}
	}


	/**
	 * This methods prints in bgr format in the standard output the coordinates of the specified variant only if it's a SCE 
	 * @param variant a {@link Variant}
	 */
	public void printSCEVariantBrgFormat(Variant variant) {
		InheritanceStateBlock variantBlock = getBlock(variant);
		if (variantBlock != null) {
			variantBlock.printSCEVariantBrgFormat(variant);
		}
	}


	/**
	 * This methods prints in vcf format in the standard output the vcf info of the specified variant only if it's a SCE 
	 * @param variant a {@link Variant}
	 * @param vcfLine line from the vcf from which the variant was extracted
	 */
	public void printSCEVariantVCFFormat(Variant variant, String vcfLine) {
		InheritanceStateBlock variantBlock = getBlock(variant);
		if (variantBlock != null) {
			variantBlock.printSCEVariantVcfFormat(variant, vcfLine);
		}
	}


	/**
	 * @param variant a Variant
	 * @return the {@link InheritanceStateBlock} containing the specified variant. Null if there is no such block
	 */
	public InheritanceStateBlock getBlock(Variant variant) {
		return getBlock(variant.getChromosome(), variant.getPosition());
	}

	
	/**
	 * @param chromosome a chromosome
	 * @param position a position on the specified chromosome
	 * @return a block on the specified chromosome containing the specified position.  Null if there is no such block
	 */
	public InheritanceStateBlock getBlock(String chromosome, int position) {
		if (ISBlockMap.containsKey(chromosome)) {
			List<InheritanceStateBlock> chromosomeBlockList = ISBlockMap.get(chromosome);
			for (InheritanceStateBlock currentBlock: chromosomeBlockList) {
				if ((position >= currentBlock.getStartPosition()) 
						&& (position < currentBlock.getStopPosition())) {
					return currentBlock;
				}
			}
		}
		return null;
	}

	
	/**
	 * Prints the statistics about the blocks
	 */
	public void printStatistics() {
		int variantCount = 0;
		int MIECount = 0;
		int SCECount = 0;
		int NICount = 0;
		System.out.println(InheritanceStateBlock.getStatisticsHeader());
		for (List<InheritanceStateBlock> currentList: ISBlockMap.values()) {
			for (InheritanceStateBlock currentBlock: currentList) {
				System.out.println(currentBlock.getStatistics());
				variantCount += currentBlock.getVariantCount();
				MIECount += currentBlock.getMIECount();
				SCECount += currentBlock.getSCECount();
				NICount += currentBlock.getNICount();				
			}
		}
		// print the genome wide stats
		String genomeWideStats = "GW\t";
		genomeWideStats += "GW\t";
		genomeWideStats += "GW\t";
		genomeWideStats += "NA\t";
		genomeWideStats += variantCount +"\t";
		genomeWideStats += MIECount +"\t";
		genomeWideStats += MIECount / (double) variantCount * 100d +"\t";
		genomeWideStats += SCECount +"\t";
		genomeWideStats += SCECount / (double) variantCount * 100d +"\t";
		genomeWideStats += NICount +"\t";
		genomeWideStats += NICount / (double) variantCount * 100d +"\t";
		System.out.println(genomeWideStats);		
	}


	/**
	 * @param intState integer representing an inheritance state as follow:
	 * -1 -> MIE block
	 * 0 -> partial block
	 * 1 -> Non-Identical block
	 * 2 -> Paternal block
	 * 3 -> Maternal block
	 * 4 -> Identical block
	 * 5 -> Not Informative block
	 * @return the {@link InheritanceState} associated to the specified integer
	 */
	private static InheritanceState intStateToBlockState(int intState) {
		switch (intState) {
		case -1:
			return InheritanceState.MIE;
		case 1:
			return InheritanceState.NON_IDENTICAL;
		case 2: 
			return InheritanceState.PATERNAL;
		case 3:
			return InheritanceState.MATERNAL;
		case 4: 
			return InheritanceState.IDENTICAL;
		case 5:
			return InheritanceState.NOT_INFORMATIVE;
		default:
			return null;			
		}
	}


	/**
	 * @param blockState an {@link InheritanceState}
	 * @return an int representing the specified inheritance state as follow:
	 * -1 -> MIE block
	 * 0 -> partial block
	 * 1 -> Non-Identical block
	 * 2 -> Paternal block
	 * 3 -> Maternal block
	 * 4 -> Identical block
	 * 5 -> Not Informative block
	 */
	private static int blockStateToInt(InheritanceState blockState) {
		switch (blockState) {
		case MIE:
			return -1;
		case NON_IDENTICAL:
			return 1;
		case PATERNAL:
			return 2;
		case MATERNAL:
			return 3;
		case IDENTICAL:
			return 4;
		case NOT_INFORMATIVE:
			return 5;
		default:
			return 0;
		}
	}


	/**
	 * @param binaryState a binary representation of an {@link InheritanceState}
	 * @return the {@link InheritanceState} of the binary state
	 */
	private static InheritanceState binaryStateToBlockState(String binaryState) {
		if (binaryState.equals("0000")) {
			return InheritanceState.IDENTICAL;
		}
		if (binaryState.equals("0001")) {
			return InheritanceState.PATERNAL;
		}
		if (binaryState.equals("0010")) {
			return InheritanceState.MATERNAL;
		}
		if (binaryState.equals("0011")) {
			return InheritanceState.NON_IDENTICAL;
		}
		return null;
	}


	/**
	 * Prints the blocks in a bgr format.
	 * The score works as follow:
	 * -1 -> MIE block
	 * 0 -> partial block
	 * 1 -> Non-Identical block
	 * 2 -> Paternal block
	 * 3 -> Maternal block
	 * 4 -> Identical block
	 * 5 -> Not Informative block
	 */
	public void printBlocksBgrFormat() {
		for (List<InheritanceStateBlock> currentList: ISBlockMap.values()) {
			for (InheritanceStateBlock currentBlock: currentList) {
				if (!currentBlock.isPartial()) {
					int score = blockStateToInt(currentBlock.getBlockState());
					System.out.println(currentBlock.getChromosome() + "\t" + currentBlock.getStartPosition() + "\t" + currentBlock.getStopPosition() + "\t " + score);
				}
			}
		}
	}
	
	
	/**
	 * @return the inheritant state blocks organized in a map sorted per chromosome
	 */
	public final Map<String, List<InheritanceStateBlock>> getBlocks() {
		return ISBlockMap;
	}
}
