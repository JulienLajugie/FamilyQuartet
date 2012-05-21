package dataStructures;

import java.util.List;
import java.util.Map;


/**
 * This class represent a list of {@link InheritanceStateBlock}
 * @author Julien Lajugie
 * @param <T> type of inheritance state block (must implement the interface {@link InheritanceState})
 */
public class InheritanceStateBlockList<T extends InheritanceState> {


	private final Map<String, List<InheritanceStateBlock<T>>> ISBlockMap; // list of blocks organised by chromosome


	/**
	 * Creates an instance of {@link InheritanceStateBlockList}
	 * @param ISBlockMap Inheritance State Blocks in a Map of association Chromosome -> List of InheritanceStateBlock
	 */
	protected InheritanceStateBlockList(Map<String, List<InheritanceStateBlock<T>>> ISBlockMap) {
		this.ISBlockMap = ISBlockMap;
	}


	/**
	 * @param variant a Variant
	 * @return the {@link InheritanceStateBlock} containing the specified variant. Null if there is no such block
	 */
	public InheritanceStateBlock<T> getBlock(Variant variant) {
		return getBlock(variant.getChromosome(), variant.getPosition());
	}

	
	/**
	 * @param chromosome a chromosome
	 * @param position a position on the specified chromosome
	 * @return a block on the specified chromosome containing the specified position.  Null if there is no such block
	 */
	public InheritanceStateBlock<T> getBlock(String chromosome, int position) {
		if (ISBlockMap.containsKey(chromosome)) {
			List<InheritanceStateBlock<T>> chromosomeBlockList = ISBlockMap.get(chromosome);
			for (InheritanceStateBlock<T> currentBlock: chromosomeBlockList) {
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
		System.out.println(InheritanceStateBlock.STAT_HEADER);
		for (List<InheritanceStateBlock<T>> currentList: ISBlockMap.values()) {
			for (InheritanceStateBlock<T> currentBlock: currentList) {
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
	 * Prints the genome wide statistics only
	 */
	public void printGenomeWideStatistics() {
		int variantCount = 0;
		int MIECount = 0;
		int SCECount = 0;
		int NICount = 0;
		System.out.println(InheritanceStateBlock.STAT_HEADER);
		for (List<InheritanceStateBlock<T>> currentList: ISBlockMap.values()) {
			for (InheritanceStateBlock<T> currentBlock: currentList) {
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
		for (List<InheritanceStateBlock<T>> currentList: ISBlockMap.values()) {
			for (InheritanceStateBlock<T> currentBlock: currentList) {
				int score = currentBlock.getBlockState().getScore();
				System.out.println(currentBlock.getChromosome() + "\t" + currentBlock.getStartPosition() + "\t" + currentBlock.getStopPosition() + "\t " + score);
			}
		}
	}
	
	
	/**
	 * @return the inheritant state blocks organized in a map sorted per chromosome
	 */
	public final Map<String, List<InheritanceStateBlock<T>>> getBlocks() {
		return ISBlockMap;
	}
}
