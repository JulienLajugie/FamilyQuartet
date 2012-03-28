package mains;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Creates inheritance block in a family trio (founder / child1 / child2)
 * The input file is a bgr file with a score of: 
 *  * 1 when the children received the same chromosome from the founder
 *  * -1 when the children received a different chromosome from the founder
 * @author Julien Lajugie
 */
public class CreateBlocksFromTrioBgr {

	/**
	 * Size of the windows used to determine the blocks
	 * When a variant is studied, the algorithm is going to look at the n previous variants
	 * and the n next variants to determine if the variant is a boundary of a block
	 */
	private final static int WINDOW_SIZE = 50;

	/**
	 * Number of differences needed between the state of the variants before and after 
	 * the studied variant to consider it as potentially at the boundary of a block
	 */
	private final static int SUM_SCORE_THRESHOLD = 8; 


	/**
	 * Usage: java CreateBlocksFromTrioBgr -f <path to the file>
	 * The input file is a bgr file with a score of: 
	 *  * 1 when the children received the same chromosome from the founder
	 *  * -1 when the children received a different chromosome from the founder
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java CreateBlocksFromTrioBgr -f <path to the file>");
			System.exit(-1);
		} else {
			try {
				createBlocksFromTrioBgr(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Creates inheritance block in a family trio (founder / child1 / child2)
	 * @param bgrFile input file bgr file with a score of: 
	 *  * 1 when the children received the same chromosome from the founder
	 *  * -1 when the children received a different chromosome from the founder
	 * @throws IOException
	 */
	private static void createBlocksFromTrioBgr(String bgrFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(bgrFile));
			String line = null;
			String previousChromo = null;
			List<Integer> positionList = new ArrayList<>();
			List<Integer> scoreList = new ArrayList<>();			
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromo = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					int score = (int) Double.parseDouble(splitLine[3].trim());
					if (previousChromo == null) {
						previousChromo = chromo;
					}
					if (!chromo.equals(previousChromo)) {
						computeCurrentChromoBlocks(previousChromo, positionList, scoreList);
						positionList.clear();
						scoreList.clear();
						previousChromo = chromo;
					} else {
						positionList.add(position);
						scoreList.add(score);
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
	 * Defines the blocks for a chromosome and print the result as a bgr
	 * @param chromo studied chromosome
	 * @param positionList list of the position of the variants on the chromosome
	 * @param scoreList list of the score of the variants 
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 */
	private static void computeCurrentChromoBlocks(String chromo, List<Integer> positionList, List<Integer> scoreList) {
		List<Integer> summitStartIndexList = new ArrayList<>();
		List<Integer> summitStopIndexList = new ArrayList<>();
		boolean inASummit = false;
		for (int i = WINDOW_SIZE; i < scoreList.size() - WINDOW_SIZE; i++) {
			int scoreBefore = 0;
			int scoreAfter = 0;
			for (int j = i - WINDOW_SIZE; j < i; j++) {
				scoreBefore += scoreList.get(j);
			}
			for (int j = i + WINDOW_SIZE; j > i; j--) {
				scoreAfter += scoreList.get(j);
			}
			if ((Math.abs(scoreBefore + scoreAfter) < SUM_SCORE_THRESHOLD) && (!inASummit)) {
				inASummit = true;				
				summitStartIndexList.add(i);
			}
			if ((Math.abs(scoreBefore + scoreAfter) >= SUM_SCORE_THRESHOLD) && (inASummit)) {
				inASummit = false;				
				summitStopIndexList.add(i);
			}
		}
		if (!summitStartIndexList.isEmpty() && !summitStopIndexList.isEmpty()) {
			int blockType = findBlockType(scoreList, 0, summitStartIndexList.get(0));
			int indexBlockStart = findFirstBlockStartIndex(scoreList, blockType);
			int indexBlockStop = 0;
			for (int i = 0; i < summitStartIndexList.size() && i < summitStopIndexList.size(); i++) {
				indexBlockStop = findIndexBlockStop(scoreList, summitStartIndexList.get(i), summitStopIndexList.get(i));
				blockType = findBlockType(scoreList, indexBlockStart, indexBlockStop);
				if (blockType == -1) {
					System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1)+ '\t' + 0.5);
				} else {
					System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1)+ '\t' + 1.5);
				}
				//System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1)+ '\t' + blockType);
				indexBlockStart = findIndexBlockStart(scoreList, summitStartIndexList.get(i), summitStopIndexList.get(i));			
			}
			blockType = findBlockType(scoreList, indexBlockStart, scoreList.size());
			indexBlockStop = findLastBlockStopIndex(scoreList, blockType);
			if (blockType == -1) {
				System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1)+ '\t' + 0.5);
			} else {
				System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1)+ '\t' + 1.5);
			}
			//System.out.println(chromo + '\t' + positionList.get(indexBlockStart) + '\t' + (positionList.get(indexBlockStop) + 1) + '\t' + blockType);
		}
	}


	/**
	 * @param scoreList list of variant scores
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @param indexSummitStart the index of the beginning of the boundary region
	 * @param indexSmmitStop the index of the end of the of the boundary region
	 * @return the index of the first variant of a block
	 */
	private static int findIndexBlockStart(List<Integer> scoreList,	int indexSummitStart, int indexSmmitStop) {
		List<Integer> summitScoreList = new ArrayList<>();
		for (int i = indexSummitStart; i <= indexSmmitStop; i++) {
			int summitScore = 0; 
			for (int j = i; j <= i + WINDOW_SIZE; j++) {
				summitScore += scoreList.get(j);
			}
			summitScoreList.add(Math.abs(summitScore));			
		}
		int indexMax = summitScoreList.size() - 1;
		for (int i = summitScoreList.size() - 2; i >= 0; i--) {
			if (summitScoreList.get(i) >= summitScoreList.get(indexMax)) {
				indexMax = i;
			}
		}		
		return indexSummitStart + indexMax;
	}


	/**
	 * @param scoreList list of variant scores
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @param indexSummitStart the index of the beginning of the boundary region
	 * @param indexSmmitStop the index of the end of the of the boundary region
	 * @return the index of the last variant of a block
	 */
	private static int findIndexBlockStop(List<Integer> scoreList, int indexSummitStart, int indexSmmitStop) {
		List<Integer> summitScoreList = new ArrayList<>();
		for (int i = indexSummitStart; i <= indexSmmitStop; i++) {
			int summitScore = 0; 
			for (int j = i - WINDOW_SIZE; j <= i; j++) {
				summitScore += scoreList.get(j);
			}
			summitScoreList.add(Math.abs(summitScore));			
		}
		int indexMax = 0;
		for (int i = 1; i < summitScoreList.size(); i++) {
			if (summitScoreList.get(i) >= summitScoreList.get(indexMax)) {
				indexMax = i;
			}
		}		
		return indexSummitStart + indexMax;
	}


	/**
	 * Finds the type of the block 
	 * @param scoreList score of the variants
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @param startIndex index of the first variant of the block
	 * @param stopIndex index of the last variant of the block
	 * @return 1 if the children received the same allele from the founder, -1 otherwise
	 */
	private static int findBlockType(List<Integer> scoreList, int startIndex, int stopIndex) {
		int score = 0;
		for (int i = startIndex; i < stopIndex; i++) {
			score += scoreList.get(i);			
		}
		if (score < 0) {
			return -1;
		} else {
			return 1;
		}
	}


	/** 
	 * @param scoreList list of the scores of the variants
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @param blockType type of the first block
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @return the index of the first variant of the first block of the chromosome
	 */
	private static int findFirstBlockStartIndex(List<Integer> scoreList, int blockType) {
		for (int i = 0; i < scoreList.size(); i++) {
			if (scoreList.get(i) == blockType) {
				return i;
			}
		}
		return -1;
	}


	/** 
	 * @param scoreList list of the scores of the variants
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @param blockType type of the last block
	 * (1 if the children received the same allele from the founder)
	 * (-1 if the children received different alleles from the founder)
	 * @return the index of the last variant of the last block of the chromosome
	 */
	private static int findLastBlockStopIndex(List<Integer> scoreList, int blockType) {
		for (int i = scoreList.size() - 1; i >= 0; i--) {
			if (scoreList.get(i) == blockType) {
				return i;
			}
		}
		return -1;
	}	
}
