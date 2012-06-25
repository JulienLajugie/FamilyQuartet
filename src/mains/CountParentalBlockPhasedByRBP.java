package mains;

import java.io.File;
import java.io.IOException;
import java.util.List;

import dataStructures.QuartetMember;
import dataStructures.SegmentalDuplication;
import dataStructures.SegmentalDuplicationList;


/**
 * Counts the number of parental blocks that have been phased by RBP and prints the location of the phased junctions
 * @author Julien Lajugie
 */
public class CountParentalBlockPhasedByRBP {

	/**
	 * Usage: java CountParentalBlockPhasedByRBP.java -b <path to the inheritance state blocks> -p <path to the read backed phasing block file> -m <quartet member (FATHER, MOTHER)>
	 * @param args -b <path to the inheritance state blocks> -p <path to the read backed phasing block file> -m <quartet member (FATHER, MOTHER)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java CountParentalBlockPhasedByRBP.java -b <path to the inheritance state blocks> -p <path to the read backed phasing block file> -m <quartet member (FATHER, MOTHER)>");
			System.exit(-1);
		} else {
			QuartetMember member = null;
			File inheritanceStateFile = null;
			File RBPFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-m")) {
					member = QuartetMember.valueOf(args[i + 1]);
					if ((member == null) || (member != QuartetMember.FATHER) && (member != QuartetMember.MOTHER)) {
						System.out.println("Usage: java CountParentalBlockPhasedByRBP.java -b <path to the inheritance state blocks> -p <path to the read backed phasing block file> -m <quartet member (FATHER, MOTHER)>");
						System.exit(-1);						
					}
				}
				if (args[i].equals("-b")) {
					inheritanceStateFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					RBPFile = new File(args[i + 1]);
				}
			}
			try {
				countParentalBlockPhasedByRBP(inheritanceStateFile, RBPFile, member);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * @param args parameters from the main function
	 * @return true if the parameters are valid
	 */
	private static boolean areParametersValid(String[] args) {
		if (args == null) {
			return false;
		}
		if (args.length != 6) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b")) {
			return false;
		}
		// case with no -p parameter
		if (!args[0].equals("-p") && !args[2].equals("-p") && !args[4].equals("-p")) {
			return false;
		}
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m") && !args[4].equals("-m")) {
			return false;
		}
		return true;
	}


	/**
	 * Counts the number of parental blocks that have been phased by RBP and prints the location of the phased junctions
	 * @param inheritanceStateFile file with the inheritance blocks of the specified member
	 * @param RBPFile bgr file with the blocks phased by RBP
	 * @param member founder member to study
	 * @throws IOException
	 */
	private static void countParentalBlockPhasedByRBP(File inheritanceStateFile, File RBPFile, QuartetMember member) throws IOException {
		SegmentalDuplicationList inheritanceBlocks = new SegmentalDuplicationList();
		inheritanceBlocks.loadBedOrBgr(inheritanceStateFile);
		SegmentalDuplicationList readBackedBlocks = new SegmentalDuplicationList();
		readBackedBlocks.loadBedOrBgr(RBPFile);
		int phasedBlockCount = 0;
		for (String currentChromosome: inheritanceBlocks.getBlocks().keySet()) {
			List<SegmentalDuplication> currentInheritanceBlockList = inheritanceBlocks.getBlocks().get(currentChromosome);
			for (int i = 0; i < currentInheritanceBlockList.size() - 1; i++) {
				int junctionStart = currentInheritanceBlockList.get(i).getStopPosition();
				int junctionStop =  currentInheritanceBlockList.get(i + 1).getStartPosition();
				SegmentalDuplication RBPBlockJunctionStart = readBackedBlocks.getBlock(currentChromosome, junctionStart);
				SegmentalDuplication RBPBlockJunctionStop = readBackedBlocks.getBlock(currentChromosome, junctionStop);
				if ((RBPBlockJunctionStart != null) && (RBPBlockJunctionStop != null) && (RBPBlockJunctionStart == RBPBlockJunctionStop)) {
					System.out.println(currentChromosome + ":" + junctionStart + "-" + junctionStop);
					phasedBlockCount++;
				}				
			}			
		}
		System.out.println("Phased Block Count: " + phasedBlockCount);
	}
}
