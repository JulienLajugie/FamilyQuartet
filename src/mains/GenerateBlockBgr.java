package mains;

import java.io.File;
import java.io.IOException;

import dataStructures.InheritanceStateBlockList;

/**
 * Prints the blocks in a bgr format
 * @author Julien Lajugie
 */
public class GenerateBlockBgr {
	
	/**
	 * Usage: java GenerateBlockBgr -b <path to the block file>
	 * @param args -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateBlockBgr.java -b <path to the block file>");
			System.exit(-1);
		} else {
			try {
				File blockFile = new File(args[1]);
				generateBlockBgr(blockFile);
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
		if (args.length != 2) {
			return false;
		}
		if (args[0].equals("-b") ) {
			return true;
		} else {
			return false;
		}
	}


	/**
	 * 
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @throws IOException if the block file is not valid
	 */
	private static void generateBlockBgr(File blockFile) throws IOException {
		InheritanceStateBlockList blockList = new InheritanceStateBlockList();
		blockList.loadFromISCAFile(blockFile);		
		blockList.printBlocksBgrFormat();
	}
}
