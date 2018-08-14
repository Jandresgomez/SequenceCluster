package main;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import ngsep.sequences.RawRead;

public class Organizer {
	private SequenceComparator sc;
	private File[] files;
	
	public static void main(String[] args) {
		
	}
	
	public Organizer() {
		sc = new SequenceComparator();
	}
	
	public void sort(String inPath, String outPath) {
		File dir = new File(inPath);
		files = dir.listFiles();
	}
}
