package main;

import java.io.File;
import java.util.Iterator;

import ngsep.sequences.DNASequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

public class Organizer {
	private SequenceComparator sc;
	private File[] files;
	
	public static void main(String[] args) {
		Organizer organizer = new Organizer();
		organizer.sort(args[0], args[1]);
	}
	
	public Organizer() {
		sc = new SequenceComparator();
	}
	
	public void sort(String inPath, String outPath) {
		File dir = new File(inPath);
		files = dir.listFiles();
		
		for(File f : files) {
			try {
				Iterator<RawRead> reader = (new FastqFileReader(f)).iterator();
				while(reader.hasNext()) {
					DNASequence seq = new DNASequence(reader.next().getSequenceString());
					System.out.println(f.getName() + "__" + seq.toString());
				}
			} catch (Exception e) { e.printStackTrace(); }
		}
	}
}
 