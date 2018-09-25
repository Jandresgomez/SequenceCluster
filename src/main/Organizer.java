package main;

import java.io.File;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;

import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

public class Organizer {
	private static final int KMER_LENGTH = 31;
	private static String SEPARATOR = "\\"; 
	private static final String RP = "report";
	private static final String TRASH = ".DS_Store";
	private static final int F_RECALL = 500000;
	private static final int SAMPLES = 1000;
	private String outPath;
	private PrintWriter reportFileWriter;
	private String inPath;
	private KMersTable table;
	private SequenceComparator sc;
	private File[] files;
	
	public Organizer(String inPath, String outPath) {
		sc = new SequenceComparator();
		this.outPath = outPath;
		this.inPath = inPath;
		table = new KMersTable();
		if(outPath.contains("/")) SEPARATOR = "/";
	}
	
	public void load() {
		File dir = new File(inPath);
		files = dir.listFiles();
	}
	
	public void report(String s) {
		reportFileWriter.println(s);
	}
	
	public void openReport(String fileName) throws Exception {
		String date = new SimpleDateFormat("dd_HHmm").format(new Date());
		reportFileWriter = new PrintWriter(new File(outPath + SEPARATOR + fileName + "_" +  date + ".txt"));
	}
	
	public void closeReport() {
		reportFileWriter.close();
	}
	
	public void createTable() throws Exception {
		load();
		int count = 0;
		
		for(File f : files) {
			if(!f.getName().equals(TRASH)) {
				System.out.println(f.getName());
				
				FastqFileReader openFile = new FastqFileReader(f);
				
				Iterator<RawRead> reader = openFile.iterator();
				while(reader.hasNext()) {
					RawRead read = reader.next();
					String s = read.getSequenceString();
					if(!s.contains("N")) table.add(new DNAShortKmer(s.substring(0, KMER_LENGTH)));
					count++;
					if(count%F_RECALL == 0) table.status();
				}  openFile.close();
				
				table.status();
				//report(f.getName());
			}
		}
		table.status();
	}

	public void sample() throws Exception {
		//openReport(RP);
		long time = System.currentTimeMillis();
		int[] freq = new int[KMER_LENGTH + 1];
		int count = 0;
		
		for(int k = 0; k < SAMPLES; k++) {
			DNAShortKmer sample = table.getRandom();
			Iterator<DNAShortKmer> it = table.getKeyIterator();
			while(it.hasNext()) {
				int i = sc.compare(sample, it.next());
				freq[i]++; count++;
				if(count%F_RECALL == 0) System.out.println(Arrays.toString(freq));
			}
		}
		
		
		report(Arrays.toString(freq));
		report((System.currentTimeMillis() - time) + "");
		closeReport();
		
	}
	public static void main(String[] args) throws Exception {
		Organizer organizer = new Organizer(args[0], args[1]);
		organizer.openReport(RP);
		organizer.createTable();
		organizer.sample();
	}
}
 